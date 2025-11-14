#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ASE (Allele-Specific Expression) pipeline:
- Parse gene lengths from GFF3 (sum of merged exon lengths per gene)
- Read allele-paired genes and haplotype-specific count matrices
- Convert counts to TPM and filter expressed genes
- Prepare DESeq2 input per treatment (hap1 vs hap2 within each sample)
- Run DESeq2 (allele effect: hap2 vs hap1) and collect results
- Identify significant ASE gene pairs per treatment

Outputs are written under OUTPUT_DIR.
"""

import os
import pandas as pd
import numpy as np
import subprocess
from collections import defaultdict
from tqdm import tqdm

# ----------------------------- Configuration -----------------------------
PAIRS_FILE   = "last.txt"                          # gene pair table (columns: hap1_gene, hap2_gene)
HAP1_COUNTS  = "hap1.txt"                          # counts matrix (rows: GeneID, columns: samples)
HAP2_COUNTS  = "hap2.txt"                          # counts matrix (rows: GeneID, columns: samples)
HAP1_GFF     = "../../01data/hap1_genome.gff3"     # GFF3 for hap1
HAP2_GFF     = "../../01data/hap2_genome.gff3"     # GFF3 for hap2
OUTPUT_DIR   = "ASE_results_test"

TPM_THRESHOLD     = 0.5   # minimum TPM threshold for “expressed” genes (applied across samples)
LOG2FC_THRESHOLD  = 0.5   # |log2FC| threshold for ASE calling
PADJ_THRESHOLD    = 0.1   # adjusted p-value threshold for ASE calling

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ----------------------------- GFF3 parsing -----------------------------
def parse_gff_lengths(gff_file: str) -> dict:
    """
    Parse a GFF3 and compute gene lengths as the merged length of all exons
    mapped from transcripts to their parent genes.

    Returns
    -------
    dict: {gene_id: merged_exon_length}
    """
    trans_to_gene = {}
    gene_exons = defaultdict(list)

    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feat_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            attrs = fields[8]

            # parse attributes into a dict
            attr_dict = {}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k] = v

            if feat_type == "mRNA":
                # map transcript ID -> parent gene
                if "ID" in attr_dict and "Parent" in attr_dict:
                    trans_to_gene[attr_dict["ID"]] = attr_dict["Parent"]
            elif feat_type == "exon":
                # record exon intervals per gene (via transcript parent)
                if "Parent" in attr_dict:
                    parents = attr_dict["Parent"].split(",")
                    for parent in parents:
                        if parent in trans_to_gene:
                            gene_id = trans_to_gene[parent]
                            gene_exons[gene_id].append((start, end))

    gene_lengths = {}
    for gene_id, intervals in gene_exons.items():
        if not intervals:
            continue
        intervals.sort()
        merged = []
        current_start, current_end = intervals[0]

        for s, e in intervals[1:]:
            if s <= current_end + 1:
                current_end = max(current_end, e)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = s, e
        merged.append((current_start, current_end))

        total_length = sum(e - s + 1 for s, e in merged)
        gene_lengths[gene_id] = total_length

    return gene_lengths


# ----------------------------- I/O helpers -----------------------------
def read_allelic_pairs(pairs_file: str) -> pd.DataFrame:
    """
    Read allele gene pairs and create a unique pair_id.

    Accepts legacy header names 'Mene_hap2'/'Mene_hap1' and renames them.
    """
    pairs = pd.read_csv(pairs_file, sep="\t")
    if pairs.columns[0] == "Mene_hap2" and pairs.columns[1] == "Mene_hap1":
        pairs.columns = ["hap2_gene", "hap1_gene"]
    # Ensure expected column order exists
    if "hap1_gene" not in pairs.columns or "hap2_gene" not in pairs.columns:
        raise ValueError("Pairs file must contain columns: 'hap1_gene' and 'hap2_gene'")
    pairs["pair_id"] = pairs["hap1_gene"] + "|" + pairs["hap2_gene"]
    return pairs


def read_counts_data(hap1_file: str, hap2_file: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read hap1/hap2 count matrices (tab-separated; 'GeneID' as row index).
    """
    hap1_df = pd.read_csv(hap1_file, sep="\t", index_col="GeneID")
    hap2_df = pd.read_csv(hap2_file, sep="\t", index_col="GeneID")
    return hap1_df, hap2_df


# ----------------------------- TPM utilities -----------------------------
def calculate_tpm(counts_df: pd.DataFrame, gene_lengths: dict) -> pd.DataFrame:
    """
    Convert raw counts to TPM using gene lengths (in bp).
    """
    tpm_df = pd.DataFrame(index=counts_df.index, columns=counts_df.columns)
    for sample in counts_df.columns:
        counts = counts_df[sample].astype(float)
        lengths = counts_df.index.map(lambda x: float(gene_lengths.get(x, 1.0)))
        lengths = np.maximum(lengths, 1.0)  # guard against zeros
        rpk = counts / (lengths / 1000.0)
        denom = rpk.sum()
        tpm = rpk / (denom / 1e6) if denom > 0 else rpk * 0.0
        tpm_df[sample] = tpm
    return tpm_df


def filter_genes_by_tpm(hap1_df: pd.DataFrame,
                        hap2_df: pd.DataFrame,
                        gene_lengths: dict,
                        tpm_threshold: float = 1.0) -> set:
    """
    Keep genes whose TPM > threshold in **all** samples (using hap1+hap2 counts).
    """
    total_counts = hap1_df.add(hap2_df, fill_value=0)
    tpm = calculate_tpm(total_counts, gene_lengths)
    expressed = tpm.index[(tpm > tpm_threshold).all(axis=1)]
    return set(expressed)


# ----------------------------- DESeq2 prep/run -----------------------------
def prepare_deseq_input(pairs: pd.DataFrame,
                        hap1_df: pd.DataFrame,
                        hap2_df: pd.DataFrame,
                        treatment: str,
                        expressed_genes: set,
                        output_dir: str):
    """
    Build DESeq2 input (counts + metadata) for a given treatment.
    Columns used for this treatment are those starting with '<treatment>_'.
    """
    samples = [c for c in hap1_df.columns if c.startswith(treatment + "_")]

    # Map pair_id -> genes
    pair_id_to_genes = pairs.set_index("pair_id")[["hap1_gene", "hap2_gene"]].to_dict("index")

    dfs = []
    for sample in samples:
        hap1_sample = hap1_df[sample]
        hap2_sample = hap2_df[sample]

        sample_df = pairs.copy()
        sample_df["sample"] = sample
        sample_df["hap1_count"] = sample_df["hap1_gene"].map(hap1_sample.to_dict())
        sample_df["hap2_count"] = sample_df["hap2_gene"].map(hap2_sample.to_dict())
        # retain only expressed genes on both haplotypes
        sample_df = sample_df[
            sample_df["hap1_gene"].isin(expressed_genes)
            & sample_df["hap2_gene"].isin(expressed_genes)
        ]
        dfs.append(sample_df)

    if len(dfs) == 0:
        raise ValueError(f"No samples found for treatment prefix '{treatment}_'")

    combined_df = pd.concat(dfs, ignore_index=True)

    # Wide table: rows = pair_id, columns = (hap1_count/hap2_count) x sample
    deseq_counts = (
        combined_df.groupby(["pair_id", "sample"])[["hap1_count", "hap2_count"]]
        .sum()
        .unstack()
    ).fillna(0)

    # Build flat column names like: "<sample>_hap1" / "<sample>_hap2"
    new_columns = []
    for count_type, sample_name in deseq_counts.columns:
        allele_type = count_type.split("_")[0]  # 'hap1' or 'hap2'
        new_columns.append(f"{sample_name}_{allele_type}")
    deseq_counts.columns = new_columns

    # Save counts
    count_file = os.path.join(output_dir, f"{treatment}_deseq_counts.tsv")
    deseq_counts.to_csv(count_file, sep="\t", index=True)

    # Build metadata: one row per counts column
    meta_rows = []
    for col in deseq_counts.columns:
        parts = col.split("_")
        sample_name = "_".join(parts[:-1])   # everything except the last token
        allele_type = parts[-1]              # 'hap1' or 'hap2'
        replicate = sample_name.split("_")[-1]  # assume replicate tag is the last piece
        meta_rows.append({"sample": col, "allele": allele_type, "replicate": replicate})

    meta_df = pd.DataFrame(meta_rows)
    meta_file = os.path.join(output_dir, f"{treatment}_deseq_meta.tsv")
    meta_df.to_csv(meta_file, sep="\t", index=False)

    return count_file, meta_file, pair_id_to_genes


def run_deseq2(count_file: str, meta_file: str, treatment: str, output_dir: str) -> pd.DataFrame:
    """
    Write and run an R script that performs DESeq2 with design: ~ replicate + allele,
    and extracts the allele effect (hap2 vs hap1).
    """
    r_script = f"""
    suppressPackageStartupMessages({{library(DESeq2)}})
    counts <- read.table("{count_file}", header=TRUE, sep="\\t", row.names=1, check.names=FALSE)
    meta   <- read.table("{meta_file}",   header=TRUE, sep="\\t", check.names=FALSE)
    rownames(meta) <- meta$sample

    # Align columns/rows
    common <- intersect(colnames(counts), rownames(meta))
    if (length(common) < ncol(counts)) {{
        message("Warning: Some columns in counts are not present in metadata (they will be dropped).")
    }}
    counts <- counts[, common, drop=FALSE]
    meta   <- meta[common, , drop=FALSE]

    # Construct DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                                  colData   = meta,
                                  design    = ~ replicate + allele)

    dds <- DESeq(dds)
    res <- results(dds, contrast=c("allele","hap2","hap1"))
    res$pair_id  <- rownames(res)
    res$treatment <- "{treatment}"
    write.table(as.data.frame(res),
                file="{output_dir}/{treatment}_deseq_results.tsv",
                sep="\\t", quote=FALSE, row.names=FALSE)
    """

    script_path = os.path.join(output_dir, f"run_deseq2_{treatment}.R")
    with open(script_path, "w") as f:
        f.write(r_script)

    try:
        subprocess.run(["Rscript", script_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running DESeq2 for {treatment}: {e}")
        # Attempt to read partial results if present
        results_file = os.path.join(output_dir, f"{treatment}_deseq_results.tsv")
        if os.path.exists(results_file):
            return pd.read_csv(results_file, sep="\t")
        raise

    results_file = os.path.join(output_dir, f"{treatment}_deseq_results.tsv")
    return pd.read_csv(results_file, sep="\t")


# ----------------------------- ASE identification -----------------------------
def identify_ase_genes_by_treatment(all_results: pd.DataFrame,
                                    pair_id_to_genes: dict,
                                    log2fc_threshold: float = 1.0,
                                    padj_threshold: float = 0.05) -> dict[str, pd.DataFrame]:
    """
    For each treatment, select significant ASE pairs using thresholds and
    restore hap1/hap2 gene IDs from pair_id.
    """
    out = {}
    for treatment in all_results["treatment"].unique():
        tr = all_results[all_results["treatment"] == treatment]
        sig = tr[
            (np.abs(tr["log2FoldChange"]) >= log2fc_threshold)
            & (tr["padj"] < padj_threshold)
        ].copy()

        def map_pair(row):
            info = pair_id_to_genes.get(row["pair_id"], {})
            return pd.Series({"hap1_gene": info.get("hap1_gene", ""), "hap2_gene": info.get("hap2_gene", "")})

        if not sig.empty:
            sig[["hap1_gene", "hap2_gene"]] = sig.apply(map_pair, axis=1)

        out[treatment] = sig[["pair_id", "hap1_gene", "hap2_gene"]].drop_duplicates()
    return out


def identify_dynamic_ase_genes(ase_results_long: pd.DataFrame):
    """
    Identify dynamic ASE pairs: those significant in >1 treatments.
    Expects a long DataFrame with columns ['pair_id', 'treatment', ...].
    """
    counts = ase_results_long.groupby("pair_id")["treatment"].nunique().reset_index()
    counts.columns = ["pair_id", "n_treatments"]
    dynamic_ids = counts[counts["n_treatments"] > 1]
    dynamic_details = ase_results_long.merge(dynamic_ids, on="pair_id")
    return dynamic_ids, dynamic_details


# ----------------------------- Main -----------------------------
def main():
    # 0) Parse gene lengths
    print("Parsing gene length information from GFF3 ...")
    hap1_lengths = parse_gff_lengths(HAP1_GFF)
    hap2_lengths = parse_gff_lengths(HAP2_GFF)

    # 1) Read gene pairs and counts
    print("Step 1/8: Reading allele gene pairs ...")
    pairs = read_allelic_pairs(PAIRS_FILE)

    print("Step 2/8: Reading count matrices ...")
    hap1_df, hap2_df = read_counts_data(HAP1_COUNTS, HAP2_COUNTS)

    # Build a joint gene-length dictionary (average if length exists in both)
    all_genes = set(hap1_df.index) | set(hap2_df.index)
    gene_lengths = {}
    for g in all_genes:
        l1 = hap1_lengths.get(g, 0)
        l2 = hap2_lengths.get(g, 0)
        gene_lengths[g] = (l1 + l2) / 2 if (l1 and l2) else (l1 or l2 or 1)

    # 2) Filter expressed genes
    print(f"Step 3/8: Filtering genes with TPM > {TPM_THRESHOLD} in all samples ...")
    expressed_genes = filter_genes_by_tpm(hap1_df, hap2_df, gene_lengths, TPM_THRESHOLD)

    # 3) Determine treatments from column prefixes (e.g., 'AL_rep1' -> 'AL')
    treatments = sorted({col.split("_")[0] for col in hap1_df.columns})

    # 4) Prepare and run DESeq2 per treatment
    print("Step 4/8–6/8: Preparing inputs and running DESeq2 per treatment ...")
    all_results = []
    pair_id_to_genes = {}

    for treatment in tqdm(treatments, desc="Treatments"):
        count_file, meta_file, pair_map = prepare_deseq_input(
            pairs, hap1_df, hap2_df, treatment, expressed_genes, OUTPUT_DIR
        )
        if not pair_id_to_genes:
            pair_id_to_genes = pair_map

        res = run_deseq2(count_file, meta_file, treatment, OUTPUT_DIR)
        res["treatment"] = treatment
        all_results.append(res)

    all_results = pd.concat(all_results, ignore_index=True)

    # 5) Identify ASE pairs per treatment
    print("Step 7/8: Identifying ASE gene pairs per treatment ...")
    ase_by_tr = identify_ase_genes_by_treatment(
        all_results, pair_id_to_genes, LOG2FC_THRESHOLD, PADJ_THRESHOLD
    )

    # Save per-treatment ASE lists
    for tr, df in ase_by_tr.items():
        out_path = os.path.join(OUTPUT_DIR, f"{tr}_ase_gene_pairs.tsv")
        if not df.empty:
            df.to_csv(out_path, sep="\t", index=False)
        else:
            print(f"Warning: No significant ASE pairs detected for treatment '{tr}'.")
            pd.DataFrame(columns=["pair_id", "hap1_gene", "hap2_gene"]).to_csv(out_path, sep="\t", index=False)

    print(f"✅ Analysis complete. Results saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
