#!/usr/bin/env bash
#########################
# Config / input paths
#########################

# Reference and query genome assemblies (e.g. two haplotypes)
REF_GENOME="Hap1.fa"
QRY_GENOME="Hap2.fa"

# Prefix for alignment outputs
PREFIX="hap1_vs_hap2"

# Threads
THREADS=20

# Gene BED files (4 columns: chr, start, end, gene_id)
HAP1_GENE_BED="hap1_genes_4col.bed"
HAP2_GENE_BED="hap2_genes_4col.bed"

# Expression matrix (rows: gene_id, columns: samples, values: TPM)
TPM_MATRIX="tpm_matrix.tsv"

# TE annotation BEDs for hap1 & hap2 (example paths)
TE_DIR="./TE"
HAP1_500BP_BED="hap1_500bp_fixed.bed"
HAP2_500BP_BED="hap2_500bp_fixed.bed"

#########################
# 0. Whole-genome alignment: nucmer + SyRI
#########################

# Align two genomes with nucmer
nucmer --maxmatch -c 100 -b 500 -l 50 -t "${THREADS}" \
  "${REF_GENOME}" "${QRY_GENOME}" -p "${PREFIX}"

# Filter alignments (identity >=90%, length >=100bp etc.)
delta-filter -m -i 90 -l 100 "${PREFIX}.delta" > "${PREFIX}.filtered.delta"

# Convert delta to tab-delimited coords
show-coords -THrd "${PREFIX}.filtered.delta" > "${PREFIX}.filtered.coords"

# Run SyRI to detect SVs and SNPs
# -c: coords file; -d: delta file; -r: reference; -q: query
SyRI \
  -c "${PREFIX}.filtered.coords" \
  -d "${PREFIX}.filtered.delta" \
  -r "${REF_GENOME}" \
  -q "${QRY_GENOME}" \
  -F B \
  -k > syri.log

# SyRI produces:
# - syri.out (SV annotations)
# - syri.vcf  or syri_output.vcf depending on version
# Adjust the VCF name below if needed.
SYRI_VCF="syri_output.vcf"

#########################
# 1. Extract SNP / InDel / SV from SyRI VCF
#########################

# SNPs (single-nucleotide polymorphisms)
bcftools view -v snps "${SYRI_VCF}" -o snp.vcf

# Indels (< 50 bp by allele-length difference)
bcftools view -v indels "${SYRI_VCF}" \
  | awk 'BEGIN{OFS="\t"}
         /^#/ {print; next}
         {
           d = length($4) - length($5);
           if (d < 0) d = -d;
           if (d < 50) print;
         }' > indel.vcf

# Structural variants (SVs: >= 50 bp by allele-length difference)
bcftools view -v indels,other "${SYRI_VCF}" \
  | awk 'BEGIN{OFS="\t"}
         /^#/ {print; next}
         {
           d = length($4) - length($5);
           if (d < 0) d = -d;
           if (d >= 50) print;
         }' > sv.vcf

#########################
# 2. Parse SyRI classes into BED for hap1 and hap2
#########################
# syri.out: typical columns:
# 1: ref_chr 2: ref_start 3: ref_end ...
# 6: qry_chr 7: qry_start 8: qry_end ...
# 10: type (SYN, INV, DUP, TRANS, NOTAL, etc.)
#########################

# DUP
grep -w "DUP" syri.out | awk '$10 == "-" {print}' | cut -f1-3 > hap1_DUP.bed
grep -w "DUP" syri.out | awk '$10 == "-" {print}' | cut -f6-8 > hap2_DUP.bed

# INV
grep -w "INV" syri.out | awk '$10 == "-" {print}' | cut -f1-3 > hap1_INV.bed
grep -w "INV" syri.out | awk '$10 == "-" {print}' | cut -f6-8 > hap2_INV.bed

# SYN (syntenic blocks)
grep -w "SYN" syri.out | awk '$10 == "-" {print}' | cut -f1-3 > hap1_SYN.bed
grep -w "SYN" syri.out | awk '$10 == "-" {print}' | cut -f6-8 > hap2_SYN.bed

# NOTAL (non-aligned regions)
grep -w "NOTAL" syri.out | awk '$6 == "-" {print}' | cut -f1-3 > hap1_NOTAL.bed
grep -w "NOTAL" syri.out | awk '$6 == "-" {print}' | cut -f6-8 > hap2_NOTAL.bed

# TRANS (translocations)
grep -w "TRANS" syri.out | awk '$10 == "-" {print}' | cut -f1-3 > hap1_TRANS.bed
grep -w "TRANS" syri.out | awk '$10 == "-" {print}' | cut -f6-8 > hap2_TRANS.bed

#########################
# 3. Generic BED sort helper (example)
#########################
# sort -k1,1 -k2,2n input.bed > sorted_input.bed

#########################
# 4. Effect of SVs on gene expression (TPM mean + two-sided Wilcoxon)
#########################
# Idea:
# (1) Assign genes to SV regions (INV + TRANS + DUP) or syntenic regions.
#     If a gene overlaps SV and SYN, classify as SV.
# (2) Compute mean TPM across all samples per gene.
# (3) Compare SV vs syntenic gene expression using two-sided Mann–Whitney U.
#########################

#########################
# 4.1 Group genes by SV vs syntenic (generic example)
#########################

# 1) Merge all SV classes (INV + TRANS + DUP) into a single BED
cat hap1_INV.bed hap1_TRANS.bed hap1_DUP.bed \
  | bedtools sort -i - \
  | bedtools merge -i - > hap1_sv_merged.bed

# For syntenic regions, you may want to use hap1_SYN.bed as syntenic blocks
bedtools sort -i hap1_SYN.bed | bedtools merge -i - > hap1_syntenic_merged.bed

# 2) Genes overlapping any SV (at least 1 bp)
bedtools intersect -u -a "${HAP1_GENE_BED}" -b hap1_sv_merged.bed > hap1_genes_in_SV.bed

# 3) Genes overlapping syntenic regions
bedtools intersect -u -a "${HAP1_GENE_BED}" -b hap1_syntenic_merged.bed > hap1_genes_in_syn.bed

# 4) Remove genes already assigned to SV from syntenic gene list
cut -f4 hap1_genes_in_SV.bed  > hap1_SV.ids
cut -f4 hap1_genes_in_syn.bed > hap1_SYN.ids

grep -v -F -f hap1_SV.ids hap1_SYN.ids > hap1_SYN_only.ids

# For hap2, same logic can be applied (hap2_sv_merged.bed, hap2_syntenic_merged.bed, etc.)

#########################
# 4.2 Mann–Whitney U test in Python
#########################

cat > stats_sv_vs_syn.py << 'PYCODE'
import pandas as pd
from scipy.stats import mannwhitneyu

# TPM matrix: rows = gene_id, columns = samples
tpm = pd.read_csv("tpm_matrix.tsv", sep="\t", index_col=0)

# Group gene lists
sv_ids  = pd.read_csv("hap1_SV.ids",  header=None)[0].tolist()
syn_ids = pd.read_csv("hap1_SYN_only.ids", header=None)[0].tolist()

# Keep only genes that exist in expression matrix
sv_ids  = [g for g in sv_ids  if g in tpm.index]
syn_ids = [g for g in syn_ids if g in tpm.index]

# Mean TPM per gene across samples
gene_mean = tpm.mean(axis=1)

sv_expr  = gene_mean.loc[sv_ids]
syn_expr = gene_mean.loc[syn_ids]

# Two-sided Mann–Whitney U test
u_stat, p_val = mannwhitneyu(sv_expr.values, syn_expr.values, alternative="two-sided")

print(f"SV genes (n={len(sv_expr)}), Syntenic genes (n={len(syn_expr)})")
print(f"Mann–Whitney U two-sided p-value: {p_val:.3e}")

# Optional: export data for plotting
out = pd.DataFrame({
    "gene_id": sv_expr.index.tolist() + syn_expr.index.tolist(),
    "group":   ["SV"]*len(sv_expr) + ["Syntenic"]*len(syn_expr),
    "mean_TPM": pd.concat([sv_expr, syn_expr]).values
})
out.to_csv("sv_vs_syn_meanTPM.tsv", sep="\t", index=False)
PYCODE

# Run the statistical test
python3 stats_sv_vs_syn.py

#########################
# 5. Classify genes by SV type per haplotype
#########################
# hap1_genes_4col.bed is created from a richer BED (e.g. with attributes in col10)
#########################

# Example to create 4-column gene BED from a GFF-like BED
# (Assuming gene ID is in the 10th field as "ID=xxx;...")
# awk 'BEGIN {OFS="\t"} NF>=10 {
#     match($10, /ID=([^;]+)/, id);
#     gene_id = (id[1] != "") ? id[1] : $10;
#     print $1, $2, $3, gene_id
# }' hap1_genes.bed > hap1_genes_4col.bed

# Intersect hap1 genes with each SV class
bedtools intersect -u -a "${HAP1_GENE_BED}" -b hap1_INV.bed   > hap1_INV_genes.bed
bedtools intersect -u -a "${HAP1_GENE_BED}" -b hap1_DUP.bed   > hap1_DUP_genes.bed
bedtools intersect -u -a "${HAP1_GENE_BED}" -b hap1_TRANS.bed > hap1_TRANS_genes.bed
bedtools intersect -u -a "${HAP1_GENE_BED}" -b hap1_SYN.bed   > hap1_SYN_genes.bed
bedtools intersect -u -a "${HAP1_GENE_BED}" -b hap1_NOTAL.bed > hap1_not_aligned_genes.bed

bedtools intersect -u -a "${HAP2_GENE_BED}" -b hap2_INV.bed   > hap2_INV_genes.bed
bedtools intersect -u -a "${HAP2_GENE_BED}" -b hap2_DUP.bed   > hap2_DUP_genes.bed
bedtools intersect -u -a "${HAP2_GENE_BED}" -b hap2_TRANS.bed > hap2_TRANS_genes.bed
bedtools intersect -u -a "${HAP2_GENE_BED}" -b hap2_SYN.bed   > hap2_SYN_genes.bed
bedtools intersect -u -a "${HAP2_GENE_BED}" -b hap2_NOTAL.bed > hap2_not_aligned_genes.bed

# Extract gene IDs for each category
cut -f4 hap1_INV_genes.bed          > hap1_INV.ids
cut -f4 hap1_DUP_genes.bed          > hap1_DUP.ids
cut -f4 hap1_TRANS_genes.bed        > hap1_TRANS.ids
cut -f4 hap1_SYN_genes.bed          > hap1_SYN.ids
cut -f4 hap1_not_aligned_genes.bed  > hap1_not_aligned.ids

cut -f4 hap2_INV_genes.bed          > hap2_INV.ids
cut -f4 hap2_DUP_genes.bed          > hap2_DUP.ids
cut -f4 hap2_TRANS_genes.bed        > hap2_TRANS.ids
cut -f4 hap2_SYN_genes.bed          > hap2_SYN.ids
cut -f4 hap2_not_aligned_genes.bed  > hap2_not_aligned.ids

#########################
# 6. TE overlap with fixed windows around breakpoints
#########################

# TE overlap for multiple TE classes (example)
bedtools intersect -a "${HAP1_500BP_BED}" -b "${TE_DIR}/hap1.Gypsy.bed"       -wo > hap1_Gypsy_overlap.txt
bedtools intersect -a "${HAP2_500BP_BED}" -b "${TE_DIR}/hap2.Gypsy.bed"       -wo > hap2_Gypsy_overlap.txt
bedtools intersect -a "${HAP1_500BP_BED}" -b "${TE_DIR}/hap1.DNA.bed"         -wo > hap1_DNA_overlap.txt
bedtools intersect -a "${HAP2_500BP_BED}" -b "${TE_DIR}/hap2.DNA.bed"         -wo > hap2_DNA_overlap.txt
bedtools intersect -a "${HAP1_500BP_BED}" -b "${TE_DIR}/hap1.Helitron.bed"    -wo > hap1_Helitron_overlap.txt
bedtools intersect -a "${HAP2_500BP_BED}" -b "${TE_DIR}/hap2.Helitron.bed"    -wo > hap2_Helitron_overlap.txt
bedtools intersect -a "${HAP1_500BP_BED}" -b "${TE_DIR}/hap1.LINE.bed"        -wo > hap1_LINE_overlap.txt
bedtools intersect -a "${HAP2_500BP_BED}" -b "${TE_DIR}/hap2.LINE.bed"        -wo > hap2_LINE_overlap.txt
bedtools intersect -a "${HAP1_500BP_BED}" -b "${TE_DIR}/hap1.LTR_unknown.bed" -wo > hap1_LTR_unknown_overlap.txt
bedtools intersect -a "${HAP2_500BP_BED}" -b "${TE_DIR}/hap2.LTR_unknown.bed" -wo > hap2_LTR_unknown_overlap.txt

#########################
# 7. Try different window sizes around breakpoints (e.g. for INV)
#########################

# hap1_INV.txt: assume columns (chr, break_start, break_end, ...)
WINDOW_SIZES=(150 500 1000 2000 5000)

for WINDOW in "${WINDOW_SIZES[@]}"; do
  awk -v window="${WINDOW}" 'BEGIN{OFS="\t"} {
      left_start  = ($2 - window > 0) ? $2 - window : 0;
      left_end    = $2 + window;
      right_start = ($3 - window > 0) ? $3 - window : 0;
      right_end   = $3 + window;
      print $1, left_start,  left_end;
      print $1, right_start, right_end;
  }' hap1_INV.txt > "hap1_${WINDOW}.bed"
done

echo "Pipeline finished (check intermediate outputs for details)."
