#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summarize associations between structural variants (SVs) and ASE gene-pair classes.

Inputs
------
- class.txt             : TSV with columns: hap1_gene, hap2_gene, classification
- overlaps_hap1.tsv     : bedtools-style overlaps for hap1 genes vs SVs
- overlaps_hap2.tsv     : (optional) same as above for hap2; skipped if absent

Outputs (in sv_assoc_out/)
--------------------------
- pair_sv_merged.tsv        : per-pair merged SV indicators and totals
- summary_overall.tsv       : descriptive stats of total SVs per class
- summary_by_svtype.tsv     : proportion of pairs with ≥1 SV by type per class
- stats.txt                 : Kruskal–Wallis test on total SVs across classes
- pairwise_MWU.tsv          : pairwise Mann–Whitney U tests (BH-adjusted)
- fisher_by_svtype.tsv      : Fisher’s exact tests per SV type (BH-adjusted)
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests

# ------------------ Input files ------------------
CLASS_FILE = "class.txt"            # hap1_gene  hap2_gene  classification
OVL_H1     = "overlaps_hap1.tsv"    # bedtools output (hap1)
OVL_H2     = "overlaps_hap2.tsv"    # optional (hap2); skipped if absent
OUT_DIR    = "sv_assoc_out"
os.makedirs(OUT_DIR, exist_ok=True)


# ------------------ Load classification ------------------
cls = pd.read_csv(CLASS_FILE, sep=r"\s+|\t", engine="python")
# Build a pair_id for easier tracking
cls["pair_id"] = cls["hap1_gene"] + "__" + cls["hap2_gene"]


# ------------------ Load overlaps -> per-gene SV summary ------------------
def load_overlap(path: str) -> pd.DataFrame:
    """Load an overlap table (bedtools-like) into a DataFrame."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "gene_chr", "gene_start", "gene_end", "GeneID",
            "sv_chr", "sv_start", "sv_end", "svtype",
        ],
    )
    return df


def per_gene_sv(df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize SVs per gene:
      - SV_total : total number of overlaps
      - One-hot (0/1) flags for having at least one of each SV type
        among: INS, DEL, DUP, INV, CNV
    """
    if df.empty:
        return pd.DataFrame(
            columns=["GeneID", "SV_total", "INS", "DEL", "DUP", "INV", "CNV"]
        )
    grp = df.groupby("GeneID", sort=False)
    out = grp.size().rename("SV_total").to_frame().reset_index()
    for t in ["INS", "DEL", "DUP", "INV", "CNV"]:
        out[t] = grp.apply(lambda g: int(any(g["svtype"] == t))).values
    return out


# hap1
h1 = per_gene_sv(load_overlap(OVL_H1)) if os.path.exists(OVL_H1) else pd.DataFrame()
h1 = h1.rename(
    columns={
        "SV_total": "hap1_SV_total",
        "INS": "hap1_INS",
        "DEL": "hap1_DEL",
        "DUP": "hap1_DUP",
        "INV": "hap1_INV",
        "CNV": "hap1_CNV",
    }
)

# hap2 (optional)
if os.path.exists(OVL_H2):
    h2 = per_gene_sv(load_overlap(OVL_H2)).rename(
        columns={
            "SV_total": "hap2_SV_total",
            "INS": "hap2_INS",
            "DEL": "hap2_DEL",
            "DUP": "hap2_DUP",
            "INV": "hap2_INV",
            "CNV": "hap2_CNV",
        }
    )
else:
    h2 = pd.DataFrame(
        columns=["GeneID", "hap2_SV_total", "hap2_INS", "hap2_DEL", "hap2_DUP", "hap2_INV", "hap2_CNV"]
    )

# ------------------ Merge to gene pairs ------------------
merged = cls.copy()

# Merge hap1 summaries
merged = merged.merge(h1, how="left", left_on="hap1_gene", right_on="GeneID")
merged = merged.drop(columns=["GeneID"])

# Merge hap2 summaries (fill zeros if not provided)
if not h2.empty:
    merged = merged.merge(h2, how="left", left_on="hap2_gene", right_on="GeneID")
    merged = merged.drop(columns=["GeneID"])
else:
    for col in ["hap2_SV_total", "hap2_INS", "hap2_DEL", "hap2_DUP", "hap2_INV", "hap2_CNV"]:
        merged[col] = 0

# Fill remaining NA with 0
merged = merged.fillna(0)

# Pair-level metrics
merged["pair_SV_total"] = merged["hap1_SV_total"].astype(int) + merged["hap2_SV_total"].astype(int)
for t in ["INS", "DEL", "DUP", "INV", "CNV"]:
    merged[f"pair_has_{t}"] = (
        (merged.get(f"hap1_{t}", 0) > 0) | (merged.get(f"hap2_{t}", 0) > 0)
    ).astype(int)

# Export merged details
merged.to_csv(os.path.join(OUT_DIR, "pair_sv_merged.tsv"), sep="\t", index=False)


# ------------------ Descriptive summaries ------------------
desc = (
    merged.groupby("classification", sort=False)
    .agg(
        mean_SV=("pair_SV_total", "mean"),
        median_SV=("pair_SV_total", "median"),
        prop_any_SV=("pair_SV_total", lambda x: np.mean(x > 0)),
    )
    .reset_index()
)

by_type = []
for t in ["INS", "DEL", "DUP", "INV", "CNV"]:
    tmp = (
        merged.groupby("classification", sort=False)[f"pair_has_{t}"]
        .agg(["mean", "sum", "count"])
        .reset_index()
    )
    tmp["svtype"] = t
    by_type.append(tmp)
by_type = pd.concat(by_type, ignore_index=True)

desc.to_csv(os.path.join(OUT_DIR, "summary_overall.tsv"), sep="\t", index=False)
by_type.to_csv(os.path.join(OUT_DIR, "summary_by_svtype.tsv"), sep="\t", index=False)


# ------------------ Between-class statistical tests ------------------
# 1) Kruskal–Wallis: compare total SVs across classes
groups = [g["pair_SV_total"].values for _, g in merged.groupby("classification", sort=False)]
kw_stat, kw_p = kruskal(*[x for x in groups if len(x) > 0])  # guard against empty groups
with open(os.path.join(OUT_DIR, "stats.txt"), "w") as fw:
    fw.write(
        f"Kruskal-Wallis (pair_SV_total ~ classification): stat={kw_stat:.3f}, p={kw_p:.3e}\n"
    )

# 2) Pairwise Mann–Whitney U with BH correction
pairs = [("Hap", "NoDiff"), ("Hap", "SubNeo"), ("SubNeo", "NoDiff")]
p_raw, labels = [], []
for a, b in pairs:
    xa = merged.loc[merged["classification"] == a, "pair_SV_total"].values
    xb = merged.loc[merged["classification"] == b, "pair_SV_total"].values
    if len(xa) > 0 and len(xb) > 0:
        stat, p = mannwhitneyu(xa, xb, alternative="two-sided")
    else:
        p = np.nan
    p_raw.append(p)
    labels.append(f"{a} vs {b}")

# Replace NaNs with 1.0 before multiple testing
p_for_adj = [1.0 if (pd.isna(p)) else p for p in p_raw]
rej, p_bh, _, _ = multipletests(p_for_adj, method="fdr_bh")

pair_df = pd.DataFrame(
    {"contrast": labels, "p_raw": p_raw, "p_bh": p_bh, "reject_at_0.05": rej}
)
pair_df.to_csv(os.path.join(OUT_DIR, "pairwise_MWU.tsv"), sep="\t", index=False)

# 3) Enrichment by SV type using Fisher’s exact test
#    Test whether the proportion of pairs with ≥1 SV of a given type differs between classes
def fisher_table(classA: str, classB: str, col: str) -> np.ndarray:
    a1 = merged.loc[merged["classification"] == classA, col].sum()
    a0 = (merged["classification"] == classA).sum() - a1
    b1 = merged.loc[merged["classification"] == classB, col].sum()
    b0 = (merged["classification"] == classB).sum() - b1
    return np.array([[a1, a0], [b1, b0]])

fis_rows = []
for t in ["INS", "DEL", "DUP", "INV", "CNV"]:
    for (a, b) in [("Hap", "NoDiff"), ("Hap", "SubNeo")]:
        table = fisher_table(a, b, f"pair_has_{t}")
        try:
            OR, p = fisher_exact(table)
        except Exception:
            OR, p = np.nan, np.nan
        fis_rows.append(
            {"svtype": t, "contrast": f"{a} vs {b}", "odds_ratio": OR, "p_raw": p}
        )

fis = pd.DataFrame(fis_rows)
# Replace NaN p-values with 1.0 before BH correction
p_for_adj = fis["p_raw"].fillna(1.0).to_numpy()
rej, p_bh, _, _ = multipletests(p_for_adj, method="fdr_bh")
fis["p_bh"] = p_bh
fis["reject_at_0.05"] = rej
fis.to_csv(os.path.join(OUT_DIR, "fisher_by_svtype.tsv"), sep="\t", index=False)

print("Done. Outputs saved in:", OUT_DIR)
