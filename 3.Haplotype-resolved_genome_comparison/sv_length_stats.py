#!/usr/bin/env python3
"""
sv_length_stats.py

Summarize length, count, and length distribution of structural variants (SVs)
from:
  1) BED files (per-class SV intervals)
  2) SyRI syri.out file

Outputs:
  - SV_statistics_report.txt   : detailed text report
  - all_sv_data.csv            : all parsed SV records
  - all_SV_SV_length_distribution.png : length distribution plots
"""

import os
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ----------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------

def read_bed_file(bed_file: str):
    """
    Read a BED file and return a list of SV records.

    Each record is a dict with keys:
      chrom, start, end, length, type

    SV type is inferred from the BED filename (basename without ".bed").
    """
    sv_data = []
    if not os.path.exists(bed_file):
        print(f"[WARN] BED file not found: {bed_file}")
        return sv_data

    sv_type = os.path.basename(bed_file).replace(".bed", "")

    with open(bed_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            length = max(end - start, 0)
            sv_data.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "length": length,
                    "type": sv_type,
                }
            )
    return sv_data


def read_syri_out(syri_file: str):
    """
    Read syri.out file and extract SV information.

    Expected SyRI columns (simplified):
      0: ref_chr
      1: ref_start
      2: ref_end
      ...
      5: qry_chr
      6: qry_start
      7: qry_end
      8: ID
      ...
      10: type (INV, DUP, INVDP, SYN, NOTAL, TRANS, etc.)

    Returns a list of dicts with keys:
      ref_chrom, ref_start, ref_end,
      qry_chrom, qry_start, qry_end,
      length, type, id
    """
    sv_data = []
    if not os.path.exists(syri_file):
        print(f"[WARN] syri.out file not found: {syri_file}")
        return sv_data

    keep_types = {"INV", "DUP", "INVDP", "SYN", "NOTAL", "TRANS"}

    with open(syri_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 11:
                continue

            sv_type = parts[10]
            # Only keep major SV classes of interest
            if sv_type not in keep_types:
                continue

            try:
                ref_start = int(parts[1])
                ref_end = int(parts[2])
                qry_start = int(parts[6])
                qry_end = int(parts[7])
            except ValueError:
                continue

            ref_chrom = parts[0]
            qry_chrom = parts[5]
            length = max(ref_end - ref_start, 0)
            sv_id = parts[8] if len(parts) > 8 else "N/A"

            sv_data.append(
                {
                    "ref_chrom": ref_chrom,
                    "ref_start": ref_start,
                    "ref_end": ref_end,
                    "qry_chrom": qry_chrom,
                    "qry_start": qry_start,
                    "qry_end": qry_end,
                    "length": length,
                    "type": sv_type,
                    "id": sv_id,
                }
            )
    return sv_data


# ----------------------------------------------------------------------
# Statistics & plotting
# ----------------------------------------------------------------------

def analyze_sv_statistics(sv_data, source_type: str):
    """
    Print basic statistics for a list of SV records and
    return a pandas DataFrame.
    """
    if not sv_data:
        print(f"[WARN] No SV records found for: {source_type}")
        return None

    df = pd.DataFrame(sv_data)

    print(f"\n=== SV statistics: {source_type} ===")
    print(f"Total SV count: {len(df)}")

    if "type" in df.columns:
        print("\nSV type distribution:")
        type_counts = df["type"].value_counts()
        print(type_counts)

        print("\nSV length per type:")
        length_stats = (
            df.groupby("type")["length"]
            .agg(["count", "sum", "mean", "median", "min", "max"])
            .round(2)
        )
        print(length_stats)
    else:
        print("\nNo 'type' column found. Only length summary:")
        print(df["length"].describe())

    return df


def plot_length_distribution(df: pd.DataFrame, title: str):
    """
    Plot multiple views of SV length distribution:
      - histogram of log10(length+1)
      - boxplot by type
      - cumulative distribution
      - barplot of length bins per type

    Saves a PNG: f"{title}_SV_length_distribution.png"
    """
    if df is None or df.empty:
        print("[WARN] No data for plotting.")
        return

    # Work on a copy to avoid modifying original DataFrame in-place.
    df_plot = df.copy()

    if "length" not in df_plot.columns:
        print("[WARN] No 'length' column in dataframe; skip plotting.")
        return

    plt.figure(figsize=(12, 8))
    sv_types = df_plot["type"].unique() if "type" in df_plot.columns else ["ALL"]

    # ------------------------------------------------------------------ #
    # 1) Histogram of log10(length + 1) per type
    # ------------------------------------------------------------------ #
    plt.subplot(2, 2, 1)
    for sv_type in sv_types:
        sub = df_plot if sv_type == "ALL" else df_plot[df_plot["type"] == sv_type]
        if sub.empty:
            continue
        vals = np.log10(sub["length"].values + 1)
        plt.hist(vals, bins=30, alpha=0.7, label=sv_type)
    plt.xlabel("SV length (log10(bp + 1))")
    plt.ylabel("Count")
    plt.title(f"{title} - SV length histogram")
    if len(sv_types) > 1:
        plt.legend()

    # ------------------------------------------------------------------ #
    # 2) Boxplot by type
    # ------------------------------------------------------------------ #
    plt.subplot(2, 2, 2)
    df_plot["log_length"] = np.log10(df_plot["length"] + 1)
    if "type" in df_plot.columns:
        sns.boxplot(data=df_plot, x="type", y="log_length")
        plt.xlabel("SV type")
    else:
        sns.boxplot(data=df_plot, y="log_length")
        plt.xlabel("")
    plt.ylabel("SV length (log10(bp + 1))")
    plt.title(f"{title} - SV length (boxplot)")
    plt.xticks(rotation=45, ha="right")

    # ------------------------------------------------------------------ #
    # 3) Cumulative distribution per type
    # ------------------------------------------------------------------ #
    plt.subplot(2, 2, 3)
    for sv_type in sv_types:
        sub = df_plot if sv_type == "ALL" else df_plot[df_plot["type"] == sv_type]
        if sub.empty:
            continue
        vals = np.sort(sub["length"].values)
        y = np.arange(1, len(vals) + 1) / len(vals)
        plt.plot(vals, y, label=sv_type, linewidth=2)
    plt.xscale("log")
    plt.xlabel("SV length (bp)")
    plt.ylabel("Cumulative fraction")
    plt.title(f"{title} - cumulative length distribution")
    if len(sv_types) > 1:
        plt.legend()

    # ------------------------------------------------------------------ #
    # 4) Binned length distribution by type
    # ------------------------------------------------------------------ #
    plt.subplot(2, 2, 4)
    length_bins = [0, 100, 500, 1000, 5000, 10000, 50000, float("inf")]
    bin_labels = ["<100", "100–500", "500–1k", "1k–5k", "5k–10k", "10k–50k", ">50k"]

    df_plot["length_bin"] = pd.cut(df_plot["length"], bins=length_bins, labels=bin_labels)
    if "type" in df_plot.columns:
        bin_counts = df_plot.groupby(["type", "length_bin"]).size().unstack(fill_value=0)
        bin_counts.plot(kind="bar", stacked=True, ax=plt.gca())
        plt.xlabel("SV type")
    else:
        bin_counts = df_plot.groupby(["length_bin"]).size()
        bin_counts.plot(kind="bar", ax=plt.gca())
        plt.xlabel("Length bin")
    plt.ylabel("Count")
    plt.title(f"{title} - counts by length bin")
    plt.xticks(rotation=45, ha="right")
    if "type" in df_plot.columns:
        plt.legend(title="Length bin (bp)")
    else:
        plt.legend().set_visible(False)

    plt.tight_layout()
    out_png = f"{title}_SV_length_distribution.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[INFO] Saved plot: {out_png}")


def generate_detailed_report(df: pd.DataFrame, output_file: str):
    """
    Generate a plain-text detailed summary report of SV statistics.
    """
    if df is None or df.empty:
        print("[WARN] No data for detailed report.")
        return

    with open(output_file, "w") as f:
        f.write("SV detailed statistics report\n")
        f.write("=" * 60 + "\n\n")

        f.write("Global statistics:\n")
        f.write(f"  Total SV count : {len(df)}\n")
        f.write(f"  Total SV length: {df['length'].sum():,} bp\n")
        f.write(f"  Mean SV length : {df['length'].mean():.2f} bp\n")
        f.write(f"  Median length  : {df['length'].median():.2f} bp\n\n")

        if "type" in df.columns:
            f.write("Per-type statistics:\n")
            type_stats = (
                df.groupby("type")["length"]
                .agg(["count", "sum", "mean", "median", "min", "max"])
                .round(2)
            )

            for sv_type in type_stats.index:
                stats = type_stats.loc[sv_type]
                f.write(f"\n{sv_type}:\n")
                f.write(f"  Count    : {stats['count']}\n")
                f.write(f"  Total bp : {stats['sum']:,} bp\n")
                f.write(f"  Mean bp  : {stats['mean']:.2f}\n")
                f.write(f"  Median bp: {stats['median']:.2f}\n")
                f.write(f"  Min bp   : {stats['min']:.2f}\n")
                f.write(f"  Max bp   : {stats['max']:.2f}\n")

        # Length-bin distribution
        f.write("\nLength-bin distribution (all SVs):\n")
        length_bins = [0, 100, 500, 1000, 5000, 10000, 50000, 100000, float("inf")]
        bin_labels = [
            "<100",
            "100–500",
            "500–1k",
            "1k–5k",
            "5k–10k",
            "10k–50k",
            "50k–100k",
            ">100k",
        ]

        df_bins = df.copy()
        df_bins["length_bin"] = pd.cut(df_bins["length"], bins=length_bins, labels=bin_labels)
        bin_distribution = df_bins["length_bin"].value_counts().reindex(bin_labels, fill_value=0)

        total = len(df_bins)
        for bin_label in bin_labels:
            count = bin_distribution[bin_label]
            percentage = (count / total * 100.0) if total > 0 else 0.0
            f.write(f"  {bin_label:>9} bp : {count:6d} ({percentage:5.2f}%)\n")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    """
    Main entry point.
    Currently uses hard-coded filenames; adjust as needed or replace with argparse.
    """

    # ------------------------------------------------------------------
    # File configuration
    # ------------------------------------------------------------------
    bed_files = {
        "INV": "INV.bed",
        "SYN": "SYN.bed",
        "TRANS": "TRANS.bed",
        "hap1_DUP": "hap1_DUP.bed",
        "hap2_DUP": "hap2_DUP.bed",
        "hap1_not_aligned": "hap1_not_aligned.bed",
        "hap2_not_aligned": "hap2_not_aligned.bed",
    }

    syri_file = "syri.out"

    all_sv_data = []

    # ------------------------------------------------------------------
    # Read BED files
    # ------------------------------------------------------------------
    print("[INFO] Reading BED files...")
    for sv_label, bed_file in bed_files.items():
        sv_data = read_bed_file(bed_file)
        all_sv_data.extend(sv_data)
        print(f"  {bed_file}: {len(sv_data)} SVs ({sv_label})")

    # ------------------------------------------------------------------
    # Read syri.out
    # ------------------------------------------------------------------
    print("\n[INFO] Reading syri.out...")
    syri_sv_data = read_syri_out(syri_file)
    all_sv_data.extend(syri_sv_data)
    print(f"  {syri_file}: {len(syri_sv_data)} SVs")

    if not all_sv_data:
        print("[ERROR] No SV data found from any source. Exiting.")
        return

    # ------------------------------------------------------------------
    # Global statistics
    # ------------------------------------------------------------------
    df_all = analyze_sv_statistics(all_sv_data, "All sources combined")
    if df_all is None:
        print("[ERROR] Failed to build DataFrame from SV data.")
        return

    # BED-only statistics
    bed_types = set(bed_files.keys())
    bed_sv_data = [sv for sv in all_sv_data if sv.get("type") in bed_types]
    if bed_sv_data:
        analyze_sv_statistics(bed_sv_data, "BED files")

    # SyRI-only statistics (by SyRI major SV classes)
    syri_types = {"INV", "DUP", "INVDP", "SYN", "NOTAL", "TRANS"}
    syri_only_data = [sv for sv in all_sv_data if sv.get("type") in syri_types]
    if syri_only_data:
        analyze_sv_statistics(syri_only_data, "SyRI (syri.out)")

    # ------------------------------------------------------------------
    # Plot distributions
    # ------------------------------------------------------------------
    print("\n[INFO] Plotting SV length distributions...")
    # Use an ASCII-safe title to avoid weird filenames
    plot_length_distribution(df_all, "all_SV")

    # ------------------------------------------------------------------
    # Detailed report
    # ------------------------------------------------------------------
    print("\n[INFO] Generating detailed report...")
    generate_detailed_report(df_all, "SV_statistics_report.txt")

    # ------------------------------------------------------------------
    # Save full table
    # ------------------------------------------------------------------
    df_all.to_csv("all_sv_data.csv", index=False)

    print(f"\n[INFO] Done. Total SVs analyzed: {len(df_all)}")
    print("Output files:")
    print("  - SV_statistics_report.txt")
    print("  - all_sv_data.csv")
    print("  - all_SV_SV_length_distribution.png")


if __name__ == "__main__":
    main()
