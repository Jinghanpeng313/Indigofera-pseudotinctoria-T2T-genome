import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path

def read_syri(syri_file):
    """Read the syri.out file and extract relevant columns."""
    cols = ["ref_chr", "ref_start", "ref_end", "qry_chr", "qry_start", "qry_end", "type"]
    df = pd.read_csv(syri_file, sep=r"\s+", header=None, names=cols)
    df["ref_start"] = pd.to_numeric(df["ref_start"], errors='coerce')
    df["ref_end"] = pd.to_numeric(df["ref_end"], errors='coerce')
    df["qry_start"] = pd.to_numeric(df["qry_start"], errors='coerce')
    df["qry_end"] = pd.to_numeric(df["qry_end"], errors='coerce')
    df["ref_len"] = df["ref_end"] - df["ref_start"] + 1
    df["qry_len"] = df["qry_end"] - df["qry_start"] + 1
    return df

def read_bed(bed_file):
    """Read a BED file and calculate the lengths of regions."""
    df = pd.read_csv(bed_file, sep=r"\s+", header=None, names=["chr", "start", "end"])
    df["length"] = df["end"] - df["start"]
    return df

def intersect_bed(syri_df, bed_df, sv_type):
    """Intersect SYRI data with a BED file and calculate length and count."""
    # Intersect the regions (you can use bedtools intersect if preferred)
    intersected = pd.merge(syri_df, bed_df, how="inner", left_on="ref_chr", right_on="chr")
    intersected = intersected[(intersected["ref_start"] < intersected["end"]) & 
                              (intersected["ref_end"] > intersected["start"])]
    intersected["intersect_len"] = np.minimum(intersected["ref_end"], intersected["end"]) - \
                                   np.maximum(intersected["ref_start"], intersected["start"])
    intersected["type"] = sv_type
    return intersected

def calculate_statistics(syri_df, bed_files, sv_types):
    """Calculate statistics for each SV type."""
    sv_stats = {}
    for sv_type, bed_file in zip(sv_types, bed_files):
        bed_df = read_bed(bed_file)
        intersected_df = intersect_bed(syri_df, bed_df, sv_type)
        total_length = intersected_df["intersect_len"].sum()
        count = len(intersected_df)
        sv_stats[sv_type] = {"total_length": total_length, "count": count}
    return sv_stats

def plot_statistics(sv_stats, output_prefix):
    """Plot statistics for each SV type."""
    types = list(sv_stats.keys())
    lengths = [sv_stats[sv]["total_length"] for sv in types]
    counts = [sv_stats[sv]["count"] for sv in types]

    fig, ax = plt.subplots(figsize=(10, 6))
    width = 0.35
    x = np.arange(len(types))
    ax.bar(x - width / 2, lengths, width, label="Total Length (bp)")
    ax.bar(x + width / 2, counts, width, label="Count")

    ax.set_xlabel("SV Type")
    ax.set_ylabel("Length / Count")
    ax.set_title("SV Type Statistics")
    ax.set_xticks(x)
    ax.set_xticklabels(types)
    ax.legend()

    plt.tight_layout()
    fig.savefig(f"{output_prefix}_sv_statistics.png")
    plt.show()

def main():
    syri_file = "syri.out"  # Path to the syri.out file
    bed_files = [
        "SYN.bed",  # Path to the SYN.bed file
        "INV.bed",  # Path to the INV.bed file
        "TRANS.bed",  # Path to the TRANS.bed file
        "hap1_DUP.bed",  # Path to hap1_DUP.bed file
        "hap2_DUP.bed",  # Path to hap2_DUP.bed file
        "hap1_not_aligned.bed",  # Path to hap1_not_aligned.bed file
        "hap2_not_aligned.bed",  # Path to hap2_not_aligned.bed file
    ]
    sv_types = ["SYN", "INV", "TRANS", "DUP", "NOTAL"]  # SV types to calculate

    syri_df = read_syri(syri_file)
    sv_stats = calculate_statistics(syri_df, bed_files, sv_types)
    plot_statistics(sv_stats, "sv_stats_output")

if __name__ == "__main__":
    main()
