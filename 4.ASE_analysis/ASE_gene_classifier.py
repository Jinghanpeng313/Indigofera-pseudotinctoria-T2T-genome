#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Classify allele-specific expression (ASE) gene pairs into three categories:

1. Hap     : Haplotype-dominant gene pairs
2. SubNeo  : Subfunctionalized / neofunctionalized pairs
3. NoDiff  : No clear expression bias

Definitions (per gene pair, across tissues):
- For each tissue, we compare hap1 vs hap2 expression (FPKM/TPM).
- A 2-fold change is used as threshold:
    - ratio = hap1 / hap2 >= 2.0  -> hap1 dominant in that tissue
    - ratio = hap1 / hap2 <= 0.5  -> hap2 dominant in that tissue
- A haplotype is considered "dominant" for the pair if it is dominant in at
  least 1/3 of all tissues (at least 1 tissue).

Classification logic:
- Hap:
    - hap1 is dominant in >= min_tissues_for_dominance AND hap2 never dominant
      OR
    - hap2 is dominant in >= min_tissues_for_dominance AND hap1 never dominant
- SubNeo:
    - hap1 dominant in at least 1 tissue AND hap2 dominant in at least 1 tissue
- NoDiff:
    - All other cases (no clear dominance pattern)
"""

import pandas as pd
import numpy as np


def classify_ase_genes(gene_pairs_file, hap1_fpkm_file, hap2_fpkm_file, output_file):
    """
    Classify ASE gene pairs into three categories:
      1. Hap     : Haplotype-dominant gene pairs
      2. SubNeo  : Subfunctionalized / neofunctionalized pairs
      3. NoDiff  : No obvious dominant pattern

    Parameters
    ----------
    gene_pairs_file : str
        TSV file containing at least two columns: 'hap1_gene' and 'hap2_gene'.
    hap1_fpkm_file : str
        TSV expression file (FPKM/TPM) for haplotype 1; must contain column 'GeneID'.
    hap2_fpkm_file : str
        TSV expression file (FPKM/TPM) for haplotype 2; must contain column 'GeneID'.
    output_file : str
        Output TSV file for classification results.

    Returns
    -------
    pandas.DataFrame or None
        DataFrame with classification results, or None if no pairs could be processed.
    """

    print("Loading input data...")

    # Read gene pair file
    pairs_df = pd.read_csv(gene_pairs_file, sep="\t")
    print(f"Gene pair table shape: {pairs_df.shape}")

    # Read FPKM/TPM data
    hap1_df = pd.read_csv(hap1_fpkm_file, sep="\t", index_col="GeneID")
    hap2_df = pd.read_csv(hap2_fpkm_file, sep="\t", index_col="GeneID")
    print(f"hap1 expression table shape: {hap1_df.shape}")
    print(f"hap2 expression table shape: {hap2_df.shape}")

    # Get tissue/sample columns (exclude 'GeneID' or 'gene_id' if present)
    tissues = [col for col in hap1_df.columns if col not in ["GeneID", "gene_id"]]
    num_tissues = len(tissues)
    # At least 1/3 of tissues must show dominance; minimum is 1 tissue
    min_tissues_for_dominance = max(1, num_tissues // 3)

    print(f"Detected {num_tissues} tissues/samples: {tissues}")
    print(f"Minimum tissues required for haplotype dominance: {min_tissues_for_dominance}")

    results = []
    skipped_genes = 0

    print(f"Analyzing {len(pairs_df)} gene pairs...")

    for idx, row in pairs_df.iterrows():
        hap1_gene = row["hap1_gene"]
        hap2_gene = row["hap2_gene"]

        # Check if both genes exist in the expression tables
        if hap1_gene not in hap1_df.index or hap2_gene not in hap2_df.index:
            skipped_genes += 1
            continue

        # Get expression vectors for the two genes
        hap1_expr = hap1_df.loc[hap1_gene, tissues].values.astype(float)
        hap2_expr = hap2_df.loc[hap2_gene, tissues].values.astype(float)

        dominant_hap1_tissues = []
        dominant_hap2_tissues = []
        ratios = []

        for i, tissue in enumerate(tissues):
            expr1 = hap1_expr[i]
            expr2 = hap2_expr[i]

            # Avoid division by zero
            if expr1 == 0 and expr2 == 0:
                ratio = 1.0  # both zero, treat as equal
            elif expr2 == 0:
                ratio = float("inf")  # hap1 >> hap2
            elif expr1 == 0:
                ratio = 0.0  # hap2 >> hap1
            else:
                ratio = expr1 / expr2

            ratios.append(ratio)

            # 2-fold change threshold for dominance in a tissue
            if ratio >= 2.0:
                dominant_hap1_tissues.append(tissue)
            elif ratio <= 0.5:
                dominant_hap2_tissues.append(tissue)

        # Classification logic
        classification = "NoDiff"
        details = ""

        hap1_dominant_count = len(dominant_hap1_tissues)
        hap2_dominant_count = len(dominant_hap2_tissues)

        # Haplotype-dominant gene pairs (Hap)
        if (
            hap1_dominant_count >= min_tissues_for_dominance
            and hap2_dominant_count == 0
        ):
            classification = "Hap"
            details = f"hap1 dominant in {hap1_dominant_count} tissues"
        elif (
            hap2_dominant_count >= min_tissues_for_dominance
            and hap1_dominant_count == 0
        ):
            classification = "Hap"
            details = f"hap2 dominant in {hap2_dominant_count} tissues"

        # Subfunctionalized / neofunctionalized pairs (SubNeo)
        elif hap1_dominant_count >= 1 and hap2_dominant_count >= 1:
            classification = "SubNeo"
            details = (
                f"hap1 dominant in {hap1_dominant_count} tissues, "
                f"hap2 dominant in {hap2_dominant_count} tissues"
            )

        # Save results for this pair
        results.append(
            {
                "pair_id": f"{hap1_gene}|{hap2_gene}",
                "hap1_gene": hap1_gene,
                "hap2_gene": hap2_gene,
                "classification": classification,
                "details": details,
                "hap1_dominant_count": hap1_dominant_count,
                "hap2_dominant_count": hap2_dominant_count,
                "hap1_dominant_tissues": ",".join(dominant_hap1_tissues)
                if dominant_hap1_tissues
                else "None",
                "hap2_dominant_tissues": ",".join(dominant_hap2_tissues)
                if dominant_hap2_tissues
                else "None",
                "total_tissues": num_tissues,
            }
        )

    print(f"Skipped {skipped_genes} gene pairs not found in expression tables.")

    if not results:
        print("ERROR: No gene pairs were successfully processed.")
        print("Please check that gene IDs match between pair and expression files.")
        return None

    # Build result DataFrame
    results_df = pd.DataFrame(results)
    print(f"Successfully processed {len(results_df)} gene pairs.")

    # Count by classification
    class_counts = results_df["classification"].value_counts()
    print("\nClassification summary:")
    for cls, count in class_counts.items():
        print(f"{cls}: {count} gene pairs")

    # Save results
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nResults written to: {output_file}")

    return results_df


# Example usage
if __name__ == "__main__":
    # Input file paths
    gene_pairs_file = "same_cds_ase_gene_piars.txt"  # gene pairs with 'hap1_gene' and 'hap2_gene'
    hap1_fpkm_file = "hap1.fpkm.txt"                 # hap1 expression table
    hap2_fpkm_file = "hap2.fpkm.txt"                 # hap2 expression table
    output_file = "ase_gene_classification_results.txt"

    # Run classification
    results = classify_ase_genes(
        gene_pairs_file, hap1_fpkm_file, hap2_fpkm_file, output_file
    )

    if results is not None:
        print("\nTop 10 classified gene pairs:")
        print(results.head(10))
