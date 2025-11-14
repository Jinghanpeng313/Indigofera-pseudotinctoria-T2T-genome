#!/usr/bin/env python3
import pandas as pd, sys, re, os
# Inputs from OrthoFinder
og = sys.argv[1]  # Orthogroups.tsv
cnt = sys.argv[2] # Orthogroups_GeneCount.tsv
species = sys.argv[3].split(",")  # 16 species IDs, order-insensitive

# strict single-copy families: each listed species has exactly 1 gene
df = pd.read_csv(og, sep="\t")
dc = pd.read_csv(cnt, sep="\t")

# Identify strict single-copy OGs across all species
cols = [c for c in dc.columns if c != "Orthogroup"]
sc = dc[(dc[cols] == 1).all(axis=1)]
print(f"Strict single-copy orthogroups: {len(sc)}")

# Species-specific (unique paralogs + unclustered)
# From Orthofinder gene count table: per species genes not assigned or private OGs
# A simple proxy: per OG, count species present
presence = (dc[cols] > 0).sum(axis=1)
ss_ogs = dc.loc[presence.eq(1),"Orthogroup"]

# Prepare a per-species summary
spec_rows=[]
for s in cols:
    # unique paralogs: OGs present only in one species and with count>1
    up = dc.loc[(presence.eq(1)) & (dc[s]>1), s].sum()
    # unclustered: OrthoFinder’s “UnassignedGenes.tsv” is better; here we read if present
    unclustered = 0
    ua = os.path.join(os.path.dirname(og), "UnassignedGenes.tsv")
    if os.path.exists(ua):
        u = pd.read_csv(ua, sep="\t")
        unclustered = (u["Species"]==s).sum()
    spec_rows.append([s, up, unclustered, up+unclustered])

out = pd.DataFrame(spec_rows, columns=["species","unique_paralogs","unclustered","species_specific"])
out.to_csv("orthofinder_species_specific_summary.tsv", sep="\t", index=False)
