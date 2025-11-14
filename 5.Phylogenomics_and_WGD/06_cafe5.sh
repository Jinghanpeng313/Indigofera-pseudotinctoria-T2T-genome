#!/usr/bin/env bash
set -euo pipefail
# Inputs:
#   Orthogroups_GeneCount.tsv (OrthoFinder)
#   species_tree_ultrametric.nwk (from step 3)
# Convert counts to CAFE5 input: remove header row/col names formatting as needed.

python - <<'PY'
import pandas as pd
df = pd.read_csv("orthofinder_out/Results_*/Orthogroups/Orthogroups_GeneCount.tsv", sep="\t")
df = df.rename(columns={"Orthogroup":"FAMILY"})
df.to_csv("cafe_families.tsv", sep="\t", index=False)
PY

cafe5 -i cafe_families.tsv -t species_tree_ultrametric.nwk -o cafe5_out \
      --pvalue 0.05 --threads 32

# Significant expansions/contractions will be in cafe5_out/
