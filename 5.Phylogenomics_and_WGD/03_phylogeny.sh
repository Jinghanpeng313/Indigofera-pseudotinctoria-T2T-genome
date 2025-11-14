#!/usr/bin/env bash
set -euo pipefail
THREADS=48
ALN_DIR="align_sc_aa"
SC_DIR="orthofinder_out/Results_*/Single_Copy_Orthologue_Sequences"
mkdir -p ${ALN_DIR}

# Align each single-copy orthogroup (amino acids), trim, and concatenate
for faa in ${SC_DIR}/*.fa; do
  base=$(basename "$faa" .fa)
  mafft --auto --thread ${THREADS} "$faa" > ${ALN_DIR}/${base}.aln.fa
  trimal -automated1 -in ${ALN_DIR}/${base}.aln.fa -out ${ALN_DIR}/${base}.trim.fa
done

# Concatenate (FASconCAT-G or AMAS)
AMAS.py concat -f fasta -d aa -i ${ALN_DIR}/*.trim.fa -t supermatrix.aa.fas -p partitions.txt

# ML tree with IQ-TREE2
iqtree2 -s supermatrix.aa.fas -p partitions.txt \
  -m MFP+MERGE -B 1000 -T ${THREADS} -pre iqtree_sc

# Time calibration: produce ultrametric tree (treePL or chronos)
# Example with 'chronos' (R/ape) using secondary calibration points
cat > chronos_calibrate.R <<'RS'
library(ape)
tr <- read.tree("iqtree_sc.treefile")
# Provide node calibration(s): example placeholders (replace with real calibrations!)
# e.g., mrca(Gmax,Vrad) ~ 10-15 Mya; mrca(Ipse_Hap1,Ccaj) ~ 21.9-40.2 Mya
# Here we just run chronos with default lambda; for publication use proper priors.
utr <- chronos(tr, lambda=1)
write.tree(utr, file="species_tree_ultrametric.nwk")
RS
Rscript chronos_calibrate.R
