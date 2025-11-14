#!/usr/bin/env bash
set -euo pipefail
GENOME=Ipse_Hap1
CDS=${GENOME}.cds.fa
OUT=dupgen_${GENOME}/results
THREADS=16

# For each duplication mode, compute Ka/Ks for gene pairs (KaKs_Calculator or PAML/codeml)
for MODE in WGD TD PD DSD; do
  PAIRS=${OUT}/${MODE}.txt   # assume two columns: geneA  geneB
  mkdir -p kaks_${MODE}
  # Prepare pairwise CDS for KaKs_Calculator 2.0 (AXT format)
  while read -r a b; do
    seqkit grep -p "^${a}$" ${CDS} > tmp_a.fa
    seqkit grep -p "^${b}$" ${CDS} > tmp_b.fa
    pair2axt tmp_a.fa tmp_b.fa ${MODE}_${a}_${b}.axt   # use your favorite converter; or write a small python align+pal2nal
    mv ${MODE}_${a}_${b}.axt kaks_${MODE}/
  done < ${PAIRS}

  KaKs_Calculator -i kaks_${MODE} -o kaks_${MODE}.results -m YN
done
