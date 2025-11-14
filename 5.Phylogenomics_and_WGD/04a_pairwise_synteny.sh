#!/usr/bin/env bash
set -euo pipefail
# For each species vs Ipse_Hap1
A=Ipse_Hap1
B=Mtru   # change per species in a loop
THREADS=32

# 1) Create BED (gene positions) from GFF3
python -m jcvi.formats.gff bed --type=mRNA --key=ID ${A}.gff3 -o ${A}.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID ${B}.gff3 -o ${B}.bed

# 2) Diamond blastp (A vs B)
diamond makedb --in ${B}.protein.fa -d ${B}
diamond blastp -d ${B} -q ${A}.protein.fa -o ${A}.${B}.blast \
  -p ${THREADS} --evalue 1e-5 --outfmt 6

# 3) Anchors + collinearity
python -m jcvi.compara.catalog ortholog ${A} ${B} --no_strip_names --cscore=0.99 --dbtype diamond

# 4) Ks for anchors (within-species A vs A for self-WGD; A vs B for speciation Ks)
python -m jcvi.compara.ks ${A}.${B}.anchors ${A}.cds.fa ${B}.cds.fa -o ${A}.${B}.ks.tsv
