#!/usr/bin/env bash

THREADS=48
ORTHO_OUT="orthofinder_out"

orthofinder -f proteomes \
  -S diamond_ultra_sens \
  -t ${THREADS} \
  -M msa -A mafft -T fasttree \
  -o ${ORTHO_OUT}
