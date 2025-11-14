#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(data.table)
})

# Inputs:
#  - target gene list file (one gene ID per line) e.g. Ipse_Hap1_species_specific.txt
#  - background (all protein-coding genes from Ipse_Hap1)
#  - gene2go mapping (GMT or two-column mapping)
#  - KEGG mapping: use KEGG API or a custom mapping if organism not in KEGG
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<4) stop("Usage: 02_enrichment_clusterProfiler.R target.txt bg.txt gene2go.tsv gene2kegg.tsv")
target <- fread(args[1], header=FALSE)$V1
bg     <- fread(args[2], header=FALSE)$V1
g2go   <- fread(args[3], header=TRUE) # columns: gene,go
g2ko   <- fread(args[4], header=TRUE) # columns: gene,ko

# GO enrichment (custom mapping)
go_list <- split(g2go$gene, g2go$go)
ego <- enricher(gene=target, universe=bg, TERM2GENE=g2go[,.(go,gene)], pAdjustMethod="BH", qvalueCutoff=0.05)
fwrite(as.data.frame(ego), "go_enrichment.tsv", sep="\t")

# KEGG enrichment (custom KO mapping, then KO->pathway via KEGG API or local table)
# expect a TERM2GENE table 'ko2pathway.tsv': columns pathway,ko
if(file.exists("ko2pathway.tsv")){
  ko2path <- fread("ko2pathway.tsv")
  # convert gene->ko to TERM2GENE by joining through KO
  t2g <- merge(g2ko, ko2path, by="ko")[,.(pathway,gene)]
  eke <- enricher(gene=target, universe=bg, TERM2GENE=t2g, pAdjustMethod="BH", qvalueCutoff=0.05)
  fwrite(as.data.frame(eke), "kegg_enrichment.tsv", sep="\t")
}
