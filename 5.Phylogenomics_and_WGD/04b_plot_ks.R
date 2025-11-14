#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop("Usage: 04b_plot_ks.R ks_list.txt out.pdf")
# ks_list.txt: lines like "Ipse_Hap1.Ipse_Hap1.ks.tsv\tI. pseudotinctoria"
x <- fread(args[1], header=FALSE, sep="\t")
plots <- list()
all <- list()
for(i in 1:nrow(x)){
  ks <- fread(x$V1[i])
  ks <- ks[is.finite(ks$Ks) & ks$Ks>0 & ks$Ks<5]
  ks$species <- x$V2[i]
  all[[i]] <- ks
}
dt <- rbindlist(all)
p <- ggplot(dt, aes(Ks, ..density.., color=species)) + 
  geom_freqpoly(binwidth=0.03, size=0.6) + theme_bw() +
  xlab("Ks") + ylab("Density")
ggsave(args[2], p, width=6, height=4)
