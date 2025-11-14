#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(stringr)})
modes <- c("WGD","TD","PD","DSD")
all <- list()
for(m in modes){
  f <- sprintf("kaks_%s.results", m)
  if(!file.exists(f)) next
  dt <- fread(f)
  # Expect columns like: Sequence, Ka, Ks, Ka/Ks
  dt$mode <- m
  all[[m]] <- dt
}
d <- rbindlist(all, use.names=TRUE, fill=TRUE)
d <- d[is.finite(Ka) & is.finite(Ks) & Ks>0]
p <- ggplot(d, aes(x=mode, y=`Ka/Ks`, fill=mode)) + geom_violin() + geom_boxplot(width=0.1, outlier.size=0.4) + theme_bw()
ggsave("kaks_by_mode.pdf", p, width=5, height=4)
