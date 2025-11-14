#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ----------------------------
# Inputs & knobs
# ----------------------------
in_bedgraph <- "hap2_20kb_pca.PC1.bedGraph"  # columns: chrom start end PC1
bin_size    <- 50000                          # for optional aggregation (bp)
use_agg     <- FALSE                          # TRUE = pick max |PC1| per bin
bar_alpha   <- 0.85                           # bar transparency
bar_width   <- bin_size                       # bar width (bp) for plotting

# ----------------------------
# Load bedGraph
# ----------------------------
df <- fread(in_bedgraph, header = FALSE, sep = "\t",
            col.names = c("chrom","start","end","pc1"))

# Midpoint for plotting and sanity checks
df <- df %>%
  mutate(position = (start + end) / 2)

# ----------------------------
# Optional: aggregate within fixed bins
# (choose one record per bin with maximum |PC1|)
# ----------------------------
if (isTRUE(use_agg)) {
  df <- df %>%
    mutate(bin_start = floor(position / bin_size) * bin_size,
           bin_end   = bin_start + bin_size) %>%
    group_by(chrom, bin_start, bin_end) %>%
    arrange(desc(abs(pc1)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      start    = bin_start,
      end      = bin_end,
      position = (start + end) / 2
    ) %>%
    select(chrom, start, end, pc1, position)
}

# ----------------------------
# A/B calls from PC1 sign
# ----------------------------
df <- df %>%
  mutate(compartment = ifelse(pc1 >= 0, "A", "B"))

# ----------------------------
# Plot per chromosome
# -------------------
