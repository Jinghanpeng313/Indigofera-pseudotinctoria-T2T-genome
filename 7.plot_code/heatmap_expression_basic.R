#!/usr/bin/env Rscript
# Basic expression heatmap with scaling.

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
  library(ComplexHeatmap)
  library(circlize)
})

expr_file <- "data.csv"     # rows = genes, cols = samples
exp <- read.csv(expr_file, row.names = 1, check.names = FALSE)

df_scaled <- t(scale(t(as.matrix(exp))))
df_scaled[is.na(df_scaled)] <- 0
df_scaled[is.nan(df_scaled)] <- 0
df_scaled[is.infinite(df_scaled)] <- 0

my_colors <- colorRampPalette(c("#91b2c6","#f9cfad","#f1a16c","#e56b2f"))(100)

ht <- Heatmap(df_scaled,
              row_order = rownames(df_scaled),
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 5),
              column_names_rot = 60,
              column_split = 6,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              show_row_names = TRUE,
              name = "Exp",
              na_col = "white",
              col = my_colors)

png("heatmap_expression_basic.png", width = 1600, height = 1200, res = 150)
draw(ht, heatmap_legend_side = "right")
dev.off()
