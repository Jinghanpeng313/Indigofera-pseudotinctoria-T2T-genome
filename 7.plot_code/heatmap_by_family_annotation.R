#!/usr/bin/env Rscript
# Expression heatmap with row splits and left-side family colors.

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
  library(ComplexHeatmap)
  library(circlize)
})

expr_file <- "data.csv"
exp <- read.csv(expr_file, row.names = 1, check.names = FALSE)

df_scaled <- t(scale(t(as.matrix(exp))))
df_scaled[is.na(df_scaled)] <- 0
df_scaled[is.nan(df_scaled)] <- 0
df_scaled[is.infinite(df_scaled)] <- 0

gene_prefix <- sub("_.*", "", rownames(df_scaled))

prefix_colors <- c(
  "4CL"="#e56b2f","ANR"="#FFB6C1","ANS"="#7AC5CD","C4H"="#BC8F8F","CHI"="#e8d37d",
  "CHS"="#91b2c6","DFR"="#8da0cb","F33H"="#66c2a5","F35H"="#fc8d62","F3H"="#e78ac3",
  "FLS"="#a6d854","FNSII"="#ffd92f","LAR"="#e5c494","LDOX"="#b3b3b3","PAL"="#1f78b4"
)

row_anno <- rowAnnotation(
  Family = gene_prefix,
  col = list(Family = prefix_colors),
  show_annotation_name = FALSE
)

ht <- Heatmap(df_scaled,
              row_order = rownames(df_scaled),
              row_split = gene_prefix,
              left_annotation = row_anno,
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 8),
              column_names_rot = 60,
              column_split = 6,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              show_row_names = TRUE,
              name = "Exp",
              na_col = "grey",
              col = colorRampPalette(c("#91b2c6","#f7f1e3","#f1a16c","#e56b2f"))(100))

png("heatmap_by_family_annotation.png", width = 1800, height = 1200, res = 150)
draw(ht, heatmap_legend_side = "right")
dev.off()
