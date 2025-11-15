#!/usr/bin/env Rscript

# ---------------------------------------------------------------------
# Coupled correlation plots:
# - Read module–module correlation (ME_cor) + p-values (ME_pvalue)
# - Read module×metabolite R and P matrices (wide), convert to long
# - Clip negative R to 0; keep metabolites with any R>=0.75
# - Bin R (rd) and P (pd) for edge aesthetics
# - Draw two quadrant-style correlation plots (lower/upper triangles)
# - Combine and export
# ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(optparse)
  library(ggnewscale)     # for new_scale()
  # The following come from the 'ggcor' ecosystem:
  #   qcorrplot(), geom_square(), geom_couple(), nice_curvature()
  # Install/load: remotes::install_github("lianchunfei/ggcor")
  library(ggcor)
  library(cowplot)
})

# ------------------------- CLI options -------------------------
opt_list <- list(
  make_option("--me_cor",      type="character", help="Module–module correlation matrix (TSV). Rows and columns are modules."),
  make_option("--me_pvalue",   type="character", help="Module–module p-value matrix (TSV). Same dimension/order as --me_cor."),
  make_option("--r_mat",       type="character", help="Module×metabolite R matrix (CSV/TSV). First column = Module, others = metabolites."),
  make_option("--p_mat",       type="character", help="Module×metabolite P matrix (CSV/TSV). First column = Module, others = metabolites."),
  make_option("--outdir",      type="character", help="Output directory."),
  make_option("--r_keep",      type="double", default=0.7, help="Keep metabolites if any R>= this threshold (after clipping negatives to 0)."),
  make_option("--rd_bins",     type="character", default="<0.2,0.2-0.4,>=0.4", help="Bins for Mantel's r (display)."),
  make_option("--pd_bins",     type="character", default="<0.01,0.01-0.05,>=0.05", help="Bins for P values (display)."),
  make_option("--label_size",  type="double", default=3.5, help="Label font size for edges."),
  make_option("--dpi",         type="integer", default=300, help="Figure DPI."),
  make_option("--width",       type="double", default=10, help="Output width (inches) for combined panel."),
  make_option("--height",      type="double", default=8,  help="Output height (inches) for combined panel."),
  make_option("--outfile",     type="character", default="mantel_combo.png", help="Output file name inside --outdir.")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$me_cor), !is.null(opt$me_pvalue), !is.null(opt$r_mat), !is.null(opt$p_mat), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------- helpers -------------------------
read_table_auto <- function(path){
  ext <- tools::file_ext(path)
  if (tolower(ext) %in% c("csv")) {
    readr::read_csv(path, show_col_types = FALSE) %>% clean_names()
  } else {
    readr::read_tsv(path, show_col_types = FALSE) %>% clean_names()
  }
}

# Convert a wide "Module × metabolite" matrix into long form:
# requires a column that identifies the module (first column after clean_names()).
wide_to_long <- function(df){
  stopifnot(ncol(df) >= 2)
  mod_col <- colnames(df)[1]
  df %>%
    rename(Module = !!mod_col) %>%
    pivot_longer(cols = -Module, names_to = "Metabolite", values_to = "value")
}

# ------------------------- 1) Read matrices -------------------------
# ME_cor / ME_pvalue (TSV with row/col names)
ME_cor <- read.table(opt$me_cor, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) %>% as.matrix()
ME_pvalue <- read.table(opt$me_pvalue, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) %>% as.matrix()

# R and P matrices (CSV/TSV): wide (Module × metabolites)
r_wide <- read_table_auto(opt$r_mat)
p_wide <- read_table_auto(opt$p_mat)

# ------------------------- 2) Melt R/P to long & merge -------------------------
r_long <- wide_to_long(r_wide) %>% rename(r = value)
p_long <- wide_to_long(p_wide) %>% rename(p = value)

mantel_data <- r_long %>% left_join(p_long, by = c("Module", "Metabolite"))

# ------------------------- 3) Clip negatives; keep metabolites -------------------------
mantel_data <- mantel_data %>%
  mutate(r_modified = pmax(0, as.numeric(r)),
         p = as.numeric(p))

keep_mets <- mantel_data %>%
  group_by(Metabolite) %>%
  summarise(any_high = any(r_modified >= opt$r_keep, na.rm = TRUE), .groups="drop") %>%
  filter(any_high) %>%
  pull(Metabolite)

mantel_filtered <- mantel_data %>% filter(Metabolite %in% keep_mets)

# ------------------------- 4) Bin r (rd) and p (pd), finalize table -------------------------
# Bins are cosmetic for edge legends. You can change cut points via CLI if desired.
rd_breaks <- c(-Inf, 0.2, 0.4, Inf)
rd_labels <- c("< 0.2", "0.2 - 0.4", ">= 0.4")
pd_breaks <- c(-Inf, 0.01, 0.05, Inf)
pd_labels <- c("< 0.01", "0.01 - 0.05", ">= 0.05")

mantel_formatted <- mantel_filtered %>%
  mutate(
    rd = cut(r_modified, breaks = rd_breaks, labels = rd_labels, right = FALSE, include.lowest = TRUE),
    pd = cut(p,           breaks = pd_breaks, labels = pd_labels, right = FALSE, include.lowest = TRUE)
  ) %>%
  rename(spec = Metabolite, env = Module) %>%
  transmute(spec, env, r = r_modified, p, rd, pd)

# Quick console summary
message("[INFO] Kept metabolites: ", paste(sort(unique(mantel_formatted$spec)), collapse = ", "))
message("[INFO] n edges (spec–env pairs): ", nrow(mantel_formatted))

# ------------------------- 5) Build two triangle plots -------------------------
# Color palette for heat tile
fill_pal <- colorRampPalette(c("#91b2c6","#f7f1e3","#f1a16c"))(100)

# Lower triangle
p_lower <- qcorrplot(
  ME_cor,
  diag = TRUE,
  type = "lower",
  grid_col = "grey",
  grid_size = 0.6,
  is_corr = TRUE
) +
  geom_square(colour = "grey", size = 0.4) +
  geom_square() +
  guides(size = "none") +
  scale_size_continuous(range = c(0, 10)) +
  new_scale("size") +
  scale_fill_gradientn(colors = fill_pal, na.value = "white") +
  geom_couple(
    aes(colour = pd, size = rd),
    data = mantel_formatted,
    label.colour = "black",
    curvature = nice_curvature(),
    label.fontface = 0,
    label.size = opt$label_size,
    drop = TRUE,
    node.colour = c("white", "white"),
    node.fill   = c("#984EA3", "#5785C1"),
    node.size   = c(4.5, 4),
    node.shape  = c(23, 21)
  ) +
  coord_cartesian(clip = "off") +
  scale_size_manual(values = c(0.5, 1.5, 1)) +
  scale_colour_manual(values = c("#984EA3", "grey","#4DAF4A")) +
  guides(
    fill   = guide_colorbar(title = "Pearson's r", order = 1),
    color  = guide_legend(title = "*P* value", order = 2),
    size   = guide_legend(title = "Mantel's r", order = 3)
  ) +
  theme(
    plot.margin       = unit(c(0.5, 0, 0.5, 0.5), "cm"),
    axis.text.y       = element_blank(),
    axis.text.x       = element_text(size = 9, color = "black"),
    axis.ticks        = element_blank(),
    panel.background  = element_blank(),
    legend.key        = element_blank(),
    legend.background = element_blank()
  )

# Upper triangle
p_upper <- qcorrplot(
  ME_cor,
  diag = TRUE,
  type = "upper",
  grid_col = "grey",
  grid_size = 0.6,
  is_corr = TRUE
) +
  geom_square(colour = "grey", size = 0.4) +
  geom_square() +
  guides(size = "none") +
  scale_size_continuous(range = c(0, 10)) +
  new_scale("size") +
  scale_fill_gradientn(colors = fill_pal, na.value = "white") +
  geom_couple(
    aes(colour = pd, size = rd),
    data = mantel_formatted,
    label.colour = "black",
    curvature = nice_curvature(),
    label.fontface = 0,
    label.size = opt$label_size,
    drop = TRUE,
    node.colour = c("white", "white"),
    node.fill   = c("#984EA3", "#5785C1"),
    node.size   = c(4.5, 4),
    node.shape  = c(23, 21)
  ) +
  coord_cartesian(clip = "off") +
  scale_size_manual(values = c(0.5, 1.5, 1)) +
  scale_colour_manual(values = c("#984EA3", "grey","#4DAF4A")) +
  guides(
    fill   = guide_colorbar(title = "Pearson's r", order = 1),
    color  = guide_legend(title = "*P* value", order = 2),
    size   = guide_legend(title = "Mantel's r", order = 3)
  ) +
  theme(
    plot.margin       = unit(c(0.5, 0, 0.5, 0.5), "cm"),
    axis.text.y       = element_blank(),
    axis.text.x       = element_text(size = 9, color = "black"),
    axis.ticks        = element_blank(),
    panel.background  = element_blank(),
    legend.key        = element_blank(),
    legend.background = element_blank()
  )

# ------------------------- 6) Combine & save -------------------------
combo <- ggdraw() +
  draw_plot(p_upper, x = 0.30, y = 0.21, scale = 0.60) +
  draw_plot(p_lower, x = -0.21, y = -0.14, scale = 0.60)

outfile <- file.path(opt$outdir, opt$outfile)
ggsave(outfile, combo, width = opt$width, height = opt$height, dpi = opt$dpi)
message("[OK] Saved: ", outfile)
