#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Transcriptome–Metabolome integration:
# - Read correlation (gene–metabolite) table + DEG + DEM tables
# - Standardize columns
# - Filter significant correlations (|r| >= r_min, q <= q_max)
# - Join with DEG/DEM
# - (Optionally) flip log2FC signs to match desired display
# - Plot 9-quadrant scatter with quadrant shading and labels
# - Extract pairs significant in both layers
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(readr)
  library(stringr)
  library(scales)
  library(janitor)
})

# ===================== File paths =====================
file_corr <- "hap2_top_hits_per_metabolite.csv"  # correlation table: must contain r/p/q; auto-detected
file_deg  <- "hap2_RNA_OI_vs_UOI_limma_all.csv"  # DEG table: GeneID, logFC, adj.P.Val / padj / q ...
file_dem  <- "MET_OI_vs_UOI_limma_all.csv"       # DEM table: metabolite/compound/name, logFC, adj.P.Val ...

# ===================== Helper: fix encoding & trim =====================
fix_chars <- function(x){
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x <- iconv(x, from = "", to = "UTF-8", sub = "")  # drop invalid bytes
    x <- trimws(x)
  }
  x
}
fix_chars_df <- function(df){
  df %>% mutate(across(where(is.character), fix_chars))
}

# ===================== Safe CSV reader (UTF-8 with fallback) =====================
safe_read_csv <- function(path){
  suppressWarnings(read_csv(path, locale = locale(encoding = "UTF-8"), show_col_types = FALSE)) %>%
    fix_chars_df()
}

corr_raw <- safe_read_csv(file_corr) %>% clean_names()
deg_raw  <- safe_read_csv(file_deg)  %>% clean_names()
dem_raw  <- safe_read_csv(file_dem)  %>% clean_names()

# ===================== Identify correlation table (must contain r/p/q) =====================
has_corr_cols <- function(df) all(c("r","p","q") %in% names(df))

if (!has_corr_cols(corr_raw) && has_corr_cols(dem_raw)) {
  message("Detected that 'file_corr' lacks r/p/q while 'file_dem' has them — swapping inputs.")
  tmp <- corr_raw; corr_raw <- dem_raw; dem_raw <- tmp
}

# ===================== Standardizers =====================
# Correlation: return (gene, metabolite, r, p, q)
standardize_corr <- function(df){
  pick <- function(nms, pat) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  gcol <- pick(names(df), "^(gene|gene_id|geneid)$")
  mcol <- pick(names(df), "^(metabolite|compound|name)$")
  if (is.null(gcol)) stop("Correlation table: gene column not found (gene/gene_id/geneid).")
  if (is.null(mcol)) stop("Correlation table: metabolite column not found (metabolite/compound/name).")
  if (!all(c("r","p","q") %in% names(df))) stop("Correlation table must contain r, p, q columns.")

  df %>%
    transmute(
      gene       = fix_chars(.data[[gcol]]),
      metabolite = fix_chars(.data[[mcol]]),
      r = as.numeric(r),
      p = as.numeric(p),
      q = as.numeric(q)
    )
}

# DEG: return (gene, log2fc_gene, padj_gene)
standardize_deg <- function(df){
  rename(df,
         gene        = any_of(c("gene","gene_id","geneid","id")),
         log2fc_gene = any_of(c("log2fc","log_fc","log2_fc","logfc","lfc")),
         padj_gene   = any_of(c("adj_p_val","adj_p_value","adj_p.value","adj_pvalue","padj","fdr","q","adj_p_val_bh"))
  ) %>%
    transmute(
      gene        = fix_chars(gene),
      log2fc_gene = as.numeric(log2fc_gene),
      padj_gene   = as.numeric(padj_gene)
    )
}

# DEM: return (metabolite, log2fc_metab, padj_metab)
standardize_dem <- function(df){
  rename(df,
         metabolite   = any_of(c("metabolite","compound","name")),
         log2fc_metab = any_of(c("log2fc","log_fc","log2_fc","logfc","lfc")),
         padj_metab   = any_of(c("adj_p_val","adj_p_value","adj_p.value","adj_pvalue","padj","fdr","q","adj_p_val_bh"))
  ) %>%
    transmute(
      metabolite    = fix_chars(metabolite),
      log2fc_metab  = as.numeric(log2fc_metab),
      padj_metab    = as.numeric(padj_metab)
    )
}

corr <- standardize_corr(corr_raw)
deg  <- standardize_deg(deg_raw)
dem  <- standardize_dem(dem_raw)

# ===================== Filters / thresholds =====================
# Correlation significance (pair-level)
r_min <- 0.6
q_max <- 0.10

# 9-quadrant thresholds and padj cutoff
x_thr <- 1.0   # |log2FC| threshold for transcriptome
y_thr <- 1.0   # |log2FC| threshold for metabolome
p_thr <- 0.05  # adjusted p-value threshold

# Optional focus metabolites (set only_focus = TRUE to enable)
focus_mets <- c("Calycosin","Butein","Sulfuretin","Chrysoeriol","Genistin")
only_focus <- FALSE

# Whether to flip log2FC signs (as in your original script)
flip_log2fc <- TRUE

# ===================== Merge by significant correlations =====================
merged <- corr %>%
  filter(!is.na(r), !is.na(q), abs(r) >= r_min, q <= q_max) %>%
  inner_join(deg, by = "gene") %>%
  inner_join(dem, by = "metabolite")

if (only_focus) {
  merged <- merged %>% filter(metabolite %in% focus_mets)
}

if (flip_log2fc) {
  merged <- merged %>%
    mutate(
      log2fc_gene  = -log2fc_gene,
      log2fc_metab = -log2fc_metab
    )
}

# ===================== 9-quadrant plotting =====================
plot_9quadrant <- function(
    df,
    x_col = "log2fc_gene",
    y_col = "log2fc_metab",
    padj_gene_col = "padj_gene",
    padj_metab_col = "padj_metab",
    label_col = "metabolite",
    x_thr = 1, y_thr = 1, p_thr = 0.05,
    label_top_n = 15,
    point_size = 2.2, alpha_pt = 0.9
){
  palette <- c(
    "up_up"   = "#e56b2f",
    "up_mid"  = "#f1a16c",
    "up_down" = "#f9cfad",
    "mid_up"  = "#91b2c6",
    "mid_mid" = "#d9d9d9",
    "mid_down"= "#9ecae1",
    "down_up" = "#fdbf6f",
    "down_mid"= "#abd9e9",
    "down_down"="#1a9850"
  )

  dat <- df %>%
    mutate(
      x = .data[[x_col]],
      y = .data[[y_col]],
      sig_gene  = (!is.na(.data[[padj_gene_col]])  & .data[[padj_gene_col]]  < p_thr & abs(x) >= x_thr),
      sig_metab = (!is.na(.data[[padj_metab_col]]) & .data[[padj_metab_col]] < p_thr & abs(y) >= y_thr),
      x_zone = case_when(x >  x_thr ~ "up",  x < -x_thr ~ "down", TRUE ~ "mid"),
      y_zone = case_when(y >  y_thr ~ "up",  y < -y_thr ~ "down", TRUE ~ "mid"),
      quad = paste0(x_zone, "_", y_zone),
      sig_class = case_when(
        sig_gene & sig_metab ~ "both_signif",
        sig_gene & !sig_metab ~ "gene_signif",
        !sig_gene & sig_metab ~ "metab_signif",
        TRUE ~ "none"
      )
    )

  # Background 3x3 tiles
  bg <- tibble(
    xmin=c(-Inf,-x_thr,x_thr,-Inf,-x_thr,x_thr,-Inf,-x_thr,x_thr),
    xmax=c(-x_thr,x_thr,Inf,-x_thr,x_thr,Inf,-x_thr,x_thr,Inf),
    ymin=c(-Inf,-Inf,-Inf,-y_thr,-y_thr,-y_thr,y_thr,y_thr,y_thr),
    ymax=c(-y_thr,-y_thr,-y_thr,y_thr,y_thr,y_thr,Inf,Inf,Inf),
    quad=c("down_down","mid_down","up_down","down_mid","mid_mid","up_mid","down_up","mid_up","up_up")
  ) %>% mutate(fill_col = alpha(palette[quad], 0.13))

  # Label priority: both significant first, then farther from origin
  lab_dat <- dat %>%
    filter(sig_class != "none") %>%
    mutate(dist0 = sqrt(x^2 + y^2)) %>%
    arrange(desc(sig_class == "both_signif"), desc(dist0)) %>%
    slice_head(n = label_top_n)

  ggplot() +
    geom_rect(data = bg, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = bg$fill_col, color = NA) +
    geom_point(data = dat, aes(x=x, y=y, color=quad, shape=sig_class),
               size = point_size, alpha = alpha_pt) +
    geom_vline(xintercept = c(0,-x_thr,x_thr),
               linetype = c("solid","dashed","dashed"), size=0.3, color="grey40") +
    geom_hline(yintercept = c(0,-y_thr,y_thr),
               linetype = c("solid","dashed","dashed"), size=0.3, color="grey40") +
    geom_label_repel(data = lab_dat, aes(x=x, y=y, label=.data[[label_col]]),
                     size=3, label.padding=unit(0.12,"lines"),
                     box.padding=unit(0.25,"lines"), max.overlaps=Inf, min.segment.length=0) +
    scale_color_manual(values = palette, name = "Quadrant") +
    scale_shape_manual(values = c("both_signif"=16,"gene_signif"=17,"metab_signif"=15,"none"=1),
                       name = "Significance") +
    labs(
      x = expression(paste("Transcriptome  ", log[2], "(FC)")),
      y = expression(paste("Metabolome  ", log[2], "(FC)")),
      title = "Transcriptome–Metabolome 9-Quadrant",
      subtitle = paste0("Corr: |r|≥", r_min, ", q≤", q_max,
                        " ; DE/DM: |log2FC|≥", x_thr, "/", y_thr, ", padj<", p_thr)
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
}

# ===================== Plot =====================
p <- plot_9quadrant(
  df = merged,
  x_col = "log2fc_gene",
  y_col = "log2fc_metab",
  padj_gene_col = "padj_gene",
  padj_metab_col = "padj_metab",
  label_col = "metabolite",
  x_thr = x_thr, y_thr = y_thr, p_thr = p_thr,
  label_top_n = 100
)
print(p)

# ===================== Extract pairs significant in both layers =====================
both_tbl <- merged %>%
  filter(
    !is.na(padj_gene),  padj_gene  < p_thr,  abs(log2fc_gene)  >= x_thr,
    !is.na(padj_metab), padj_metab < p_thr,  abs(log2fc_metab) >= y_thr
  ) %>%
  mutate(
    quadrant = dplyr::case_when(
      log2fc_gene  >  x_thr & log2fc_metab >  y_thr ~ "up_up",
      log2fc_gene  >  x_thr & log2fc_metab < -y_thr ~ "up_down",
      log2fc_gene  < -x_thr & log2fc_metab >  y_thr ~ "down_up",
      log2fc_gene  < -x_thr & log2fc_metab < -y_thr ~ "down_down",
      TRUE ~ "central"
    )
  ) %>%
  arrange(dplyr::desc(abs(log2fc_gene) + abs(log2fc_metab))) %>%
  select(
    gene, metabolite,
    log2fc_gene,  padj_gene,
    log2fc_metab, padj_metab,
    r, q, quadrant
  )

# Optional: nicer printing (rounded)
both_tbl_print <- both_tbl %>%
  mutate(
    log2fc_gene  = round(log2fc_gene, 3),
    log2fc_metab = round(log2fc_metab, 3),
    padj_gene    = signif(padj_gene, 3),
    padj_metab   = signif(padj_metab, 3),
    r            = ifelse(is.na(r), NA, round(r, 3)),
    q            = ifelse(is.na(q), NA, signif(q, 3))
  )

print(both_tbl_print, n = min(100, nrow(both_tbl_print)))

# ===================== Summary =====================
cat("\n[Summary]\n")
both_counts <- both_tbl %>% count(quadrant) %>% arrange(desc(n))
print(both_counts)

cat("\nTotal pairs:", nrow(both_tbl),
    " | unique genes:", both_tbl %>% distinct(gene) %>% nrow(),
    " | unique metabolites:", both_tbl %>% distinct(metabolite) %>% nrow(), "\n")

both_genes <- both_tbl %>% distinct(gene) %>% arrange(gene)
both_mets  <- both_tbl %>% distinct(metabolite) %>% arrange(metabolite)
cat("\nFirst 20 genes:\n"); print(head(both_genes, 20), n = 20)
cat("\nFirst 20 metabolites:\n"); print(head(both_mets, 20), n = 20)

# ===================== Export (optional) =====================
readr::write_csv(both_tbl, "both_signif_pairs.csv")
# readr::write_csv(both_genes, "both_signif_genes_unique.csv")
# readr::write_csv(both_mets,  "both_signif_metabolites_unique.csv")
