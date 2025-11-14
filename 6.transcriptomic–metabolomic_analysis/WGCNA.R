#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(tidyverse)
  library(readxl)
  library(optparse)
})

# ------------------------------------------------------------------
# CLI: all parameters are provided at runtime (nothing hard-coded)
# ------------------------------------------------------------------
opt_list <- list(
  make_option("--fpkm",        type="character", help="Gene expression matrix (TSV). Rows=GeneID, columns=samples."),
  make_option("--metabolites", type="character", help="Metabolite abundance table (XLSX). Rows=metabolites, columns=samples."),
  make_option("--metabolite-id-col", type="character", default=NA,
              help="Name of metabolite ID column in the XLSX. If omitted, the first column is used."),
  make_option("--targets",     type="character", help="Plain-text file with one metabolite name per line to keep."),
  make_option("--outdir",      type="character", help="Output directory."),
  make_option("--minClusterSize", type="integer", default=50, help="Minimum module size."),
  make_option("--mergeCutHeight", type="double", default=0.25, help="Module eigengene merge threshold."),
  make_option("--topExprQuantile", type="double", default=0.30,
              help="Keep genes with mean expression >= this quantile (0..1)."),
  make_option("--qc-pattern",  type="character", default="QC",
              help="Regex to remove QC columns from metabolite data."),
  make_option("--softPower",   type="integer", default=NA,
              help="If provided, use this soft-threshold; otherwise auto-pick via pickSoftThreshold."),
  make_option("--log_expr_offset", type="double", default=1.0, help="Offset for log2 transform of expression."),
  make_option("--log_metab_offset", type="double", default=1e-6, help="Offset for log2 transform of metabolites.")
)

opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$fpkm) || is.null(opt$metabolites) || is.null(opt$targets) || is.null(opt$outdir)) {
  stop("Required: --fpkm --metabolites --targets --outdir", call.=FALSE)
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# 1) Load data (expression matrix and metabolite table)
# ------------------------------------------------------------------
fpkm <- fread(opt$fpkm, data.table = FALSE)
stopifnot("GeneID" %in% colnames(fpkm))
rownames(fpkm) <- fpkm$GeneID
fpkm <- fpkm[, setdiff(colnames(fpkm), "GeneID"), drop = FALSE]

expr_samples <- colnames(fpkm)

metab_raw <- as.data.frame(read_excel(opt$metabolites))
metab_colnames <- colnames(metab_raw)
id_col <- if (!is.na(opt$metabolite_id_col)) {
  if (!opt$metabolite_id_col %in% metab_colnames)
    stop("--metabolite-id-col not found in metabolite table.")
  opt$metabolite_id_col
} else {
  metab_colnames[1]
}

# match sample columns by name overlap
sample_regex <- paste0("(^", paste(gsub("([\\W])", "\\\\\\1", expr_samples), collapse="|"), "$)")
matched_samples <- metab_colnames[grepl(sample_regex, metab_colnames)]
if (length(matched_samples) == 0) stop("No matching sample columns between metabolite and expression matrices.")

# ------------------------------------------------------------------
# 2) Keep only targeted metabolites (provided as a file)
# ------------------------------------------------------------------
target_metabolites <- readr::read_lines(opt$targets) %>% setdiff(c("", NA))
if (!(id_col %in% colnames(metab_raw))) stop("Metabolite ID column not found in the XLSX.")
metab_raw_filtered <- metab_raw[metab_raw[[id_col]] %in% target_metabolites, , drop = FALSE]
if (nrow(metab_raw_filtered) == 0) stop("None of the target metabolites were found in the table.")

# ------------------------------------------------------------------
# 3) Build metabolite matrix aligned to expression samples
# ------------------------------------------------------------------
metab_mat <- metab_raw_filtered[, matched_samples, drop = FALSE]
rownames(metab_mat) <- metab_raw_filtered[[id_col]]

# remove QC columns
qc_cols <- grep(opt$qc_pattern, colnames(metab_mat), value = TRUE, ignore.case = TRUE)
if (length(qc_cols) > 0) metab_mat <- metab_mat[, setdiff(colnames(metab_mat), qc_cols), drop = FALSE]

# align order and names to expression samples
# match ignoring case, then reorder both to common set
metab_names_lower <- tolower(colnames(metab_mat))
expr_names_lower  <- tolower(expr_samples)
map <- match(expr_names_lower, metab_names_lower)
keep <- which(!is.na(map))
if (length(keep) == 0) stop("No common samples after QC removal.")

fpkm     <- fpkm[, expr_samples[keep], drop = FALSE]
metab_mat <- metab_mat[, map[keep], drop = FALSE]
colnames(metab_mat) <- colnames(fpkm)  # exact same sample names

# sample x metabolite, log2 and z-score
metab_for_traits <- t(as.matrix(metab_mat))
metab_log <- log2(metab_for_traits + opt$log_metab_offset)
metab_z   <- scale(metab_log)

# ------------------------------------------------------------------
# 4) Process expression matrix
# ------------------------------------------------------------------
expr <- log2(as.matrix(fpkm) + opt$log_expr_offset)

# filter by mean expression quantile (e.g., keep >= 30th percentile)
meanFPKM <- rowMeans(pmax(2^expr - opt$log_expr_offset, 0), na.rm=TRUE)
threshold <- as.numeric(quantile(meanFPKM, opt$topExprQuantile, na.rm = TRUE))
keepGenes <- meanFPKM >= threshold
expr_filt <- expr[keepGenes, , drop = FALSE]

datExpr <- t(expr_filt)  # samples x genes
gsg <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# ------------------------------------------------------------------
# 5) Network construction and module–trait association
# ------------------------------------------------------------------
allowWGCNAThreads()
# pick soft-threshold automatically if not supplied
if (is.na(opt$softPower)) {
  powers <- c(1:10, seq(12, 30, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, networkType = "signed")
  # choose smallest power with scale-free R^2 >= 0.8, otherwise the max of tested powers
  r2 <- sft$fitIndices[, "SFT.R.sq"]
  cand <- which(r2 >= 0.8)
  softPower <- if (length(cand) > 0) sft$fitIndices[min(cand), "Power"] else max(powers)
} else {
  softPower <- opt$softPower
}

adjacency <- adjacency(datExpr, power = softPower, type = "signed")
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 3,
  pamRespectsDendro = FALSE,
  minClusterSize = opt$minClusterSize
)
dynamicColors <- labels2colors(dynamicMods)

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
merged <- mergeCloseModules(datExpr, dynamicColors, cutHeight = opt$mergeCutHeight, verbose = 0)
mergedColors <- merged$colors
mergedMEs <- orderMEs(merged$newMEs)

moduleTraitCor <- cor(mergedMEs, metab_z, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# ------------------------------------------------------------------
# 6) Write outputs
# ------------------------------------------------------------------
write.table(round(moduleTraitCor, 4),
            file = file.path(opt$outdir, "moduleTraitCor.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

write.table(round(moduleTraitP, 6),
            file = file.path(opt$outdir, "moduleTraitPvalues.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

write.table(
  data.frame(Gene = colnames(datExpr), Module = mergedColors),
  file = file.path(opt$outdir, "gene_module_colors.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

png(file.path(opt$outdir, "module_trait_heatmap.png"), width = 6000, height = 5500, res = 150)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(metab_z),
  yLabels = colnames(mergedMEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("#91b2c6", "#f9cfad", "#f1a16c"))(100),
  textMatrix = signif(moduleTraitCor, 2),
  main = "Module–trait relationships"
)
dev.off()

png(file.path(opt$outdir, "gene_dendrogram_modules.png"), width = 2200, height = 1400, res = 150)
plotDendroAndColors(
  geneTree,
  mergedColors,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()
