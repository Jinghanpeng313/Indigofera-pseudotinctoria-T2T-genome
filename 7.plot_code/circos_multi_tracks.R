#!/usr/bin/env Rscript
# Circos plot with multiple tracks: chromosome backbone, gene density,
# TE density, SNP density, INDEL density, and SV density.
# Inputs:
#   - chrom_length.xlsx: two columns [chrom, end]; 0-based start will be added
#   - hap1.genedensity.txt: 4 cols [Chr, Start, End, Density]
#   - TEdensity.txt:       4 cols [Chr, Start, End, Density]
#   - INS_density.bedGraph: 4 cols [Chr, Start, End, Density]  (rename if needed)
#   - SNP_density.bedGraph: 4 cols [Chr, Start, End, Density]
#   - DEL_density.bedGraph: 4 cols [Chr, Start, End, Density]

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(circlize)
  library(RColorBrewer)
})

# ---------- 1) Read chromosome lengths ----------
chr.len <- read_excel("chrom_length.xlsx", col_names = c("chrom","end")) |>
  mutate(start = 0)

# ---------- 2) Circos init ----------
circos.clear()
n_sec <- nrow(chr.len)
gaps <- rep(1.1, n_sec); if (n_sec >= 1) gaps[n_sec] <- 10

circos.par(
  start.degree = 85,
  track.height = 0.07,
  track.margin = c(0, 0),
  gap.degree = gaps
)

circos.initialize(
  factors = chr.len$chrom,
  xlim = cbind(chr.len$start, chr.len$end)
)

# ---------- 3) Chromosome backbone (labels + ticks) ----------
circos.track(
  factors = chr.len$chrom,
  ylim = c(0, 1),
  track.height = .05,
  bg.border = NA,
  bg.col = "#8AAECD",
  panel.fun = function(x, y) {
    circos.text(
      x = mean(CELL_META$xlim),
      y = mean(CELL_META$ylim) + 3,
      labels = CELL_META$sector.index,
      facing = "outside", niceFacing = TRUE, cex = 1, font = 2
    )
    circos.axis(
      h = "top",
      labels.cex = .6, labels.font = 2,
      major.tick = TRUE,
      major.at = seq(0, 1e8, by = 1e7),
      labels   = seq(0, 100, by = 10),
      major.tick.length = .3,
      minor.ticks = 1,
      labels.facing = "inside"
    )
  }
)

# Helper: harmonize chr names like "chr1" -> "chr01" (if needed)
pad_chr <- function(v) ifelse(nchar(v)==4, sub("chr","chr0",v), v)

# Generic drawer for a dense rectangle track
draw_track_rect <- function(df, palette_name = "Purples", log_transform = FALSE,
                            q_low = 0, q_high = 1, track_height = 0.085) {
  df <- df %>% mutate(Chr = pad_chr(Chr))
  z <- if (log_transform) log10(df$Density + 1e-5) else df$Density
  if (q_low > 0 || q_high < 1) {
    lo <- quantile(z, q_low, na.rm = TRUE); hi <- quantile(z, q_high, na.rm = TRUE)
  } else {
    lo <- min(z, na.rm = TRUE); hi <- max(z, na.rm = TRUE)
  }
  cols <- brewer.pal(8, palette_name)
  col_fun <- circlize::colorRamp2(seq(lo, hi, length.out = 8), cols)
  df$Color <- col_fun(z)

  circos.track(
    factors = chr.len$chrom,
    ylim = c(0, 10),
    track.height = track_height,
    bg.col = NA, bg.border = NA,
    panel.fun = function(x, y) {
      current_chr <- CELL_META$sector.index
      chr_data <- df[df$Chr == current_chr, ]
      if (nrow(chr_data) > 0) {
        circos.rect(
          xleft = chr_data$Start, xright = chr_data$End,
          ybottom = 0, ytop = 10,
          col = chr_data$Color, border = NA
        )
      }
    }
  )
}

# ---------- 4) Gene density track ----------
genedensity <- read.table("hap1.genedensity.txt", sep = "\t", header = FALSE,
                          col.names = c("Chr","Start","End","Density"))
draw_track_rect(genedensity, palette_name = "Purples", log_transform = FALSE)

# ---------- 5) TE density track (light gradient, log-scaled, 5–95% robust range) ----------
TEdensity <- read.table("TEdensity.txt", sep = "\t", header = FALSE,
                        col.names = c("Chr","Start","End","Density"))
draw_track_rect(TEdensity, palette_name = "YlOrBr",
                log_transform = TRUE, q_low = 0.05, q_high = 0.95)

# ---------- 6) SNP density ----------
snp <- read.table("SNP_density.bedGraph", sep = "\t", header = FALSE,
                  col.names = c("Chr","Start","End","Density"))
draw_track_rect(snp, palette_name = "BuGn")

# ---------- 7) INDEL density ----------
indel <- read.table("INS_density.bedGraph", sep = "\t", header = FALSE,
                    col.names = c("Chr","Start","End","Density"))
draw_track_rect(indel, palette_name = "Blues")

# ---------- 8) SV density (e.g., DEL) ----------
sv <- read.table("DEL_density.bedGraph", sep = "\t", header = FALSE,
                 col.names = c("Chr","Start","End","Density"))
draw_track_rect(sv, palette_name = "Oranges")

message("Circos figure drawn. Save from the graphics device if needed.")
