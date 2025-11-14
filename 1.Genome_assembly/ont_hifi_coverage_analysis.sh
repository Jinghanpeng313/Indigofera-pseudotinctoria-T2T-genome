#!/bin/bash
# -------------------------------------------------------------------
# 1. Align ONT Data to the Genome using Minimap2
# -------------------------------------------------------------------
minimap2 -t 8 -ax map-ont maji.genome_gap_filled.fa pass.ul.fq.gz > mj_ont.sam

# -------------------------------------------------------------------
# 2. Convert SAM to BAM format using Samtools (ONT Data)
# -------------------------------------------------------------------
samtools view -@ 8 -Sb mj_ont.sam > mj_ont.bam

# -------------------------------------------------------------------
# 3. Sort BAM file (ONT Data)
# -------------------------------------------------------------------
samtools sort -@ 8 mj_ont.bam -o mj_ont_sorted.bam

# -------------------------------------------------------------------
# 4. Calculate Coverage using Bedtools (ONT Data)
# -------------------------------------------------------------------
bedtools genomecov -ibam mj_ont_sorted.bam -bg > mj_ont_cov2.txt

# Extract the first four columns to simplify the coverage data for HiFi
cut -f1-4 window_cov.txt > hifi_cov.txt

# -------------------------------------------------------------------
# 5. Align HiFi Reads to the hap1 Genome using Minimap2
# -------------------------------------------------------------------
minimap2 -t 8 -ax map-ont hap1.fa pass.ul.fq.gz > hifi.sam

# -------------------------------------------------------------------
# 6. Convert SAM to BAM format using Samtools (HiFi Data)
# -------------------------------------------------------------------
samtools view -@ 8 -Sb hifi.sam > hifi.bam

# -------------------------------------------------------------------
# 7. Sort BAM file (HiFi Data)
# -------------------------------------------------------------------
samtools sort -@ 8 hifi.bam -o hifi_sorted.bam

# -------------------------------------------------------------------
# 8. Calculate Coverage using Bedtools (HiFi Data)
# -------------------------------------------------------------------
bedtools genomecov -ibam hifi_sorted.bam -bg > hifi_cov2.txt

# -------------------------------------------------------------------
# 9. R Plotting for Visualization (karyoploteR)
# -------------------------------------------------------------------

Rscript -e "
# Load necessary library
library(karyoploteR)

# 9.1 Read and prepare karyotype data for hap1 genome
karyotype <- read.table('A.chr.len.txt', sep = '\t', header = TRUE)
kar <- toGRanges(karyotype)

# Read HiFi coverage data and format it
hifi <- read.table('hifi_cov.txt', sep = '\t', header = FALSE)
colnames(hifi) <- c('chr', 'start', 'end', 'cov')
hifi_cov <- toGRanges(hifi[, 1:4])

# Create karyotype plot for hap1 genome
kp_1 <- plotKaryotype(genome = kar, plot.type = 2, chromosomes = 'all', ideogram.plotter = NULL)

# Limit HiFi coverage to a max of 300
hifi$cov <- pmin(hifi$cov, 300)

# Add HiFi coverage to the plot
kpPlotRegions(kp_1, data = kar, data.panel = 'ideogram', col = '#eab883')  # Background color
kpBars(kp_1, hifi_cov, x0 = hifi$start, x1 = hifi$end, y1 = hifi$cov, data.panel = 1, r0 = 0.1, r1 = 0.55, ymin = 0, ymax = 300, border = NA, col = c('#7895c1'))

# 9.2 Read and prepare ONT coverage data
ont_cov <- read.table('ont_cov.txt', sep = '\t', header = FALSE)
colnames(ont_cov) <- c('chr', 'start', 'end', 'cov')
ont_cov_gr <- toGRanges(ont_cov[, 1:4])

# Limit ONT coverage to a max of 300
ont_cov$cov <- pmin(ont_cov$cov, 300)

# Add ONT coverage to the plot
kpBars(kp_1, ont_cov_gr, x0 = ont_cov$start, x1 = ont_cov$end, y1 = ont_cov$cov, data.panel = 2, r0 = 0, r1 = 0.45, ymin = 0, ymax = 300, border = NA, col = c('#84ba42'))

# 9.3 Add additional genomic features like telomeres, centromeres, and gaps
tels <- read.table('mj.tel.txt', sep = '\t', header = TRUE)
tel40 <- toGRanges(tels)
kpPlotRegions(kp_1, data = tel40, data.panel = 'ideogram', col = 'black')

cen <- read.table('mj.cen.txt', sep = '\t', header = TRUE)
cen11 <- toGRanges(cen)
kpPlotRegions(kp_1, data = cen11, data.panel = 'ideogram', col = '#26445e')

gap <- read.table('mj.gap.txt', sep = '\t', header = TRUE)
gap11 <- toGRanges(gap)
kpPlotRegions(kp_1, data = gap11, data.panel = 'ideogram', col = '#d4562e')
"
