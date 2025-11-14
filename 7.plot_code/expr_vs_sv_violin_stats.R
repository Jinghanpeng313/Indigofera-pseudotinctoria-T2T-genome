#!/usr/bin/env Rscript
# Gene expression (log10(TPM+1)) vs SV classes: violin + box + mean±SE + pairwise Wilcoxon.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

df <- read_tsv("hap1_genes.txt",
               col_names = c("gene_id","TPM","SV"),
               show_col_types = FALSE) |>
  mutate(logTPM = log10(TPM + 1),
         SV = recode(SV,
                     "SYN"="Synteny","INV"="Inversion","TRANS"="Translocation",
                     "DUP"="Duplication","NOTAL"="Not aligned",
                     .default = SV),
         SV = factor(SV, levels = c("Synteny","Inversion","Translocation","Duplication","Not aligned")))

groups <- levels(df$SV)
comparisons <- combn(groups, 2, simplify = FALSE)
y_max <- max(df$logTPM, na.rm = TRUE)

p <- ggplot(df, aes(x = SV, y = logTPM, fill = SV)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", width = 0.95) +
  geom_boxplot(width = 0.18, outlier.shape = 21, outlier.fill = "white",
               outlier.size = 1.8, fatten = 1, alpha = 0.9, color = "black") +
  stat_summary(fun = mean,
               fun.min = \(z) mean(z) - sd(z)/sqrt(length(z)),
               fun.max = \(z) mean(z) + sd(z)/sqrt(length(z)),
               geom = "pointrange", shape = 22, size = 0.35, fill = "white",
               color = "black", position = position_dodge(width = 0.9)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE,
                     tip.length = 0.01, size = 4, vjust = 0.6) +
  stat_compare_means(method = "kruskal.test", label.y = y_max * 1.12, label.x = 3, size = 4.2) +
  labs(x = NULL, y = expression(paste("Gene expression (", log[10], "(TPM + 1))"))) +
  coord_cartesian(ylim = c(0, y_max * 1.25)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, margin = margin(r = 8)))
ggsave("expr_vs_sv_violin_stats.png", p, width = 8.6, height = 4.8, dpi = 300); print(p)
