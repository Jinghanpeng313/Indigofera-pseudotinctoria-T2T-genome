#!/usr/bin/env Rscript
# Same as above but styled with ggsci::scale_fill_npg.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ggsci)
})

df <- read_tsv("hap1_genes.txt",
               col_names = c("gene_id","TPM","SV"),
               show_col_types = FALSE) |>
  mutate(logTPM = log10(TPM + 1),
         SV = recode(SV,
                     "SYN"="Synteny","INV"="Inversion","TRANS"="Translocation",
                     "DUP"="Duplication","NOTAL"="Not aligned", .default = SV),
         SV = factor(SV, levels = c("Synteny","Inversion","Translocation","Duplication","Not aligned")))

comparisons <- combn(levels(df$SV), 2, simplify = FALSE)
y_max <- max(df$logTPM, na.rm = TRUE)

p <- ggplot(df, aes(SV, logTPM, fill = SV)) +
  geom_violin(trim = FALSE, alpha = 0.85, width = 0.95, color = "gray20", linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.shape = 21, outlier.fill = "white",
               outlier.size = 1.6, fatten = 1, alpha = 0.95,
               color = "gray10", linewidth = 0.45) +
  stat_summary(fun = mean,
               fun.min = \(z) mean(z)-sd(z)/sqrt(length(z)),
               fun.max = \(z) mean(z)+sd(z)/sqrt(length(z)),
               geom = "pointrange", shape = 22, size = 0.35, fill = "white", color = "gray10") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                     label = "p.signif", hide.ns = FALSE, tip.length = 0.01,
                     size = 4.2, step.increase = 0.09) +
  stat_compare_means(method = "kruskal.test", label.y = y_max * 1.18, size = 4.2) +
  scale_fill_npg() +
  labs(x = NULL, y = expression(paste("Gene expression (", log[10], "(TPM + 1))"))) +
  coord_cartesian(ylim = c(0, y_max * 1.25), clip = "off") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold", margin = margin(t = 6)),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)))
ggsave("expr_vs_sv_violin_npg.png", p, width = 8.6, height = 4.8, dpi = 300); print(p)
