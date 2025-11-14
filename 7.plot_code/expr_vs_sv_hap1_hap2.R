#!/usr/bin/env Rscript
# Combine hap1 and hap2 gene tables; compare expression per SV class, colored by haplotype.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

hap1 <- read.table("hap1_genes.txt", header = FALSE, sep = "\t",
                   col.names = c("Gene","Expression","SV_Type"), stringsAsFactors = FALSE) |>
  mutate(Subgenome = "Hap1")
hap2 <- read.table("hap2_genes.txt", header = FALSE, sep = "\t",
                   col.names = c("Gene","Expression","SV_Type"), stringsAsFactors = FALSE) |>
  mutate(Subgenome = "Hap2")
df <- bind_rows(hap1, hap2) |>
  mutate(logTPM = log10(Expression + 1),
         SV_Type = recode(SV_Type,
                          "SYN"="Synteny","INV"="Inversion","TRANS"="Translocation",
                          "DUP"="Duplication","NOTAL"="Not aligned", .default = SV_Type),
         SV_Type = factor(SV_Type, levels = c("Synteny","Inversion","Translocation","Duplication","Not aligned")),
         Subgenome = factor(Subgenome, levels = c("Hap1","Hap2")))

p <- ggplot(df, aes(x = SV_Type, y = logTPM, fill = Subgenome)) +
  geom_violin(position = position_dodge(0.9), alpha = 0.75, color = NA) +
  geom_boxplot(width = 0.18, position = position_dodge(0.9), alpha = 0.3) +
  scale_fill_manual(values = c("Hap1"="#f1a16c","Hap2"="#91b2c6")) +
  labs(y = "Gene expression (log10(TPM + 1))", x = "SV type", fill = "Haplotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))
# Optional global test per facet could be added if desired.

ggsave("expr_vs_sv_hap1_hap2.png", p, width = 8.6, height = 4.8, dpi = 300); print(p)
