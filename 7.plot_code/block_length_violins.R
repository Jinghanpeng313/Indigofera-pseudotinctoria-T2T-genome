#!/usr/bin/env Rscript
# Violin + inner box for log10 block lengths per SV/syntenic type.

suppressPackageStartupMessages(library(tidyverse))

df <- read_tsv("last.txt", col_names = c("type","length_bp"), show_col_types = FALSE) |>
  filter(length_bp > 0) |>
  mutate(type = case_when(
    type %in% c("TRANS","TRANSAL") ~ "Translocation",
    type %in% c("DUP","DUPAL")     ~ "Duplication",
    type %in% c("INV","INVDP","INVTR","INVAL","INVDPAL","INVTRAL") ~ "Inversion",
    type == "SYN"   ~ "Syntenic region",
    type == "NOTAL" ~ "Not aligned",
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(type)) |>
  mutate(type = factor(type, levels = c("Translocation","Duplication","Inversion",
                                        "Syntenic region","Not aligned")))

pal <- c("#f1a16c","#e8d37d","#91b2c6","#69a893","#8da0cb"); names(pal) <- levels(df$type)

p <- ggplot(df, aes(x = log10(length_bp), y = type, fill = type)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.95, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = NA, color = "black", alpha = 0.6) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = expression(log[10](length~"(bp)")), y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11))
ggsave("block_length_violins.png", p, width = 7.5, height = 4.8, dpi = 300); print(p)
