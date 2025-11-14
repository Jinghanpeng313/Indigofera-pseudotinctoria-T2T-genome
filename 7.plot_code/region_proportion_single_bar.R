#!/usr/bin/env Rscript
# Single stacked bar showing region composition.

suppressPackageStartupMessages({
  library(tidyverse)
})

df <- tibble::tibble(
  region = factor(c("Down 2k","Exon","Intron","Up 2k"),
                  levels = c("Down 2k","Exon","Intron","Up 2k")),
  pct = c(0.2116, 0.1067, 0.3509, 0.3308)
)

df_pos <- df %>%
  mutate(ymax = cumsum(pct),
         ymin = lag(ymax, default = 0),
         ymid = (ymin + ymax)/2)

cols <- c("Down 2k"="#f1a16c","Exon"="#e8d37d","Intron"="#91b2c6","Up 2k"="#69a893")

p <- ggplot() +
  geom_col(data = df, aes(x = " ", y = pct, fill = region),
           width = 0.7, color = "black", linewidth = 0.5) +
  geom_text(data = df_pos, aes(x = " ", y = ymid, label = region),
            angle = 90, vjust = 0.5, size = 5) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25),
                     labels = scales::percent_format(accuracy = 1),
                     expand = c(0,0), position = "right") +
  scale_fill_manual(values = cols) +
  coord_cartesian(clip = "off") +
  guides(fill = "none") +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.y.right = element_text(size = 12, face = "bold"),
        plot.margin = margin(10, 30, 10, 10),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

ggsave("region_proportion_bar.png", p, width = 4, height = 6, dpi = 300)
print(p)
