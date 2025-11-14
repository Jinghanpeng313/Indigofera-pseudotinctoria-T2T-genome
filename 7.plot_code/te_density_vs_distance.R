#!/usr/bin/env Rscript
# Boxplots of TE density across distance bins to inversions, with linear fit on raw points
# and on group means.

suppressPackageStartupMessages(library(tidyverse))

dat <- read.table("box_hap2.txt", header = TRUE, sep = "\t")
dat$Site <- factor(dat$Site, levels = c("center","10K","100K","500K","1M","5M"))

# Regression on all points
model <- lm(Te_density ~ Site, data = dat)

p1 <- ggplot(dat, aes(Site, Te_density)) +
  geom_boxplot(aes(fill = Site), outlier.shape = 16, outlier.color = "black", outlier.size = 2) +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white") +
  geom_smooth(method = "lm", aes(group = 1), color = "blue", se = FALSE) +
  theme_minimal() +
  labs(title = "TE density vs. distance to inversion",
       subtitle = sprintf("All points: R² = %.2f; p = %.3g",
                          summary(model)$r.squared,
                          coef(summary(model))[2,4]),
       x = NULL, y = "TE density")

# Regression on group means
data_avg <- dat |>
  group_by(Site) |>
  summarise(Avg_Te_density = mean(Te_density, na.rm = TRUE), .groups = "drop")
model_avg <- lm(Avg_Te_density ~ Site, data = data_avg)

p2 <- ggplot(dat, aes(Site, Te_density)) +
  geom_boxplot(aes(fill = Site), outlier.shape = 16, outlier.color = "black", outlier.size = 2) +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white") +
  geom_smooth(data = data_avg, aes(x = Site, y = Avg_Te_density, group = 1),
              method = "lm", color = "blue", se = FALSE) +
  theme_minimal() +
  labs(title = "TE density vs. distance (fit on group means)",
       subtitle = sprintf("Group means: R² = %.2f; p = %.3g",
                          summary(model_avg)$r.squared,
                          coef(summary(model))[2,4]),
       x = NULL, y = "TE density")

ggsave("te_density_vs_distance_points.png", p1, width = 7.2, height = 4.6, dpi = 300)
ggsave("te_density_vs_distance_means.png",  p2, width = 7.2, height = 4.6, dpi = 300)
print(p1); print(p2)
