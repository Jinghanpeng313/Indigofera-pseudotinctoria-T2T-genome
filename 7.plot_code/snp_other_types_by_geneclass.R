#!/usr/bin/env Rscript
# Bar chart for other variants_* (non-impact) types.

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

infile <- "hap1_SNP_ASE.txt"
dat <- read.delim(infile, check.names = FALSE) %>% clean_names()
names(dat)[1] <- "GeneType"

all_variant_cols <- grep("^variants_", names(dat), value = TRUE)
impact_cols <- grep("^variants_impact_(HIGH|LOW|MODERATE|MODIFIER)$",
                    names(dat), value = TRUE)
other_cols <- setdiff(all_variant_cols, impact_cols)
if (length(other_cols) == 0) {
  stop("No non-impact variants_* columns found.")
}

sum_by_genetype <- function(df, keep_cols) {
  df %>%
    select(all_of(c("GeneType", keep_cols))) %>%
    pivot_longer(-GeneType, names_to = "Feature", values_to = "Value") %>%
    group_by(GeneType, Feature) %>%
    summarise(Total = sum(as.numeric(Value), na.rm = TRUE), .groups = "drop")
}

other_long <- sum_by_genetype(dat, other_cols) %>%
  mutate(OtherType = sub("^variants_", "", Feature)) %>%
  select(GeneType, OtherType, Total) %>%
  mutate(Total_log = log10(Total + 1))

ord <- other_long %>%
  group_by(OtherType) %>%
  summarise(m = mean(Total_log), .groups = "drop") %>%
  arrange(desc(m)) %>%
  pull(OtherType)
other_long$OtherType <- factor(other_long$OtherType, levels = ord)

gene_colors <- c("NoDiff"="#91b2c6","SubNeo"="#f1a16c","Hap"="#69a893")

p <- ggplot(other_long,
            aes(x = OtherType, y = Total_log, fill = GeneType)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = gene_colors) +
  labs(title = "Other SNP types (non-impact) per gene class",
       x = "SNP type", y = "log10(Total + 1)", fill = "Gene class") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("snp_other_types_by_geneclass.png", p, width = 9, height = 5, dpi = 300)
print(p)
