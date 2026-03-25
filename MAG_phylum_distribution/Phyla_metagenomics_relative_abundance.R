library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
library(patchwork)

deep_samples_data <- read.csv("oil_plot_data.csv", stringsAsFactors = FALSE)

dominant_df <- deep_samples_data %>%
  filter(Phylum != "Others") %>%
  filter(!is.na(value)) %>%
  group_by(sample) %>%
  slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sample, dominant_phylum = Phylum, dominant_abundance = value)

phylum_order <- dominant_df %>%
  count(dominant_phylum, name = "n_samples") %>%
  arrange(desc(n_samples)) %>%
  pull(dominant_phylum)

sample_order <- dominant_df %>%
  mutate(dominant_phylum = factor(dominant_phylum, levels = phylum_order)) %>%
  arrange(dominant_phylum, desc(dominant_abundance)) %>%
  pull(sample) %>%
  unique()

deep_samples_data$sample <- factor(deep_samples_data$sample, levels = sample_order)

deep_samples_data$Phylum <- factor(
  deep_samples_data$Phylum,
  levels = c("Pseudomonadota","Halobacteriota","Bacillota","Patescibacteria",
             "Actinomycetota","Campylobacterota","Bacteroidota",
             "Desulfobacterota","Omnitrophota","Methanobacteriota","Others")
)

phylum_colors <- c(
  "Pseudomonadota"    = "#D89000",
  "Halobacteriota"    = "#8DA0CB",
  "Bacillota"         = "#A3A500",
  "Patescibacteria"   = "#9590FF",
  "Actinomycetota"    = "#F8766D",
  "Campylobacterota"  = "#00B0F6",
  "Bacteroidota"      = "#00BFC4",
  "Desulfobacterota"  = "#39B600",
  "Omnitrophota"      = "#B8F1ED",
  "Methanobacteriota" = "#E78AC3",
  "Others"            = "#555555"
)

plot_one_group <- function(df, title, phylum_colors) {
  
  phylum_order <- df %>%
    group_by(Phylum) %>%
    summarise(mean_abundance = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_abundance) %>%
    pull(Phylum)
  
  df$Phylum <- factor(df$Phylum, levels = phylum_order)
  
  ggplot(df, aes(x = sample, y = value, fill = Phylum)) +
    geom_col(width = 1, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = phylum_colors) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 10),
      plot.title   = element_blank(),
      legend.title = element_blank()
    )
}

group_list <- list(
  Pseudomonadota   = c("Pseudomonadota"),
  Patescibacteria  = c("Patescibacteria"),
  Halobacteriota   = c("Halobacteriota"),
  Bacillota        = c("Bacillota"),
  Actinomycetota   = c("Actinomycetota"),
  Omnitrophota     = c("Omnitrophota"),
  Desulfobacterota = c("Desulfobacterota"),
  Methanobacteriota = c("Methanobacteriota"),
  Campylobacterota = c("Campylobacterota")
)

plot_list  <- list()
width_list <- c()

for (g in names(group_list)) {
  
  doms <- group_list[[g]]
  
  df_g <- deep_samples_data %>%
    filter(dominant_phylum %in% doms)
  
  n_sample <- length(unique(df_g$sample))
  width_list <- c(width_list, n_sample)
  
  plot_list[[g]] <- plot_one_group(
    df = df_g,
    title = g,
    phylum_colors = phylum_colors
  )
}

final_plot_oil <- wrap_plots(
  plot_list,
  nrow = 1,
  widths = width_list,
  guides = "collect"
) &
  theme(legend.position = "right")

final_plot_oil

ggsave(
  "oil_distribution.pdf",
  plot = final_plot_oil,
  width = 8,
  height = 4
)

# Note: Minor aesthetic adjustments (e.g., font size, label placement, and spacing) were performed in Adobe Illustrator to improve figure clarity.
