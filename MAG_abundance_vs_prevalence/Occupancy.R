library(ggplot2)
library(ape)
library(vegan)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggrepel)
library(dplyr)

# Data

otu_table <- read.delim("all_MAG_covem_oil.tsv", row.names = 1,
                        check.names = FALSE)
dim(otu_table)

taxa <- read.csv("anonotation.csv")

# Calculation 

occupancy <- rowSums(otu_table > 0) / ncol(otu_table) * 100

otu_relative_abundance <- sweep(otu_table, 2, colSums(otu_table), FUN = "/") * 100

otu_mean_relative_abundance <- rowMeans(otu_relative_abundance, na.rm = TRUE)

otu_summary <- data.frame(genome = rownames(otu_table),
                          Occupancy = occupancy,
                          Mean_Relative_Abundance = otu_mean_relative_abundance)

write.csv(otu_summary, "MAG_summary_oil.csv", row.names = FALSE)

otu_taxa_summary <- left_join(
  otu_summary,
  taxa,
  by = c("genome" = "MAG")
)

write.csv(otu_taxa_summary, "otu_taxa_summary_oil.csv", row.names = FALSE)

names(otu_taxa_summary)

otu_taxa_summary <- otu_taxa_summary %>%
  mutate(Phylum_Color = ifelse(Occupancy > 40, 
                               Phylum, "Other"))

top10_phyla <- c(
  "Pseudomonadota", "Bacillota", "Patescibacteria",
  "Bacteroidota", "Desulfobacterota",
  "Omnitrophota", "Methanobacteriota","Halobacteriota","Actinomycetota","Campylobacterota"
)

top10_colors <- c(
  '#D89000','#A3A500','#9590FF','#00BFC4',
  '#39B600','#B8F1ED','#E78AC3',"#8DA0CB","#F8766D","#00B0F6"
)

otu_taxa_summary <- otu_taxa_summary %>%
  mutate(
    Phylum_Color = ifelse(
      Phylum %in% top10_phyla,
      Phylum,
      "Other"
    )
  )

otu_taxa_summary$Phylum_Color <- factor(
  otu_taxa_summary$Phylum_Color,
  levels = c(top10_phyla, "Other")
)

color_map <- c(
  setNames(top10_colors, top10_phyla),
  Other = "grey80"
)

(p_oil<- ggplot(
  otu_taxa_summary,
  aes(x = Mean_Relative_Abundance, y = Occupancy)
) +
    geom_point(
      data = subset(otu_taxa_summary, Occupancy <= 40),
      aes(color = Phylum_Color),
      size = 1.5,
      alpha = 0.6
    ) +
    
    geom_point(
      data = subset(otu_taxa_summary, Occupancy > 40),
      aes(color = Phylum_Color),
      size = 2.5
    ) +
    
    geom_hline(
      yintercept = 40,
      linetype = "dashed",
      color = "grey50"
    ) +
    
    geom_text_repel(
      data = subset(otu_taxa_summary, Occupancy > 40),
      aes(label = Genus, color = Phylum_Color),
      size = 2.8,
      box.padding = 0.3,
      point.padding = 0.4,
      max.overlaps = 30,
      segment.color = "grey70"
    ) +
    
    scale_color_manual(
      values = color_map,
      name = "Phylum (Top10 by abundance)"
    ) +
    
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    
    labs(
      x = "Mean relative abundance (%)",
      y = "Occupancy (%)"
    ) +
    
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text  = element_text(size = 11),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10)
    ))

ggsave(
  "Occupancy.pdf",
  plot = p_oil,
  width = 4,
  height = 4.5
)
