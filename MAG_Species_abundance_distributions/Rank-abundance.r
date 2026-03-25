library(ggplot2)
library(BiodiversityR)
library(tidyverse)

otu_table <- read.delim(
  "all_MAG_covem_oil.txt",
  row.names = 1,
  check.names = FALSE
)

mapping <- read.csv("mapping_oil.csv", stringsAsFactors = FALSE)

head(colnames(otu_table))
head(mapping$Sample)

group_colors <- c(
  oil        = "#e67e22"
)

groups <- names(group_colors)

rank_all <- data.frame()

for (g in groups) {
  
  samples_g <- mapping %>%
    filter(group == g) %>%
    pull(Sample)
  
  samples_g <- intersect(samples_g, colnames(otu_table))
  if (length(samples_g) == 0) next
  
  otu_g <- otu_table[, samples_g, drop = FALSE]
  
  otu_rel <- sweep(otu_g, 2, colSums(otu_g), FUN = "/")
  
  mean_abundance <- rowMeans(otu_rel, na.rm = TRUE)
  
  mean_abundance <- mean_abundance[mean_abundance > 0]
  
  rank_df <- as.data.frame(
    rankabundance(
      as.data.frame(t(mean_abundance)),
      digits = 6
    )
  )
  
  rank_df <- rank_df[, c("rank", "abundance")]
  
  rank_df$group <- g
  
  rank_all <- bind_rows(rank_all, rank_df)
}

head(rank_all)

p_rank <- ggplot(
  rank_all,
  aes(x = rank,
      y = log10(abundance),
      color = group)
) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = group_colors,
    name = "Group"
  ) +
  scale_y_continuous(
    breaks = 0:-6,
    labels = c("100", "10", "1", "0.1", "0.01", "0.001", "0.0001")
  ) +
  labs(
    x = "MAG rank",
    y = "Relative abundance (%)"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 10)
  )

p_rank

ggsave(
  "rank_abundance_oil.pdf",
  p_rank,
  width = 5.5,
  height = 3
)
