rm(list = ls())

library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# pie charts 

process <- read.csv("mag_tpm_ProcessImportance_EachGroup.csv")

process <- process %>%
  select(-c("Method", "GroupBasedOn"))

process <- process[-c(5:10),]

process_long <- process %>%
  pivot_longer(
    cols = -Group,
    names_to = "Process",
    values_to = "Prop"
  )

pro <- process_long %>%
  mutate(
    Process = factor(
      Process,
      levels = c("HeS", "HoS", "DL", "HD", "DR")
    )
  )

Process_colors <- c(
  "HeS" = "#e51a1b",
  "HoS" = "#377eb8",
  "DL" = "#984ea2",
  "HD" = "#4caf4a",
  "DR" = "#ff7f00"
)

plot <- ggplot(pro, aes(x = "", y = Prop * 100, fill = Process)) +
  geom_bar(
    width = 1,
    stat = "identity",
    color = "white",
    linewidth = 0.5
  ) +
  coord_polar("y", start = 0) +
  
  geom_text(
    aes(label = paste0(round(Prop * 100, 2), "%")),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 3,
    fontface = "bold"
  ) +
  
  facet_wrap(~ Group, nrow = 1) +
  
  scale_fill_manual(values = Process_colors) +
  theme_void() +
  theme(
    strip.text = element_text(
      size = 12,
      face = "bold",
      margin = margin(0, 0, 10, 0)
    ),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.spacing = unit(1, "cm")
  ) +
  
  labs(fill = "Process") +
  guides(fill = guide_legend(override.aes = list(color = NA)))

plot

ggsave("bar_plot.pdf", plot, width = 12, height = 6)

# stacked bar plot

OTUBinClass <- read.table("mag_tpm_Taxon_Bin.csv", header = TRUE, sep = ",")

tax_Abu <- OTUBinClass

# Keep top 10 phyla by relative abundance and collapse the rest as Others
phylum_abundance <- tax_Abu %>%
  group_by(Phylum) %>%
  summarise(sum = sum(TaxonRelativeAbundance), .groups = "drop")

top10_phyla <- phylum_abundance %>%
  arrange(desc(sum)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

tax_Abu_modified <- tax_Abu %>%
  mutate(
    Phylum = if_else(Phylum %in% top10_phyla, Phylum, "Others")
  ) %>%
  mutate(
    Phylum = factor(Phylum, levels = c(top10_phyla, "Others"))
  )

ecology <- read.csv("EachBin_process.csv")

oil_field_processed <- ecology %>%
  filter(Group == "oil_field") %>%
  select(-Group) %>%
  pivot_longer(
    cols = starts_with("bin"),
    names_to = "Bin",
    values_to = "Value"
  ) %>%
  mutate(
    Value = ifelse(is.na(Value), 0, Value)
  )

AbuGene1 <- tax_Abu_modified %>%
  select(ID, Bin)

oil_field_data <- AbuGene1 %>%
  left_join(oil_field_processed, by = "Bin")

write.csv(oil_field_data, "oil_field_data.csv", row.names = FALSE)

anno <- read.csv("anonotation.csv")
anno2 <- anno %>%
  rename(ID = MAG)

run_icamp_plot <- function(df_data, tag = "oilfield") {
  
  df <- df_data %>%
    left_join(anno2, by = "ID")
  
  df_sum <- df %>%
    group_by(Phylum, Index) %>%
    summarise(SumValue = sum(Value, na.rm = TRUE), .groups = "drop")
  
  df_pct <- df_sum %>%
    mutate(
      Percent = SumValue / sum(SumValue, na.rm = TRUE) * 100
    )
  
  top10_phylum <- df_pct %>%
    group_by(Phylum) %>%
    summarise(TotalPercent = sum(Percent), .groups = "drop") %>%
    arrange(desc(TotalPercent)) %>%
    slice(1:10) %>%
    pull(Phylum)
  
  df_pct2 <- df_pct %>%
    mutate(
      Phylum = ifelse(Phylum %in% top10_phylum, Phylum, "Other")
    ) %>%
    group_by(Phylum, Index) %>%
    summarise(
      Percent = sum(Percent),
      .groups = "drop"
    ) %>%
    group_by(Phylum) %>%
    mutate(Total = sum(Percent)) %>%
    ungroup()
  
  write.csv(
    df_pct2,
    paste0("icamp_", tag, "_phylum.csv"),
    row.names = FALSE
  )
  
  p <- ggplot(df_pct2, aes(x = reorder(Phylum, Total), y = Percent, fill = Index)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = Process_colors) +
    coord_flip() +
    theme_bw() +
    labs(
      x = "Phylum",
      y = "Contribution (%)",
      fill = "Ecological Process"
    )
  
  ggsave(
    paste0("icamp_", tag, "_phylum.pdf"),
    p,
    width = 4,
    height = 3
  )
  
  return(list(data = df_pct2, plot = p))
}

# Run for oilfield only
oil_result <- run_icamp_plot(oil_field_data, "oilfield")

oil_result$plot