library(tidyverse)
library(patchwork)
library(dplyr)
library(tidyr)
library(forcats)

# taxonomic novelty

ar_tax <- read.csv("ar_taxa.txt", sep = "\t", header = TRUE) 
bac_tax <- read.csv("bac_taxa.txt", sep = "\t", header = TRUE)
tax <- rbind(ar_tax, bac_tax)

df_ratio <- tax %>%
  mutate(across(Phylum:Species, ~ na_if(., ""))) %>%
  pivot_longer(
    cols = Phylum:Species,
    names_to = "level",
    values_to = "name"
  ) %>%
  group_by(Domain, level) %>%
  summarise(
    known = sum(!is.na(name)),
    novel = sum(is.na(name)),
    total = known + novel,
    novel_pct = novel / total * 100,
    .groups = "drop"
  )

df_ratio <- df_ratio %>%
  mutate(
    level = factor(
      level,
      levels = c("Class", "Order", "Family", "Genus", "Species")
    )
  )

df_ratio <- df_ratio %>%
  filter(level %in% c("Order", "Family", "Genus", "Species"))

df_plot <- df_ratio %>%
  select(Domain, level, known, novel, total) %>%
  pivot_longer(
    cols = c(known, novel),
    names_to = "type",
    values_to = "count"
  ) %>%
  group_by(Domain, level) %>%
  mutate(
    pct = count / sum(count) * 100
  ) %>%
  ungroup()

df_plot <- df_plot %>%
  mutate(
    level = factor(level, levels = c("Order", "Family", "Genus", "Species"))
  )

plot_novel_ratio <- function(df, domain_name, tag = NULL) {
  
  df_sub <- df %>% filter(Domain == domain_name)
  
  ggplot(df_sub, aes(x = pct, y = level, fill = type)) +
    geom_col(
      width = 0.6,
      color = "black",
      size = 0.3
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = c(
        novel = "grey60",
        known = "white"
      ),
      labels = c(
        novel = "Novel taxa",
        known = "Known taxa"
      )
    ) +
    labs(
      x = "Percentage of MAGs (%)",
      y = NULL,
      fill = "type"
    ) +
    theme_classic() +
    scale_y_discrete(limits = c("Order", "Family", "Genus", "Species")) +
    annotate(
      "text",
      x = 0,
      y = Inf,
      label = tag,
      hjust = -0.5,
      vjust = 1.2,
      fontface = "bold",
      size = 5
    ) +
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      plot.margin = unit(c(1, 1, 1, 1), "lines"),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.background = element_blank(),
      axis.text.y.right = element_text(color = "#666666"),
      axis.title.y.right = element_text(color = "#666666")
    )
}

p_bac <- plot_novel_ratio(df_plot, "Bacteria", tag = "b")
p_arc <- plot_novel_ratio(df_plot, "Archaea", tag = "c")

(p <-p_arc + p_bac)

ggsave("MAG_novelty.pdf", p, width = 10, height = 4)

# Distribution of novel and known MAGs

ar_tax <- read.csv("ar_taxa.txt", sep = "\t", header = TRUE) 
bac_tax <- read.csv("bac_taxa.txt", sep = "\t", header = TRUE)
tax<-rbind(ar_tax,bac_tax)
df <-tax[,1:2]
df <- tax %>%
  filter(Phylum != "Unassigned") %>%
  mutate(
    Phylum = sub("Methanobacteriota_[A-B]", "Methanobacteriota", Phylum),
    Phylum = sub("Bacillota_[A-I]", "Bacillota", Phylum),
    Phylum = sub("Myxococcota_A", "Myxococcota", Phylum),
    Phylum = sub("Desulfobacterota_[A-I]", "Desulfobacterota", Phylum),
    Species = sub("^$", NA, Species)
  ) %>%
  group_by(Domain, Phylum) %>%
  summarise(
    total = n(),
    Novel = sum(is.na(Species))
  ) %>%
  ungroup()
data <- df %>%
  mutate(label = Phylum) %>%
  mutate(label = fct_reorder(fct_lump_n(label, 10, w=total), total, .desc = TRUE)) %>%
  group_by(label) %>% 
  summarise(
    total = sum(total), 
    Novel = sum(Novel),
    fracNovel = Novel / total,
    Known = total - Novel
  ) %>%
  pivot_longer(c(Novel, Known)) %>%
  mutate(name = fct_relevel(name, c("Novel", "Known")))
levels(data$label)
pcolor <- c("Pseudomonadota"='#D89000',
            "Bacillota"='#A3A500',
            "Actinomycetota"="#F8766D",
            "Bacteroidota"='#00BFC4',
            "Patescibacteria"='#9590FF',
            "Desulfobacterota"='#39B600',
            "Chloroflexota"='#FF62BC',
            "Halobacteriota"="#8DA0CB",
            "Planctomycetota"='#9590FF',
            "Campylobacterota"="#00B0F6",
            "Unassigned"="#DBAA77",
            "Others"="#696969") 

data_plot <- df %>%
  mutate(label = fct_other(Phylum, keep = names(pcolor)) %>% 
           fct_relevel(names(pcolor))) %>%
  group_by(label) %>%
  summarise(total=sum(total), 
            Novel = sum(Novel), 
            fracNovel = Novel / total,
            Known = total - Novel) %>%
  pivot_longer(cols = c(Novel, Known), names_to="name", values_to="value") %>%
  mutate(name = fct_relevel(name, c("Novel", "Known")))

(p<-ggplot(data_plot, aes(x=label, y=value, fill=label, alpha=name)) +
  geom_col(col="black") +
  geom_line(aes(y=fracNovel*max(total), group=1), col="#666666") +
  geom_point(aes(y=fracNovel*max(total)), col="#666666") +
  scale_y_continuous(name = "ASVs",
                     sec.axis = sec_axis(~./max(data_plot$total), name = "Novel ASV %")) +
  labs(x="", alpha="") +
  scale_fill_manual(values = pcolor) +
  scale_alpha_manual(values = c(1, 0.3), labels=c("Novel species", "Known species")) +
  guides(fill="none") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=30, hjust=1, vjust=1),
    plot.margin = unit(c(1,1,1,1), "lines"),
    panel.grid = element_blank(),
    legend.position = c(0.7, 0.45),
    legend.background = element_blank(),
    axis.text.y.right = element_text(color="#666666"),
    axis.title.y.right = element_text(color="#666666")
  ))

ggsave("MAG_Novel_per_Abundphylum.pdf", width = 5, height = 4)

# Note: Minor aesthetic adjustments (e.g., font size, label placement, and spacing) were performed in Adobe Illustrator to improve figure clarity.
