library(phyloseq)
library(devtools)
library(amplicon)
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
library(linkET)
library(dplyr)
library(vegan)
library(readr)
library(tibble)

env_common <- read.csv(
  "env_common.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mag_common <- read.csv(
  "mag_common.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

env_vars <- c("Depth", "tem", "ph", "SO42-", "NO3-", "Cl-", "Na", "K", "Mg", "Ca")

rownames(mag_common) <- mag_common[, 1]
mag_common <- mag_common[, -1, drop = FALSE]

mag_common <- as.data.frame(t(as.matrix(mag_common)), check.names = FALSE)
mag_common$Acession <- rownames(mag_common)
rownames(mag_common) <- NULL

mag_common <- mag_common %>%
  select(Acession, everything())

common_ids <- intersect(env_common$Acession, mag_common$Acession)

env_common <- env_common %>%
  filter(Acession %in% common_ids)

mag_common <- mag_common %>%
  filter(Acession %in% common_ids)

env_common <- env_common[match(common_ids, env_common$Acession), ]
mag_common <- mag_common[match(common_ids, mag_common$Acession), ]

stopifnot(all(env_common$Acession == mag_common$Acession))

run_mantel_single_all <- function(env_df, mag_df, var_name, min_n = 3, permutations = 999) {
  
  keep_ids <- env_df %>%
    filter(!is.na(.data[[var_name]])) %>%
    pull(Acession)
  
  if (length(keep_ids) < min_n) {
    return(data.frame(
      Variable = var_name,
      N = length(keep_ids),
      Mantel_r = NA,
      P_value = NA
    ))
  }
  
  env_use <- env_df %>%
    filter(Acession %in% keep_ids) %>%
    arrange(match(Acession, keep_ids))
  
  mag_use <- mag_df %>%
    filter(Acession %in% keep_ids) %>%
    arrange(match(Acession, keep_ids))
  
  mag_matrix <- mag_use[, !(colnames(mag_use) %in% "Acession"), drop = FALSE]
  mag_matrix <- as.matrix(mag_matrix)
  
  mag_matrix <- mag_matrix[, colSums(mag_matrix, na.rm = TRUE) > 0, drop = FALSE]
  
  if (ncol(mag_matrix) == 0) {
    return(data.frame(
      Variable = var_name,
      N = nrow(mag_use),
      Mantel_r = NA,
      P_value = NA
    ))
  }
  
  env_vector <- env_use[[var_name]]
  env_matrix <- matrix(env_vector, ncol = 1)
  
  comm_dist <- vegdist(mag_matrix, method = "bray")
  env_dist  <- dist(env_matrix, method = "euclidean")
  
  mt <- mantel(comm_dist, env_dist, permutations = permutations, method = "spearman")
  
  data.frame(
    Variable = var_name,
    N = nrow(env_use),
    Mantel_r = unname(mt$statistic),
    P_value = mt$signif
  )
}

result_list <- list()

for (v in env_vars) {
  tmp <- run_mantel_single_all(
    env_df = env_common,
    mag_df = mag_common,
    var_name = v,
    min_n = 3,
    permutations = 999
  )
  result_list[[v]] <- tmp
}

mantel_results_all <- do.call(rbind, result_list)

mantel_results_plot <- mantel_results_all %>%
  transmute(
    spec = "All samples", 
    env  = Variable,
    N    = N,
    r    = Mantel_r,
    p    = P_value,
    
    r_abs = abs(Mantel_r),
    
    rd = case_when(
      is.na(Mantel_r) ~ NA_character_,
      abs(Mantel_r) < 0.1 ~ "< 0.1",
      abs(Mantel_r) >= 0.1 & abs(Mantel_r) < 0.2 ~ "0.1 - 0.2",
      abs(Mantel_r) >= 0.2 & abs(Mantel_r) < 0.4 ~ "0.2 - 0.4",
      abs(Mantel_r) >= 0.4 ~ ">= 0.4"
    ),
    
    line = case_when(
      is.na(P_value) ~ NA_character_,
      P_value < 0.01 ~ "< 0.01",
      P_value >= 0.01 & P_value < 0.05 ~ "0.01 - 0.05",
      P_value >= 0.05 ~ ">= 0.05"
    ),
    
    col = case_when(
      is.na(Mantel_r) ~ NA_character_,
      Mantel_r > 0 ~ "Positive",
      Mantel_r < 0 ~ "Negative",
      Mantel_r == 0 ~ "Neutral"
    )
  )

write.csv(
  mantel_results_plot,
  "Mantel_results_all_samples_for_plot.csv",
  row.names = FALSE
)

data <- read.csv(
  "env_common.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(data)

rownames(data) <- data$Acession

env <- data[, c("Depth", "tem", "ph", "SO42-", "NO3-", "Cl-", "Na", "K", "Mg", "Ca"), drop = FALSE]
env[] <- lapply(env, function(x) as.numeric(as.character(x)))

env2 <- env[, sapply(env, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

linkdata <- read.csv(
  "Mantel_results_all_samples_for_plot.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(linkdata)

linkdata_plot <- linkdata %>%
  filter(!is.na(r), !is.na(p), !is.na(rd), !is.na(line), !is.na(col))

linkdata_plot <- linkdata_plot %>%
  filter(env %in% colnames(env2))

linkdata_plot$rd <- factor(
  linkdata_plot$rd,
  levels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4"),
  ordered = TRUE
)

linkdata_plot$line <- factor(
  linkdata_plot$line,
  levels = c("< 0.01", "0.01 - 0.05", ">= 0.05"),
  ordered = TRUE
)

linkdata_plot$col <- factor(
  linkdata_plot$col,
  levels = c("Positive", "Negative", "Neutral"),
  ordered = TRUE
)

linkdata_plot$spec <- factor(
  linkdata_plot$spec,
  levels = c("All samples"),
  ordered = TRUE
)

(p <- qcorrplot(
  correlate(
    env2,
    method = "spearman",
    use = "pairwise.complete.obs"
  ),
  type = "upper",
  diag = FALSE
) +
  geom_square(size = 0.25) +
  geom_couple(
    aes(
      colour = line,
      size = rd,
      from = spec,
      to = env
    ),
    data = linkdata_plot,
    curvature = 0.25
  ) +
  scale_size_manual(
    values = c(
      "< 0.1" = 0.15,
      "0.1 - 0.2" = 0.3,
      "0.2 - 0.4" = 0.6,
      ">= 0.4" = 1
    )
  ) +
  scale_colour_manual(
    values = c(
      "< 0.01" = "#c74546",
      "0.01 - 0.05" = "#93cc82",
      ">= 0.05" = "#DDDDDD"
    )
  ) +
  guides(
    size = guide_legend(
      title = "Mantel's r",
      override.aes = list(colour = "grey35"),
      order = 2
    ),
    colour = guide_legend(
      title = "Mantel's p",
      override.aes = list(size = 3),
      order = 1
    ),
    fill = guide_colorbar(
      title = "Spearman's r",
      order = 3
    )
  ) +
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7)
  ))

ggsave(
  "Mantel_qcorrplot.pdf",
  plot = p,
  width = 4,
  height = 4.5
)

# Note: Minor aesthetic adjustments (e.g., font size, label placement, and spacing) were performed in Adobe Illustrator to improve figure clarity.
