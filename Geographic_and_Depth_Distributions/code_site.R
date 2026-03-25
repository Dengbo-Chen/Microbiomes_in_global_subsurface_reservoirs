library(ggplot2)
library(maps)
library(dplyr)
library(patchwork)
library(raster)
library(geodata)
library(readr)
library(ggthemes)
library(RColorBrewer)
library(ggridges)
library(gghalves)
library(tidyr)

#MAP

data1 <- read.csv("sample_site.csv")
names(data1)
Worldmap <- map_data("world")

alt <- worldclim_global(var = "elev", res = 2.5, path = "./wc")
alt <- raster(alt) 

ext_pts <- extent(
  min(data1$Longitude) - 10, max(data1$Longitude) + 10,
  min(data1$Latitude)  - 10, max(data1$Latitude)  + 10
)
alt_crop <- crop(alt, ext_pts)


alt_df <- as.data.frame(alt, xy = TRUE, na.rm = TRUE)
colnames(alt_df) <- c("x", "y", "elev_m")


data1$elev_m <- raster::extract(alt, data1[, c("Longitude", "Latitude")])


(all_sample_map<- ggplot() +
  geom_raster(
    data = alt_df,
    aes(x = x, y = y, fill = elev_m)
  ) +
  scale_fill_gradient(
    name = "Elevation (m)",
    low  = "#edf8e9",
    high = "#8c510a"
  ) +
  geom_point(
    data = data1,
    aes(x = Longitude, y = Latitude, color = group),
    size = 0.5
  ) +
  scale_color_manual(values = c(
    "Deep aquifer water"        = "#BC4749",
    "mine fluids"              = "#247BA0",
    "oil field produced water" = "#e67e22",
    "Shale-produced fluids"    = "#993B90"
  )) +
  scale_x_continuous(
    name = "Longitude(°)",
    breaks = seq(-180, 180, by = 60)
  ) +
  scale_y_continuous(
    name = "Latitude(°)",
    breaks = seq(-90, 90, by = 30)
  ) +
  coord_fixed(
    ratio = 1,
    xlim = c(-180, 180),
    ylim = c(-90, 90),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  ))

ggsave("world_map_group.pdf",all_sample_map, width = 6, height = 4)

#Depth_Distributions

depth_data  <- read.csv("depth_data.csv", check.names = FALSE)
mapping_cdb <- read.csv("mapping.csv", check.names = FALSE)

data_plot <- depth_data %>%
  left_join(mapping_cdb %>% dplyr::select(Sample, group), by = "Sample") %>%
  filter(!is.na(group)) %>%
  rename(depth = `depth`)

data_plot <- data_plot %>%
  mutate(
    group = factor(group, levels = c(
      "deep_water","mine","oil","shale" 
    ))
  )

box_stats <- data_plot %>%
  group_by(group) %>%
  summarise(
    Q1 = quantile(depth, 0.25, na.rm = TRUE),
    Q3 = quantile(depth, 0.75, na.rm = TRUE),
    median = median(depth, na.rm = TRUE),
    IQR = IQR(depth, na.rm = TRUE),
    lower_whisker = max(min(depth, na.rm = TRUE), Q1 - 1.5 * IQR),
    upper_whisker = min(max(depth, na.rm = TRUE), Q3 + 1.5 * IQR),
    outliers = list(depth[depth < lower_whisker | depth > upper_whisker]),
    .groups = "drop"
  ) %>%
  mutate(
    width = 0.25,
    offset = 0.15   
  )

(fig_depth_ridge <- ggplot(data_plot, aes(x = depth, y = group, fill = group)) +
    geom_density_ridges(
      scale = 0.75,
      rel_min_height = 0.01,
      alpha = 0.8,
      color = "black",
      size = 0.2,
      panel_scaling = FALSE
    ) +
    geom_rect(
      data = box_stats,
      aes(
        xmin = Q1, xmax = Q3,
        ymin = as.numeric(group) - width / 2 - offset,
        ymax = as.numeric(group) + width / 2 - offset,
        fill = group
      ),
      color = "black",
      size = 0.2,
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = box_stats,
      aes(
        x = median, xend = median,
        y = as.numeric(group) - width / 2 - offset,
        yend = as.numeric(group) + width / 2 - offset
      ),
      color = "black",
      size = 0.3,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = box_stats,
      aes(
        x = lower_whisker, xend = Q1,
        y = as.numeric(group) - offset,
        yend = as.numeric(group) - offset
      ),
      color = "black",
      linetype = "dashed",
      size = 0.3,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = box_stats,
      aes(
        x = Q3, xend = upper_whisker,
        y = as.numeric(group) - offset,
        yend = as.numeric(group) - offset
      ),
      color = "black",
      linetype = "dashed",
      size = 0.3,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = box_stats %>% unnest(outliers),
      aes(
        x = outliers,
        y = as.numeric(group) - offset
      ),
      color = "gray30",
      size = 0.8,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = c(
        "deep_water" = "#e74c3c",
        "mine"       = "#0000FF",
        "oil"        = "#e67e22",
        "shale"      = "#993B90"
      )
    ) +
    scale_x_reverse() + 
    scale_y_discrete(position = "right") +
    labs(
      x = "Depth (m)",
      y = "Group"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 13, face = "bold"),
      axis.title.y = element_text(size = 13, face = "bold"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    coord_flip())

ggsave(paste0("fig_depth_ridge.pdf"),fig_depth_ridge, width = 4, height = 3)

# Note: Minor aesthetic adjustments (e.g., font size, label placement, and spacing) were performed in Adobe Illustrator to improve figure clarity.
