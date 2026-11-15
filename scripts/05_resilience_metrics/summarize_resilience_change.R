# Libraries ---------------------------------------------------------------
library(pacman)
p_load(glue, ggplot2, qs, sf, tidyverse, magick)

# Clean environment -------------------------------------------------------
gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Paths -------------------------------------------------------------------
in_path <- "./outputs/ResilienceMetrics_4km"
out_dir <- "./figs/redundancy_delta/2026"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load data ---------------------------------------------------------------
data <- qs::qread(glue('{in_path}/resilience_metrics_all.qs'))

data_long <- data %>% 
  pivot_longer(cols = c(cross_scale_red, within_scale_red, cross_scale_div),
               names_to = 'metric',
               values_to = 'value')

baseline <- data_long %>%
  filter(yr == "2011") %>%
  select(pixelID, gcm, metric, value_2011 = value)

delta_px <- data_long %>%
  filter(yr %in% c("2031", "2091")) %>%
  left_join(baseline, by = c("pixelID", "gcm", "metric")) %>%
  mutate(
    delta = value - value_2011,
    factoryr = case_when(
      yr == "2031" ~ "2011 - 2031",
      yr == "2091" ~ "2011 - 2091",
      TRUE ~ yr
    )
  )

# save if useful
qs::qsave(delta_px, "./outputs/ResilienceMetrics_4km/delta_px_from_2011.qs")

# ---------------------------------------------------------
# Mean delta and SE across GCMs per pixel
# ---------------------------------------------------------
mean_delta_px <- delta_px %>%
  group_by(pixelID, lat, long, metric, factoryr) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    sd_delta   = sd(delta, na.rm = TRUE),
    k          = sum(!is.na(delta)),
    se_delta   = sd_delta / sqrt(pmax(k, 1)),
    .groups = "drop"
  )

qs::qsave(mean_delta_px, "./outputs/ResilienceMetrics_4km/mean_delta_px_across_gcms.qs")

# ---------------------------------------------------------
# Load map information
# ---------------------------------------------------------
targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

limt <- sf::st_read("inputs/NT1_BCR6/NT1_BCR6_poly.shp", quiet = TRUE)
ecrg <- sf::st_read("inputs/ecoregions/NWT_ecoregions_5Regions.shp", quiet = TRUE)

limt <- sf::st_transform(limt, crs = targetCRS)
ecrg <- sf::st_transform(ecrg, crs = targetCRS)

# ---------------------------------------------------------
# Metric names
# ---------------------------------------------------------
metric_titles <- c(
  cross_scale_red  = "Cross-scale redundancy",
  within_scale_red = "Within-scale redundancy",
  cross_scale_div  = "Cross-scale diversity"
)

vars <- c("cross_scale_red", "within_scale_red", "cross_scale_div")

# ---------------------------------------------------------
# Palettes
# mean delta = diverging
# se delta   = sequential
# ---------------------------------------------------------

color_pal <- c("#543005", "#8c510a", "#bf812d", 
               "white",
               "#35978f", "#01665e", "#003c30")


se_pal <- c("white", "#dadaeb", "#bcbddc", "#9e9ac8", "#756bb1", "#54278f")

# ---------------------------------------------------------
# Function to make mean delta maps
# ---------------------------------------------------------
make_mean_delta_maps <- function(data, vars) {
  
  ggplot_list <- lapply(vars, function(fill_var) {
    
    title <- metric_titles[[fill_var]]
    
    var_data <- data %>%
      filter(metric == fill_var)
    
    max_abs <- max(abs(var_data$mean_delta), na.rm = TRUE)
    min_abs <- min(abs(var_data$mean_delta), na.rm = TRUE)
    
    ggplot() +
      geom_tile(
        data = var_data,
        aes(x = long, y = lat, fill = mean_delta),
        width = 4000,
        height = 4000
      ) +
      scale_fill_gradientn(colors = color_pal,
                           values = scales::rescale(c(-max_abs, 0, max_abs)),
                           limits = c(-max_abs, max_abs),
                           na.value = "#D8E2E7",  # light blue-grey
                           guide = guide_colorbar(barwidth = 10, barheight = 1)) +
      geom_sf(data = limt, fill = NA, col = "#36454F", linewidth = 0.3) +
      geom_sf(data = ecrg, fill = NA, col = "#36454F", linewidth = 0.2) +
      ggtitle(title) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(2, "line"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 14)
      ) +
      labs(
        x = "Longitude",
        y = "Latitude",
        fill = "Mean Δ"
      ) +
      coord_sf(lims_method = "geometry_bbox") +
      facet_wrap(~ factoryr)
  })
  
  names(ggplot_list) <- vars
  ggplot_list
}

# ---------------------------------------------------------
# Function to make SE maps
# ---------------------------------------------------------
make_se_delta_maps <- function(data, vars) {
  
  ggplot_list <- lapply(vars, function(fill_var) {
    
    title <- paste0(metric_titles[[fill_var]], " (SE)")
    
    var_data <- data %>%
      filter(metric == fill_var)
    
    max_se <- max(var_data$se_delta, na.rm = TRUE)
    
    ggplot() +
      geom_tile(
        data = var_data,
        aes(x = long, y = lat, fill = se_delta),
        width = 4000,
        height = 4000
      ) +
      scale_fill_gradientn(
        colours = se_pal,
        limits = c(0, max_se),
        na.value = "#D8E2E7",  # light blue-grey
        oob = scales::squish
      ) +
      geom_sf(data = limt, fill = NA, col = "#36454F", linewidth = 0.3) +
      geom_sf(data = ecrg, fill = NA, col = "#36454F", linewidth = 0.2) +
      ggtitle(title) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(2, "line"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 14)
      ) +
      labs(
        x = "Longitude",
        y = "Latitude",
        fill = "SE Δ"
      ) +
      coord_sf(lims_method = "geometry_bbox") +
      facet_wrap(~ factoryr)
  })
  
  names(ggplot_list) <- vars
  ggplot_list
}

# ---------------------------------------------------------
# Make plots
# ---------------------------------------------------------
mean_plots <- make_mean_delta_maps(mean_delta_px, vars)
se_plots   <- make_se_delta_maps(mean_delta_px, vars)

# ---------------------------------------------------------
# Save plots
# ---------------------------------------------------------
out_dir <- "./figs/redundancy_delta/2026"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# mean delta
lapply(seq_along(mean_plots), function(i) {
  plot <- mean_plots[[i]]
  var_name <- names(mean_plots)[i]
  
  ggsave(
    filename = glue("{out_dir}/{var_name}_mean_delta.png"),
    plot = plot,
    units = "in",
    width = 5,
    height = 5,
    dpi = 300
  )
})

# SE delta
lapply(seq_along(se_plots), function(i) {
  plot <- se_plots[[i]]
  var_name <- names(se_plots)[i]
  
  ggsave(
    filename = glue("{out_dir}/{var_name}_se_delta.png"),
    plot = plot,
    units = "in",
    width = 5,
    height = 5,
    dpi = 300
  )
})

# ---------------------------------------------------------
# Optional: print one example
# ---------------------------------------------------------
print(mean_plots[["cross_scale_red"]])
print(se_plots[["cross_scale_red"]])


## put mean and se maps together one per metric 
library(ggpubr)
library(glue)

out_dir <- "./figs/redundancy_delta/2026"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

plot_6 <- list(
  mean_plots[["cross_scale_red"]],
  mean_plots[["within_scale_red"]],
  mean_plots[["cross_scale_div"]],
  se_plots[["cross_scale_red"]],
  se_plots[["within_scale_red"]],
  se_plots[["cross_scale_div"]]
)

# Check all are ggplots
sapply(plot_6, inherits, what = "ggplot")

combined_all <- ggpubr::ggarrange(
  plotlist = plot_6,
  ncol = 3,
  nrow = 2,
  labels = c("a)", "b)", "c)", "d)", "e)", "f)"),
  font.label = list(size = 18, face = "bold")
)

combined_all_final <- ggpubr::annotate_figure(combined_all)

ggpubr::ggexport(
  combined_all_final,
  filename = "./figs/redundancy_delta/2026/combined_all_metrics_mean_se.png",
  width = 4500,
  height = 3000,
  res = 300
)
