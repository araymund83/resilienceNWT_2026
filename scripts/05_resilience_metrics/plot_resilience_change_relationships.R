### =========================================================
  # Resilience metrics plots
  # =========================================================

require(pacman)
pacman::p_load(
  dplyr, purrr, glue, qs, tibble, tidyr, sf
)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# ---------------------------------------------------------
# INPUTS
# ---------------------------------------------------------

in_dir <- "./outputs/ResilienceMetrics_4km"
data <- qs::qread(glue("{in_dir}/resilience_metrics_all.qs"))


# Load map / region polygons ---------------------------------------------
targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

ecrg <- sf::st_read("inputs/ecoregions/NWT_ecoregions_5Regions.shp", quiet = TRUE)

limt <- sf::st_read('inputs/NT1_BCR6/NT1_BCR6_poly.shp') 

# Extract by mask for the ecoregions ---------------------------------------
plot(st_geometry(ecrg))
limt <- sf::st_transform(x = limt, crs = targetCRS)
ecrg <- sf::st_transform(x = ecrg, crs = targetCRS)
plot(st_geometry(ecrg))
plot(st_geometry(limt))

# ---------------------------------------------------------
# Add region to each pixel
# ---------------------------------------------------------
pts <- data %>%
  distinct(pixelID, lat, long) %>%
  sf::st_as_sf(coords = c("long", "lat"), crs = targetCRS, remove = FALSE)

pts_region <- sf::st_join(pts, ecrg[, c("ECO3_NAM_1")], left = TRUE) %>%
  sf::st_drop_geometry()

data2 <- data %>%
  left_join(
    pts_region %>% select(pixelID, ECO3_NAM_1),
    by = "pixelID"
  )

# ---------------------------------------------------------
# Pivot resilience metrics to long format
# ---------------------------------------------------------
data_long <- data2 %>%
  pivot_longer(
    cols = c(cross_scale_red, within_scale_red, cross_scale_div),
    names_to = "metric",
    values_to = "value"
  )

# ---------------------------------------------------------
# A) Yearly summaries per GCM × Region × Metric
# ---------------------------------------------------------
yearly <- data_long %>%
  filter(!is.na(ECO3_NAM_1)) %>%
  group_by(gcm, yr, ECO3_NAM_1, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = sum(!is.na(value)),
    se   = sd / sqrt(pmax(n, 1)),
    .groups = "drop"
  )

# ---------------------------------------------------------
# B) Helper: delta change + SE by self-join
# Delta = mean(year2) - mean(year1)
# ---------------------------------------------------------
pair_delta <- function(df, y1 = "2011", y2 = "2031") {
  left  <- df %>% filter(yr == y1) %>% rename(mean1 = mean, se1 = se)
  right <- df %>% filter(yr == y2) %>% rename(mean2 = mean, se2 = se)
  
  both <- inner_join(
    left %>% select(gcm, ECO3_NAM_1, metric, mean1, se1),
    right %>% select(gcm, ECO3_NAM_1, metric, mean2, se2),
    by = c("gcm", "ECO3_NAM_1", "metric")
  )
  
  if (nrow(both) == 0) {
    return(tibble(
      gcm = character(), ECO3_NAM_1 = character(), metric = character(),
      yr_pair = character(), delta_change = numeric(), delta_change_se = numeric()
    ))
  }
  
  delta_change <- both$mean2 - both$mean1
  delta_se     <- sqrt((both$se1^2) + (both$se2^2))
  
  tibble(
    gcm = both$gcm,
    ECO3_NAM_1 = both$ECO3_NAM_1,
    metric = both$metric,
    yr_pair = glue("{y2} - {y1}"),
    delta_change = delta_change,
    delta_change_se = delta_se
  )
}

per_gcm_changes <- bind_rows(
  pair_delta(yearly, "2011", "2031"),
  pair_delta(yearly, "2011", "2091")
)

# ---------------------------------------------------------
# C) Average across GCMs
# mean delta across GCMs
# SE = sd / sqrt(k)
# ---------------------------------------------------------
avg_changes <- per_gcm_changes %>%
  group_by(ECO3_NAM_1, metric, yr_pair) %>%
  summarise(
    vals = list(delta_change[is.finite(delta_change)]),
    k    = lengths(vals),
    .groups = "drop"
  ) %>%
  mutate(
    delta_change = ifelse(k > 0, sapply(vals, mean), NA_real_),
    delta_change_se = ifelse(k > 1, sapply(vals, function(v) sd(v) / sqrt(length(v))), NA_real_)
  ) %>%
  select(-vals)


# ---------------------------------------------------------
# D) Panel A data: CSR vs WSR
# ---------------------------------------------------------
x_tbl <- avg_changes %>%
  filter(metric == "cross_scale_red") %>%
  transmute(ECO3_NAM_1, yr_pair, x = delta_change, x_se = delta_change_se)

y_tbl <- avg_changes %>%
  filter(metric == "within_scale_red") %>%
  transmute(ECO3_NAM_1, yr_pair, y = delta_change, y_se = delta_change_se)

wide_data <- inner_join(x_tbl, y_tbl, by = c("ECO3_NAM_1", "yr_pair")) %>%
  filter(!is.na(ECO3_NAM_1)) %>%
  mutate(
    ECO3_NAM_1 = factor(
      ECO3_NAM_1,
      levels = c("Mid-Boreal", "High Boreal", "Low Subarctic", "High Subarctic", "Low Arctic north")
    )
  )

# ---------------------------------------------------------
# E) Panel A plot
# ---------------------------------------------------------
plot_delta_changes <- ggplot(
  wide_data,
  aes(x = x, y = y, color = yr_pair, group = interaction(ECO3_NAM_1, yr_pair))
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_errorbarh(
    aes(xmin = x - x_se, xmax = x + x_se),
    height = 0.05, alpha = 0.6, orientation = 'y'
  ) +
  geom_errorbar(
    aes(ymin = y - y_se, ymax = y + y_se),
    width = 0.05, alpha = 0.6
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 0.6) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'right',
    legend.key.width = unit(2, 'line'),
    plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 14, margin = margin(3, 0, 3, 0))
  ) +
  scale_color_manual(
    NULL,
    values = c("2031 - 2011" = "deeppink1", "2091 - 2011" = "turquoise3")
  ) +
  labs(
    x = "Delta cross-scale redundancy",
    y = "Delta within-scale redundancy",
    title = "Delta change: Cross-scale vs. Within-scale Redundancy"
  ) +
  facet_wrap(~ ECO3_NAM_1)

# ---------------------------------------------------------
# F) Panel B plot: Delta cross-scale diversity by region
# ---------------------------------------------------------
csd_data <- avg_changes %>%
  filter(metric == "cross_scale_div") %>%
  filter(!is.na(ECO3_NAM_1)) %>%
  mutate(
    ECO3_NAM_1 = factor(
      ECO3_NAM_1,
      levels = c("Mid-Boreal", "High Boreal", "Low Subarctic", "High Subarctic", "Low Arctic north")
    )
  )

plot_delta_div <- ggplot(
  csd_data,
  aes(x = ECO3_NAM_1, y = delta_change, color = yr_pair, group = yr_pair)
) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(
    aes(ymin = delta_change - delta_change_se, ymax = delta_change + delta_change_se),
    width = 0.15,
    position = position_dodge(width = 0.3)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 0.6) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'right',
    legend.key.width = unit(2, 'line'),
    plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  scale_color_manual(
    NULL,
    values = c("2031 - 2011" = "deeppink1", "2091 - 2011" = "turquoise3")
  ) +
  labs(
    x = "Region",
    y = "Delta cross-scale diversity",
    title = "Delta change in cross-scale diversity by region"
  )

# ---------------------------------------------------------
# Save plots
# ---------------------------------------------------------
out_dir <- "./figs/Redundancy/2026/delta_change"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(
  filename = glue("{out_dir}/delta_change_cross_vs_within_region.png"),
  plot = plot_delta_changes,
  width = 10,
  height = 8
)

ggsave(
  filename = glue("{out_dir}/delta_change_cross_scale_div_region.png"),
  plot = plot_delta_div,
  width = 7,
  height = 6
) 

library(ggpubr)
library(glue)

out_dir <- "./figs/Redundancy/2026/delta_change"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


p1 <- plot_delta_changes +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

p2 <- plot_delta_div +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )

leg <- cowplot::get_legend(p1)

p1_clean <- p1 + theme(legend.position = "none")
p2_clean <- p2 + theme(legend.position = "none")

plots_stacked <- cowplot::plot_grid(
  p1_clean,
  p2_clean,
  ncol = 1,
  labels = c("a)", "b)"),
  label_size = 18,
  label_fontface = "bold",
  align = "v",
  axis = "l"
)

combined_delta <- cowplot::plot_grid(
  plots_stacked,
  leg,
  ncol = 2,
  rel_widths = c(1, 0.22)
)

ggsave(
  filename = glue("{out_dir}/delta_change_combined_region_vertical.png"),
  plot = combined_delta,
  width = 12,
  height = 10,
  units = "in",
  dpi = 400,
  bg = "white"
)

#horizontal 

# 2. Extract legend from original plot
leg <- cowplot::get_legend(
  plot_delta_changes + theme(legend.position = "right")
)

# 3. Combine plots horizontally
plots_side <- cowplot::plot_grid(
  p1_clean,
  p2_clean,
  ncol = 2,
  labels = c("a)", "b)"),
  label_size = 14,
  label_fontface = "bold",
  align = "h",
  axis = "tb"
)

# 4. Add legend on the right
combined_delta <- cowplot::plot_grid(
  plots_side,
  leg,
  ncol = 2,
  rel_widths = c(1, 0.18)
)

# 5. Save
ggsave(
  filename = glue("{out_dir}/delta_change_combined_region_horizontal.png"),
  plot = combined_delta,
  width = 14,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

# put  png figures together 
library(magick)
library(glue)

out_dir <- "./figs/Redundancy/2026/delta_change"

img_a <- image_read(glue("{out_dir}/delta_change_cross_vs_within_region.png"))
img_b <- image_read(glue("{out_dir}/delta_change_cross_scale_div_region.png"))

# Make same height
target_h <- max(image_info(img_a)$height, image_info(img_b)$height)

img_a <- image_resize(img_a, paste0("x", target_h))
img_b <- image_resize(img_b, paste0("x", target_h))

# Add panel labels
img_a <- image_annotate(
  img_a, "a)", 
  size = 60, 
  gravity = "northwest", 
  location = "+35+30",
  weight = 600
)

img_b <- image_annotate(
  img_b, "b)", 
  size = 60, 
  gravity = "northwest", 
  location = "+35+30",
  weight = 700
)

# Optional: add spacing between panels
spacer <- image_blank(width = 40, height = target_h, color = "white")

combined_vertical <- image_append(c(img_a, img_b), stack = TRUE)

image_write(
  combined_vertical,
  path = glue("{out_dir}/delta_change_combined_region_magick_vertical.png"),
  format = "png"
)

#horizontal 
# Combine horizontally
combined <- image_append(c(img_a, spacer, img_b))  # horizontal

# Save
image_write(
  combined,
  path = glue("{out_dir}/delta_change_combined_region_horizontal.png"),
  format = "png"
)



