# =========================================================
# Edehzhie resilience metric change plots
# =========================================================

require(pacman)

pacman::p_load(
  dplyr, purrr, glue, qs, tibble, tidyr,
  ggplot2, readr
)

options(scipen = 999)

# ---------------------------------------------------------
# INPUTS
# ---------------------------------------------------------
in_dir <- "./outputs/resilience_metrics_edehzhie"

data <- qs::qread(
  glue("{in_dir}/resilience_metrics_Edehzhie.qs")
)

out_dir <- "./figs/resilience_metrics_edehzhie"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ---------------------------------------------------------
# Clean data
# ---------------------------------------------------------
data <- data %>%
  mutate(
    gcm = tidyr::replace_na(gcm, "INM-CM4"),
    yr = as.character(yr)
  )

# ---------------------------------------------------------
# Pivot resilience metrics to long format
# ---------------------------------------------------------
data_long <- data %>%
  pivot_longer(
    cols = c(
      cross_scale_red,
      within_scale_red,
      cross_scale_div
    ),
    names_to = "metric",
    values_to = "value"
  )

# ---------------------------------------------------------
# A) Mean per GCM x year x metric
# ---------------------------------------------------------
yearly <- data_long %>%
  group_by(gcm, yr, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    n    = sum(!is.na(value)),
    se   = sd / sqrt(pmax(n, 1)),
    .groups = "drop"
  )

# ---------------------------------------------------------
# B) Helper function: delta change by GCM
# Delta = mean(year2) - mean(year1)
# ---------------------------------------------------------
pair_delta <- function(df, y1 = "2011", y2 = "2031") {
  
  left <- df %>%
    filter(yr == y1) %>%
    rename(mean1 = mean, se1 = se)
  
  right <- df %>%
    filter(yr == y2) %>%
    rename(mean2 = mean, se2 = se)
  
  both <- inner_join(
    left %>% select(gcm, metric, mean1, se1),
    right %>% select(gcm, metric, mean2, se2),
    by = c("gcm", "metric")
  )
  
  if (nrow(both) == 0) {
    return(
      tibble(
        gcm = character(),
        metric = character(),
        yr_pair = character(),
        delta_change = numeric(),
        delta_change_se = numeric()
      )
    )
  }
  
  tibble(
    gcm = both$gcm,
    metric = both$metric,
    yr_pair = glue("{y2} - {y1}"),
    delta_change = both$mean2 - both$mean1,
    delta_change_se = sqrt((both$se1^2) + (both$se2^2))
  )
}

# ---------------------------------------------------------
# C) Delta changes per GCM
# ---------------------------------------------------------
per_gcm_changes <- bind_rows(
  pair_delta(yearly, "2011", "2031"),
  pair_delta(yearly, "2011", "2091")
)

# ---------------------------------------------------------
# D) Average delta changes across GCMs
# mean delta across GCMs
# SE = sd / sqrt(k)
# ---------------------------------------------------------
avg_changes <- per_gcm_changes %>%
  group_by(metric, yr_pair) %>%
  summarise(
    k = sum(is.finite(delta_change)),
    delta_change = mean(delta_change, na.rm = TRUE),
    delta_change_se = ifelse(
      k > 1,
      sd(delta_change, na.rm = TRUE) / sqrt(k),
      NA_real_
    ),
    .groups = "drop"
  )

print(avg_changes)

qs::qsave(
  avg_changes,
  glue("{in_dir}/mean_delta_resilience_metrics_Edehzhie.qs")
)

write_csv(
  avg_changes,
  glue("{in_dir}/mean_delta_resilience_metrics_Edehzhie.csv")
)

# ---------------------------------------------------------
# E) Wide data for CSR vs WSR plot
# ---------------------------------------------------------
x_tbl <- avg_changes %>%
  filter(metric == "cross_scale_red") %>%
  transmute(
    yr_pair,
    x = delta_change,
    x_se = delta_change_se
  )

y_tbl <- avg_changes %>%
  filter(metric == "within_scale_red") %>%
  transmute(
    yr_pair,
    y = delta_change,
    y_se = delta_change_se
  )

wide_data <- inner_join(
  x_tbl,
  y_tbl,
  by = "yr_pair"
)

# ---------------------------------------------------------
# F) Plot 1: Delta cross-scale vs within-scale redundancy
# ---------------------------------------------------------
plot_delta_changes <- ggplot(
  wide_data,
  aes(
    x = x,
    y = y,
    color = yr_pair,
    group = yr_pair
  )
) +
  geom_point(size = 4, alpha = 0.9) +
  geom_errorbarh(
    aes(
      xmin = x - x_se,
      xmax = x + x_se
    ),
    height = 0.02,
    alpha = 0.7,
    orientation = "y"
  ) +
  geom_errorbar(
    aes(
      ymin = y - y_se,
      ymax = y + y_se
    ),
    width = 0.02,
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "darkgray",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "darkgray",
    linewidth = 0.6
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  scale_color_manual(
    NULL,
    values = c(
      "2031 - 2011" = "deeppink1",
      "2091 - 2011" = "turquoise3"
    )
  ) +
  labs(
    x = "Delta cross-scale redundancy",
    y = "Delta within-scale redundancy",
    title = "Edehzhie: Delta change in redundancy metrics"
  )

print(plot_delta_changes)

ggsave(
  filename = glue("{out_dir}/Edehzhie_delta_CSR_WSR.png"),
  plot = plot_delta_changes,
  width = 7,
  height = 5,
  dpi = 300
)

# ---------------------------------------------------------
# G) Plot 2: Delta cross-scale diversity
# ---------------------------------------------------------
csd_data <- avg_changes %>%
  filter(metric == "cross_scale_div")

plot_delta_div <- ggplot(
  csd_data,
  aes(
    x = yr_pair,
    y = delta_change,
    color = yr_pair
  )
) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(
      ymin = delta_change - delta_change_se,
      ymax = delta_change + delta_change_se
    ),
    width = 0.12,
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "darkgray",
    linewidth = 0.6
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  scale_color_manual(
    values = c(
      "2031 - 2011" = "deeppink1",
      "2091 - 2011" = "turquoise3"
    )
  ) +
  labs(
    x = "Time period",
    y = "Delta cross-scale diversity",
    title = "Edehzhie: Delta change in cross-scale diversity"
  )

print(plot_delta_div)

ggsave(
  filename = glue("{out_dir}/Edehzhie_delta_cross_scale_diversity.png"),
  plot = plot_delta_div,
  width = 6,
  height = 5,
  dpi = 300
)

message("Done.")
