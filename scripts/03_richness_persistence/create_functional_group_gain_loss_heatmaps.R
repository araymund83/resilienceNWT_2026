library(pacman)
p_load(glue, ggplot2,ggalluvial, glue, qs, tidyverse)

# Session/options
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Create output directory
dir.create("./outputs/fg_presence_by_region", showWarnings = FALSE)

# Define parameters
regions <- c("High_Boreal", "Mid_Boreal", "Low_Subarctic", 
             "High_Subarctic", "Low_Arctic_north")
gcms <- c("CanESM2", "CCSM4", "INM-CM4")
yrs <- c(2011, 2031, 2091)
path <- "./tables/fg_matrices_allgcms_by_scale"
files <- list.files(path, pattern = '4km', full.names = TRUE)

pix_eco <- qs::qread('./inputs/pixel_lookup_eco_4km.qs')

pix_eco <- pix_eco %>%
  mutate(
    ECO3_NAM_1 = case_when(
      ECO3_NAM_1 %in% c("Peel River Plateau", "Fort MacPherson Plain") ~ "Low Subarctic",
      ECO3_NAM_1 == "British-Richardson Mountains" ~ "High Subarctic",
      TRUE ~ ECO3_NAM_1
    )
  )

fg_data <- lapply(files, qs::qread)

names(fg_data) <- gcms

fg_data2 <- bind_rows(fg_data)


# Identify FG columns
fg_cols <- names(fg_data2) %>%
  str_subset("^(GRA|INS|OMN)-")

# Long format
fg_long <- fg_data2 %>%
  mutate(yr = as.integer(yr)) %>%
  pivot_longer(
    cols = all_of(fg_cols),
    names_to = "FG",
    values_to = "presence"
  ) %>%
  mutate(
    presence = if_else(presence > 0, 1L, 0L)
  )


fg_region <- fg_long %>%
  mutate(pixelID = as.character(pixelID)) %>%
  left_join(
    pix_eco %>% mutate(pixelID = as.character(pixelID)),
    by = c("pixelID", 'x', 'y')
  )

fg_region <- fg_region %>%
  mutate(
    region = ECO3_NAM_1
  ) %>%
  filter(!is.na(region))
  
# Diverging palette like TBI maps -----------------------------------------
tbi_cols <- c(
  "#a6611a", "#dfc27d", "white", "#80cdc1", "#018571"
)

region_order <- c("High_Boreal", "Mid_Boreal", "Low_Subarctic", 
                  "High_Subarctic", "Low_Arctic_north")
fg_levels <- fg_region %>%
  distinct(FG) %>%
  pull(FG) %>%
  sort()

# Function to summarize gains/losses -------------------------------------
make_gain_loss_heatmap_df <- function(transitions_df, timeframe_label) {
  
  transitions_df %>%
    mutate(
      change_value = case_when(
        transition == "appear" ~ 1L,
        transition == "disappear" ~ -1L,
        TRUE ~ 0L
      )
    ) %>%
    group_by(region, gcm, FG) %>%
    summarise(
      gains = sum(transition == "appear", na.rm = TRUE),
      losses = sum(transition == "disappear", na.rm = TRUE),
      net_change = gains - losses,
      .groups = "drop"
    ) %>%
    mutate(
      timeframe = timeframe_label,
      region = case_when(
        region == "Mid-Boreal"       ~ "Mid_Boreal",
        region == "High Boreal"      ~ "High_Boreal",
        region == "Low Subarctic"    ~ "Low_Subarctic",
        region == "High Subarctic"   ~ "High_Subarctic",
        region == "Low Arctic north" ~ "Low_Arctic_north",
        TRUE ~ region
      ),
      region = factor(region, levels = region_order),
      gcm = factor(gcm, levels = gcms),
      FG = factor(FG, levels = rev(fg_levels))
    )
}

transitions_31 <- fg_region %>%
  filter(yr %in% c(2011, 2031)) %>%
  select(pixelID, region, gcm, yr, FG, presence) %>%
  pivot_wider(
    names_from = yr,
    values_from = presence,
    names_prefix = "yr_",
    values_fill = 0
  ) %>%
  mutate(
    transition = case_when(
      yr_2011 == 0 & yr_2031 == 1 ~ "appear",
      yr_2011 == 1 & yr_2031 == 0 ~ "disappear",
      yr_2011 == 1 & yr_2031 == 1 ~ "persist",
      TRUE ~ "absent"
    )
  )

# -------------------------------------------------------------------------
# Create transitions: 2011 -> 2091
# -------------------------------------------------------------------------

transitions_91 <- fg_region %>%
  filter(yr %in% c(2011, 2091)) %>%
  select(pixelID, region, gcm, yr, FG, presence) %>%
  pivot_wider(
    names_from = yr,
    values_from = presence,
    names_prefix = "yr_",
    values_fill = 0
  ) %>%
  mutate(
    transition = case_when(
      yr_2011 == 0 & yr_2091 == 1 ~ "appear",
      yr_2011 == 1 & yr_2091 == 0 ~ "disappear",
      yr_2011 == 1 & yr_2091 == 1 ~ "persist",
      TRUE ~ "absent"
    )
  )

heat_31 <- make_gain_loss_heatmap_df(transitions_31, "2011–2031")
heat_91 <- make_gain_loss_heatmap_df(transitions_91, "2011–2091")

plot_gain_loss_heatmap <- function(df, title_text) {
  
  max_abs <- max(abs(df$net_change), na.rm = TRUE)
  
  ggplot(df, aes(x = region, y = FG, fill = net_change)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(net_change == 0, "", net_change)),
      size = 3
    ) +
    facet_wrap(~ gcm, nrow = 1) +
    scale_fill_gradientn(
      colours = tbi_cols,
      limits = c(-max_abs, max_abs),
      values = scales::rescale(c(-max_abs, -max_abs / 2, 0, max_abs / 2, max_abs)),
      name = "Net change\n(gains - losses)"
    ) +
    scale_x_discrete(
      labels = c(
        "Mid_Boreal" = "Mid Boreal",
        "High_Boreal" = "High Boreal",
        "Low_Subarctic" = "Low Subarctic",
        "High_Subarctic" = "High Subarctic",
        "Low_Arctic_north" = "Low Arctic north"
      )
    ) +
    labs(
      title = title_text,
      x = "Region",
      y = "Functional group"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

fig_gain_loss_31 <- plot_gain_loss_heatmap(
  heat_31,
  "Functional group gains and losses, 2011–2031"
)

fig_gain_loss_91 <- plot_gain_loss_heatmap(
  heat_91,
  "Functional group gains and losses, 2011–2091"
)

fig_gain_loss_31
fig_gain_loss_91

ggsave(
  "./figs/persistance/FG_gain_loss_heatmap_2011_2031.png",
  fig_gain_loss_31,
  width = 11,
  height = 5.5,
  dpi = 300
)

ggsave(
  "./figs/persistance/FG_gain_loss_heatmap_2011_2091.png",
  fig_gain_loss_91,
  width = 11,
  height = 5.5,
  dpi = 300
)

### proportions 
make_gain_loss_heatmap_df <- function(transitions_df, timeframe_label) {
  
  transitions_df %>%
    group_by(region, gcm, FG) %>%
    summarise(
      n_total = n(),
      n_gains = sum(transition == "appear", na.rm = TRUE),
      n_losses = sum(transition == "disappear", na.rm = TRUE),
      
      gain_prop = n_gains / n_total,
      loss_prop = n_losses / n_total,
      
      signed_prop = gain_prop - loss_prop,
      
      .groups = "drop"
    ) %>%
    mutate(
      timeframe = timeframe_label,
      
      region = case_when(
        region == "Mid-Boreal"       ~ "Mid_Boreal",
        region == "High Boreal"      ~ "High_Boreal",
        region == "Low Subarctic"    ~ "Low_Subarctic",
        region == "High Subarctic"   ~ "High_Subarctic",
        region == "Low Arctic north" ~ "Low_Arctic_north",
        TRUE ~ region
      ),
      
      region = factor(region, levels = region_order),
      gcm = factor(gcm, levels = gcms),
      FG = factor(FG, levels = rev(fg_levels))
    )
}


heat_31 <- make_gain_loss_heatmap_df(
  transitions_31,
  "2011–2031"
)

heat_91 <- make_gain_loss_heatmap_df(
  transitions_91,
  "2011–2091"
)


plot_gain_loss_heatmap <- function(df, subtitle_text) {
  
  max_abs <- max(abs(df$signed_prop), na.rm = TRUE)
  
  ggplot(df, aes(x = region, y = FG, fill = signed_prop)) +
    
    geom_tile(color = "white", linewidth = 0.3) +
        geom_text(
      aes(label = round(signed_prop, 2)),
      size = 3
    ) +
    facet_wrap(~ gcm, nrow = 1) +
    scale_fill_gradient2(
      low = "#a6611a",
      mid = "white",
      high = "#018571",
      midpoint = 0,
      limits = c(-max_abs, max_abs),
      breaks = c(-max_abs, 0, max_abs),
      labels = c("Loss", "Neutral", "Gain"),
      name = "Proportion"
    ) +
    labs(
      title = "Proportion of pixels showing functional group gains and losses",
      subtitle = subtitle_text,
      x = "Region",
      y = "Functional group"
    ) +
    scale_x_discrete(
      labels = c(
        "Mid_Boreal" = "Mid-Boreal",
        "High_Boreal" = "High Boreal",
        "Low_Subarctic" = "Low Subarctic",
        "High_Subarctic" = "High Subarctic",
        "Low_Arctic_north" = "Low Arctic north"
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = 'bold'),
      plot.subtitle = element_text(size = 12, face = 'bold'),
      legend.key.width = unit(2, 'line'),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = 'bold')
    )
}


fig_gain_loss_31 <- plot_gain_loss_heatmap(
  heat_31,
  "(2011–2031)"
)

fig_gain_loss_91 <- plot_gain_loss_heatmap(
  heat_91,
  "(2011–2091)"
)

fig_gain_loss_31
fig_gain_loss_91

ggsave(
  "./figs/persistance/FG_gain_loss_heatmap_proportion_2011_2031.png",
  fig_gain_loss_31,
  width = 11,
  height = 5.5,
  dpi = 300
)

ggsave(
  "./figs/persistance/FG_gain_loss_heatmap_proportion_2011_2091.png",
  fig_gain_loss_91,
  width = 11,
  height = 5.5,
  dpi = 300
)

## put figures together 
library(ggpubr)

combined_fg_heatmaps <- ggarrange(
  fig_gain_loss_31,
  fig_gain_loss_91,
  ncol = 1,
  nrow = 2,
  align = "v",
  labels = c("a)", "b)"),
  common.legend = TRUE,
  legend = "right"
)

combined_fg_heatmaps


ggsave(
  "./figs/persistance/FG_gain_loss_heatmaps_prop_combined.png",
  combined_fg_heatmaps,
  width = 12,
  height = 10,
  dpi = 300
)
