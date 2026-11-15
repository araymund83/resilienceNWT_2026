# Load libraries ----------------------------------------------------------
library(pacman)
p_load(dplyr, future, future.apply, ggplot2, glue, qs, tidyverse)
# Clean environment -------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())

# Source functions -------------------------------------------------------
source('./R/functions_discontinuities.R')


pwrtab <- read.csv('./R/GRI.constant.power.table.csv')

discont_list <- purrr::pmap(
  combos,
  function(gcm, yr) {
    nm <- glue("{gcm}_{yr}")
    message(crayon::blue(glue("Running discontinuities for {nm}")))
    
    Discontinuity_Barichievy(
      Data         = bodymass_list[[nm]],
      start.column = 4,
      end.column   = ncol(bodymass_list[[nm]]),
      IS_NA        = TRUE,
      Resolution   = 4000,
      Sample.N     = 1000,
      gcm          = gcm,
      yr           = yr,
      parallel     = FALSE
    )
  }
)

names(discont_list) <- glue("{combos$gcm}_{combos$yr}")

# ---------------------------------------------------------
# Run aggregations
# ---------------------------------------------------------
aggs_list <- purrr::pmap(
  combos,
  function(gcm, yr) {
    nm <- glue("{gcm}_{yr}")
    message(crayon::blue(glue("Assigning aggregations for {nm}")))
    
    Aggregations(
      GRI.output     = discont_list[[nm]],
      SetThresh      = FALSE,
      constantThresh = 0.8,
      yr             = yr,
      gcm            = gcm
    )
  }
)

names(aggs_list) <- glue("{combos$gcm}_{combos$yr}")

# ---------------------------------------------------------
# Flatten aggregation output to data frames
# ---------------------------------------------------------
flatten_aggs <- function(x) {
  x2 <- x[purrr::map_lgl(x, ~ is.data.frame(.x) && nrow(.x) > 0)]
  if (length(x2) == 0) return(NULL)
  bind_rows(x2)
}

aggs_df_list <- purrr::map(aggs_list, flatten_aggs)

# ---------------------------------------------------------
# Save outputs
# ---------------------------------------------------------
out_aggs <- "./outputs/Aggregations_4km_filtered"
if (!dir.exists(out_aggs)) dir.create(out_aggs, recursive = TRUE)

purrr::walk2(
  names(aggs_df_list),
  aggs_df_list,
  ~ qs::qsave(.y, glue("{out_aggs}/Aggregations_{.x}_filtered.qs"))
)

# ---------------------------------------------------------
# Optional summary
# ---------------------------------------------------------
summary_tab <- combos %>%
  mutate(
    n_pixels = purrr::map2_int(gcm, yr, ~ nrow(filter(data, gcm == .x, yr == .y)))
  )

print(summary_tab)
qs::qsave(summary_tab, "./outputs/Aggregations_4km_filtered/aggregation_summary.qs")


agg_count_list <- purrr::imap(aggs_df_list, function(df, nm) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>%
    dplyr::group_by(pixelID, lat, long, gcm, yr) %>%
    dplyr::summarise(
      n_aggs = dplyr::n_distinct(aggs),
      n_spp  = dplyr::n(),
      .groups = "drop"
    )
})

agg_count_df <- dplyr::bind_rows(agg_count_list)

str(agg_count_df)
head(agg_count_df)

qs::qsave(
  agg_count_df,
  "./outputs/Aggregations_4km_filtered/aggregation_count_per_pixel.qs"
)

agg_count_df <- qs::qread("./outputs/Aggregations_4km_filtered/aggregation_count_per_pixel.qs")
#get richness by pixel
rich_df <- data %>%
  dplyr::select(pixelID, yr, gcm, richness) %>%
  dplyr::mutate(pixelID = as.character(pixelID))
#join richness and aggregations
scatter_df <- agg_count_df %>%
  dplyr::mutate(pixelID = as.character(pixelID)) %>%
  dplyr::left_join(rich_df, by = c("pixelID", "yr", "gcm"))

p <- ggplot(scatter_df, aes(x = n_aggs, y = richness)) +
  geom_point(alpha = 0.3, size = 0.5) +
  
  # smooth trend
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  
  facet_grid(yr ~ gcm) +
  
  labs(
    title = "Species richness vs body-mass aggregations",
    x = "Number of aggregations",
    y = "Species Richness"
  ) +
  
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

#get species per aggregation
agg_sizes_df <- purrr::imap_dfr(aggs_df_list, function(df, nm) {
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df %>%
    dplyr::group_by(pixelID, lat, long, gcm, yr, aggs) %>%
    dplyr::summarise(
      n_spp_in_agg = dplyr::n(),
      .groups = "drop"
    )
})
#get pixels that have at least one singleton aggregation
singleton_pixel_df <- agg_sizes_df %>%
  dplyr::group_by(pixelID, lat, long, gcm, yr) %>%
  dplyr::summarise(
    n_aggs = dplyr::n(),
    n_singleton_aggs = sum(n_spp_in_agg == 1),
    has_singleton = n_singleton_aggs > 0,
    .groups = "drop"
  )

#how many pixels have at least one 1-species aggregation?
singleton_summary <- singleton_pixel_df %>%
  dplyr::group_by(gcm, yr) %>%
  dplyr::summarise(
    n_pixels = dplyr::n(),
    n_pixels_with_singleton = sum(has_singleton),
    pct_pixels_with_singleton = 100 * n_pixels_with_singleton / n_pixels,
    .groups = "drop"
  )

singleton_summary

singleton_distribution <- singleton_pixel_df %>%
  dplyr::group_by(gcm, yr, n_singleton_aggs) %>%
  dplyr::summarise(
    n_pixels = dplyr::n(),
    .groups = "drop"
  )

singleton_distribution

library(ggplot2)

p_singleton <- ggplot(singleton_summary,
                      aes(x = yr, y = pct_pixels_with_singleton, fill = gcm)) +
  geom_col(position = "dodge") +
  labs(
    title = "Pixels with at least one 1-species aggregation",
    x = "Year",
    y = "% of pixels"
  ) +
  theme_bw()

print(p_singleton)
