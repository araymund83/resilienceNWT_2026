# =========================================================
# Resilience metrics from aggs_df_list
# aggs_df_list = named list of flattened aggregation outputs
# one element per gcm_year, each a data.frame with:
# pixelID, species, aggs, gcm, yr, lat, long, ...
# =========================================================

require(pacman)
pacman::p_load(
  dplyr, purrr, glue, qs, tibble, tidyr
)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# ---------------------------------------------------------
# INPUTS
# ---------------------------------------------------------
# ---------------------------------------------------------
# Directory with aggregation outputs
# ---------------------------------------------------------
agg_dir <- "./outputs/Aggregations_4km_filtered"

files <- list.files(
  agg_dir,
  pattern = "Aggregations_.*\\.qs$",
  full.names = TRUE
)

# ---------------------------------------------------------
# Extract names (gcm + year)
# ---------------------------------------------------------
meta <- tibble(file = files) %>%
  mutate(
    gcm = stringr::str_match(file, "Aggregations_(.*)_(\\d{4})")[,2],
    yr  = stringr::str_match(file, "Aggregations_(.*)_(\\d{4})")[,3],
    name = glue("{gcm}_{yr}")
  )

# ---------------------------------------------------------
# Read into list
# ---------------------------------------------------------
aggs_df_list <- purrr::map(files, qs::qread)

# assign names like "CanESM2_2011"
names(aggs_df_list) <- meta$name


# ---------------------------------------------------------
# Check
# ---------------------------------------------------------
names(aggs_df_list)

str(aggs_df_list[[1]])
traitsFG <- qs::qread('./inputs/traitsMassFG_2025.qs') 
# ---------------------------------------------------------

# ---------------------------------------------------------
# Check FG lookup
# ---------------------------------------------------------
if (!"FG" %in% names(traitsFG)) {
  stop("traitsFG must contain a column named 'FG'")
}
if (is.null(rownames(traitsFG))) {
  stop("traitsFG must have species as rownames")
}

fg_lookup <- traitsFG$FG
names(fg_lookup) <- rownames(traitsFG)

# ---------------------------------------------------------
# Helper: add FG to one pixel data frame
# ---------------------------------------------------------
add_fg <- function(df_pixel, fg_lookup) {
  df_pixel %>%
    mutate(FG = fg_lookup[species]) %>%
    filter(!is.na(FG))
}

# ---------------------------------------------------------
# Metric 1: Cross-scale redundancy
# mean number of aggregations occupied per FG
# ---------------------------------------------------------
calc_cross_scale_red_pixel <- function(df_pixel, fg_lookup) {
  df_pixel <- add_fg(df_pixel, fg_lookup)
  if (nrow(df_pixel) == 0) return(NA_real_)
  
  phi_values <- xtabs(~ FG + aggs, data = df_pixel)
  
  # IMPORTANT:
  # use presence/absence across aggregations, not raw species counts
  fg_scales <- rowSums(phi_values > 0)
  
  mean(fg_scales)
}

# ---------------------------------------------------------
# Metric 2: Within-scale redundancy
# for each FG: average number of species per aggregation
# then mean across FGs
# ---------------------------------------------------------
calc_within_scale_red_pixel <- function(df_pixel, fg_lookup) {
  df_pixel <- add_fg(df_pixel, fg_lookup)
  if (nrow(df_pixel) == 0) return(NA_real_)
  
  tmp <- df_pixel %>%
    group_by(FG, aggs) %>%
    summarise(num_sp = n_distinct(species), .groups = "drop") %>%
    group_by(FG) %>%
    summarise(avg_num_sp = mean(num_sp), .groups = "drop")
  
  mean(tmp$avg_num_sp)
}

# ---------------------------------------------------------
# Metric 3: Cross-scale diversity
# Fallback Shannon version:
# diversity of FG occupancy across aggregations
#
# If you already trust get_cross_scale_diversity2() from your sourced
# functions file, you can swap this function out.
# ---------------------------------------------------------
calc_cross_scale_div_pixel <- function(df_pixel, fg_lookup) {
  df_pixel <- add_fg(df_pixel, fg_lookup)
  if (nrow(df_pixel) == 0) return(NA_real_)
  
  phi_values <- xtabs(~ FG + aggs, data = df_pixel)
  fg_scales <- rowSums(phi_values > 0)
  
  p <- fg_scales / sum(fg_scales)
  -sum(p * log(p))
}

# ---------------------------------------------------------
# Calculate all metrics for one flattened aggregation data frame
# ---------------------------------------------------------
calc_metrics_one_df <- function(df, fg_lookup) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  
  required_cols <- c("pixelID", "species", "aggs", "gcm", "yr", "lat", "long")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(glue("Missing required columns: {paste(missing_cols, collapse = ', ')}"))
  }
  
  pixel_list <- split(df, df$pixelID)
  
  metrics <- purrr::map_dfr(pixel_list, function(px) {
    tibble(
      pixelID = as.character(px$pixelID[1]),
      lat = px$lat[1],
      long = px$long[1],
      gcm = px$gcm[1],
      yr = px$yr[1],
      cross_scale_red = calc_cross_scale_red_pixel(px, fg_lookup),
      within_scale_red = calc_within_scale_red_pixel(px, fg_lookup),
      cross_scale_div = calc_cross_scale_div_pixel(px, fg_lookup),
      n_species = n_distinct(px$species),
      n_aggs = n_distinct(px$aggs)
    )
  })
  
  metrics
}

# ---------------------------------------------------------
# Run across all gcm_year elements in aggs_df_list
# ---------------------------------------------------------
resilience_metrics_list <- purrr::imap(
  aggs_df_list,
  ~ {
    message(glue("Calculating metrics for {.y}"))
    calc_metrics_one_df(.x, fg_lookup)
  }
)

# Combined table
resilience_metrics_df <- bind_rows(resilience_metrics_list, .id = "source_name")

# ---------------------------------------------------------
# Save combined table
# ---------------------------------------------------------
out_dir <- "./outputs/ResilienceMetrics_4km"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

qs::qsave(
  resilience_metrics_df,
  glue("{out_dir}/resilience_metrics_all.qs")
)

# ---------------------------------------------------------
# Save one file per gcm x yr
# ---------------------------------------------------------
split_out <- resilience_metrics_df %>%
  group_split(gcm, yr)

split_keys <- resilience_metrics_df %>%
  distinct(gcm, yr)

purrr::walk2(
  split_out,
  seq_len(nrow(split_keys)),
  function(df_sub, i) {
    gcm_i <- split_keys$gcm[i]
    yr_i  <- split_keys$yr[i]
    
    qs::qsave(
      df_sub,
      glue("{out_dir}/resilience_metrics_{gcm_i}_{yr_i}.qs")
    )
  }
)

# ---------------------------------------------------------
# Quick summary
# ---------------------------------------------------------
summary_tab <- resilience_metrics_df %>%
  group_by(gcm, yr) %>%
  summarise(
    n_pixels = n(),
    mean_cross_scale_red = mean(cross_scale_red, na.rm = TRUE),
    mean_within_scale_red = mean(within_scale_red, na.rm = TRUE),
    mean_cross_scale_div = mean(cross_scale_div, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tab)

qs::qsave(
  summary_tab,
  glue("{out_dir}/resilience_metrics_summary.qs")
)