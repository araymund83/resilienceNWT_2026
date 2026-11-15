# Load libraries ----------------------------------------------------------
library(pacman)
p_load(
  colorspace, dplyr, future, furrr, ggnewscale, ggplot2,
  glue, qs, scales, tidyverse
)

# Clean environment -------------------------------------------------------
g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Set parallel plan -------------------------------------------------------
plan(multisession, workers = 10)

# Parameters --------------------------------------------------------------
yrs <- c("2011", "2031", "2091")
gcms <- c("CanESM2", "CCSM4", "INM-CM4")

agg_dir <- "./outputs/Aggregations_4km_filtered"
out_root <- "./figs/aggregations2026/4km_filtered"

# Load FG data ------------------------------------------------------------
traitsFG <- qs::qread("./inputs/traitsMassFG_2025.qs") %>%
  dplyr::select(FG) %>%
  tibble::rownames_to_column("species")

# FG palette --------------------------------------------------------------
fg_levels <- sort(unique(traitsFG$FG))

fg_palette <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#A6CEE3", "#FB9A99", "#B2DF8A", "#CAB2D6"
)

if (length(fg_levels) > length(fg_palette)) {
  fg_palette <- qualitative_hcl(length(fg_levels), palette = "Dynamic")
}

fg_palette <- setNames(fg_palette[seq_along(fg_levels)], fg_levels)

# Aggregation line colors -------------------------------------------------
agg_colors <- c(
  "#A65041", "#E7CDC2", "#80A0C7", "#394165", "#A890A8",
  "#B1934A", "#DCA258", "#100F14", "#8B9DAF", "#EEDA9D", "#90A878"
)

# Read files --------------------------------------------------------------
files <- list.files(
  agg_dir,
  pattern = "^Aggregations_.*\\.qs$",
  full.names = TRUE
)

# Helper: make enough aggregation colors dynamically ----------------------
make_agg_palette <- function(n) {
  if (n <= length(agg_colors)) {
    setNames(agg_colors[seq_len(n)], as.character(seq_len(n)))
  } else {
    pal <- qualitative_hcl(n, palette = "Dark 3")
    setNames(pal, as.character(seq_len(n)))
  }
}

# Function to plot one pixel ----------------------------------------------
plot_one_pixel <- function(df_pixel, gcm, yr, traitsFG, out_dir) {
  
  if (is.null(df_pixel) || nrow(df_pixel) == 0) return(NULL)
  if (!all(c("species", "log.mass", "aggs", "pixelID") %in% names(df_pixel))) return(NULL)
  
  df <- df_pixel %>%
    left_join(traitsFG, by = "species") %>%
    arrange(log.mass)
  
  if (anyNA(df$FG)) {
    warning(glue("Some species did not match FG in pixel {df$pixelID[1]}"))
  }
  
  df <- df %>%
    mutate(spRank = row_number())
  
  thresholds_df <- df %>%
    group_by(aggs) %>%
    summarise(
      aggThresholds = dplyr::last(log.mass),
      .groups = "drop"
    ) %>%
    mutate(aggs = as.character(aggs))
  
  agg_pal <- make_agg_palette(nrow(thresholds_df))
  thresholds_df$aggs <- factor(thresholds_df$aggs, levels = thresholds_df$aggs)
  
  p <- ggplot(df, aes(x = log.mass, y = spRank)) +
    geom_point(aes(color = FG), size = 4.2) +
    scale_color_manual(name = "Functional Groups", values = fg_palette) +
    ggnewscale::new_scale_color() +
    geom_vline(
      data = thresholds_df,
      aes(xintercept = aggThresholds, color = aggs),
      linetype = "dashed",
      linewidth = 0.9
    ) +
    scale_color_manual(name = "Aggregations", values = agg_pal) +
    labs(
      title = glue("Discontinuity breaks for {gcm} - {yr}"),
      subtitle = glue("Pixel {df$pixelID[1]}"),
      x = "log10 body mass",
      y = "Ranked body mass"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold")
    )
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  file_out <- glue("{out_dir}/aggs_pixel_{df$pixelID[1]}_{gcm}_{yr}.png")
  ggsave(file_out, plot = p, width = 8, height = 5, units = "in", dpi = 300)
  
  invisible(NULL)
}

# Function to process one gcm-year file -----------------------------------
make_plot_discontinuities_new <- function(file, traitsFG, out_root) {
  
  nm <- basename(file)
  m <- stringr::str_match(nm, "^Aggregations_(.*)_(\\d{4}).*\\.qs$")
  
  gcm <- m[, 2]
  yr  <- m[, 3]
  
  if (is.na(gcm) || is.na(yr)) {
    warning(glue("Could not parse filename: {nm}"))
    return(NULL)
  }
  
  message(glue("Processing {gcm} - {yr}"))
  
  df_all <- qs::qread(file)
  
  if (is.null(df_all) || !is.data.frame(df_all) || nrow(df_all) == 0) {
    warning(glue("Empty or invalid file: {nm}"))
    return(NULL)
  }
  
  # If old files still use 'sp', rename to species
  if ("sp" %in% names(df_all) && !"species" %in% names(df_all)) {
    df_all <- df_all %>% rename(species = sp)
  }
  
  required_cols <- c("pixelID", "species", "log.mass", "aggs")
  if (!all(required_cols %in% names(df_all))) {
    warning(glue("Missing required columns in {nm}: {paste(setdiff(required_cols, names(df_all)), collapse=', ')}"))
    return(NULL)
  }
  
  pixel_list <- split(df_all, df_all$pixelID)
  
  out_dir <- glue("{out_root}/{gcm}/{yr}")
  
  future_walk(
    pixel_list,
    ~ plot_one_pixel(
      df_pixel = .x,
      gcm = gcm,
      yr = yr,
      traitsFG = traitsFG,
      out_dir = out_dir
    ),
    .options = furrr_options(seed = TRUE)
  )
  
  invisible(NULL)
}

# Run over all files ------------------------------------------------------
future_walk(
  files,
  ~ make_plot_discontinuities_new(
    file = .x,
    traitsFG = traitsFG,
    out_root = out_root
  ),
  .options = furrr_options(seed = TRUE)
)
