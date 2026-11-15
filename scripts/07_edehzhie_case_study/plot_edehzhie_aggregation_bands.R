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

agg_dir <- "./outputs/aggs_edehzhie/Aggregations_Edehzhie"
out_root <- "./figs/aggregations2026/edeh"

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

if (length(fg_levels) > length(fg_palette)) {
  fg_palette <- qualitative_hcl(length(fg_levels), palette = "Dynamic")
}

fg_palette <- setNames(fg_palette[seq_along(fg_levels)], fg_levels)

# Read aggregation files --------------------------------------------------
files <- list.files(
  agg_dir,
  pattern = "^Aggregations_.*\\.qs$",
  full.names = TRUE
)

# Plot function -----------------------------------------------------------
make_plot_aggregation_bands_edeh <- function(gcm, yr, files, traitsFG,
                                             fg_palette, out_root) {
  
  matched_files <- files %>%
    grep(gcm, ., value = TRUE) %>%
    grep(yr, ., value = TRUE)
  
  if (length(matched_files) == 0) {
    warning(glue("No files found for {gcm}-{yr}"))
    return(invisible(NULL))
  }
  
  invisible(lapply(seq_along(matched_files), function(i) {
    
    df <- qs::qread(matched_files[i]) %>%
      rename(sp = species) %>%
      left_join(traitsFG, by = c("sp" = "species")) %>%
      arrange(log.mass) %>%
      mutate(
        spRank = row_number(),
        aggs = factor(aggs)
      )
    
    agg_bounds <- df %>%
      group_by(aggs) %>%
      summarise(
        xmin = min(log.mass, na.rm = TRUE),
        xmax = max(log.mass, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(xmin) 
    
    p <- ggplot(df, aes(x = log.mass, y = spRank)) +
      
      # Shade aggregation intervals only; gaps remain unshaded
      geom_rect(
        data = agg_bounds,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = -Inf,
          ymax = Inf),
        inherit.aes = FALSE,
        fill = 'grey78',
        alpha = 0.45,
        show.legend = FALSE
      ) +
        geom_point(
        aes(color = FG),
        size = 1
      ) +
      
      scale_color_manual(
        name = "Functional groups",
        values = fg_palette
      ) +
      
      ggtitle(
        glue("Body-mass aggregations for {gcm} - {yr}"),
        subtitle = "Edehzhie"
      ) +
      
      labs(
        x = "log10 body mass",
        y = "Ranked body mass"
      ) +
      
      theme_bw() +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 16),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold")
      )
    
    save_dir <- file.path(out_root, "aggregation_bands_gaps_unshaded", gcm)
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    
    filename <- file.path(
      save_dir,
      paste0("2edehzhie_aggregation_bands_gaps_unshaded_", gcm, "_", yr, "_", i, ".png")
    )
    
    ggsave(
      filename = filename,
      plot = p,
      width = 10,
      height = 5,
      units = "in",
      dpi = 300
    )
  }))
  
  invisible(NULL)
}

# Run all combinations ----------------------------------------------------
combinations <- expand.grid(
  gcm = gcms,
  yr = yrs,
  stringsAsFactors = FALSE
)

future_pmap(
  combinations,
  function(gcm, yr) {
    make_plot_aggregation_bands_edeh(
      gcm = gcm,
      yr = yr,
      files = files,
      traitsFG = traitsFG,
      fg_palette = fg_palette,
      out_root = out_root
    )
  }
)

plan(sequential)
