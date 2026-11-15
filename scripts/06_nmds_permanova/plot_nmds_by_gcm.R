library(pacman)
p_load(
  qs, glue, dplyr, ggplot2, ggrepel, ggpubr, cowplot, vegan
)

get_region_from_file <- function(file, gcm) {
  sub(
    paste0("^nmds_(.*)_", gcm, "_10percent_allyrs\\.qs$"),
    "\\1",
    basename(file)
  )
}

path <- "./outputs/nmds_results_2025"
out_dir <- "./figs/nmds_sampling_2026"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

region_order_files <- c(
  "Mid_Boreal",
  "High_Boreal",
  "Low_Subarctic",
  "High_Subarctic",
  "Low_Arctic_north"
)

get_region_from_file <- function(file, gcm) {
  sub(
    paste0("^nmds_(.*)_", gcm, "_10percent_allyrs\\.qs$"),
    "\\1",
    basename(file)
  )
}

clean_region <- function(region) {
  dplyr::case_when(
    region == "Mid_Boreal" ~ "Mid-Boreal",
    region == "High_Boreal" ~ "High Boreal",
    region == "Low_Subarctic" ~ "Low Subarctic",
    region == "High_Subarctic" ~ "High Subarctic",
    region == "Low_Arctic_north" ~ "Low Arctic north",
    TRUE ~ gsub("_", " ", region)
  )
}

get_shared_limits <- function(path, gcm, region_order_files, arrow_scale = 0.45, top_n_species = 8) {
  
  all_x <- c()
  all_y <- c()
  
  for (region in region_order_files) {
    
    nmds_file <- glue("{path}/nmds_{region}_{gcm}_10percent_allyrs.qs")
    if (!file.exists(nmds_file)) next
    
    nmds_out <- qs::qread(nmds_file)
    
    site_scores <- as.data.frame(nmds_out$points)
    names(site_scores) <- sub("^MDS1$", "NMDS1", names(site_scores))
    names(site_scores) <- sub("^MDS2$", "NMDS2", names(site_scores))
    
    species_scores <- as.data.frame(scores(nmds_out, "species"))
    species_scores$species <- rownames(species_scores)
    
    species_scores <- species_scores %>%
      mutate(vec_len = sqrt(NMDS1^2 + NMDS2^2)) %>%
      arrange(desc(vec_len)) %>%
      slice_head(n = top_n_species) %>%
      mutate(
        NMDS1 = NMDS1 * arrow_scale,
        NMDS2 = NMDS2 * arrow_scale
      )
    
    all_x <- c(all_x, site_scores$NMDS1, species_scores$NMDS1, 0)
    all_y <- c(all_y, site_scores$NMDS2, species_scores$NMDS2, 0)
  }
  
  x_lim <- max(abs(all_x), na.rm = TRUE) * 1.08
  y_lim <- max(abs(all_y), na.rm = TRUE) * 1.08
  
  list(x_lim = x_lim, y_lim = y_lim)
}

make_single_nmds_plot <- function(path, region, gcm, x_lim, y_lim,
                                  arrow_scale = 0.45,
                                  top_n_species = 8,
                                  show_legend = FALSE) {
  
  region_title <- clean_region(region)
  
  sites_file <- glue("{path}/sites_{region}_{gcm}_10percent_allyrs.qs")
  nmds_file  <- glue("{path}/nmds_{region}_{gcm}_10percent_allyrs.qs")
  
  sites <- qs::qread(sites_file)
  nmds_out <- qs::qread(nmds_file)
  
  NMDSoutpoints <- as.data.frame(nmds_out$points)
  names(NMDSoutpoints) <- sub("^MDS1$", "NMDS1", names(NMDSoutpoints))
  names(NMDSoutpoints) <- sub("^MDS2$", "NMDS2", names(NMDSoutpoints))
  
  NMDSoutpoints$yr <- factor(sites$yr, levels = c("2011", "2031", "2091"))
  NMDSoutpoints$pixelID <- sites$pixelID
  
  species_scores <- as.data.frame(scores(nmds_out, "species"))
  species_scores$species <- rownames(species_scores)
  
  species_scores <- species_scores %>%
    mutate(vec_len = sqrt(NMDS1^2 + NMDS2^2)) %>%
    arrange(desc(vec_len)) %>%
    slice_head(n = top_n_species) %>%
    mutate(
      NMDS1 = NMDS1 * arrow_scale,
      NMDS2 = NMDS2 * arrow_scale
    )
  
  stress <- format(round(nmds_out$stress, 3), nsmall = 3)
  
  hull_df <- NMDSoutpoints %>%
    group_by(yr) %>%
    filter(n() >= 3) %>%
    slice(chull(NMDS1, NMDS2)) %>%
    ungroup()
  
  p <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
    
    geom_polygon(
      data = hull_df,
      aes(x = NMDS1, y = NMDS2, fill = yr, group = yr),
      alpha = 0.16,
      colour = NA
    ) +
    
    geom_point(
      data = NMDSoutpoints,
      aes(x = NMDS1, y = NMDS2, colour = yr),
      size = 0.5, alpha = 0.1
    ) +
    
    geom_segment(
      data = species_scores,
      aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
      arrow = arrow(length = grid::unit(0.08, "inches")),
      colour = "grey25",
      linewidth = 0.45
    ) +
    
    geom_text_repel(
      data = species_scores,
      aes(x = NMDS1, y = NMDS2, label = species),
      size = 3,
      box.padding = 0.25,
      point.padding = 0.2,
      max.overlaps = 50
    ) +
    
    scale_colour_manual(
      values = c("2011" = "#FFB00D", "2031" = "#0072B2", "2091" = "#AA3377"),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("2011" = "#FFB00D", "2031" = "#0072B2", "2091" = "#AA3377"),
      name = NULL
    ) +
    
    coord_fixed(
      xlim = c(-x_lim, x_lim),
      ylim = c(-y_lim, y_lim),
      expand = FALSE
    ) +
    
    labs(
      title = glue("NMDS for {region_title}"),
      subtitle = glue("{gcm} | Stress = {stress}"),
      x = "NMDS1",
      y = "NMDS2"
    ) +
    
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "#4D4D4D", fill = NA, linewidth = 0.6),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = 'none',
      legend.key = element_blank(),
      legend.text = element_text(size = 10),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5)
    )
  
  p
}

make_legend_panel <- function() {
  
  legend_df <- data.frame(
    x = c(1, 1, 1),
    y = c(0.8, 0.7, 0.6),
    lab = c("2011", "2031", "2091")
  )
  mid <- mean(legend_df$y)
  
  cols <- c("#FFB00D", "#0072B2", "#AA3377")
  

  ggplot() +
    geom_point(
      data = legend_df,
      aes(x, y),
      shape = 22,
      size = 8,
      fill = cols,
      alpha = 0.18,
      colour = cols,
      stroke = 0.8
    ) +
    geom_point(
      data = legend_df,
      aes(x, y),
      colour = cols,
      size = 3.2,
      alpha = 0.65
    ) +
    geom_text(
      data = legend_df,
      aes(x = 1.38, y = y, label = lab),
      hjust = 0,
      size = 4.8
    ) +
    coord_cartesian(
      xlim = c(0.8, 2),
      ylim = c(mid - 0.35, mid + 0.35),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      plot.margin = margin(5, 5, 5, 5)
    )
}

make_combined_nmds_figure <- function(path, gcm, out_dir,
                                      region_order_files,
                                      arrow_scale = 0.45,
                                      top_n_species = 8) {
  
  lims <- get_shared_limits(
    path = path,
    gcm = gcm,
    region_order_files = region_order_files,
    arrow_scale = arrow_scale,
    top_n_species = top_n_species
  )
  
  plot_list <- lapply(
    region_order_files,
    function(region) {
      make_single_nmds_plot(
        path = path,
        region = region,
        gcm = gcm,
        x_lim = lims$x_lim,
        y_lim = lims$y_lim,
        arrow_scale = arrow_scale,
        top_n_species = top_n_species,
        show_legend = FALSE
      )
    }
  )
  
  
  legend_panel <- make_legend_panel()

  print(legend_panel)
  
  names(plot_list) <- region_order_files
  
  
  combined <- ggpubr::ggarrange(
    plotlist = list(
      plot_list[["Mid_Boreal"]],
      plot_list[["High_Boreal"]],
      plot_list[["Low_Subarctic"]],
      plot_list[["High_Subarctic"]],
      plot_list[["Low_Arctic_north"]],
      legend_panel
    ),
    common.legend = TRUE,
    legend = 'none',
    ncol = 3,
    nrow = 2,
    align = "hv"
  )
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ggsave(
    filename = glue("{out_dir}/NMDS_combined_{gcm}_10percent_allyrs.png"),
    plot = combined,
    width = 16,
    height = 10,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  combined
}



fig_CanESM2 <- make_combined_nmds_figure(
  path = path,
  gcm = "CanESM2",
  out_dir = out_dir,
  region_order_files = region_order_files
)

fig_CCSM4 <- make_combined_nmds_figure(
  path = path,
  gcm = "CCSM4",
  out_dir = out_dir,
  region_order_files = region_order_files
)

fig_INMCM4 <- make_combined_nmds_figure(
  path = "./outputs/nmds_results_2026",
  gcm = "INM-CM4",
  out_dir = out_dir,
  region_order_files = region_order_files
)
