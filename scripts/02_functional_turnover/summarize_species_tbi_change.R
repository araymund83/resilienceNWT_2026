# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(
  dplyr, glue, purrr, qs, stringr, sf, terra,
  tidyverse, ggplot2, magick
)

# Clean environment -------------------------------------------------------
gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Options -----------------------------------------------------------------
na_col <- "#D8E2E7"   # same NA colour used in mean-delta maps
detect_NA <- TRUE

# Load TBI tables ---------------------------------------------------------
data_dir <- "./tables/tbi72spp_2026_multiscale"

tbi_df <- qs::qread(
  glue("{data_dir}/tbi72spp_2026_allgcms_allyrs_allscales.qs")
)

tbi_mean_input <- tbi_df %>%
  filter(valid_tbi) %>%
  mutate(
    yrs = glue("{yr1}-{yr2}"),
    facet_label = yrs
  )

mean_tbiChange <- tbi_mean_input %>%
  group_by(pixelID, x, y, scale, yrs, facet_label) %>%
  summarise(
    mean_tbiChange = ifelse(
      all(is.na(change3)),
      NA_real_,
      mean(change3, na.rm = TRUE)
    ),
    .groups = "drop"
  )

tbi_sig_count <- tbi_mean_input %>%
  mutate(
    significant = if_else(!is.na(p_TBI) & p_TBI < 0.05, 1L, 0L)
  ) %>%
  group_by(pixelID, x, y, scale, yrs, facet_label) %>%
  summarise(
    significant_count = ifelse(
      all(is.na(p_TBI)),
      NA_integer_,
      sum(significant, na.rm = TRUE)
    ),
    .groups = "drop"
  )

scales <- unique(tbi_mean_input$scale)

# Load shapefiles and CRS -------------------------------------------------
targetCRS <- paste(
  "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
  "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
)

limt <- sf::st_read("inputs/NT1_BCR6/NT1_BCR6_poly.shp", quiet = TRUE)
ecrg <- sf::st_read("inputs/ecoregions/NWT_ecoregions_5Regions.shp", quiet = TRUE)

limt <- sf::st_transform(limt, crs = targetCRS)
ecrg <- sf::st_transform(ecrg, crs = targetCRS)

# Limits by scale ---------------------------------------------------------
limits_by_scale <- tbi_df %>%
  filter(valid_tbi) %>%
  group_by(scale) %>%
  summarise(
    max_abs = max(abs(change3), na.rm = TRUE),
    .groups = "drop"
  )

# Palettes ----------------------------------------------------------------
pal_div <- c(
  "#a6611a",
  "#dfc27d",
  "white",
  "#80cdc1",
  "#018571"
)

sig_pal <- c(
  "0" = "white",
  "1" = "#BFD3E6",
  "2" = "#88419D",
  "3" = "#4D004B"
)

# Create complete grid for each scale -------------------------------------
make_grid_for_scale <- function(sc) {
  
  cell_size <- as.numeric(str_remove(sc, "km")) * 1000
  
  st_make_grid(
    limt,
    cellsize = cell_size,
    square = TRUE
  ) %>%
    st_as_sf() %>%
    st_filter(st_union(limt), .predicate = st_intersects) %>%
    mutate(
      grid_x = as.numeric(st_coordinates(st_centroid(.))[, 1]),
      grid_y = as.numeric(st_coordinates(st_centroid(.))[, 2]),
      x_join = round(grid_x / cell_size) * cell_size,
      y_join = round(grid_y / cell_size) * cell_size
    ) %>%
    st_drop_geometry()
}

grid_list <- lapply(scales, make_grid_for_scale)
names(grid_list) <- scales

# Prepare mean TBI data with explicit NA cells ----------------------------
prepare_mean_grid <- function(sc) {
  
  cell_size <- as.numeric(str_remove(sc, "km")) * 1000
  
  mean_sc_raw <- mean_tbiChange %>%
    filter(scale == sc) %>%
    mutate(
      x_join = round(x / cell_size) * cell_size,
      y_join = round(y / cell_size) * cell_size
    )
  
  mean_sc <- grid_list[[sc]] %>%
    crossing(facet_label = unique(mean_sc_raw$facet_label)) %>%
    left_join(
      mean_sc_raw,
      by = c("x_join", "y_join", "facet_label")
    )
  
  if (detect_NA) {
    message("--------------------------------------------------")
    message("Mean TBI scale: ", sc)
    message("Total cells plotted: ", nrow(mean_sc))
    message("NA mean_tbiChange cells: ", sum(is.na(mean_sc$mean_tbiChange)))
    message("Cells with no joined scale: ", sum(is.na(mean_sc$scale)))
  }
  
  mean_sc
}

# Prepare significance data with explicit NA cells ------------------------
prepare_sig_grid <- function(sc) {
  
  cell_size <- as.numeric(str_remove(sc, "km")) * 1000
  
  sig_sc_raw <- tbi_sig_count %>%
    filter(scale == sc) %>%
    mutate(
      x_join = round(x / cell_size) * cell_size,
      y_join = round(y / cell_size) * cell_size
    )
  
  sig_sc <- grid_list[[sc]] %>%
    crossing(facet_label = unique(sig_sc_raw$facet_label)) %>%
    left_join(
      sig_sc_raw,
      by = c("x_join", "y_join", "facet_label")
    )
  
  if (detect_NA) {
    message("--------------------------------------------------")
    message("TBI significance scale: ", sc)
    message("Total cells plotted: ", nrow(sig_sc))
    message("NA significant_count cells: ", sum(is.na(sig_sc$significant_count)))
    message("Cells with no joined scale: ", sum(is.na(sig_sc$scale)))
  }
  
  sig_sc
}

# Plot and save figures ---------------------------------------------------
invisible(lapply(scales, function(sc) {
  
  message("Processing scale: ", sc)
  
  mean_sc <- prepare_mean_grid(sc)
  sig_sc  <- prepare_sig_grid(sc)
  
  max_abs_scale <- limits_by_scale %>%
    filter(scale == sc) %>%
    pull(max_abs)
  
  tbi_breaks <- c(
    -max_abs_scale,
    -max_abs_scale / 2,
    0,
    max_abs_scale / 2,
    max_abs_scale
  )
  
  out <- glue("./figs/TBI_72spp_changemean2026/{sc}")
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  
  # --------------------------------------------------------------
  # Mean signed TBI
  # --------------------------------------------------------------
  ggTBI <- ggplot() +
    geom_tile(
      data = mean_sc,
      aes(x = grid_x, y = grid_y, fill = mean_tbiChange)
    ) +
    scale_fill_gradientn(
      colours = pal_div,
      limits = c(-max_abs_scale, max_abs_scale),
      values = scales::rescale(tbi_breaks),
      breaks = tbi_breaks,
      labels = round(tbi_breaks, 1),
      na.value = na_col,
      oob = scales::squish
    ) +
    geom_sf(data = limt, fill = NA, col = "#36454F") +
    geom_sf(data = ecrg, fill = NA, col = "#36454F") +
    ggtitle(
      "Mean TBI for species") +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(2, "line"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 14)
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "Mean TBI"
    ) +
    coord_sf(lims_method = "geometry_bbox") +
    facet_wrap(~ facet_label) 
  
  # ggsave(
  #   plot = ggTBI,
  #   filename = glue("{out}/TBI_mean_Change_72spp_{sc}.png"),
  #   units = "in",
  #   width = 4,
  #   height = 6,
  #   dpi = 300
  # )
  # 
  # --------------------------------------------------------------
  # Significant GCM count
  # --------------------------------------------------------------
  ggSig <- ggplot() +
    geom_tile(
      data = sig_sc,
      aes(
        x = grid_x,
        y = grid_y,
        fill = factor(
          significant_count,
          levels = c(0, 1, 2, 3)
        )
      )
    ) +
    scale_fill_manual(
      values = sig_pal,
      name = "Significant GCMs",
      labels = c("0", "1", "2", "3"),
      na.value = na_col,
      drop = FALSE
    ) +
    geom_sf(data = limt, fill = NA, color = "#36454F") +
    geom_sf(data = ecrg, fill = NA, color = "#36454F") +
    ggtitle(
      "Number of GCMs with significant TBI",
      subtitle = "(p < 0.05)"
    ) +
    facet_wrap(~ facet_label) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14),
      legend.key = element_rect(
        colour = "grey50",
        linewidth = 0.3),
      legend.key.width = unit(2, "line"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "Number of GCMs"
    ) +
    coord_sf(lims_method = "geometry_bbox") 
  
  if (sc == "4km") {
    p_tbi_mean <<- ggTBI
    p_tbi_sig  <<- ggSig
  }
  
  # --------------------------------------------------------------
  # Save figures
  # --------------------------------------------------------------
  
  ggsave(
    plot = ggTBI,
    filename = glue("{out}/TBI_mean_Change_72spp_{sc}.png"),
    units = "in",
    width = 6,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    filename = glue("{out}/72sppTBI_Significance_Count_Bin_{sc}.png"),
    plot = ggSig,
    width = 6,
    height = 5,
    dpi = 300
  )
}))

# Combine both figures ----------------------------------------------------
invisible(lapply(scales, function(sc) {
  
  out <- glue("./figs/TBI_72spp_changemean2026/{sc}")
  
  path1 <- glue("{out}/TBI_mean_Change_72spp_{sc}.png")
  path2 <- glue("{out}/72sppTBI_Significance_Count_Bin_{sc}.png")
  
  img1 <- magick::image_read(path1)
  img2 <- magick::image_read(path2)
  
  img1_annot <- magick::image_annotate(
    img1,
    "a)",
    size = 60,
    gravity = "northwest",
    location = "+5+5",
    color = "black",
    weight = 700
  )
  
  img2_annot <- magick::image_annotate(
    img2,
    "b)",
    size = 60,
    gravity = "northwest",
    location = "+5+5",
    color = "black",
    weight = 700
  )
  
  combined <- magick::image_append(
    c(img1_annot, img2_annot),
    stack = TRUE
  )
  
  magick::image_write(
    combined,
    path = glue("{out}/TBI_72spp_mean_and_significance_{sc}.png"),
    density = "300x300"
  )
}))


library(png)
library(grid)
library(ggpubr)
library(ggplotify)

# ---------------------------------------------------
# Read images as grobs
# ---------------------------------------------------

p_tbi_mean <- ggplotify::as.ggplot(
  grid::rasterGrob(
    png::readPNG("./figs/TBI_72spp_changemean2026/4km/TBI_mean_Change_72spp_4km.png"),
    interpolate = TRUE
  )
)

p_tbi_sig <- ggplotify::as.ggplot(
  grid::rasterGrob(
    png::readPNG("./figs/TBI_72spp_changemean2026/4km/72sppTBI_Significance_Count_Bin_4km.png"),
    interpolate = TRUE
  )
)

p_richness <- ggplotify::as.ggplot(
  grid::rasterGrob(
    png::readPNG("./figs/richness_maps_multiscale/mean_sppRichness/mean_species_richness_4km_2011.png"),
    interpolate = TRUE
  )
)

# ---------------------------------------------------
# Arrange TBI panels in left column
# ---------------------------------------------------

left_col <- ggarrange(
  p_tbi_mean,
  p_tbi_sig,
  ncol = 1,
  labels = c("a)", "b)"),
  font.label = list(size = 18, face = "bold"),
  hjust = -0.5,
  vjust = 1.5
)


# ---------------------------------------------------
# Arrange richness in right column
# ---------------------------------------------------

right_col <- ggarrange(
  p_richness,
  ncol = 1,
  labels = "c)",
  font.label = list(size = 18, face = "bold"),
  hjust = -0.5,
  vjust = 1.5
)

# ---------------------------------------------------
# Combine columns
# ---------------------------------------------------

combined <- ggpubr::ggarrange(
  left_col,
  right_col,
  ncol = 2,
  widths = c(1.25, 1)
)

# ---------------------------------------------------
# Save
# ---------------------------------------------------

ggpubr::ggexport(
  combined,
  filename = "./figs/TBI_significance_richness_2column_ggpubr_4km.png",
  width = 4200,
  height = 3000,
  res = 300
)

combined
