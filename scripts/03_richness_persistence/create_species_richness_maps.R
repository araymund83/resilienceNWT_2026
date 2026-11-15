require(pacman)
pacman::p_load(
  dplyr, fs, glue, qs, terra, ggplot2, purrr, stringr, tidyr, sf
)

rm(list = ls())
options(scipen = 999)

# ------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------
in_dir  <- "./tables/community_matrices"
out_dir <- "./figs/richness_maps_multiscale"

if (!dir_exists(out_dir)) dir_create(out_dir, recurse = TRUE)

yrs    <- c("2011", "2031", "2091")
gcms   <- c("CanESM2", "CCSM4", "INM-CM4")
scales <- c("4km", "10km", "20km")


targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

limt <- sf::st_read("inputs/TaigaPlainsStudyArea/TaigaPlainsStudyAreaPoly.shp", quiet = TRUE) %>%
  sf::st_transform(crs = targetCRS)
ecodis <- sf::st_read("./inputs/EcoDistrict_NWT_gov/ecodistrictsBCR6NWT.shp", quiet = TRUE) %>%
  sf::st_transform(crs = targetCRS)

# ------------------------------------------------------------------
# Helper: read one community matrix
# ------------------------------------------------------------------
read_comm <- function(scale_sel, yr_sel, gcm_sel, in_dir) {
  
  f <- glue("{in_dir}/comm_{scale_sel}_{yr_sel}_{gcm_sel}.qs")
  
  if (!file.exists(f)) {
    stop("Missing community matrix: ", f)
  }
  
  qs::qread(f)
}

# ------------------------------------------------------------------
#Richness from community matrix
# ------------------------------------------------------------------
calc_richness <- function(mat) {
  
  mat <- as.data.frame(mat)
  browser()
  non_spp <- c("pixelID", "x", "y")
  spp_cols <- setdiff(names(mat), non_spp)
  
  rowSums(mat[, spp_cols, drop = FALSE], na.rm = TRUE)
}

# ------------------------------------------------------------------
# Get global limits by scale
# ------------------------------------------------------------------
get_rich_range <- function(scale_sel, yr_sel, gcm_sel) {
  
  mat <- read_comm(scale_sel, yr_sel, gcm_sel, in_dir)
  rich <- calc_richness(mat)
  
  range(rich, na.rm = TRUE)
}

range_tbl <- tidyr::crossing(
  scale_sel = scales,
  yr_sel    = yrs,
  gcm_sel   = gcms
) %>%
  mutate(
    rng = purrr::pmap(
      list(scale_sel, yr_sel, gcm_sel),
      get_rich_range
    )
  )

global_lims <- range_tbl %>%
  mutate(
    min_r = purrr::map_dbl(rng, 1),
    max_r = purrr::map_dbl(rng, 2)
  ) %>%
  group_by(scale_sel) %>%
  summarise(
    global_min = min(min_r, na.rm = TRUE),
    global_max = max(max_r, na.rm = TRUE),
    .groups = "drop"
  )

print(global_lims)

# ------------------------------------------------------------------
# Make one richness map
# ------------------------------------------------------------------
make_richness_map <- function(scale_sel, yr_sel, gcm_sel) {
  
  message("--------------------------------------------------")
  message(glue("Mapping richness: {scale_sel} | {yr_sel} | {gcm_sel}"))
  
  mat <- read_comm(scale_sel, yr_sel, gcm_sel, in_dir)
  mat <- as.data.frame(mat)
  
  rich <- calc_richness(mat)
  
  df_plot <- tibble(
    pixelID  = mat$pixelID,
    x        = mat$x,
    y        = mat$y,
    richness = rich
  ) %>%
    filter(!is.na(x), !is.na(y))
  
  # save richness table too
  qs::qsave(
    df_plot,
    glue("{out_dir}/richness_{scale_sel}_{yr_sel}_{gcm_sel}.qs")
  )
  
  # choose scale-specific legend limits
  lims <- global_lims %>%
    filter(scale_sel == !!scale_sel)

  # color palette
  color_pal <- c(
    "#f7f7f7", "#96685C", "#FFBE32",
    "#71928D", "#36585D", "#7b3294"
  )
  
  p <- ggplot()+
    geom_tile(data = df_plot, aes(x = x, y = y, fill = richness)) +
    scale_fill_gradientn(
      colours = color_pal,
      limits = c(lims$global_min, lims$global_max),
      oob = scales::squish,
      na.value = "white",
      breaks = pretty(c(lims$global_min, lims$global_max), n = 5),
      name = "Richness"
    ) +
    theme_bw() +
    geom_sf(data = limt, fill = NA, color = "#36454F", linewidth = 0.3) +
    geom_sf(data = ecodis, fill = NA, color = "#36454F", linewidth = 0.2) +
    coord_sf(crs = targetCRS, datum = NA) +
    labs(
      title = glue("Species richness ({scale_sel})"),
      subtitle = glue("{gcm_sel} — {yr_sel}"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12)
    ) + geom_sf()
  
  
  ggsave(
    filename = glue("{out_dir}/richness_{scale_sel}_{yr_sel}_{gcm_sel}.png"),
    plot = p,
    width = 6,
    height = 8,
    dpi = 300
  )
  
  invisible(TRUE)
}

# ------------------------------------------------------------------
# Run all combinations
# ------------------------------------------------------------------
combos <- tidyr::crossing(
  scale_sel = scales,
  yr_sel    = yrs,
  gcm_sel   = gcms
)

purrr::pwalk(
  combos,
  make_richness_map
)

message("Done.")
#__________________________________________________________
##put figures together ____________________________________
#_________________________________________________________
require(pacman)
pacman::p_load(magick, fs, glue)

yrs    <- c("2011", "2031", "2091")
gcms   <- c("CanESM2", "CCSM4", "INM-CM4")
scales <- c("4km", "10km", "20km")

fig_dir  <- "./figs/richness_maps_multiscale"
comb_dir <- "./figs/richness_maps_multiscale/combined"

if (!dir_exists(comb_dir)) dir_create(comb_dir, recurse = TRUE)

# ----------------------------------------------------------------------
# Combine maps by scale and GCM
# ----------------------------------------------------------------------
for (scale_sel in scales) {
  
  for (gcm in gcms) {
    
    paths <- glue("{fig_dir}/richness_{scale_sel}_{yrs}_{gcm}.png")
    
    if (any(!file.exists(paths))) {
      warning(
        "Missing one or more files for ", scale_sel, " / ", gcm, ":\n",
        paste(paths[!file.exists(paths)], collapse = "\n")
      )
      next
    }
    
    imgs <- lapply(paths, image_read)
    
    # make sure all same height
    max_h <- max(sapply(imgs, \(im) image_info(im)$height))
    imgs  <- lapply(imgs, \(im) image_scale(im, geometry = paste0("x", max_h)))
    
    combined <- image_append(do.call(c, imgs), stack = FALSE)
    
    image_write(
      combined,
      path = glue("{comb_dir}/richness_{scale_sel}_{gcm}_2011_2031_2091.png")
    )
    
    message("Saved combined figure: ",
            glue("{comb_dir}/richness_{scale_sel}_{gcm}_2011_2031_2091.png"))
  }
}

message("Done.")
