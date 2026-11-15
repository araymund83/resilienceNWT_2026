#Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(dplyr,glue, purrr, qs, RColorBrewer, stringr, sf, terra, tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Load data ---------------------------------------------------------------
#load files for extracting TBI
in_dir <- './outputs/tbi_72spp_multiscale'

files <- list.files(in_dir, pattern = 'table', full.names = TRUE)
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[2]
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
year_range <- glue::glue("{yr1}-{yr2}")
#gcm <- gcms[1]

#load shapefiles and CRS for maps. 
targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#limt <- sf::st_read("inputs/TaigaPlainsStudyArea/TaigaPlainsStudyAreaPoly.shp")
limt <- sf::st_read('inputs/NT1_BCR6/NT1_BCR6_poly.shp')
ecrg <- sf::st_read('inputs/ecoregions/NWT_ecoregions_5Regions.shp')

# Extract by mask for the ecoregions 
plot(st_geometry(ecrg))
limt <- sf::st_transform(x = limt, crs = targetCRS)
ecrg <- sf::st_transform(x = ecrg, crs = targetCRS)
plot(st_geometry(ecrg))
plot(st_geometry(limt))

# ---- Load TBI results ----
tbi_all <- purrr::map_dfr(files, \(f) {
  qs::qread(f) %>%
    mutate(source_file = basename(f))
})


#function to extract the TBI, signficance and if the change was + or - element for each gcm
# ---- Combine all GCMs into one tibble ----
tbi_all <- tbi_all %>%
  mutate(
    pValue2 = case_when(
      !is.na(p_TBI) & p_TBI < 0.05 ~ "significant",
      !is.na(p_TBI) & p_TBI >= 0.05 ~ "no_sig",
      TRUE ~ NA_character_
    ),
    change = stringr::str_trim(change),
    change2 = case_when(
      change == "+" ~ "positive",
      change == "-" ~ "negative",
      change == "0" ~ "no_change",
      TRUE ~ NA_character_
    ),
    change3 = case_when(
      change == "-" ~ -TBI,
      change == "+" ~ TBI,
      change == "0" ~ 0,
      TRUE ~ NA_real_
    )
  )
# ---- Save combined results ----
output_dir <- './tables/tbi72spp_2026_multiscale'

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
qs::qsave(tbi_all, file = glue('{output_dir}/tbi72sp_2026_allgcms_allyrs_allscales.qs'))

tbi_df <- qs::qread(glue('{output_dir}/tbi72sp_2026_allgcms_allyrs_allscales.qs'))

limits_by_scale <- tbi_df %>%
   group_by(scale) %>%
  summarise(
    min_abs = min(abs(change3), na.rm = TRUE),
    max_abs = max(abs(change3), na.rm = TRUE)
  )

make_TBI_72SppBin_Ras <- function(tbi_df, yr1, yr2, var, scale) {
  max_abs_scale <- limits_by_scale %>%
    filter(scale == !!scale) %>%
    pull(max_abs)
#browser()
  rasDF <- tbi_df %>%
    filter(yr1 == !!yr1, yr2 == !!yr2, scale == !!scale, valid_tbi == TRUE) %>%
    group_by(gcm) %>% 
    mutate(
      yr_label = glue("{yr1}-{yr2}")) %>% 
    rename(
      tbi = TBI,
      p_val = p_TBI,
      Longitude = x,
      Latitude = y,
    ) %>% 
    ungroup() 
  # Save output
  out <- glue('./tables/tbi72sppBin2026')
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  qs::qsave(rasDF, glue('{out}/TBIpval72sppBin_allgcms_{yr1}_{yr2}_{scale}.qs'))
  
  # Plot
  message(glue("Creating {var} map for {yr1}-{yr2}..."))
  
  if (var == 'pValue2') {
    color_2vals <- c('#b3b3b3', '#4d004b')
    ggVar <- ggplot() +
      geom_tile(data = rasDF, aes(x = Longitude, y = Latitude, fill = pValue2)) +
      scale_fill_manual(
        values = c(
          "no_sig" = "white",
          "significant" = "#4d004b"
        ),
        labels = c(
          "no_sig" = "No significant",
          "significant" = "Significant"
        ),
        na.value = "transparent"
      ) +
      geom_sf(data = limt, fill = NA, col = '#36454F') +
      geom_sf(data = ecrg, fill = NA, col = '#36454F') +
      ggtitle('pValue TBI for species', subtitle = glue('{yr1}-{yr2} | {scale}')) +
      theme_bw() +
      theme(panel.background = element_rect(fill = "white", colour = NA),
            legend.position = 'bottom',
            legend.key.width = unit(2, 'line'),
            plot.title = element_text(size = 16, face = 'bold'),
            plot.subtitle = element_text(size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 12, face = 'bold'),
            strip.text = element_text(size = 14)) +
      labs(x = 'Longitude', y = 'Latitude', fill = 'pValue') +
      coord_sf() +
      facet_wrap(~ factor(gcm, c('INM-CM4', 'CCSM4', 'CanESM2')))
    
    out_plot <- './figs/tbi_72spp_multiscale_maps'
    if (!dir.exists(out_plot)) dir.create(out_plot, recursive = TRUE)
    ggsave(plot = ggVar, filename = glue('{out_plot}/{var}_72spp_{yr1}_{yr2}_{scale}.png'),
           units = 'in', width = 10, height = 7, dpi = 300)
  } else if (var == 'tbi') {

   
     ggTBI <- ggplot() +
      geom_tile(data = rasDF, aes(x = Longitude, y = Latitude, fill = change3)) +
      scale_fill_gradient2(
        low = "#a6611a",
        mid = "white",
        high = "#018571",
        midpoint = 0,
        limits = c(-max_abs_scale, max_abs_scale),
        na.value = "transparent"
      ) +
      geom_sf(data = limt, fill = NA, col = '#36454F') +
      geom_sf(data = ecrg, fill = NA, col = '#36454F') +
      ggtitle('TBI index for species', subtitle = glue('{yr1}-{yr2} | {scale}')) +
      theme_bw() +
      theme(panel.background = element_rect(fill = "white", colour = NA),
            legend.position = 'bottom',
            legend.key.width = unit(2, 'line'),
            plot.title = element_text(size = 16, face = 'bold'),
            plot.subtitle = element_text(size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 12, face = 'bold'),
            strip.text = element_text(size = 14)) +
      labs(x = 'Longitude', y = 'Latitude', fill = 'TBI') +
      coord_sf(lims_method = 'geometry_bbox') +
      facet_wrap(~ factor(gcm, c('INM-CM4', 'CCSM4', 'CanESM2')))
    
    out_plot <-  './figs/tbi_72spp_multiscale_maps'
    if (!dir.exists(out_plot)) dir.create(out_plot, recursive = TRUE)
    ggsave(plot = ggTBI, filename = glue('{out_plot}/{var}_72spp_{yr1}_{yr2}_{scale}.png'),
           units = 'in', width = 10, height = 7, dpi = 300)
  }
}

# ---- Run plotting function ----
# make_TBI_72SppBin_Ras(
#   tbi_df = tbi_df,
#   yr1 = yr1,
#   yr2 = yr2,
#   var = 'tbi',##var options are 'tbi' and 'pValue2'
#   scale = scale
# )

scales <- c("4km", "10km", "20km")

out_list <- lapply(scales, function(sc) {
  make_TBI_72SppBin_Ras(
    tbi_df = tbi_all,
    yr1 = "2011",
    yr2 = "2031",
    var = "pValue2",##var options are 'tbi' and 'pValue2'
    scale = sc
  )
})

