# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(fs, glue, qs, stringr, tidyverse)

rm(list = ls())
options(scipen = 999)
# Load data ---------------------------------------------------------------
in_dir  <- "./tables/comm_matrices_allgcms_by_scale"
out_dir <- "./tables/fg_pa_by_gcm_scale"
if (!dir_exists(out_dir)) dir_create(out_dir, recurse = TRUE)

# Inputs ------------------------------------------------------------------
files <- list.files(in_dir, pattern = "\\.qs$", full.names = TRUE)
gcms  <- c("CanESM2", "CCSM4", "INM-CM4")

fg_lookup <- qs::qread('./inputs/traitsMassFG_2025.qs')

make_fg_matrix_one_scale <- function(f, fg_lookup, out_dir) {
  
  fname <- basename(f)
  scale <- str_extract(fname, "4km|10km|20km")
  
  message("Reading: ", fname)
  dat <- qs::qread(f)
  
  gcm_list <- split(dat, dat$gcm)
  
  for (gcm in names(gcm_list)) {
    
    message("Processing: ", scale, " | ", gcm)
    d <- gcm_list[[gcm]]
    
    id_cols  <- c("pixelID", "x", "y", "yr", "gcm")
    spp_cols <- setdiff(names(dat), id_cols)
   
    d_long <- d %>%
      pivot_longer(
        cols = all_of(spp_cols),
        names_to = "species",
        values_to = "pa"
      ) %>%
      mutate(pa = if_else(is.na(pa), 0L, as.integer(pa > 0))) %>%
      left_join(fg_lookup, by = "species") %>%
      filter(!is.na(fg))
 
    fg_wide <- d_long %>%
      group_by(pixelID, x, y, yr, gcm, fg) %>%
      summarise(pa_fg = as.integer(any(pa == 1L)), .groups = "drop") %>%
      pivot_wider(
        names_from = fg,
        values_from = pa_fg,
        values_fill = 0
      ) %>%
      arrange(yr, pixelID)
    
    out_file <- glue("{out_dir}/fgPA_{gcm}_allyrs_{scale}.qs")
    qs::qsave(fg_wide, out_file)
    
    message("Saved: ", basename(out_file))
  }
}

for (f in files) {
  make_fg_matrix_one_scale(f, fg_lookup, out_dir)
}
