# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(
  crayon, data.table, dplyr, glue, progress, qs, stringr,
  tidyverse, vegan
)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999)

# Parameters --------------------------------------------------------------
fg_path <- "./tables/comm_matrices_filtered_4km/fgCommMatrices_filtered_4km.qs"
dir_samples <- "./tables/sampling"

yrs <- c("2011", "2031", "2091")
gcms <- c("CanESM2", "CCSM4", "INM-CM4")

# Read sampled region files ----------------------------------------------
regions_fls <- list.files(
  dir_samples,
  pattern = "10percent.*\\.qs$",
  full.names = TRUE
)

file_data <- lapply(regions_fls, qs::qread)
names(file_data) <- basename(regions_fls)


fg_all <- qs::qread(fg_path)



## FG columns --------------------------------------------------------------
fg_cols <- setdiff(names(fg_all), c("pixelID", "x", "y", "yr", "gcm"))


# Main NMDS function ------------------------------------------------------

run_nmds <- function(region_file, pix_sample, gcm_sel, yrs, fg_all, fg_cols) {
  
  if (is.null(pix_sample) || nrow(pix_sample) == 0) {
    cat("Skipping empty sample for:", basename(region_file), "\n")
    return(NULL)
  }
  
  region_name <- basename(region_file) |>
    str_remove("_10percent_sampling\\.qs$")
  
  cat("Processing region:", region_name, "| GCM:", gcm_sel, "\n")
  
  pix_sample <- pix_sample %>%
    mutate(pixelID = as.numeric(pixelID))
  
  fg_df_yrs <- fg_all %>%
    filter(
      pixelID %in% pix_sample$pixelID,
      gcm == gcm_sel,
      yr %in% yrs
    ) %>%
    left_join(
      pix_sample %>% select(pixelID, ECO2_NAM_1, ECO3_NAM_1, ECO4_NAM_1),
      by = "pixelID"
    )
  
  if (nrow(fg_df_yrs) < 5) {
    cat("Skipping region", region_name, "- insufficient rows\n")
    return(NULL)
  }
  
  fg_df_yrs <- fg_df_yrs %>%
    mutate(total_sum = rowSums(across(all_of(fg_cols)), na.rm = TRUE)) %>%
    filter(total_sum != 0)
  
  if (nrow(fg_df_yrs) < 5) {
    cat("Skipping region", region_name, "- insufficient non-empty rows\n")
    return(NULL)
  }
  
  fg_df_sp <- fg_df_yrs %>%
    select(all_of(fg_cols))
  
  colnames(fg_df_sp) <- gsub("-", "_", colnames(fg_df_sp))
  
  fg_df_sites <- fg_df_yrs %>%
    select(pixelID, x, y, ECO2_NAM_1, ECO3_NAM_1, ECO4_NAM_1,yr, gcm)
  
  fit <- adonis2(
    fg_df_sp ~ yr,
    data = fg_df_sites,
    permutations = 999,
    method = "bray"
  )
  
  distances_data <- vegdist(fg_df_sp, method = "bray")
  homo_test <- anova(betadisper(distances_data, fg_df_sites$yr))
  
  nmds_time <- system.time({
    nmds_out <- metaMDS(
      fg_df_sp,
      distance = "bray",
      k = 2,
      trymax = 100,
      trace = FALSE
    )
  })
  
  cat("Processing time:", nmds_time[[3]], "seconds\n")
  
  dir.create("./outputs/nmds_results_2026", recursive = TRUE, showWarnings = FALSE)
  dir.create("./outputs/bootstrapping_2026", recursive = TRUE, showWarnings = FALSE)
  dir.create("./outputs/permanova_2026", recursive = TRUE, showWarnings = FALSE)
  
  qs::qsave(
    nmds_out,
    glue("./outputs/nmds_results_2026/nmds_{region_name}_{gcm_sel}_10percent_allyrs.qs")
  )
  
  qs::qsave(
    fg_df_sites,
    glue("./outputs/nmds_results_2026/sites_{region_name}_{gcm_sel}_10percent_allyrs.qs")
  )
  
  qs::qsave(
    fit,
    glue("./outputs/bootstrapping_2026/{region_name}_{gcm_sel}_10percent_allyrs.qs")
  )
  
  qs::qsave(
    homo_test,
    glue("./outputs/permanova_2026/{region_name}_{gcm_sel}_10percent_allyrs.qs")
  )
  
  list(
    region = region_name,
    gcm = gcm_sel,
    n_rows = nrow(fg_df_yrs),
    n_fg = length(fg_cols),
    stress = nmds_out$stress,
    permanova = fit,
    dispersion = homo_test,
    nmds = nmds_out
  )
}

gcm_sel <- gcms[3]

nmds_results <- lapply(seq_along(file_data), function(i) {
  run_nmds(
    region_file = regions_fls[i],
    pix_sample = file_data[[i]],
    gcm_sel = gcm_sel,
    yrs = yrs,
    fg_all = fg_all,
    fg_cols = fg_cols
  )
})

names(nmds_results) <- basename(regions_fls)
