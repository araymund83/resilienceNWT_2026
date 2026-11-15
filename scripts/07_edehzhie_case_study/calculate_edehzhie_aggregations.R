# ------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------
library(pacman)

p_load(
  dplyr, tidyr, purrr, glue, qs, crayon, readr
)

# ------------------------------------------------------------------
# Source functions
# ------------------------------------------------------------------
source("./R/functions_discontinuities.R")

pwrtab <- read.csv("./R/GRI.constant.power.table.csv")

# ------------------------------------------------------------------
# Clean GCM values
# ------------------------------------------------------------------
comm_matrix_ede <- qs::qread( "./tables/community_matrices/comm_matrix_Edehzhie.qs")

comm_matrix_ede <- comm_matrix_ede %>%
  mutate(
    gcm = tidyr::replace_na(gcm, "INM-CM4"),
    yr  = as.character(yr)
  )

# ------------------------------------------------------------------
# Get species columns
# ------------------------------------------------------------------
species_cols <- setdiff(
  names(comm_matrix_ede),
  c("polyID", "gcm", "yr")
)

# ------------------------------------------------------------------
# Make sure FG lookup has needed species
# ------------------------------------------------------------------

fg_lookup <- qs::qread("./inputs/traitsMassFG_2025.qs")
fg_lookup <- fg_lookup %>%
  mutate(species = as.character(species))

missing_mass <- setdiff(species_cols, fg_lookup$species)

if (length(missing_mass) > 0) {
  stop(
    "Missing log10Mass values for species: ",
    paste(missing_mass, collapse = ", ")
  )
}

mass_vec <- fg_lookup %>%
  filter(species %in% species_cols) %>%
  select(species, log10Mass) %>%
  deframe()

# ------------------------------------------------------------------
# Convert presence/absence matrix to body-mass matrix
# present species = log10Mass
# absent species = NA
# ------------------------------------------------------------------
bodymass_data <- comm_matrix_ede %>%
  mutate(
    pixelID = polyID,
    lat = NA_real_,
    long = NA_real_
  ) %>%
  select(pixelID, lat, long, gcm, yr, all_of(species_cols))

for (sp in species_cols) {
  bodymass_data[[sp]] <- ifelse(
    bodymass_data[[sp]] == 1,
    mass_vec[[sp]],
    NA_real_
  )
}

# ------------------------------------------------------------------
# Create combos
# ------------------------------------------------------------------
combos <- bodymass_data %>%
  distinct(gcm, yr) %>%
  arrange(gcm, yr)

print(combos)

# ------------------------------------------------------------------
# Create bodymass list by GCM/year
# IMPORTANT:
# first 3 columns must be pixelID, lat, long
# species columns start at column 4
# ------------------------------------------------------------------
bodymass_list <- purrr::pmap(
  combos,
  function(gcm, yr) {
    
    nm <- glue("{gcm}_{yr}")
    message(crayon::blue(glue("Preparing bodymass matrix for {nm}")))
    
    bodymass_data %>%
      filter(gcm == !!gcm, yr == !!yr) %>%
      select(pixelID, lat, long, all_of(species_cols))
  }
)

names(bodymass_list) <- glue("{combos$gcm}_{combos$yr}")

# ------------------------------------------------------------------
# Run discontinuities
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# Run aggregations
# ------------------------------------------------------------------
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

# ------------------------------------------------------------------
# Flatten aggregation output
# ------------------------------------------------------------------
flatten_aggs <- function(x) {
  x2 <- x[purrr::map_lgl(x, ~ is.data.frame(.x) && nrow(.x) > 0)]
  if (length(x2) == 0) return(NULL)
  bind_rows(x2)
}

aggs_df_list <- purrr::map(aggs_list, flatten_aggs)

# ------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------
out_aggs <- "./outputs/aggs_edehzhie/Aggregations_Edehzhie"

if (!dir.exists(out_aggs)) {
  dir.create(out_aggs, recursive = TRUE)
}

purrr::walk2(
  names(aggs_df_list),
  aggs_df_list,
  ~ qs::qsave(.y, glue("{out_aggs}/Aggregations_{.x}_filtered.qs"))
)

qs::qsave(
  discont_list,
  "./outputs/aggs_edehzhie/discontinuities_Edehzhie.qs"
)

qs::qsave(
  aggs_list,
  "./outputs/aggs_edehzhie/aggs_list_Edehzhie.qs"
)

qs::qsave(
  aggs_df_list,
  "./outputs/aggs_edehzhie/aggs_df_list_Edehzhie.qs"
)

# ------------------------------------------------------------------
# Aggregation count per polygon
# ------------------------------------------------------------------
agg_count_list <- purrr::imap(
  aggs_df_list,
  function(df, nm) {
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df %>%
      group_by(pixelID, gcm, yr) %>%
      summarise(
        n_aggs = n_distinct(aggs),
        n_spp  = n(),
        .groups = "drop"
      )
  }
)

agg_count_df <- bind_rows(agg_count_list)

qs::qsave(
  agg_count_df,
  "./outputs/aggs_edehzhie/aggregation_count_per_polygon.qs"
)

write_csv(
  agg_count_df,
  "./outputs/aggs_edehzhie/aggregation_count_per_polygon.csv"
)

message("Done.")
