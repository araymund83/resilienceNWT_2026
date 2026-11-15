# Load libraries ----------------------------
require(pacman)
pacman::p_load(dplyr,glue,qs,tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())

# Load data -----------------------------------------
in_dir <- './tables/comm_matrices_filtered_4km'
yrs <- c('2011', '2031','2091')
yr1 <- yrs[1]
yr2 <- yrs[3]
yr <- yr2
gcms <- c('CanESM2', 'CCSM4', 'INM-CM4')
log10Mass <-qs::qread('./inputs/traitsMassFG_2025.qs') 
log10Mass<- log10Mass %>% select(log10Mass)

data <- qs::qread(glue('{in_dir}/sppCommMatrices_filtered_4km.qs'))

# ---------------------------------------------------------
# Species columns
# ---------------------------------------------------------
meta_cols <- c("pixelID", "x", "y", "yr", "gcm", "all_na", "richness")
species <- setdiff(names(data), meta_cols)

if (length(species) == 0) {
  stop("No species columns detected.")
}

# ---------------------------------------------------------
# Check species in mass table
# ---------------------------------------------------------
if (!all(species %in% rownames(log10Mass))) {
  missing_spp <- setdiff(species, rownames(log10Mass))
  stop(glue("Missing species in log10Mass: {paste(missing_spp, collapse = ', ')}"))
}

# if log10Mass has one column, this works
mass_vec <- log10Mass[species, , drop = FALSE][, 1]
names(mass_vec) <- species

# ---------------------------------------------------------
# Build body-mass matrix for one gcm x yr subset
# presence -> log10 mass
# absence  -> NA
# also rename x/y to long/lat because DD uses those names
# ---------------------------------------------------------
make_bodymass_matrix_for_DD <- function(df, mass_vec, gcm, yr,
                                        out_dir = "./tables/logMassMat_4km_filtered") {
  
  df <- as.data.frame(df)
  
  # keep species in known order
  spp_df <- df[, names(mass_vec), drop = FALSE]
  
  # convert PA -> log mass matrix
  bm <- sweep(as.matrix(spp_df), 2, mass_vec, `*`)
  bm[bm == 0] <- NA_real_
  bm <- as.data.frame(bm)
  
  out <- bind_cols(
    df %>%
      transmute(
        pixelID = as.character(pixelID),
        lat = y,
        long = x
      ),
    bm
  )
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  qs::qsave(out, glue("{out_dir}/logMassMatrix_{gcm}_{yr}_filtered.qs"))
  
  return(out)
}

# ---------------------------------------------------------
# Build list of body-mass matrices by gcm and year
# ---------------------------------------------------------
combos <- tidyr::crossing(
  gcm = sort(unique(data$gcm)),
  yr  = sort(unique(data$yr))
)

bodymass_list <- purrr::pmap(
  combos,
  function(gcm, yr) {
    message(crayon::blue(glue("Building body-mass matrix for {gcm} {yr}")))
    
    df_sub <- data %>%
      filter(gcm == !!gcm, yr == !!yr)
    
    make_bodymass_matrix_for_DD(
      df = df_sub,
      mass_vec = mass_vec,
      gcm = gcm,
      yr = yr
    )
  }
)

names(bodymass_list) <- glue("{combos$gcm}_{combos$yr}")

