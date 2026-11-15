require(pacman)
pacman::p_load(
  terra, dplyr, fs, glue, stringr, purrr, qs, readr, tidyr
)

rm(list = ls())
options(scipen = 999)

# ------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------
in_dir_10km <- "./outputs/pres_abs_10km"
in_dir_20km <- "./outputs/pres_abs_20km"

out_dir <- "./tables/community_matrices"
if (!fs::dir_exists(out_dir)) fs::dir_create(out_dir, recurse = TRUE)

yrs  <- c("2011", "2031", "2091")
gcms <- c("CanESM2", "CCSM4", "INM-CM4")
scales <- c("10km", "20km")

target_species <- c(
  "ALFL","AMCR","AMRE","AMRO","ATSP","ATTW","BAWW","BBWA","BBWO","BCCH",
  "BHCO","BHVI","BLPW","BOCH","BRBL","BRCR","BTNW","CAWA","CCSP","CHSP",
  "CONW","CORA","COYE","DEJU","EAKI","EVGR","FOSP","GCKI","GCTH","GRAJ",
  "HAFL","HETH","HOLA","LCSP","LEFL","LEYE","LISP","MAWA","NOWA","OCWA",
  "OSFL","OVEN","PAWA","PHVI","PISI","PUFI","RBGR","RBNU","RCKI","REVI",
  "RUBL","RUGR","RWBL","SAVS","SOSP","SWSP","SWTH","TEWA","TRES","VATH",
  "WAVI","WCSP","WETA","WEWP","WIWA","WIWR","WTSP","WWCR","YBFL","YBSA",
  "YEWA","YRWA"
)

# ------------------------------------------------------------------
# Helper to choose input folder
# ------------------------------------------------------------------
get_input_dir <- function(scale_sel) {
  if (scale_sel == "10km") return(in_dir_10km)
  if (scale_sel == "20km") return(in_dir_20km)
  stop("Unknown scale: ", scale_sel)
}

# ------------------------------------------------------------------
# Build metadata from existing files
# ------------------------------------------------------------------
build_meta <- function(scale_sel) {
  
  in_dir <- get_input_dir(scale_sel)
  
  files <- list.files(in_dir, pattern = "\\.tif$", full.names = TRUE)
  
  tibble(
    file = files,
    fname = basename(files)
  ) %>%
    mutate(
      scale   = scale_sel,
      species = str_extract(fname, "(?<=pa_)[A-Z]{4}(?=_)"),
      yr      = str_extract(fname, "(?<=_)(2011|2031|2091)(?=_)"),
      gcm     = str_extract(fname, "CanESM2|CCSM4|INM-CM4")
    ) %>%
    filter(
      species %in% target_species,
      !is.na(yr),
      !is.na(gcm)
    )
}

meta_all <- bind_rows(
  build_meta("10km"),
  build_meta("20km")
)

# quick check
print(meta_all %>% count(scale, yr, gcm))

# ------------------------------------------------------------------
# Function to build one community matrix
# ------------------------------------------------------------------
make_comm_matrix <- function(scale_sel, yr_sel, gcm_sel, meta_all) {
  
  message("--------------------------------------------------")
  message(glue("Building matrix: {scale_sel} | {yr_sel} | {gcm_sel}"))
  
  sub <- meta_all %>%
    filter(
      scale == scale_sel,
      yr    == yr_sel,
      gcm   == gcm_sel
    ) %>%
    arrange(species)
  
  if (nrow(sub) == 0) {
    warning(glue("No files found for: {scale_sel} | {yr_sel} | {gcm_sel}"))
    return(NULL)
  }
  
  # stack rasters
  r <- terra::rast(sub$file)
  names(r) <- sub$species
  
  # convert to data frame
  df <- terra::as.data.frame(
    r,
    xy = TRUE,
    cells = TRUE,
    na.rm = FALSE
  ) %>%
    rename(pixelID = cell)
  
  # save
  stub <- glue("comm_{scale_sel}_{yr_sel}_{gcm_sel}")
  
  qs::qsave(
    df,
    glue("{out_dir}/{stub}.qs")
  )
  
   message(glue("Saved: {stub}"))
  message(glue("Rows: {nrow(df)} | Species: {nrow(sub)}"))
  
  return(df)
}

# ------------------------------------------------------------------
# Run all combinations
# ------------------------------------------------------------------
combos <- tidyr::crossing(
  scale_sel = scales,
  yr_sel    = yrs,
  gcm_sel   = gcms
)

comm_list <- purrr::pmap(
  combos,
  make_comm_matrix,
  meta_all = meta_all
)

message("Done.")