require(pacman)
pacman::p_load(
  terra, dplyr, fs, glue, stringr, tibble, readr, purrr
)

rm(list = ls())
options(scipen = 999)

# ------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------
pred_dir <- "E:/GitHub/birdsNWT/outputs"

out_dir_250m <- "./outputs/pres_abs_250m"
out_dir_10km <- "./outputs/pres_abs_10km"
out_dir_20km <- "./outputs/pres_abs_20km"

if (!dir_exists(out_dir_250m)) dir_create(out_dir_250m, recurse = TRUE)
if (!dir_exists(out_dir_10km)) dir_create(out_dir_10km, recurse = TRUE)
if (!dir_exists(out_dir_20km)) dir_create(out_dir_20km, recurse = TRUE)

thrs <- read_csv("./inputs/prevOcc3.csv", show_col_types = FALSE)

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
# Find files
# ------------------------------------------------------------------
species_dirs <- dir_ls(pred_dir, type = "directory")
sp_dirs <- species_dirs[path_file(species_dirs) %in% target_species]

ras_files <- dir_ls(
  sp_dirs,
  recurse = TRUE,
  glob = "*.tif"
)

meta <- tibble(file = ras_files) %>%
  mutate(
    fname   = basename(file),
    species = path_file(path_dir(file)),
    yr      = str_extract(fname, "(?<=_)\\d{4}(?=_)"),
    gcm     = str_extract(fname, "CanESM2|CCSM4|INM-CM4")
  ) %>%
  filter(
    species %in% target_species,
    yr %in% c("2011", "2031", "2091"),
    !is.na(gcm),
    str_detect(fname, "^mean_.*\\.tif$")
  ) %>%
  arrange(species, gcm, yr)

# ------------------------------------------------------------------
# Join thresholds
# ------------------------------------------------------------------
meta <- meta %>%
  left_join(
    thrs %>% select(spec, meandens),
    by = c("species" = "spec")
  )

missing_thrs <- meta %>%
  filter(is.na(meandens)) %>%
  distinct(species)

if (nrow(missing_thrs) > 0) {
  stop("Missing thresholds for species: ",
       paste(missing_thrs$species, collapse = ", "))
}

# ------------------------------------------------------------------
# Aggregation factors
# 250 m -> 10 km = 10000 / 250 = 40
# 250 m -> 20 km = 20000 / 250 = 80
# ------------------------------------------------------------------
fact_10km <- 40
fact_20km <- 80

# ------------------------------------------------------------------
# Function
# ------------------------------------------------------------------
make_pa_rasters <- function(file, species, yr, gcm, meandens) {
  
  message("--------------------------------------------------")
  message("Processing: ", species, " | ", gcm, " | ", yr)
  message("Threshold: ", signif(meandens, 4))
  
  r <- rast(file)
  
  # 1) threshold at 250 m
  # presence = 1 if density >= threshold, else 0
  r_pa_250m <- ifel(r >= meandens, 1, 0)
  names(r_pa_250m) <- glue("pa_{species}_{yr}_{gcm}_250m")
  
  # 2) aggregate to 10 km using max
  r_pa_10km <- aggregate(
    r_pa_250m,
    fact = fact_10km,
    fun = max,
    na.rm = TRUE
  )
  names(r_pa_10km) <- glue("pa_{species}_{yr}_{gcm}_10km")
  
  # 3) aggregate to 20 km using max
  r_pa_20km <- aggregate(
    r_pa_250m,
    fact = fact_20km,
    fun = max,
    na.rm = TRUE
  )
  names(r_pa_20km) <- glue("pa_{species}_{yr}_{gcm}_20km")
  
  # save
  f250 <- file.path(out_dir_250m, glue("pa_{species}_{yr}_{gcm}_250m.tif"))
  f10  <- file.path(out_dir_10km, glue("pa_{species}_{yr}_{gcm}_10km.tif"))
  f20  <- file.path(out_dir_20km, glue("pa_{species}_{yr}_{gcm}_20km.tif"))
  
  writeRaster(r_pa_250m, f250, overwrite = TRUE)
  writeRaster(r_pa_10km, f10, overwrite = TRUE)
  writeRaster(r_pa_20km, f20, overwrite = TRUE)
  
  message("Saved:")
  message("  ", f250)
  message("  ", f10)
  message("  ", f20)
  
  invisible(TRUE)
}

# ------------------------------------------------------------------
# Run
# ------------------------------------------------------------------
pwalk(
  meta %>% select(file, species, yr, gcm, meandens),
  make_pa_rasters
)

message("Done.")
