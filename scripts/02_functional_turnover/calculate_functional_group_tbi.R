# ======================================================================
# TBI analysis for FG presence-absence matrices
# Comparisons: 2011-2031 and 2011-2091
# ======================================================================

require(pacman)
pacman::p_load(
  dplyr, tidyr, purrr, qs, fs, glue, stringr, tibble, adespatial
)

rm(list = ls())
options(scipen = 999)

# ----------------------------------------------------------------------
# Directories
# ----------------------------------------------------------------------
in_dir  <- "./tables/fg_matrices_allgcms_by_scale"
out_dir <- "./outputs/tbi_fg_multiscale"
if (!fs::dir_exists(out_dir)) fs::dir_create(out_dir, recurse = TRUE)

# ----------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------
files <- list.files(in_dir, pattern = "\\.qs$", full.names = TRUE)

comparisons <- list(
  c("2011", "2031"),
  c("2011", "2091")
)

# ----------------------------------------------------------------------
# Function
# ----------------------------------------------------------------------
run_tbi_one_file <- function(f, yr1, yr2, out_dir) {
  
  fname <- basename(f)
  gcm   <- stringr::str_extract(fname, "CanESM2|CCSM4|INM-CM4")
  scale <- stringr::str_extract(fname, "4km|10km|20km")
  
  message("Running TBI: ", gcm, " | ", scale, " | ", yr1, "-", yr2)
  
  dat <- qs::qread(f)
  
  # keep only the two years
  d_sub <- dat %>%
    filter(yr %in% c(yr1, yr2))
  
  if (nrow(d_sub) == 0) {
    warning("No rows for ", gcm, " | ", scale, " | ", yr1, "-", yr2)
    return(NULL)
  }
  
  id_cols <- c("pixelID", "x", "y", "yr", "gcm")
  fg_cols <- setdiff(names(d_sub), id_cols)
  
  if (length(fg_cols) == 0) {
    warning("No FG columns for ", gcm, " | ", scale)
    return(NULL)
  }
  
  # split by year
  d1 <- d_sub %>%
    filter(yr == yr1) %>%
    arrange(pixelID)
  
  d2 <- d_sub %>%
    filter(yr == yr2) %>%
    arrange(pixelID)

  if (nrow(d1) == 0 || nrow(d2) == 0) {
    warning("Missing rows for ", fname, " | ", yr1, "-", yr2)
    return(NULL)
  }
  
  if (nrow(d1) != nrow(d2)) {
    stop("Different number of rows for ", fname, " | ", yr1, "-", yr2)
  }
  
  if (!identical(d1$pixelID, d2$pixelID)) {
    stop("pixelID order does not match for ", fname, " | ", yr1, "-", yr2)
  }
  
  # matrices
  mat1 <- d1 %>%
    select(pixelID, all_of(fg_cols)) %>%
    column_to_rownames("pixelID") %>%
   as.matrix()
  
  mat2 <- d2 %>%
    select(pixelID, all_of(fg_cols)) %>%
    column_to_rownames("pixelID") %>%
    as.matrix()
 
  # identify empty rows
  empty1 <- rowSums(mat1, na.rm = TRUE) == 0
  empty2 <- rowSums(mat2, na.rm = TRUE) == 0
  keep   <- !empty1 & !empty2
  
  # full pixel table for mapping
  full_tbl <- tibble(
    pixelID = d1$pixelID,
    x       = d1$x,
    y       = d1$y,
    gcm     = gcm,
    scale   = scale,
    yr1     = yr1,
    yr2     = yr2,
    empty_yr1 = empty1,
    empty_yr2 = empty2,
    valid_tbi = keep
  )
  
  # subset only valid pixels for TBI
  mat1_keep <- mat1[keep, , drop = FALSE]
  mat2_keep <- mat2[keep, , drop = FALSE]
  
  d1_keep <- d1[keep, ]
  
    # Get Temporal Beta diversity Index (TBI) ---------------------------------
  tbi_res <- adespatial::TBI(mat1_keep, mat2_keep, method = 'sorensen', BCD = TRUE, nperm = 999,
                                   test.BC = TRUE, save.BC = TRUE, seed. = 1234,
                                   test.t.perm = TRUE, 
                                   clock = TRUE)

  tbi_tbl <- tibble(
    pixelID = d1_keep$pixelID,
    TBI     = tbi_res$BCD.mat[, 1],
    loss    = tbi_res$BCD.mat[, 2],
    gain    = tbi_res$BCD.mat[, 3],
    change  = tbi_res$BCD.mat[, 4],
    p_TBI   = tbi_res$p.TBI,
    p_loss  = tbi_res$p.adj
    )
  
  # join back to full table
  out_tbl <- full_tbl %>%
    left_join(tbi_tbl, by = "pixelID")
  out_stub <- glue("TBI_fg_{gcm}_{scale}_{yr1}_{yr2}")
  
  qs::qsave(out_tbl, glue("{out_dir}/{out_stub}_table.qs"))
  qs::qsave(tbi_res, glue("{out_dir}/{out_stub}_raw.qs"))
  
  message("Saved: ", out_stub)
  
  return(out_tbl)
}

# ----------------------------------------------------------------------
# Run all files and comparisons
# ----------------------------------------------------------------------
results <- list()

for (f in files) {
  for (cmp in comparisons) {
    yr1 <- cmp[1]
    yr2 <- cmp[2]
    
    results[[paste(basename(f), yr1, yr2, sep = "_")]] <- run_tbi_one_file(
      f = f,
      yr1 = yr1,
      yr2 = yr2,
      out_dir = out_dir
    )
  }
}

