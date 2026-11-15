library(pacman)
p_load(glue, fs, tydiverse)


# Folder containing the 72 *varimp8.csv files
cov_dir <- "E:/GitHub/birdsNWT/covariates"

keep_species <- c(
  "ALFL","AMCR","AMRE","AMRO","ATSP","ATTW","BAWW","BBWA","BBWO",
  "BCCH","BHCO","BHVI","BLPW","BOCH","BRBL","BRCR","BTNW","CAWA",
  "CCSP","CHSP","CONW","CORA","COYE","DEJU","EAKI","EVGR","FOSP",
  "GCKI","GCTH","GRAJ","HAFL","HETH","HOLA","LCSP","LEFL","LEYE",
  "LISP","MAWA","NOWA","OCWA","OSFL","OVEN","PAWA","PHVI","PISI",
  "PUFI","RBGR","RBNU","RCKI","REVI","RUBL","RUGR","RWBL","SAVS",
  "SOSP","SWSP","SWTH","TEWA","TRES","VATH","WAVI","WCSP","WETA",
  "WEWP","WIWA","WIWR","WTSP","WWCR","YBFL","YBSA","YEWA","YRWA"
)

# Read and combine all variable importance files
varimp_all <- list.files(
  cov_dir,
  pattern = "varimp8\\.csv$",
  full.names = TRUE
) %>%
  map_dfr(function(f) {
    read_csv(f, show_col_types = FALSE) %>%
      mutate(
        species = str_extract(basename(f), "^[A-Z]{4}")
      )
  })

varimp_all <- varimp_all %>%
  filter(species %in% keep_species)

# Order predictors by median relative influence
var_order <- varimp_all %>%
  group_by(var) %>%
  summarise(
    median_rel = median(rel.inf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_rel) %>%
  pull(var)

varimp_all <- varimp_all %>%
  mutate(
    var = factor(var, levels = var_order)
  )

# Median values for labels
median_vals <- varimp_all %>%
  group_by(var) %>%
  summarise(
    median_rel = median(rel.inf, na.rm = TRUE),
    .groups = "drop"
  )

# Make figure
p <- ggplot(varimp_all, aes(x = rel.inf, y = var)) +
  geom_boxplot(
    fill = "white",
    colour = "black",
    outlier.alpha = 0.5
  ) +
  geom_text(
    data = median_vals,
    aes(
      x = 62,
      y = var,
      label = sprintf("%.1f", median_rel)
    ),
    hjust = 0,
    size = 3
  ) +
  annotate(
    "text",
    x = 62,
    y = nlevels(varimp_all$var) + 1.5,
    label = "Median (%)",
    fontface = "bold",
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Relative influence (%)",
    y = "Predictor"
  ) +
  coord_cartesian(
    xlim = c(0, 72),
    clip = "off"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 80, b = 10, l = 10)
  )

p

# Save figure
ggsave(
  filename = file.path(cov_dir, "Figure_S1_relative_influence_boxplots_median.png"),
  plot = p,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Create output folder if it doesn't exist
dir_create("figs/BRT_predictors")

if (!dir.exists("figs/BRT_predictors")) {
  dir.create(
    "figs/BRT_predictors",
    recursive = TRUE
  )
}

ggsave(
  filename = "figs/BRT_predictors/Figure_S1_relative_influence_boxplots_median.png",
  plot = p,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)


# Top 5 predictors per species
table_S2_top5 <- varimp_all %>%
  group_by(species) %>%
  arrange(desc(rel.inf), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Wide format: one row per species
table_S2_top5_wide <- table_S2_top5 %>%
  select(species, rank, var, rel.inf) %>%
  pivot_wider(
    names_from = rank,
    values_from = c(var, rel.inf),
    names_glue = "{.value}_{rank}"
  ) %>%
  arrange(species)

# Save
write_csv(
  table_S2_top5_wide,
  "tables/BRT_predictors/Table_S2_top5_predictors_per_species.csv")

  
  # Create species × predictor matrix
table_S2 <- varimp_all %>%
    select(species, var, rel.inf) %>%
    pivot_wider(
      names_from = var,
      values_from = rel.inf,
      values_fill = 0
    ) %>%
    arrange(species)
  
table_S2 <- varimp_all %>%
  select(species, var, rel.inf) %>%
  pivot_wider(
    names_from = var,
    values_from = rel.inf,
    values_fill = 0
  ) %>%
  arrange(species) %>%
  mutate(
    across(
      -species,
      ~ ifelse(.x == 0, NA, round(.x, 3))
    )
  )

write_csv(
  table_S2,
  "tables/BRT_predictors/Table_S2_species_variable_importance_matrix.csv",
  na = ""
)
