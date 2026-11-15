

## ResilienceNWT_2026
This repository contains the R scripts and derived outputs supporting the manuscript entitled 
*Northern landbird communities show geographically constrained resilience to projected 
climate and forest changes*. 

The workflow analyzes projected boreal landbird communities in the Taiga Plains, Northwest Territories, 
across three GCMs and the 2011, 2031, and 2091 time periods.

The analyses cover presence–absence and community matrices, taxonomic and 
functional-group temporal beta diversity (TBI), richness and persistence,
body-mass aggregations, cross-scale resilience metrics, and NMDS/PERMANOVA. 

The Edéhzhíe workflow supports the body-mass aggregation illustration in Figure 2.


## Repository structure

```text
resilienceNWT_2026/
├── inputs/     # Traits, spatial boundaries, and other input data
├── tables/     # Processed matrices, summaries, and manuscript tables
├── outputs/    # Intermediate TBI, aggregation, resilience, and NMDS results
├── figs/       # Main-text and supplementary figures
├── R/          # Shared functions and supporting files
│   ├── balanceSampling.R
│   ├── functions_discontinuities.R
│   └── GRI.constant.power.table.csv
└── scripts/    # R scripts used in the workflow
    ├── 01_prepare_communities/  # Community preparation
    │   ├── create_presence_absence_rasters.R
    │   ├── build_species_community_matrices.R
    │   └── build_functional_group_matrices.R
    ├── 02_functional_turnover/  # Functional turnover
    │   ├── calculate_species_tbi.R
    │   ├── calculate_functional_group_tbi.R
    │   ├── create_species_tbi_rasters.R
    │   ├── create_functional_group_tbi_rasters.R
    │   ├── summarize_species_tbi_change.R
    │   └── summarize_functional_group_tbi_change.R
    ├── 03_richness_persistence/  # Richness and persistence
    │   ├── create_species_richness_maps.R
    │   ├── create_functional_group_richness_maps.R
    │   └── create_functional_group_gain_loss_heatmaps.R
    ├── 04_body_mass_aggregations/  # Body-mass aggregations
    │   ├── build_log_mass_matrices.R
    │   ├── calculate_aggregations_and_discontinuities.R
    │   └── plot_aggregations.R
    ├── 05_resilience_metrics/  # Resilience metrics
    │   ├── calculate_resilience_metrics.R
    │   ├── summarize_resilience_change.R
    │   ├── plot_resilience_change_relationships.R
    │   └── map_baseline_resilience_metrics.R
    ├── 06_nmds_permanova/  # NMDS and PERMANOVA
    │   ├── create_balanced_nmds_sample.R
    │   ├── run_nmds_by_region.R
    │   ├── summarize_permanova.R
    │   └── plot_nmds_by_gcm.R
    ├── 07_edehzhie_case_study/  # Edéhzhíe case study
    │   ├── extract_edehzhie_communities.R
    │   ├── build_edehzhie_community_matrices.R
    │   ├── build_edehzhie_log_mass_matrices.R
    │   ├── calculate_edehzhie_aggregations.R
    │   ├── plot_edehzhie_aggregations.R
    │   ├── plot_edehzhie_resilience_changes.R
    │   └── plot_edehzhie_aggregation_bands.R
    └── 08_manuscript_support/  # Manuscript support
        ├── plot_brt_predictor_importance.R
        └── create_table_s2_species_predictors.R
```

## Step-by-Step Workflow

Paths are relative to the repository root. `{gcm}`, `{yr}`, `{yr1}`, `{yr2}`, `{scale}`, `{region}`, and `{var}` represent values that vary by model, year, scale, region, or metric.

| Step | Script name | Purpose | Inputs | Outputs |
|---:|---|---|---|---|
| 1 | `create_presence_absence_rasters.R` | Apply the presence threshold to 250 m density rasters and aggregate presence–absence data to analysis scales. | 250 m density rasters **(available upon request).** | `outputs/pres_abs_250m/`; `outputs/pres_abs_4km/`; `outputs/pres_abs_10km/`; `outputs/pres_abs_20km/` |
| 2 | `build_species_community_matrices.R` | Build species community matrices from the presence–absence rasters at each scale. | `outputs/pres_abs_250m/`, `outputs/pres_abs_4km/`, `outputs/pres_abs_10km/`, `outputs/pres_abs_20km/` | `tables/community_matrices/` |
| 3 | `build_functional_group_matrices.R` | Convert species community matrices to functional-group presence–absence matrices. | `tables/comm_matrices_allgcms_by_scale/`; `inputs/traitsMassFG_2025.qs` | `tables/fg_pa_by_gcm_scale/` |
| 4 | `calculate_species_tbi.R` | Calculate taxonomic TBI for the 72 species for 2011–2031 and 2011–2091. | Species community matrices in `tables/community_matrices/` `tables/fg_matrices_allgcms_by_scale/` | `outputs/tbi_72spp_multiscale/` |
| 5 | `calculate_functional_group_tbi.R` | Calculate TBI for functional groups for 2011–2031 and 2011–2091. | `tables/fg_matrices_allgcms_by_scale/` | `outputs/tbi_fg_multiscale/` |
| 6 | `create_species_tbi_rasters.R` | Combine species TBI, p-value, and change results; attach coordinates and create raster/figure outputs. |  `outputs/tbi_fg_multiscale/TBI_fg_{gcm}_{scale}_{yr1}_{yr2}_table.qs`,  `outputs/tbi_72spp_multiscale/`. | `tables/tbi72sppBin2026/TBIpval72sppBin_allgcms_{yr1}_{yr2}_{scale}.qs`; `figs/tbi_72spp_multiscale_maps/` |
| 7 | `create_functional_group_tbi_rasters.R` | Combine functional-group TBI, p-value, and change results; attach coordinates and create raster/figure outputs. | `outputs/tbi_fg_multiscale/`; `TBI_fg_{gcm}_{scale}_{yr1}_{yr2}_table.qs` | `tables/tbifg_2026_multiscale/tbiFG_2026_allgcms_allyrs_allscales.qs`; `figs/tbi_fg_multiscales_maps/` |
| 8 | `summarize_species_tbi_change.R` | Summarize mean species TBI change and the number of GCMs with significant change. | `tables/tbi72sppBin2026/tbi72spp_2026_allgcms_allyrs_allscales.qs` | `figs/TBI_72spp_changemean2026/{scale}/TBI_mean_Change_72spp_{scale}.png`; `figs/TBI_72spp_changemean2026/{scale}/72sppTBI_Significance_Count_Bin_{scale}.png` |
| 9 | `summarize_functional_group_tbi_change.R` | Summarize mean functional-group TBI change and significant-change counts. | Combined functional-group TBI results from Step 7  | `figs/TBI_fg_changemean2026/{scale}/TBI_mean_Change_FG_{scale}.png`; `figs/TBI_fg_changemean2026/{scale}/FGTBI_Significance_Count_Bin_{scale}.png` |
| 10 | `create_species_richness_maps.R` | Create species richness maps at 4, 10, and 20 km for each year and GCM, plus combined maps. | `tables/community_matrices/` | `figs/richness_maps_multiscale/`; `figs/richness_maps_multiscale/combined/` |
| 11 | `create_functional_group_richness_maps.R` | Map functional-group richness by year and GCM. | `tables/fg_matrices_allgcms_by_scale/` | `figs/FG_richness_maps/FG_richness_ALLGCMs4.png` |
| 12 | `create_functional_group_gain_loss_heatmaps.R` | Summarize functional-group presence/gain/loss by region and plot persistence heatmaps. | `tables/fg_matrices_allgcms_by_scale/`; `inputs/traitsMassFG_2025.qs` | `outputs/fg_presence_by_region/`; `figs/persistence/FG_gain_loss_heatmap_2011_2031.png`; `figs/persistence/FG_gain_loss_heatmap_2011_2091.png`; `figs/persistence/FG_gain_loss_heatmaps_prop_combined.png` |
| 13 | `build_log_mass_matrices.R` | Create log body-mass matrices for 4 km community matrices. | `tables/comm_matrices_filtered_4km/ppCommMatrices_filtered_4km.qs` | `tables/logMassMat_4km_filtered/logMassMatrix_{gcm}_{yr}_filtered.qs` |
| 14 | `calculate_aggregations_and_discontinuities.R` | Calculate body-mass aggregations and discontinuities. | `R/functions_discontinuities.R`; `R/GRI.constant.power.table.csv`; log-mass matrices  `tables/logMassMat/allLogMassXY_recoded.qs`; reconcile with Step 13). | `outputs/Aggregations_4km_filtered/`; `outputs/Aggregations_4km_filtered/aggregation_count_per_pixel.qs` |
| 15 | `plot_aggregations.R` | Plot aggregation results. | `outputs/Aggregations_4km_filtered/` | `figs/aggregations2026/4km_filtered/` |
| 16 | `calculate_resilience_metrics.R` | Calculate cross-scale redundancy, within-scale redundancy, and cross-scale diversity for three GCMs and three years. Preserve `pixelID` in saved results. | `outputs/Aggregations_4km_filtered/` | `outputs/ResilienceMetrics_4km/resilience_metrics_all.qs` |
| 17 | `summarize_resilience_change.R` | Calculate per-pixel absolute change from 2011 and mean change across GCMs; create mean-change maps. | `outputs/ResilienceMetrics_4km/resilience_metrics_all.qs` | `outputs/ResilienceMetrics_4km/delta_px_from_2011.qs`; `outputs/ResilienceMetrics_4km/mean_delta_px_across_gcms.qs`; `figs/redundancy_delta/2026/{var}_mean_delta.png` |
| 18 | `plot_resilience_change_relationships.R` | Plot regional relationships among changes in resilience metrics. | `outputs/ResilienceMetrics_4km/resilience_metrics_all.qs` | `figs/Redundancy/2026/delta_change/delta_change_cross_vs_within_region.png`; `figs/Redundancy/2026/delta_change/delta_change_cross_scale_div_region.png` |
| 19 | `map_baseline_resilience_metrics.R` | Create maps of baseline 2011 resilience metrics. | `outputs/ResilienceMetrics_4km/resilience_metrics_all.qs` | `figs/redundancy/baseline/2011/facet_metrics_mean_2011.png` |
| 20 | `create_balanced_nmds_sample.R` | Select balanced sampling units for NMDS. | `R/balanceSampling.R`; `inputs/pix_xy_alleco_4km_2026.qs` | `tables/sampling/` |
| 21 | `run_nmds_by_region.R` | Run NMDS for five regions, three years, and selected GCMs; run associated resampling and PERMANOVA analyses. | `tables/comm_matrices_filtered_4km/fgCommMatrices_filtered_4km.qs`; `tables/sampling/` | `outputs/nmds_results_2026/nmds_{region}_{gcm}_10percent_allyrs.qs`; `outputs/nmds_results_2026/sites_{region}_{gcm}_10percent_allyrs.qs`; `outputs/bootstrapping_2026/`; `outputs/permanova_2026/` |
| 22 | `summarize_permanova.R` | Extract and format PERMANOVA results for the manuscript. | `outputs/permanova_2026/` | `tables/Table4_PERMANOVA_summary_{gcm}.csv` |
| 23 | `plot_nmds_by_gcm.R` | Plot NMDS ordinations by GCM. | `outputs/nmds_results_2026/` | `figs/nmds_sampling_2026/NMDS_combined_{gcm}_10percent_allyrs.png` |
| 24 | `extract_edehzhie_communities.R` | Extract or aggregate presence–absence data to the Edéhzhíe boundary. | `inputs/Edehzhie/Edehzhie_Boundary.shp` | `outputs/pres_abs_by_polygon/pres_abs_Edehzhie.gpkg` |
| 25 | `build_edehzhie_community_matrices.R` | Build community matrices for Edéhzhíe. | `outputs/pres_abs_by_polygon/pres_abs_Edehzhie.gpkg` | `tables/community_matrices/comm_matrix_Edehzhie.qs` |
| 26 | `build_edehzhie_log_mass_matrices.R` | Build Edéhzhíe log-mass matrices. | `tables/decho_edeh_communities/sppCommMatrices_dechoEdeh_4km.qs` | `tables/logMassMat_4km_decho/logMassMatrix_{gcm}_{yr}_filtered.qs` |
| 27 | `calculate_edehzhie_aggregations.R` | Calculate Edéhzhíe body-mass aggregations and discontinuities. | `R/functions_discontinuities.R`; `R/GRI.constant.power.table.csv`; `inputs/traitsMassFG_2025.qs`; `tables/community_matrices/comm_matrix_Edehzhie.qs` | `outputs/aggs_edehzhie/Aggregations_Edehzhie.qs` |
| 28 | `plot_edehzhie_aggregations.R` | Plot Edéhzhíe body-mass aggregation results for Figure 2. | `outputs/aggs_edehzhie/Aggregations_Edehzhie/` | `figs/aggregations2026/edeh/` |
| 29 | `plot_edehzhie_resilience_changes.R` | Plot Edéhzhíe changes in redundancy and cross-scale diversity. | `outputs/resilience_metrics_edehzhie/resilience_metrics_Edehzhie.qs` | `figs/resilience_metrics_edehzhie/Edehzhie_delta_CSR_WSR.png`; `figs/resilience_metrics_edehzhie/Edehzhie_delta_cross_scale_diversity.png` |
| 30 | `plot_edehzhie_aggregation_bands.R` | Plot Edéhzhíe aggregation bands. | `outputs/aggs_edehzhie/Aggregations_Edehzhie/` | `outputs/aggs_edehzhie/Aggregations_Edehzhie/` |
| 31 | `plot_brt_predictor_importance.R` | Plot relative influence of BRT predictors across species (Supplementary Figure S2). | `birdsNWT/covariates/` and species-level BRT results  | `figs/BRT_predictors/Figure_S2_relative_influence_boxplots_median.png` |
| 32 | `create_table_s2_species_predictors.R` | Format the species-level predictor importance results for the manuscript/supplement. | `tables/BRT_predictors/Table_S2_top5_predictors_per_species.csv` | `tables/Table_S2_species_variable_importance_matrix.csv` |

