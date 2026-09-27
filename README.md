# What Moves Climate Action? Structural Uncertainty and Intervention in Climate Attitude Networks

Code for a project examining how climate-related attitudes, beliefs, and context variables are
connected, which of those connections have a data-supported causal direction, and what a working
causal model predicts about intervening on them to change later climate behavior. It also examines
how those predictions change when uncertain relationships are given different directions or
attributed to unmeasured common causes. Uses five waves of a U.S. survey (N = 1,987 for the
attitude network, N = 870 for the behavior outcome).

The analysis combines a Gaussian graphical model (descriptive network), constraint-based causal
discovery (FCI, with PC-stable as a sensitivity check) with bootstrap stability analysis, a fully
directed working structural causal model fit in `lavaan`, and several layers of sensitivity analysis
(directional-completion, latent-confounding-compatible, and a PAG-compatible LV-IDA analysis that
avoids committing to one directed completion at all).

## Repository structure

- `clean_pipeline/` -- the confirmatory analysis, pulled out of the exploratory `.qmd` into scripts
  that run top to bottom in a fresh R session. Start here. See its own README for run order and
  what each script does.
- `r_patches/` -- figure scripts and additional sensitivity/diagnostic analyses, numbered roughly in
  the order they were written. Not every script feeds the paper. Superseded figure variants (for
  example the `10_figure7_*` scripts) and one-off check scripts are kept as a record. The scripts
  behind the paper's figures and tables are listed [below](#manuscript-figures-and-tables).
- `lv-ida/` -- LV-IDA implementation (Malinsky and Spirtes, 2016), used for the PAG-compatible
  effect analysis.
- `figures/`, `tables/` -- rendered figures and output tables.
- `pipeline_outputs/` -- intermediate and final output from `clean_pipeline/` (bootstrap arrays,
  fitted model objects, CSVs).
- `clean_pipeline_alt_will/`, `pipeline_outputs_alt_will/` -- a side-by-side sensitivity check using
  `cvcc4_will` (a descriptive-norm survey item) instead of `cvcc4_should` for the social-norms node.
- `climate_analysis_avg_v2_altweather.qmd` -- the exploratory analysis notebook (measurement
  reduction, clustering, EFA, construct-validity decisions). Kept as the record of those one-time
  decisions rather than turned into a script.
- `analysis_decisions_log.md` -- running log of analysis decisions and why they were made.

## Data

The raw participant-level survey data is not included in this repository.

Get in touch for access to the raw data if you want to reproduce the pipeline end to end.

## Reproducing the analysis

From the repository root, in a fresh R session:

```r
source("clean_pipeline/run_all.R")
```

This runs `clean_pipeline/01` through `14` in order and freezes the result (with a manifest and a
baseline sanity check) into `pipeline_outputs/reference_runs/`. Scripts `15`-`19` are follow-on
sensitivity and extension analyses. They are not needed for the core network, causal-discovery, and
working-model results, but several paper figures use their output (see the table below), so run them
too if you want to rebuild every figure. `20_stability_matrix_figures.R` regenerates the
supplementary stability-matrix figure straight from `03`'s bootstrap output and doesn't need the
rest of the chain. See `clean_pipeline/README.md` for details and dependencies.

Figure and supplementary-analysis scripts in `r_patches/` are run separately once
`clean_pipeline/`'s output is in place. `r_patches/41_regenerate_all_figures.R` re-exports
Figures 1-5 and most supplementary figures in one go. The two measurement-reduction figures
(`20_supplementary_measurement_figures.R`) need objects from the `.qmd` and are run on their own.

Several `r_patches/` scripts start with a `setwd()` pointing to the author's local copy of the
repository. Remove or edit that line and run from the repository root.

## Manuscript figures and tables

| Paper item | Output file | Script | Also needs |
|---|---|---|---|
| Figure 1, conditional-dependence networks | `figures/fig2_ggm_main_fixed.pdf`, `figures/fig3_ggm_extended_fixed.pdf` | `r_patches/36_revision_figures_1_2A.R` | `clean_pipeline/01`-`03` |
| Figure 2, causal-discovery stability and PAG | `figures/fig_causal_combined.pdf` | `r_patches/38_figure_causal_combined.R` | `clean_pipeline/03`, `16` |
| Figure 3, working structural causal model | `figures/fig5_scm_hierarchical.pdf` | `r_patches/07_figure5_scm_hierarchical_v3.R` | `clean_pipeline/04` |
| Figure 4, predicted intervention effects | `figures/fig4_workingscm.pdf` | `r_patches/34_figure_bootstrap_lvida_ridge.R` | `clean_pipeline/07`, `09`, `18` |
| Figure 5 (NCC version), sensitivity to relationships compatible with unmeasured common causes | `figures/fig_confound_sensitivity.pdf` | `r_patches/51_figure_confound_sensitivity.R` | `clean_pipeline/15` |
| Figure 5 (longer version) and Supplementary figure (NCC version), dependence of the two leading effects on specific relationships | `figures/fig_pathway_dependence.pdf` | `r_patches/44_figure_pathway_dependence.R` | `clean_pipeline/15`, `r_patches/37` |
| Optional NCC Figure 4 (not yet used), Figures 4 and 5 combined as panels A and B | `figures/fig_sensitivity_combined.pdf` | `r_patches/52_figure_sensitivity_combined.R` | `clean_pipeline/07`, `09`, `15` |
| Supplementary, construct correlation heatmap | `figures/fig1_heatmap_construct.pdf` | `r_patches/04_figure1_heatmap.R` | `clean_pipeline/01` (`df_main` in session) |
| Supplementary, item dendrogram and PCA scree plot | `figures/figS_item_dendrogram.pdf`, `figures/figS_pca_scree.pdf` | `r_patches/20_supplementary_measurement_figures.R` | Sections 3.5-3.8 of the `.qmd` run in the same session |
| Supplementary, FCI and PC-stable stability matrices | `figures/figS_stability_matrix_detailed_ext.pdf`, `figures/figS_stability_pc_ext.pdf` | `clean_pipeline/20_stability_matrix_figures.R` | `clean_pipeline/03` |
| Supplementary, all pair and triple interventions | `figures/figS_all_combo_interventions.pdf` | `r_patches/08b_supp_figure_all_combo_interventions.R` | `clean_pipeline/08` |
| Supplementary, feedback (cyclic) specifications | `figures/fig_supp_cyclic_feedback.pdf` | `r_patches/35_figure_cyclic_feedback_equilibrium.R` | `clean_pipeline/07`, `19` |
| Supplementary, bootstrap LV-IDA effects | `figures/figS_bootstrap_lvida_pointcloud.pdf` | `r_patches/34_figure_bootstrap_lvida_ridge.R` | `clean_pipeline/18` |
| Supplementary table, network connectedness versus predicted effect | `pipeline_outputs/tables/ggm_connectedness_vs_scm_effect.csv` | `r_patches/46_ggm_connectedness_vs_effect.R` | `clean_pipeline/01`, `02`, `15` |

Node labels and colours for all figures are set in `r_patches/03_figure_style.R`. The internal
variable names (for example `belief_concern`) are unchanged from the analysis code. Only the
displayed labels (for example "Belief and concern") differ.

## Requirements

R packages: `dplyr`, `tidyr`, `readr`, `purrr`, `huge`, `qgraph`, `igraph`, `graph`, `pcalg`,
`lavaan`, `ggplot2`, `patchwork`, `tibble`, `scales`, `psych`, `bootnet`, `arrow`, `furrr`,
`future`, and `pheatmap` (for `r_patches/20_supplementary_measurement_figures.R`). `graph` and
`RBGL`, which `pcalg` depends on, are installed from Bioconductor
(`BiocManager::install(c("graph", "RBGL"))`). The sensitivity analyses in `r_patches/14`-`16`
additionally need `RCIT`/`RCoT` and `micd`. The CCI skeleton check in `r_patches/10` additionally
needs `CCI.KP` (or its predecessor `CCI`).
