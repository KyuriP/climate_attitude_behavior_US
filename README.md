# What Moves Climate Action? Structural Uncertainty and Intervention in Climate Attitude Networks

Code for a project examining how climate-related attitudes, beliefs, and context variables are
connected, which of those connections have a data-supported causal direction, and what that means
for where an intervention would actually have leverage on later climate-related behavior. Uses five
waves of a U.S. survey (N = 1,987 for the attitude network, N = 870 for the behavior outcome).

The analysis combines a Gaussian graphical model (descriptive network), constraint-based causal
discovery (FCI, with PC-stable as a sensitivity check) with bootstrap stability analysis, a fully
directed working structural causal model fit in `lavaan`, and several layers of sensitivity analysis
(directional-completion, latent-confounding-compatible, and a PAG-compatible LV-IDA analysis that
avoids committing to one directed completion at all).

## Repository structure

- `clean_pipeline/` -- the confirmatory analysis, pulled out of the exploratory `.qmd` into scripts
  that run top to bottom in a fresh R session. Start here; see its own README for run order and
  what each script does.
- `r_patches/` -- figure scripts and additional sensitivity/diagnostic analyses.
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

The raw participant-level survey data is not included in this repository. `clean_pipeline/01_data_prep.R`
expects it in a `data_henry/` folder (parquet files); that folder is gitignored since it's
individual-level survey data with demographic fields, not something to publish alongside the code.
Get in touch for access to the raw data if you want to reproduce the pipeline end to end.

## Reproducing the analysis

From the repository root, in a fresh R session:

```r
source("clean_pipeline/run_all.R")
```

This runs `clean_pipeline/01` through `14` in order and freezes the result (with a manifest and a
baseline sanity check) into `pipeline_outputs/reference_runs/`. Scripts `15`-`18` are optional
follow-on sensitivity analyses; see `clean_pipeline/README.md` for details and dependencies. Figure
and supplementary-analysis scripts in `r_patches/` are run separately once `clean_pipeline/`'s output
is in place.

## Requirements

R packages: `dplyr`, `tidyr`, `readr`, `huge`, `qgraph`, `pcalg`, `lavaan`, `ggplot2`, and (for the
sensitivity analyses in `r_patches/14`-`16`) `RCIT`/`RCoT` and `micd`.
