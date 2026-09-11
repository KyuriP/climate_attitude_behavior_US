# config -- n_boot, alphas, node lists, paths. everything downstream sources
# this so nothing hardcodes these values separately again.

# causal-discovery bootstrap count -- this is what's actually behind table 2 /
# fig 2-3 / the 16-edge scm right now. don't bump this without rerunning
# 03_bootstrap_causal_discovery.R + 04_scm_finalize.R and rechecking table 2 /
# the edge set / orientation sensitivity against the new output.
N_BOOT  <- 1000
ALPHAS  <- c("0.05" = 0.05, "0.01" = 0.01)
SEED    <- 42

# ---- node sets ----
NODE_ORDER_MAIN <- c(
  "belief_concern", "harm_present", "harm_future",
  "policy_support", "trust_science", "social_norms",
  "politics", "weather_risk_prep"
)
NODE_ORDER_EXT <- c(NODE_ORDER_MAIN, "climate_behavior")

# ---- node abbreviations, used by the figure scripts (r_patches/) for axis labels ----
ABBR <- c(
  belief_concern    = "BC",
  harm_present      = "HP",
  harm_future       = "HF",
  policy_support    = "PS",
  trust_science     = "TS",
  social_norms      = "SN",
  politics          = "POL",
  weather_risk_prep = "WR"
)
ABBR_EXT <- c(ABBR, climate_behavior = "CB")

# ---- scm edge-inclusion rule ----
EXISTENCE_MIN <- 0.60

# an edge only counts as bootstrap-resolved if the pooled FCI orientation
# asymmetry clears this band in magnitude (and is positive, i.e. in the
# asserted direction). Below the band -- or negative -- the direction is
# either theory-completed on substantive grounds or, for harm_future->
# harm_present specifically, flagged uncertain because its per-alpha
# asymmetry also flips sign across the two bootstrap alpha levels (see
# r_patches/21_orientation_crossalpha_table.R and
# r_patches/27_diagnose_hf_hp_edge.R). Its pooled asymmetry (+.076) already
# falls under this .10 band on its own, so the two criteria agree for the
# current bootstrap array -- adopted 2026-09-11, replacing the old .01
# sign-only threshold (see analysis_decisions_log.md).
ORIENTATION_ASYMMETRY_EPS <- 0.10

# ---- input data ----
RAW_DATA_DIR <- "data_henry"

# ---- paths ----
OUTPUT_DIR <- "pipeline_outputs"
MANIFEST_PATH <- file.path(OUTPUT_DIR, "bootstrap_manifest.csv")
LATEST_POINTER_PATH <- file.path(OUTPUT_DIR, "LATEST_bootstrap_run.csv")

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# separate bootstrap from N_BOOT above -- resamples participants (N=870) and
# refits the scm path coefficients, for the intervention-effect cis (results:
# "1,000 participant bootstrap resamples", rank frequencies). don't confuse
# with N_BOOT, they're unrelated numbers.
N_BOOT_INTERVENTION <- 1000
BOOT_SEED_INTERVENTION <- 2026
