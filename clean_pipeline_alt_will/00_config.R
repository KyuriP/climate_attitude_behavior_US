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

# ---- scm edge-inclusion rule ----
EXISTENCE_MIN <- 0.60
ORIENTATION_ASYMMETRY_EPS <- 0.01

# ---- input data ----
RAW_DATA_DIR <- "data_henry"

# ---- paths ----
OUTPUT_DIR <- "pipeline_outputs_alt_will"  # keeps this fully separate from the real pipeline_outputs
MANIFEST_PATH <- file.path(OUTPUT_DIR, "bootstrap_manifest.csv")
LATEST_POINTER_PATH <- file.path(OUTPUT_DIR, "LATEST_bootstrap_run.csv")

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# separate bootstrap from N_BOOT above -- resamples participants (N=870) and
# refits the scm path coefficients, for the intervention-effect cis (results:
# "1,000 participant bootstrap resamples", rank frequencies). don't confuse
# with N_BOOT, they're unrelated numbers.
N_BOOT_INTERVENTION <- 1000
BOOT_SEED_INTERVENTION <- 2026
