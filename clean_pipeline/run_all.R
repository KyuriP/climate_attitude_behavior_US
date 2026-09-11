# runs 01 through 14 end to end, in one fresh R session, and freezes the
# result as a dated reference run -- exact output files this run produced
# (copied by mtime, not a hand-typed filename list, so nothing can be missed
# or stale-pulled from an older session), sessionInfo(), the exact
# scm_edges_finalized.csv this run generated, and a check of the 8 single-
# node baseline ATEs against the values locked in on 2026-09-08
# (.2022/.2157/.1399/.0921/.1080/.0481/.0115/.0077 for belief_concern/
# harm_present/harm_future/weather_risk_prep/politics/social_norms/
# trust_science/policy_support).
#
# why this exists: we've now hit three separate stale-artifact bugs from
# mixing outputs across different pipeline states or sessions (the S9
# interaction numbers, the "24 of 28 edges" prose, and 06's own hardcoded
# baseline comment) -- run_all.R is the fix. one script, one uninterrupted
# session, one frozen bundle with its own provenance record, so "which run
# did this number actually come from" stops being a question we have to
# reconstruct after the fact.
#
# run this in a FRESH R session -- restart R first. don't source it into a
# session that's already run 15/16/17 or any other script; the whole point
# is that this bundle isn't contaminated by leftover state from something
# else. run from the repo root (belief_network_US/), same as every other
# script here.

if (exists("confound_candidates") || exists("fci_ext_05") || exists("classified")) {
  stop("this session already has objects from 15/16/17 (confound_candidates/",
       "fci_ext_05/classified) -- restart R to a fresh session before running ",
       "run_all.R, so the frozen bundle can't inherit any leftover state.")
}

if (!file.exists("clean_pipeline/00_config.R")) {
  stop("clean_pipeline/00_config.R not found from the current working directory -- ",
       "run this from the belief_network_US/ repo root, same as every other script here.")
}

source("clean_pipeline/00_config.R")

run_start_time <- Sys.time()
cat("=== clean_pipeline reference run starting",
    format(run_start_time, "%Y-%m-%d %H:%M:%S %Z"), "===\n")

scripts_in_order <- c(
  "01_data_prep.R",
  "02_ggm.R",
  "03_bootstrap_causal_discovery.R",
  "04_scm_finalize.R",
  "06_intervention_ates_singlenode.R",
  "07_intervention_ates_8node.R",
  "08_intervention_ates_combo.R",
  "09_intervention_bootstrap_ci.R",
  "10_orientation_enumeration_fit.R",
  # "11_interaction_moderation.R" -- PERMANENTLY REMOVED 2026-09-11. The
  # belief_concern x politics interaction/moderation analysis (and Supp.
  # Section S9) was dropped from the paper entirely: politics is no longer
  # exogenous to belief_concern after the edge reversal, and policy_support
  # is a sink node so the interaction had no downstream route to
  # climate_behavior anyway. See analysis_decisions_log.md Section 38. The
  # script itself is now a stub explaining this; do not re-add it here.
  "12_wave5_attrition_check.R",
  "13_ipw_attrition_sensitivity.R",
  "14_behavior_outcome_sensitivity.R"
)
# 05 is sourced automatically by 06 onward, not listed here -- matches the
# documented run order in clean_pipeline/README.md.

for (s in scripts_in_order) {
  cat("\n--- sourcing", s, "---\n")
  source(file.path("clean_pipeline", s))
}

run_end_time <- Sys.time()
cat("\n=== reference run finished", format(run_end_time, "%Y-%m-%d %H:%M:%S %Z"),
    "( elapsed:", format(run_end_time - run_start_time), ") ===\n")

# ---- freeze: dated snapshot dir, only files this run actually produced ----

run_tag <- format(run_start_time, "%Y%m%d_%H%M%S")
freeze_dir <- file.path(OUTPUT_DIR, "reference_runs", run_tag)
dir.create(file.path(freeze_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

freeze_by_mtime <- function(dir_path, dest_path, since) {
  files <- list.files(dir_path, full.names = TRUE)
  info <- file.info(files)
  fresh <- files[!info$isdir & info$mtime >= since]
  if (length(fresh) > 0) file.copy(fresh, dest_path, overwrite = TRUE)
  basename(fresh)
}

top_level_frozen <- freeze_by_mtime(OUTPUT_DIR, freeze_dir, run_start_time)
tables_frozen <- freeze_by_mtime(file.path(OUTPUT_DIR, "tables"), file.path(freeze_dir, "tables"), run_start_time)

cat("\nFroze", length(top_level_frozen), "top-level output files and",
    length(tables_frozen), "table files into", freeze_dir, "\n")

# ---- verify the 8 baseline single-node ATEs against the locked values -----

# corrected 2026-09-08: the .2022/.2157/etc this originally had came from a
# session that turned out to be contaminated (stale df_extended, not a real
# spec difference) -- a genuinely clean run_all.R execution reproduces these
# values instead, confirmed two independent ways in the same run: 06's exact
# mean-propagation baseline AND 09's participant-bootstrap point estimate
# agree to 4 decimals, and 09's bootstrapped CI (.1832-.2278 / .1613-.2288)
# matches 09's own original sanity-check comment, which nobody had touched.
# UPDATED 2026-09-11 after the politics/belief_concern and policy_support/
# social_norms edge reversal (see analysis_decisions_log.md Section 34/35) --
# these replace the pre-reversal 2026-09-08 values, which are now stale.
# Cross-validated two independent ways on the post-reversal model: 07's
# exact mean-propagation combo_0 output and 09's participant-bootstrap point
# estimate agree to 4+ decimals for all 8 nodes (see chat/log). Note
# politics drops from .1109 to exactly 0 -- a real, expected consequence of
# the reversal: politics is no longer upstream of belief_concern, and its
# one remaining outgoing edge (politics -> policy_support) dead-ends at
# policy_support, which itself has no outgoing edges in this 16-edge SCM.
expected_baseline <- c(
  belief_concern = 0.201272438461644, harm_present = 0.186965335962417,
  harm_future = 0.129295277380915, weather_risk_prep = 0.0930759991784261,
  politics = 0, social_norms = 0.0633652634284836,
  trust_science = 0.0436997792836055, policy_support = 0
)

baseline_8node_path <- file.path(OUTPUT_DIR, "tables", "orientation_enumeration_ate_deterministic_8node.csv")
stopifnot(file.exists(baseline_8node_path))
baseline_8node <- utils::read.csv(baseline_8node_path, stringsAsFactors = FALSE)
actual_row <- baseline_8node[baseline_8node$scenario == "combo_0", c("node", "ate_climate_behavior")]
actual_full <- setNames(actual_row$ate_climate_behavior, actual_row$node)[names(expected_baseline)]

baseline_compare <- data.frame(
  node = names(expected_baseline),
  expected = unname(expected_baseline),
  actual = unname(round(actual_full, 4)),
  actual_full_precision = unname(actual_full),
  match = unname(abs(actual_full - expected_baseline) < 5e-5)
)

cat("\n=== baseline reference-run check: 8 single-node ATEs vs the 2026-09-08 locked values ===\n")
print(baseline_compare, row.names = FALSE)

all_match <- all(baseline_compare$match)
if (!all_match) {
  warning("REFERENCE RUN MISMATCH: one or more baseline ATEs did NOT match the locked ",
          "2026-09-08 values. Do NOT treat this run as frozen/authoritative -- see ",
          "baseline_compare above and ", file.path(freeze_dir, "reference_run_manifest.txt"),
          " -- figure out why before running 15/16/17 on top of it.")
} else {
  cat("\nAll 8 baseline ATEs match the locked 2026-09-08 values (within 5e-5). ",
      "This run is clean -- treat ", freeze_dir, " as the frozen reference bundle.\n", sep = "")
}

# ---- manifest: timestamp, sessionInfo(), the finalized-edges file used ----

scm_path <- file.path(OUTPUT_DIR, "scm_edges_finalized.csv")
scm_info <- if (file.exists(scm_path)) {
  c(paste("  path:", scm_path),
    paste("  mtime:", format(file.info(scm_path)$mtime)),
    paste("  md5:", unname(tools::md5sum(scm_path))))
} else {
  "  NOT FOUND -- 04_scm_finalize.R did not produce it this run."
}

manifest_lines <- c(
  "clean_pipeline reference run",
  paste("run started: ", format(run_start_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste("run finished:", format(run_end_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste("elapsed:     ", format(run_end_time - run_start_time)),
  "",
  "finalized 16-edge SCM file used this run:",
  scm_info,
  "",
  "baseline check (8 single-node ATEs vs 2026-09-08 locked values):",
  utils::capture.output(print(baseline_compare, row.names = FALSE)),
  "",
  paste("ALL 8 MATCH:", all_match),
  "",
  paste("files frozen (", OUTPUT_DIR, "/):", sep = ""),
  paste(" -", top_level_frozen),
  "",
  paste("files frozen (", OUTPUT_DIR, "/tables/):", sep = ""),
  paste(" -", tables_frozen),
  "",
  "sessionInfo():",
  utils::capture.output(sessionInfo())
)
manifest_path <- file.path(freeze_dir, "reference_run_manifest.txt")
writeLines(manifest_lines, manifest_path)

cat("\nWrote", manifest_path, "\n")
cat("Frozen reference bundle:", freeze_dir, "\n")
cat("(15/16/17 read pipeline_outputs/... and pipeline_outputs/tables/... directly, not\n",
    "the reference_runs/ copy -- so as long as nothing is re-run in ", OUTPUT_DIR,
    " between now and running 15/16/17, they're operating on this exact frozen state.\n",
    "the reference_runs/ copy is the audit trail, not the live input path.)\n", sep = "")
