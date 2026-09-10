# resume helper for run_all.R -- for exactly the situation where 01-09
# already ran successfully in this live session (maybe across an interrupt
# and a manual re-source of one script, like 09's parallel-hang fix on
# 2026-09-08) and it's wasteful to redo them just to get run_all.R's
# freeze/manifest step at the end. source 10 through 14 yourself first, in
# order, then source this -- it's the exact tail of run_all.R (the freeze +
# baseline check + manifest), unchanged, just split out so it doesn't
# require re-running everything above it.
#
# needs run_start_time to already exist (set at the top of run_all.R, and
# untouched by an interrupt/resume since interrupting only stops the
# in-progress script, it doesn't clear the environment).

if (!exists("run_start_time")) {
  stop("run_start_time not found -- this expects you to have started with ",
       "run_all.R (which sets it before sourcing 01), not a from-scratch ",
       "session. if you really are starting fresh, just run run_all.R itself ",
       "instead of this file.")
}

run_end_time <- Sys.time()
cat("\n=== reference run finished", format(run_end_time, "%Y-%m-%d %H:%M:%S %Z"),
    "( elapsed since run_all.R started:", format(run_end_time - run_start_time), ") ===\n")

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
# corrected 2026-09-08 -- see run_all.R's own comment on this constant for
# why these specific numbers, not the earlier .2022/.2157 set.

expected_baseline <- c(
  belief_concern = 0.2059, harm_present = 0.1951, harm_future = 0.1322,
  weather_risk_prep = 0.0928, politics = 0.1109, social_norms = 0.0571,
  trust_science = 0.0342, policy_support = 0.0162
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
  "note: 09 was interrupted once mid-run (future::multisession hang on macOS)",
  "and re-sourced standalone after being fixed to run sequentially, then",
  "10-14 were sourced individually and this freeze step run manually via",
  "finish_reference_freeze.R rather than run_all.R end to end. same session",
  "throughout (no restart), so df_extended/base_edges etc. are consistent",
  "across every step -- just not literally one unbroken source() call.",
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
