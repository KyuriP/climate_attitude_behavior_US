# diffs pipeline_outputs (cvcc4_should) against pipeline_outputs_alt_will
# (cvcc4_will), read-only, nothing written to disk. run from belief_network_US/
# after running as much of clean_pipeline_alt_will/01..14 as you want:
#
#   source("clean_pipeline_alt_will/compare_to_primary.R")
#
# joins each pair of files on their key columns, rounds to ROUND_DIGITS, and
# flags anything that moves past TOL, or any categorical column
# (orientation_label etc) that changes at all. rounding first because a plain
# sort+diff on these csvs makes everything look different just from
# floating-point print precision, not real changes.

PRIMARY_DIR <- "pipeline_outputs"
ALT_DIR     <- "pipeline_outputs_alt_will"
ROUND_DIGITS <- 4
TOL <- 5 * 10^(-ROUND_DIGITS)   # 0.0005 at 4 digits -- one rounding-unit of slack

if (!dir.exists(PRIMARY_DIR)) stop("Can't find '", PRIMARY_DIR, "' -- run this from belief_network_US/.")
if (!dir.exists(ALT_DIR)) stop(
  "Can't find '", ALT_DIR, "' -- the alt-will pipeline hasn't produced any output yet. ",
  "Run clean_pipeline_alt_will/01_data_prep.R onward first."
)

.hr <- function(ch = "-") cat(strrep(ch, 78), "\n")

# reads rel_path from both dirs, joins on key_cols, rounds numeric_cols and
# flags diffs past TOL, flags any change at all in categorical_cols.
# numeric_cols=NULL auto-detects (any col that's numeric in both files).
compare_csv <- function(rel_path, key_cols, numeric_cols = NULL, categorical_cols = NULL) {
  .hr("=")
  cat(rel_path, "\n")
  .hr("=")

  p_path <- file.path(PRIMARY_DIR, rel_path)
  a_path <- file.path(ALT_DIR, rel_path)

  if (!file.exists(p_path)) { cat("  [skip] not found in primary:", p_path, "\n\n"); return(invisible(NULL)) }
  if (!file.exists(a_path)) { cat("  [skip] not found in alt-will yet:", a_path, "\n\n"); return(invisible(NULL)) }

  df_p <- utils::read.csv(p_path, stringsAsFactors = FALSE)
  df_a <- utils::read.csv(a_path, stringsAsFactors = FALSE)

  missing_key_p <- setdiff(key_cols, names(df_p))
  missing_key_a <- setdiff(key_cols, names(df_a))
  if (length(missing_key_p) || length(missing_key_a)) {
    cat("  [skip] key column(s) missing -- schema changed. primary missing:",
        paste(missing_key_p, collapse = ","), " alt missing:",
        paste(missing_key_a, collapse = ","), "\n\n")
    return(invisible(NULL))
  }

  if (is.null(numeric_cols)) {
    candidate <- setdiff(intersect(names(df_p), names(df_a)), key_cols)
    numeric_cols <- candidate[vapply(df_p[candidate], is.numeric, logical(1)) &
                                 vapply(df_a[candidate], is.numeric, logical(1))]
  }
  numeric_cols <- intersect(numeric_cols, intersect(names(df_p), names(df_a)))
  categorical_cols <- intersect(categorical_cols, intersect(names(df_p), names(df_a)))

  # if a column I listed as numeric turns out not to be (like current_evidence,
  # which is a text label) fall back to comparing it as categorical instead
  # of crashing round()
  is_actually_numeric <- vapply(df_p[numeric_cols], is.numeric, logical(1)) &
    vapply(df_a[numeric_cols], is.numeric, logical(1))
  not_numeric <- numeric_cols[!is_actually_numeric]
  if (length(not_numeric)) {
    cat("  [note] not numeric in both files -- comparing as categorical instead:",
        paste(not_numeric, collapse = ", "), "\n")
    categorical_cols <- union(categorical_cols, not_numeric)
  }
  numeric_cols <- numeric_cols[is_actually_numeric]

  merged <- merge(df_p, df_a, by = key_cols, all = TRUE, suffixes = c("_primary", "_alt"))

  only_p <- merged[Reduce(`|`, lapply(numeric_cols, function(c) is.na(merged[[paste0(c, "_alt")]]))), , drop = FALSE]
  only_a <- merged[Reduce(`|`, lapply(numeric_cols, function(c) is.na(merged[[paste0(c, "_primary")]]))), , drop = FALSE]

  if (nrow(only_p) > 0) {
    cat("  ROWS ONLY IN PRIMARY (", nrow(only_p), "):\n", sep = "")
    print(only_p[, key_cols, drop = FALSE], row.names = FALSE)
  }
  if (nrow(only_a) > 0) {
    cat("  ROWS ONLY IN ALT-WILL (", nrow(only_a), "):\n", sep = "")
    print(only_a[, key_cols, drop = FALSE], row.names = FALSE)
  }

  flagged_any <- FALSE
  n_compared <- 0

  for (col in numeric_cols) {
    cp <- merged[[paste0(col, "_primary")]]
    ca <- merged[[paste0(col, "_alt")]]
    ok <- !is.na(cp) & !is.na(ca)
    n_compared <- max(n_compared, sum(ok))
    if (!any(ok)) next
    diff <- round(ca[ok], ROUND_DIGITS) - round(cp[ok], ROUND_DIGITS)
    flag <- abs(diff) > TOL
    if (any(flag)) {
      flagged_any <- TRUE
      cat("\n  [", col, "] ", sum(flag), " row(s) differ by more than ", TOL, ":\n", sep = "")
      out <- data.frame(
        merged[ok, key_cols, drop = FALSE][flag, , drop = FALSE],
        primary = round(cp[ok][flag], ROUND_DIGITS),
        alt_will = round(ca[ok][flag], ROUND_DIGITS),
        diff = round(diff[flag], ROUND_DIGITS)
      )
      print(out, row.names = FALSE)
    }
  }

  for (col in categorical_cols) {
    cp <- merged[[paste0(col, "_primary")]]
    ca <- merged[[paste0(col, "_alt")]]
    ok <- !is.na(cp) & !is.na(ca)
    flag <- ok & (cp != ca)
    if (any(flag)) {
      flagged_any <- TRUE
      cat("\n  [", col, "] STRUCTURAL CHANGE in ", sum(flag), " row(s):\n", sep = "")
      out <- data.frame(
        merged[flag, key_cols, drop = FALSE],
        primary = cp[flag],
        alt_will = ca[flag]
      )
      print(out, row.names = FALSE)
    }
  }

  if (!flagged_any) {
    cat("  OK -- ", n_compared, " rows compared on ", length(numeric_cols),
        " numeric column(s); no difference exceeds ", TOL,
        if (length(categorical_cols)) paste0("; no change in ", paste(categorical_cols, collapse = ", ")) else "",
        ".\n", sep = "")
  }
  cat("\n")
  invisible(merged)
}

# ============================================================================
# The comparisons
# ============================================================================

cat("Comparing", PRIMARY_DIR, "(cvcc4_should) vs", ALT_DIR, "(cvcc4_will)\n")
cat("Rounding to", ROUND_DIGITS, "decimals; flagging diffs >", TOL, "\n\n")

# --- SCM structure: edges, existence/asymmetry, and orientation/tier labels ---
# This is the first and most important check -- see README's "what to actually
# look at". A flagged orientation_label/final_tier change means the SCM's
# causal structure itself would differ under cvcc4_will, which matters far
# more than any downstream ATE shifting a little.
compare_csv(
  "scm_edges_finalized.csv",
  key_cols = c("from", "to"),
  numeric_cols = c("p_adjacent", "avg_n", "asymmetry"),
  categorical_cols = c("current_evidence", "existence_band", "orientation_label", "final_tier")
)

# --- Single-node / 8-node deterministic ATEs (Figure 7 data) ---
compare_csv(
  "tables/orientation_enumeration_ate_deterministic_8node.csv",
  key_cols = c("scenario", "node"),
  numeric_cols = "ate_climate_behavior"
)

# --- Combo (pair/triple) deterministic ATEs ---
compare_csv(
  "tables/orientation_enumeration_combo_ate_deterministic.csv",
  key_cols = c("scenario", "target_label"),
  numeric_cols = "ate_climate_behavior"
)

# --- Bootstrap CIs on the intervention ATEs ---
compare_csv(
  "intervention_bootstrap_ci.csv",
  key_cols = "target",
  numeric_cols = c("point_estimate", "boot_mean", "ci_lower_95", "ci_upper_95")
)
compare_csv(
  "intervention_bootstrap_pair_rank1_freq.csv",
  key_cols = "target",
  numeric_cols = c("times_ranked_first", "pct_ranked_first")
)
compare_csv(
  "intervention_bootstrap_triple_rank1_freq.csv",
  key_cols = "target",
  numeric_cols = c("times_ranked_first", "pct_ranked_first")
)

# --- Interaction/moderation (Study 1's S9 numbers) ---
compare_csv(
  "tables/interaction_coefficient.csv",
  key_cols = "version",
  numeric_cols = "beta"
)
compare_csv(
  "tables/interaction_shift_results_exact.csv",
  key_cols = "scenario",
  numeric_cols = c("belief_concern_mean", "policy_support_mean", "climate_behavior_mean")
)
compare_csv(
  "tables/interaction_shift_results_montecarlo.csv",
  key_cols = "scenario",
  numeric_cols = c("belief_concern_mean", "policy_support_mean", "climate_behavior_mean")
)

# --- Attrition / IPW sensitivity ---
compare_csv(
  "attrition_waves1to4_comparison.csv",
  key_cols = "variable",
  numeric_cols = c("mean_retained", "mean_attrited", "smd")
)
compare_csv(
  "ipw_balance_check.csv",
  key_cols = "variable",
  numeric_cols = c("smd_unweighted", "smd_ipw_weighted")
)
compare_csv(
  "ipw_single_node_ate_compare.csv",
  key_cols = "node",
  numeric_cols = c("ate_unweighted", "ate_ipw_weighted", "diff")
)

# --- Behavior-outcome sensitivity (mitig4 replacement check) ---
compare_csv(
  "behavior_sensitivity_ate.csv",
  key_cols = c("node", "outcome"),
  numeric_cols = "ate"
)

.hr("=")
cat("Done. Any section above printing rows (rather than 'OK') is where cvcc4_will\n")
cat("moves something beyond rounding noise -- that's what to bring to Sara.\n")
.hr("=")
