# =============================================================================
# 22_intervention_bootstrap_ci.R
#
# Purpose (methodological review, item #3 / Blocker 3): quantify sampling
# uncertainty in the fitted path coefficients and propagate it into the
# intervention-effect estimates, kept explicitly separate from the
# structural-orientation uncertainty already handled by script 02.
#
# Design, per team decision:
#   - Participant bootstrap only, baseline orientation only (the 16-edge
#     working SCM as specified, no edge flips). Do NOT nest the 8 orientation
#     specifications inside every resample -- these are two different
#     uncertainties (sampling/parameter vs. structural) and are kept apart.
#   - B = 1000 resamples of participants (N = 870), refit the SCM in each,
#     recompute ALL intervention targets: 7 singles + 21 pairs + 35 triples
#     (63 combo targets over the 7 non-political attitude nodes) PLUS a
#     standalone political-orientation single-node target (not part of any
#     pair/triple, since politics is a background covariate rather than a
#     plausible intervention target) -- 64 targets total, giving bootstrap
#     CIs for all 8 attitude/context nodes' single-node effects.
#     REVISED 2026-09-07: extended from the original 6-node combo search
#     (41 targets) to 7 nodes (harm_future added back in) per author
#     request, so single-node bootstrap CIs now cover all 8 nodes and the
#     pair/triple search is consistent with Figure 4/Table 2 showing all 8
#     single-node effects. best_pair_label/best_triple_label below are now
#     chosen dynamically from the point estimates rather than hardcoded,
#     since harm_future's addition can change which pair/triple wins.
#   - Because the model is linear and recursive, the mean intervention effect
#     is propagated DETERMINISTICALLY through the fitted structural
#     equations for each resample (no within-resample Monte Carlo draws --
#     the exogenous/residual noise is mean-zero by construction, so a
#     20,000-draw simulation and the deterministic mean agree in expectation
#     and only the deterministic version is free of simulation noise). The
#     existing Monte Carlo simulate_scm_generic() pipeline (script 02) is
#     untouched and can still serve as a validation check against this.
#   - Report 95% percentile CIs for all eight single-node ATEs (seven
#     combo-eligible nodes plus politics) and for the best pair and best
#     triple, now identified dynamically from the point estimates on the
#     real sample (see step below) rather than hardcoded to the old 6-node
#     winners. Because "best pair"/"best triple" are selected from the same
#     data, also report how often each candidate pair/triple ranks first
#     across the bootstrap resamples -- a more honest summary than attaching
#     an ordinary CI to a winner chosen once.
#
# HOW TO USE: run after the .qmd has df_extended (N = 870 Wave-5 subsample)
# in the session, same precondition as script 02. Self-contained otherwise
# (own copy of base_edges / helper functions).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lavaan)
  library(igraph)
  library(purrr)
  library(tibble)
  library(furrr)
  library(future)
})

stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

set.seed(2026)
N_BOOT <- 1000

# --- Baseline 16-edge working SCM (unchanged; identical to script 02's
# base_edges, kept as its own copy so this script has no dependency on
# script 02 having been sourced first). ---------------------------------------
base_edges <- tibble::tribble(
  ~from,               ~to,
  "politics",          "belief_concern",
  "belief_concern",    "harm_future",
  "belief_concern",    "harm_present",
  "harm_future",       "harm_present",
  "belief_concern",    "trust_science",
  "harm_future",       "trust_science",
  "belief_concern",    "policy_support",
  "trust_science",     "policy_support",
  "politics",          "policy_support",
  "policy_support",    "social_norms",
  "trust_science",     "social_norms",
  "harm_present",      "weather_risk_prep",
  "belief_concern",    "weather_risk_prep",
  "harm_present",      "climate_behavior",
  "social_norms",      "climate_behavior",
  "weather_risk_prep", "climate_behavior"
)
stopifnot(nrow(base_edges) == 16)

all_nodes <- unique(c(base_edges$from, base_edges$to))
# REVISED 2026-09-07: harm_future added back into the combo-eligible node
# set (7 nodes now), excluding only politics (a background ideological
# covariate, not a plausible intervention target -- see below for its
# separate standalone single-node CI).
intervene_nodes <- c("belief_concern", "harm_present", "harm_future",
                      "weather_risk_prep", "social_norms", "trust_science",
                      "policy_support")

build_lavaan_syntax <- function(edges) {
  edges |>
    dplyr::group_by(to) |>
    dplyr::summarise(rhs = paste(from, collapse = " + "), .groups = "drop") |>
    dplyr::mutate(line = paste(to, "~", rhs)) |>
    dplyr::pull(line) |>
    paste(collapse = "\n")
}

topo_order <- function(edges) {
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  igraph::topo_sort(g, mode = "out") |> names()
}

build_intervene_targets <- function(nodes, max_size = 3) {
  targets <- list()
  for (k in seq_len(max_size)) {
    for (combo in utils::combn(nodes, k, simplify = FALSE)) {
      targets[[paste(combo, collapse = "+")]] <- combo
    }
  }
  targets
}
intervene_targets <- build_intervene_targets(intervene_nodes, max_size = 3)
stopifnot(length(intervene_targets) == choose(7, 1) + choose(7, 2) + choose(7, 3))  # 7+21+35 = 63

# Politics is added AFTER combo generation, as its own length-1 target only
# -- it is never combined with anything into a pair/triple, but still gets
# the same deterministic-propagation + bootstrap CI treatment as every other
# single-node target, via the exact same target_ate()/boot_one() machinery.
intervene_targets[["politics"]] <- "politics"
stopifnot(length(intervene_targets) == 64)

# --- Deterministic mean propagation (replaces Monte Carlo for this script
# only). For a linear-Gaussian recursive SCM, E[node] under a given
# intervention assignment is exactly sum(beta_parent * E[parent]) for
# unintervened nodes (residual noise is mean-zero by construction) and the
# fixed value for intervened nodes; exogenous nodes have standardized mean 0.
scm_mean_propagate <- function(fit, edges, all_nodes, intervene = list()) {
  ord <- union(topo_order(edges), all_nodes)
  std <- lavaan::standardizedSolution(fit) |> dplyr::filter(op == "~")
  means <- setNames(numeric(length(ord)), ord)
  for (node in ord) {
    if (!is.null(intervene[[node]])) {
      means[[node]] <- intervene[[node]]
      next
    }
    parents <- edges$from[edges$to == node]
    if (length(parents) == 0) {
      means[[node]] <- 0
      next
    }
    coefs <- std |> dplyr::filter(lhs == node, rhs %in% parents)
    beta <- setNames(coefs$est.std, coefs$rhs)
    means[[node]] <- sum(vapply(parents, function(p) beta[[p]] * means[[p]], numeric(1)))
  }
  means
}

target_ate <- function(fit, edges, all_nodes, target_combo) {
  baseline_cb <- scm_mean_propagate(fit, edges, all_nodes, intervene = list())[["climate_behavior"]]
  do_cb <- scm_mean_propagate(
    fit, edges, all_nodes,
    intervene = setNames(as.list(rep(0.5, length(target_combo))), target_combo)
  )[["climate_behavior"]]
  do_cb - baseline_cb
}

# --- Point estimate on the real (unresampled) sample, for reference --------
model_syntax_pt <- build_lavaan_syntax(base_edges)
fit_pt <- lavaan::sem(model_syntax_pt, data = df_extended, estimator = "MLR", fixed.x = FALSE)
stopifnot(lavaan::lavInspect(fit_pt, "converged"))
point_ate <- purrr::imap_dbl(intervene_targets, ~ target_ate(fit_pt, base_edges, all_nodes, .x))
cat("Point-estimate ATEs (deterministic propagation, real sample) for all eight singles:\n")
print(round(point_ate[c(intervene_nodes, "politics")], 4))

# --- Bootstrap loop (parallelized across resamples with furrr) -------------
# Each resample sets its own seed (2026 + b) inside the worker, independent of
# execution order or worker count, so results are identical whether this runs
# sequentially or in parallel and however many workers are used.
n_row <- nrow(df_extended)

boot_one <- function(b) {
  set.seed(2026 + b)
  idx <- sample.int(n_row, n_row, replace = TRUE)
  df_b <- df_extended[idx, ]

  fit_b <- tryCatch(
    lavaan::sem(model_syntax_pt, data = df_b, estimator = "MLR", fixed.x = FALSE),
    error = function(e) NULL
  )
  converged <- !is.null(fit_b) &&
    isTRUE(tryCatch(lavaan::lavInspect(fit_b, "converged"), error = function(e) FALSE))
  if (!converged) {
    return(setNames(rep(NA_real_, length(intervene_targets)), names(intervene_targets)))
  }
  purrr::map_dbl(intervene_targets, ~ target_ate(fit_b, base_edges, all_nodes, .x))
}

plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))
boot_list <- furrr::future_map(
  seq_len(N_BOOT), boot_one,
  .options = furrr::furrr_options(seed = TRUE),
  .progress = TRUE
)
plan(sequential)

boot_ate <- do.call(rbind, boot_list)
rownames(boot_ate) <- NULL
n_converged <- sum(!is.na(boot_ate[, 1]))

cat("\nBootstrap complete:", n_converged, "/", N_BOOT, "resamples converged.\n")
if (n_converged < N_BOOT) {
  cat("(", N_BOOT - n_converged, "resamples were skipped for non-convergence/fit failure ",
      "and are NA in boot_ate -- CIs below use only converged rows.)\n", sep = "")
}

# --- 95% percentile CIs for all eight singles + best pair + best triple ----
# REVISED 2026-09-07: best_pair_label/best_triple_label are now chosen
# dynamically from the real-sample point estimates (point_ate), rather than
# hardcoded to the old 6-node winners -- with harm_future now combo-eligible,
# a pair/triple involving it may outrank the previous best_pair/best_triple.
pair_labels_pt   <- names(intervene_targets)[purrr::map_int(intervene_targets, length) == 2]
triple_labels_pt <- names(intervene_targets)[purrr::map_int(intervene_targets, length) == 3]
best_pair_label   <- pair_labels_pt[which.max(point_ate[pair_labels_pt])]
best_triple_label <- triple_labels_pt[which.max(point_ate[triple_labels_pt])]
cat("\nBest pair by real-sample point estimate:", best_pair_label,
    "(", round(point_ate[[best_pair_label]], 4), ")\n")
cat("Best triple by real-sample point estimate:", best_triple_label,
    "(", round(point_ate[[best_triple_label]], 4), ")\n")

ci_targets <- c(intervene_nodes, "politics", best_pair_label, best_triple_label)

ci_table <- purrr::map_dfr(ci_targets, function(tgt) {
  vals <- boot_ate[, tgt]
  vals <- vals[!is.na(vals)]
  tibble::tibble(
    target = tgt,
    point_estimate = round(point_ate[[tgt]], 4),
    boot_mean      = round(mean(vals), 4),
    ci_lower_95    = round(quantile(vals, .025, names = FALSE), 4),
    ci_upper_95    = round(quantile(vals, .975, names = FALSE), 4),
    n_boot         = length(vals)
  )
})

cat("\n=== 95% bootstrap CIs (participant resampling, baseline orientation only) ===\n")
print(as.data.frame(ci_table), row.names = FALSE)

# --- How often each pair / triple ranks first across bootstrap resamples ---
pair_labels   <- names(intervene_targets)[purrr::map_int(intervene_targets, length) == 2]
triple_labels <- names(intervene_targets)[purrr::map_int(intervene_targets, length) == 3]

rank1_frequency <- function(labels) {
  sub <- boot_ate[, labels, drop = FALSE]
  ok  <- stats::complete.cases(sub)
  sub <- sub[ok, , drop = FALSE]
  winner <- labels[apply(sub, 1, which.max)]
  tibble::tibble(target = labels) |>
    dplyr::left_join(
      tibble::tibble(target = winner) |> dplyr::count(target, name = "times_ranked_first"),
      by = "target"
    ) |>
    dplyr::mutate(
      times_ranked_first = tidyr::replace_na(times_ranked_first, 0L),
      pct_ranked_first    = round(100 * times_ranked_first / nrow(sub), 1)
    ) |>
    dplyr::arrange(dplyr::desc(pct_ranked_first))
}

pair_rank_table   <- rank1_frequency(pair_labels)
triple_rank_table <- rank1_frequency(triple_labels)

cat("\n=== How often each pair ranks first across bootstrap resamples ===\n")
print(as.data.frame(pair_rank_table), row.names = FALSE)
cat("\n=== How often each triple ranks first across bootstrap resamples ===\n")
print(as.data.frame(triple_rank_table), row.names = FALSE)

# --- Write outputs -----------------------------------------------------------
dir.create("pipeline_outputs", showWarnings = FALSE)
write.csv(ci_table,          "pipeline_outputs/intervention_bootstrap_ci.csv", row.names = FALSE)
write.csv(pair_rank_table,   "pipeline_outputs/intervention_bootstrap_pair_rank1_freq.csv", row.names = FALSE)
write.csv(triple_rank_table, "pipeline_outputs/intervention_bootstrap_triple_rank1_freq.csv", row.names = FALSE)
saveRDS(boot_ate, "pipeline_outputs/intervention_bootstrap_ate_matrix.rds")

cat("\nWrote pipeline_outputs/intervention_bootstrap_ci.csv, ",
    "intervention_bootstrap_pair_rank1_freq.csv, ",
    "intervention_bootstrap_triple_rank1_freq.csv, and the full ",
    "boot_ate matrix (intervention_bootstrap_ate_matrix.rds) for any further ",
    "inspection.\n", sep = "")
