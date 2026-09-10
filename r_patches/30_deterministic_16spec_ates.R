# =============================================================================
# Patch 30: DETERMINISTIC re-computation of the 16-specification single-node
#           ATEs, replacing script 02 v4's Monte Carlo simulation
#           (simulate_scm_generic(), n=20000 draws) with exact mean
#           propagation (scm_mean_propagate() / single_node_ates(), the same
#           deterministic helpers already used in scripts 22/23/26/27).
#
# WHY THIS SCRIPT EXISTS
# -----------------------
# Script 02 v4 enumerated all 2^4=16 orientation combinations for the four
# directionally-unresolved edges and computed six single-node ATEs per
# scenario via Monte Carlo simulation (rnorm()-based sampling, n=20000).
# Because the SCM is linear-recursive, E[node] = sum(beta_parent *
# E[parent]) is EXACT and needs no simulation at all -- mean propagation
# gives the same answer as the Monte Carlo estimate up to sampling noise
# (~1/sqrt(20000) ~ .007 SD). That noise floor is invisible for nodes with
# a real nonzero ATE, but it produced a spurious NEGATIVE lower bound for
# policy_support (-.00707) in scenarios where policy_support's only
# outgoing edge (policy_support -> social_norms) is flipped, removing all
# downstream reach and making its true ATE on climate_behavior structurally
# exactly 0. Decided to rerun the 16-spec ATEs deterministically rather than
# just footnoting the negative lower bound as Monte Carlo noise -- this
# script is that rerun.
#
# WHAT THIS SCRIPT DOES NOT CHANGE
# ----------------------------------
# base_edges, flip_candidates, flip_edges(), is_acyclic(), build_lavaan_
# syntax(), topo_order(), and the scenario_list/scenario_flip_labels
# generation are copied VERBATIM from 02_full_orientation_enumeration_v4.R
# -- the 16 fitted lavaan models are IDENTICAL to v4's (same data, same
# model syntax per scenario), so CFI/TLI/RMSEA/SRMR/AIC/BIC will not change
# by a single digit versus the already-reported v4 fit table. Only the ATE
# COMPUTATION METHOD changes (exact mean propagation vs. Monte Carlo). This
# script therefore does NOT redo the model-fit reporting (tables/
# orientation_enumeration_fit.csv / fit_supplement_table.csv, already
# correct from v4) or the 41-target combo-ATE table (not implicated by the
# policy_support issue, and comparably exact under mean propagation, but
# out of scope for this targeted fix -- worth doing that
# rerun deterministically too, later).
#
# OUTPUT
# ------
# tables/orientation_uncertainty_band_full_deterministic.csv -- the
# corrected 6-row single-node ATE range table, SAME shape/columns as v4's
# tables/orientation_uncertainty_band_full.csv, to be compared side by side
# against it. Expect all six nodes' ate_baseline/ate_min/ate_max to match
# the Monte Carlo version to within ~.01, EXCEPT policy_support's ate_min,
# which should now read exactly 0 (or a value indistinguishable from 0 to
# reporting precision) instead of -.00707.
#
# HOW TO USE
# ----------
# Run AFTER script 02 v4 (or independently -- this script is self-contained
# and refits its own 16 lavaan models; it does not depend on any object
# left behind by 02 v4). Requires df_extended (N=870, pre-Wave-5 corrected
# data) already in scope, exactly as v4 required.
# =============================================================================

library(dplyr)
library(tidyr)
library(lavaan)
library(igraph)
library(purrr)
library(tibble)

dir.create("tables", showWarnings = FALSE)

# --- 0: sanity check ---------------------------------------------------------
stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

# --- 1: baseline 16-edge working-SCM skeleton -- VERBATIM from 02 v4 --------
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

flip_candidates <- tibble::tribble(
  ~edge_label,                          ~from,            ~to,
  "politics -> belief_concern",          "politics",       "belief_concern",
  "policy_support -> social_norms",      "policy_support",  "social_norms",
  "social_norms -> climate_behavior",    "social_norms",    "climate_behavior",
  "harm_future -> harm_present",         "harm_future",     "harm_present"
)

n_flip <- nrow(flip_candidates)
stopifnot(n_flip == 4)

# --- 2: helpers --------------------------------------------------------------
# flip_edges / is_acyclic / build_lavaan_syntax / topo_order: VERBATIM from
# 02 v4. scm_mean_propagate / single_node_ates: VERBATIM from
# 27_diagnose_hf_hp_edge.R (already used and verified there for the single
# baseline-vs-reversed HF/HP comparison; this script just applies the same
# two functions across all 16 scenarios instead of 2).

flip_edges <- function(edges, flip_set_idx) {
  e <- edges
  if (length(flip_set_idx) > 0) {
    for (i in flip_set_idx) {
      row <- flip_candidates[i, ]
      match_idx <- which(e$from == row$from & e$to == row$to)
      if (length(match_idx) == 1) {
        tmp <- e$from[match_idx]
        e$from[match_idx] <- e$to[match_idx]
        e$to[match_idx]   <- tmp
      }
    }
  }
  e
}

is_acyclic <- function(edges) {
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  igraph::is_dag(g)
}

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

scm_mean_propagate <- function(fit, edges, all_nodes, intervene = list()) {
  ord <- union(topo_order(edges), all_nodes)
  std <- lavaan::standardizedSolution(fit) |> dplyr::filter(op == "~")
  means <- setNames(numeric(length(ord)), ord)
  for (node in ord) {
    if (!is.null(intervene[[node]])) { means[[node]] <- intervene[[node]]; next }
    parents <- edges$from[edges$to == node]
    if (length(parents) == 0) { means[[node]] <- 0; next }
    coefs <- std |> dplyr::filter(lhs == node, rhs %in% parents)
    beta <- setNames(coefs$est.std, coefs$rhs)
    means[[node]] <- sum(vapply(parents, function(p) beta[[p]] * means[[p]], numeric(1)))
  }
  means
}

single_node_ates <- function(fit, edges, all_nodes, nodes) {
  base_cb <- scm_mean_propagate(fit, edges, all_nodes, list())[["climate_behavior"]]
  purrr::map_dbl(nodes, function(nd) {
    do_cb <- scm_mean_propagate(fit, edges, all_nodes, setNames(list(0.5), nd))[["climate_behavior"]]
    do_cb - base_cb
  }) |> setNames(nodes)
}

# --- 3: enumerate ALL 16 combinations -- VERBATIM generation logic from v4 --
scenario_list <- setNames(
  lapply(0:(2^n_flip - 1), function(k) which(as.logical(intToBits(k)[1:n_flip]))),
  paste0("combo_", 0:(2^n_flip - 1))
)

scenario_flip_labels <- purrr::imap_chr(scenario_list, function(idx, nm) {
  if (length(idx) == 0) "(none -- baseline)" else paste(flip_candidates$edge_label[idx], collapse = "; ")
})

all_nodes <- unique(c(base_edges$from, base_edges$to))
intervene_nodes <- c("belief_concern", "harm_present", "weather_risk_prep",
                      "social_norms", "trust_science", "policy_support")

# --- 4: main loop: refit + DETERMINISTIC ATEs for every acyclic scenario ----
# Same acyclicity/convergence gating as v4 (status "ok" / "cyclic_skipped" /
# "fit_failed"), so the same scenarios that were valid in v4 will be valid
# here -- this script does not change which scenarios count, only how their
# ATEs are computed.

det_out <- scenario_list |>
  purrr::imap(function(flip_idx, scen_name) {

    e <- flip_edges(base_edges, flip_idx)

    if (!is_acyclic(e)) {
      message(scen_name, " (", scenario_flip_labels[[scen_name]], "): cyclic, skipped.")
      return(list(ate = NULL, status = "cyclic_skipped"))
    }

    model_syntax <- build_lavaan_syntax(e)
    fit <- tryCatch(
      lavaan::sem(model_syntax, data = df_extended, estimator = "MLR", fixed.x = FALSE),
      error = function(err) NULL
    )
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      message(scen_name, " (", scenario_flip_labels[[scen_name]], "): fit failed/non-converged, skipped")
      return(list(ate = NULL, status = "fit_failed"))
    }

    ates <- single_node_ates(fit, e, all_nodes, intervene_nodes)
    ate_rows <- tibble::tibble(
      scenario = scen_name,
      node = intervene_nodes,
      ate_climate_behavior = unname(ates[intervene_nodes])
    )

    list(ate = ate_rows, status = "ok")
  })

det_status <- purrr::map_chr(det_out, "status")
det_results <- det_out |> purrr::map("ate") |> purrr::compact() |> dplyr::bind_rows()

cat("\n--- Scenarios attempted:", length(scenario_list),
    "| acyclic+converged:", sum(det_status == "ok"),
    "| cyclic (excluded):", sum(det_status == "cyclic_skipped"),
    "| fit failed:", sum(det_status == "fit_failed"), "---\n")
cat("(This should read identically to script 02 v4's own status line -- same",
    "16 scenarios, same acyclicity/convergence gating, only the ATE",
    "computation method differs. If the counts here do NOT match v4's",
    "reported counts, stop and report the discrepancy before trusting",
    "anything below.)\n")

readr_ok <- requireNamespace("readr", quietly = TRUE)

# --- 5: headline uncertainty band, deterministic version --------------------
det_uncertainty_band <- det_results |>
  dplyr::group_by(node) |>
  dplyr::summarise(
    ate_baseline = ate_climate_behavior[scenario == "combo_0"],
    ate_min = min(ate_climate_behavior),
    ate_max = max(ate_climate_behavior),
    range   = ate_max - ate_min,
    n_specs = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(range))

cat("\n--- Six single-node intervention ranges, DETERMINISTIC (exact mean",
    "propagation, no Monte Carlo noise) -- compare directly against",
    "tables/orientation_uncertainty_band_full.csv (v4, Monte Carlo). Expect",
    "all values to match to within ~.01 except policy_support's ate_min,",
    "which should now read exactly/near-exactly 0 instead of -.00707. ---\n")
print(det_uncertainty_band)

if (readr_ok) {
  readr::write_csv(det_uncertainty_band, "tables/orientation_uncertainty_band_full_deterministic.csv")
  readr::write_csv(det_results, "tables/orientation_enumeration_ate_deterministic.csv")
} else {
  write.csv(det_uncertainty_band, "tables/orientation_uncertainty_band_full_deterministic.csv", row.names = FALSE)
  write.csv(det_results, "tables/orientation_enumeration_ate_deterministic.csv", row.names = FALSE)
}

# --- 6: explicit side-by-side vs. the Monte Carlo version, if available -----
mc_path <- "tables/orientation_uncertainty_band_full.csv"
if (file.exists(mc_path)) {
  mc_band <- if (readr_ok) readr::read_csv(mc_path, show_col_types = FALSE) else read.csv(mc_path)
  compare_tbl <- det_uncertainty_band |>
    dplyr::select(node, det_baseline = ate_baseline, det_min = ate_min, det_max = ate_max) |>
    dplyr::inner_join(
      mc_band |> dplyr::select(node, mc_baseline = ate_baseline, mc_min = ate_min, mc_max = ate_max),
      by = "node"
    ) |>
    dplyr::mutate(
      diff_baseline = round(det_baseline - mc_baseline, 4),
      diff_min = round(det_min - mc_min, 4),
      diff_max = round(det_max - mc_max, 4)
    )
  cat("\n--- Deterministic vs. Monte Carlo (v4) side by side -- diff columns",
      "should all be small (~.01 or less), confirming this rerun is a pure",
      "noise-removal exercise and not a substantive change. policy_support's",
      "diff_min is the one expected to be meaningfully nonzero (removing the",
      "spurious -.00707). ---\n")
  print(as.data.frame(compare_tbl), row.names = FALSE)
  write.csv(compare_tbl, "tables/orientation_ate_deterministic_vs_montecarlo.csv", row.names = FALSE)
} else {
  cat("\n(", mc_path, "not found in this session -- skipping side-by-side comparison table;",
      "det_uncertainty_band above is still the corrected headline result.)\n")
}

cat("\n--- Done. Send back: (1) the 'Scenarios attempted' status line above,",
    "(2) the six-row det_uncertainty_band printout, (3) the deterministic-vs-",
    "Monte-Carlo comparison table if it printed. ---\n")
