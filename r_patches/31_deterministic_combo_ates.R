# =============================================================================
# Patch 31: DETERMINISTIC re-computation of all 63 combo-intervention ATEs
#           (7 singles + 21 pairs + 35 triples) x 16 structural specifications,
#           REVISED 2026-09-07: extended from 6 to 7 combo-eligible nodes
#           (harm_future added back in, politics still excluded -- see
#           r_patches/22_intervention_bootstrap_ci.R for politics' separate
#           standalone single-node treatment) per author request, so the
#           pair/triple search is consistent with Figure 4/Table 2 now
#           reporting single-node effects for all 8 attitude/context nodes.
#           replacing 02 v4's Monte Carlo simulation (simulate_scm_generic(),
#           n=20000 draws per cell) with exact mean propagation
#           (scm_mean_propagate(), same deterministic helper script 30 already
#           used for the 6 single-node targets).
#
# Script 30 fixed the single-node ATEs (script 02 v4's Monte Carlo estimate
# had produced a spurious negative lower bound for policy_support, since its
# true ATE in some scenarios is structurally exactly 0 once its only outgoing
# edge is flipped away). That same Monte Carlo noise floor
# (~1/sqrt(20000) ~ .007 SD) is sitting underneath all 35 pair/triple targets
# too -- script 30's own header flagged this as "out of scope for this
# targeted fix -- worth doing that rerun deterministically too, later." This
# is that rerun, same method, extended from 6 targets to all 41.
#
# what stays the same: base_edges, flip_candidates, flip_edges(), is_acyclic(), build_lavaan_
# syntax(), topo_order(), scm_mean_propagate(), and the scenario_list
# generation are copied VERBATIM from script 30 (which copied the first six
# from 02 v4 verbatim in turn) -- the 16 fitted lavaan models are IDENTICAL
# to v4's and script 30's, so nothing about the model fit changes here, only
# the ATE computation method for the 35 combo targets script 30 didn't touch.
# intervene_targets (the 41-target enumeration) is copied verbatim from
# 02 v4's build_intervene_targets().
#
# output:
# tables/orientation_uncertainty_band_combos_full_deterministic.csv -- the
# corrected 41-row combo ATE range table, same shape/columns as v4's
# tables/orientation_uncertainty_band_combos_full.csv, to be compared side by
# side against it.
# tables/orientation_enumeration_combo_ate_deterministic.csv -- the full
# per-scenario x per-target detail table (up to 41 x 16 = 656 rows), same
# shape as v4's tables/orientation_enumeration_combo_ate.csv.
# Expect every target's ate_baseline/ate_min/ate_max to match the Monte Carlo
# version to within ~.01, except for any target whose range currently
# includes a spurious small-negative value from a scenario where one of its
# member nodes has no remaining downstream path (the same policy_support
# pattern script 30 already found and fixed for the single-node case).
#
# run in the same R session as 02 v4 / script 30 (or independently -- this
# script is self-contained and refits its own 16 lavaan models). Requires
# df_extended (N=870, pre-Wave-5 corrected data) already in scope, exactly as
# 02 v4 and script 30 required.
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

# --- 1: baseline 16-edge working-SCM skeleton -- VERBATIM from script 30 ----
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

# --- 2: helpers -- VERBATIM from script 30 (which took them from 02 v4 / -----
# 27_diagnose_hf_hp_edge.R) ---------------------------------------------------

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

# scm_mean_propagate() already takes intervene as a named list, so setting
# more than one node at once (the whole reason pairs/triples exist) needs no
# changes to the helper itself -- just pass more entries.

# --- 3: enumerate ALL 16 orientation combinations -- VERBATIM from script 30

scenario_list <- setNames(
  lapply(0:(2^n_flip - 1), function(k) which(as.logical(intToBits(k)[1:n_flip]))),
  paste0("combo_", 0:(2^n_flip - 1))
)

scenario_flip_labels <- purrr::imap_chr(scenario_list, function(idx, nm) {
  if (length(idx) == 0) "(none -- baseline)" else paste(flip_candidates$edge_label[idx], collapse = "; ")
})

all_nodes <- unique(c(base_edges$from, base_edges$to))
# REVISED 2026-09-07: harm_future added back into the combo-eligible node
# set (7 nodes now); politics remains excluded (background covariate, not a
# plausible intervention target).
intervene_nodes <- c("belief_concern", "harm_present", "harm_future",
                      "weather_risk_prep", "social_norms", "trust_science",
                      "policy_support")

# --- 3b: enumerate all 63 combo targets -------------------------------------
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
stopifnot(length(intervene_targets) ==
            choose(7, 1) + choose(7, 2) + choose(7, 3))  # 7+21+35 = 63

# --- 4: main loop: refit + DETERMINISTIC combo ATEs for every acyclic scenario
# Same acyclicity/convergence gating as v4 and script 30, so the same
# scenarios that were valid before are valid here -- this script does not
# change which scenarios count, only how their ATEs are computed, and it
# reaches every one of the 63 targets (7 singles + 21 pairs + 35 triples).

det_combo_out <- scenario_list |>
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

    base_means <- scm_mean_propagate(fit, e, all_nodes, list())
    base_cb    <- base_means[["climate_behavior"]]

    combo_rows <- purrr::imap(intervene_targets, function(combo, label) {
      do_means <- scm_mean_propagate(
        fit, e, all_nodes,
        setNames(as.list(rep(0.5, length(combo))), combo)
      )
      tibble::tibble(
        scenario = scen_name,
        target_label = label,
        target_size = length(combo),
        ate_climate_behavior = do_means[["climate_behavior"]] - base_cb
      )
    }) |> dplyr::bind_rows()

    list(ate = combo_rows, status = "ok")
  })

det_combo_status  <- purrr::map_chr(det_combo_out, "status")
det_combo_results <- det_combo_out |> purrr::map("ate") |> purrr::compact() |> dplyr::bind_rows()

cat("\n--- Scenarios attempted:", length(scenario_list),
    "| acyclic+converged:", sum(det_combo_status == "ok"),
    "| cyclic (excluded):", sum(det_combo_status == "cyclic_skipped"),
    "| fit failed:", sum(det_combo_status == "fit_failed"), "---\n")
cat("(This should read identically to 02 v4's and script 30's own status",
    "lines -- same 16 scenarios, same acyclicity/convergence gating, only",
    "the ATE computation method differs. If the counts here do NOT match,",
    "stop and report the discrepancy before trusting anything below.)\n")

readr_ok <- requireNamespace("readr", quietly = TRUE)

# --- 5: headline uncertainty band, deterministic version, all 63 targets ---
det_combo_uncertainty_band <- det_combo_results |>
  dplyr::group_by(target_label, target_size) |>
  dplyr::summarise(
    ate_baseline = ate_climate_behavior[scenario == "combo_0"],
    ate_min = min(ate_climate_behavior),
    ate_max = max(ate_climate_behavior),
    range   = ate_max - ate_min,
    n_specs = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(ate_baseline))

cat("\n--- All 63 intervention targets, DETERMINISTIC (exact mean",
    "propagation, no Monte Carlo noise) -- compare directly against",
    "tables/orientation_uncertainty_band_combos_full.csv (v4, Monte Carlo).",
    "Expect all values to match to within ~.01. ---\n")
print(as.data.frame(det_combo_uncertainty_band), row.names = FALSE)

if (readr_ok) {
  readr::write_csv(det_combo_uncertainty_band, "tables/orientation_uncertainty_band_combos_full_deterministic.csv")
  readr::write_csv(det_combo_results, "tables/orientation_enumeration_combo_ate_deterministic.csv")
} else {
  write.csv(det_combo_uncertainty_band, "tables/orientation_uncertainty_band_combos_full_deterministic.csv", row.names = FALSE)
  write.csv(det_combo_results, "tables/orientation_enumeration_combo_ate_deterministic.csv", row.names = FALSE)
}

# --- 6: explicit side-by-side vs. the Monte Carlo version, if available ----
mc_path <- "tables/orientation_uncertainty_band_combos_full.csv"
if (file.exists(mc_path)) {
  mc_band <- if (readr_ok) readr::read_csv(mc_path, show_col_types = FALSE) else read.csv(mc_path)
  compare_tbl <- det_combo_uncertainty_band |>
    dplyr::select(target_label, target_size, det_baseline = ate_baseline, det_min = ate_min, det_max = ate_max) |>
    dplyr::inner_join(
      mc_band |> dplyr::select(target_label, mc_baseline = ate_baseline, mc_min = ate_min, mc_max = ate_max),
      by = "target_label"
    ) |>
    dplyr::mutate(
      diff_baseline = round(det_baseline - mc_baseline, 4),
      diff_min = round(det_min - mc_min, 4),
      diff_max = round(det_max - mc_max, 4)
    ) |>
    dplyr::arrange(dplyr::desc(det_baseline))
  cat("\n--- Deterministic vs. Monte Carlo (v4) side by side (41 of the 63 targets have a",
    "Monte Carlo counterpart; the 22 new harm_future combos do not) --",
      "diff columns should all be small (~.01 or less), confirming this",
      "rerun is a pure noise-removal exercise and not a substantive change.",
      "Flag any target where diff_min/diff_max looks like more than noise",
      "(the policy_support single-node pattern script 30 found: a spurious",
      "small-negative Monte Carlo value where the true ATE is exactly 0). ---\n")
  print(as.data.frame(compare_tbl), row.names = FALSE)
  write.csv(compare_tbl, "tables/orientation_combo_ate_deterministic_vs_montecarlo.csv", row.names = FALSE)

  cat("\n--- Top pair and top triple, deterministic version (for the",
      "manuscript's Table 3 / Figure 6): ---\n")
  print(as.data.frame(det_combo_uncertainty_band |> dplyr::filter(target_size == 2) |> dplyr::slice(1)), row.names = FALSE)
  print(as.data.frame(det_combo_uncertainty_band |> dplyr::filter(target_size == 3) |> dplyr::slice(1)), row.names = FALSE)
} else {
  cat("\n(", mc_path, "not found in this session -- skipping side-by-side comparison table;",
      "det_combo_uncertainty_band above is still the corrected headline result.)\n")
}

cat("\n--- Done. Send back: (1) the 'Scenarios attempted' status line above,",
    "(2) the full 63-row det_combo_uncertainty_band printout (or at minimum",
    "the top pair and top triple lines) so I can identify whichever pair and",
    "triple now come out on top over the 7-node search (harm_future is",
    "combo-eligible for the first time in this rerun, so the previous",
    "6-node winners -- belief_concern+weather_risk_prep and",
    "harm_present+weather_risk_prep+social_norms -- may no longer be best),",
    "(3) the deterministic-vs-Monte-Carlo comparison table if it printed, so",
    "I can update Table 2 / Figure 4 / the Supplement all-combo figure with",
    "the corrected numbers. ---\n")
