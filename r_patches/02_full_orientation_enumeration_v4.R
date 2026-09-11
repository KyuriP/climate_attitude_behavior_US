# =============================================================================
# Patch 02 (v4): FULL enumeration of the 2^4 = 16 orientation combinations for
#           the four edges current data leaves genuinely direction-unclear,
#           with model-fit indices attached to every specification.
#
# Script 27 (r_patches/27_diagnose_hf_hp_edge.R) ran the 5-step diagnostic on
# harm_future -> harm_present, the single largest coefficient in the model
# (beta=.636). Findings: (1) adjacency is 100% at BOTH alpha=.05 and
# alpha=.01 -- never in question; (2) orientation asymmetry has OPPOSITE
# SIGNS at the two alpha levels (-.108 at .05, +.261 at .01), and this sign
# pattern replicated in ALL 8 independent seed replicates (1,000 resamples
# each) at both thresholds -- not a single-bootstrap-realization artifact;
# (3) reversing the edge produces IDENTICAL model fit (CFI/TLI/RMSEA/SRMR/
# AIC/BIC match to reported precision -- consistent with Markov equivalence,
# since harm_future and harm_present share belief_concern as a common
# parent) and changes only one of the six policy-relevant single-node ATEs,
# by .0046 SD (negligible next to its bootstrap CI).
#
# Orientation rule used throughout this project: an edge is directionally
# supported only if the bootstrap asymmetry sign agrees at BOTH alpha=.05
# and alpha=.01; if the signs disagree, the direction is unresolved even if
# the POOLED (across-alpha) asymmetry happens to be positive. That is
# exactly what happens here: harm_future->harm_present's pooled asymmetry is
# a weak +.077 (see pipeline_outputs/scm_edges_finalized.csv), which is why
# script 18's audit calls it "bootstrap-resolved" -- pooling masks a real
# per-alpha disagreement (see r_patches/21_orientation_crossalpha_table.R,
# which is what first surfaced this edge as one of only two -- alongside
# social_norms->climate_behavior -- where the sign flips across alpha). Per
# this rule, harm_future -> harm_present joins the flip set below as its 4th
# edge. No other edge in the 16-edge skeleton meets the sign-disagreement
# criterion (script 21's table has been checked for all 16; only these two
# qualify), so the flip set is 4 edges and the full enumeration is 16
# specifications. Mechanically this leaves the rest of the script alone:
# scenario_list / scenario_flip_labels / the acyclicity+convergence loop /
# the combo-ATE loop are all written generically in terms of
# nrow(flip_candidates), so only the flip_candidates table itself (one new
# row) needed to change. Nothing here assumes all 16 specifications turn out
# acyclic and converged -- the explicit per-scenario is_acyclic() and
# convergence checks below are what actually decide that for each of the
# 16, and their results should be read off the console/full_fit table, not
# assumed.
#
# base_edges below is the same 16-edge skeleton used elsewhere (including
# the fixed edge harm_future->trust_science, which is not one of the four
# flippable edges). Output file paths are unchanged (tables/orientation_
# enumeration_fit.csv, _ate.csv, _combo_ate.csv, orientation_uncertainty_
# band_full.csv, orientation_uncertainty_band_combos_full.csv, fit_
# supplement_table.csv) so that 07_figure5_scm_hierarchical_v3.R and
# 10_figure7_uncertainty_pub.R / 08b_supp_figure_all_combo_interventions.R
# need no code changes to pick up the 16-specification results -- only their
# "8 specifications" comments/captions needed a text update.
#
# Run this AFTER script 01 (or at minimum after Part A's data-correction
# step, so that df_extended in the session is already the pre-Wave-5
# primary N = 870 dataset), AND after the SCM model spec (in the .qmd,
# Section 7.5/7.6) already includes all 16 base edges -- no NEW base edge is
# being added here, only a new edge is being made FLIPPABLE. This script is
# self-contained: it redefines its own copies of base_edges,
# flip_candidates, and the helper functions.
#
# output: 16 rows/specifications. Every one is EXPLICITLY checked for
# acyclicity and lavaan convergence -- status will be "ok", "cyclic_skipped",
# or "fit_failed" per scenario; do not assume 16/16 will be "ok". Read the
# actual counts off the console output and full_fit table below before
# trusting anything downstream.
# ===============================================================================

library(dplyr)
library(tidyr)
library(lavaan)
library(igraph)
library(ggplot2)
library(purrr)
library(tibble)

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# --- 0: sanity check that the corrected data is in scope --------------------
stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

# --- 1: baseline 16-edge working-SCM skeleton -- UNCHANGED from v3 ----------
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

# The 4 candidate flip edges under the cross-alpha sign-consistency
# criterion (see header) -- 3 UNCHANGED from v2/v3, PLUS harm_future ->
# harm_present, NEW in v4 per script 27's findings.
flip_candidates <- tibble::tribble(
  ~edge_label,                          ~from,            ~to,
  "politics -> belief_concern",          "politics",       "belief_concern",
  "policy_support -> social_norms",      "policy_support",  "social_norms",
  "social_norms -> climate_behavior",    "social_norms",    "climate_behavior",
  "harm_future -> harm_present",         "harm_future",     "harm_present"
)

n_flip <- nrow(flip_candidates)  # 4, up from 3 in v2/v3
stopifnot(n_flip == 4)

# --- 2: helpers (self-contained copies, unchanged from v1/v2/v3) -----------

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

simulate_scm_generic <- function(fit, edges, all_nodes, intervene = list(),
                                  n = 20000, seed = 42) {
  set.seed(seed)
  ord <- topo_order(edges)
  ord <- union(ord, all_nodes)

  std <- lavaan::standardizedSolution(fit) |> dplyr::filter(op == "~")
  resid_var <- lavaan::standardizedSolution(fit) |>
    dplyr::filter(op == "~~", lhs == rhs) |>
    dplyr::select(node = lhs, resid_var = est.std)

  dat <- as.data.frame(matrix(NA_real_, nrow = n, ncol = length(ord)))
  names(dat) <- ord

  for (node in ord) {
    if (!is.null(intervene[[node]])) {
      dat[[node]] <- intervene[[node]]
      next
    }
    parents <- edges$from[edges$to == node]
    if (length(parents) == 0) {
      dat[[node]] <- rnorm(n, 0, 1)
      next
    }
    coefs <- std |> dplyr::filter(lhs == node, rhs %in% parents)
    beta <- setNames(coefs$est.std, coefs$rhs)
    rv <- resid_var$resid_var[resid_var$node == node]
    if (length(rv) == 0) rv <- 1 - sum(beta^2)
    lin_pred <- Reduce(`+`, lapply(parents, function(p) beta[[p]] * dat[[p]]))
    dat[[node]] <- lin_pred + rnorm(n, 0, sqrt(max(rv, 0.01)))
  }
  dat
}

# --- 3: enumerate ALL 16 combinations ---------------------------------------
# combo_0 = baseline (no flips) through combo_15 = all four flipped.
# Bit k (0-indexed) of the combination index corresponds to flip_candidates
# row k+1 (politics->belief_concern = bit 0, policy_support->social_norms =
# bit 1, social_norms->climate_behavior = bit 2, harm_future->harm_present =
# bit 3, NEW in v4). This generation logic itself is unchanged from v2/v3 --
# it was already written generically in terms of n_flip.

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

# Combo-target generation (6 singles + 15 pairs + 20 triples = 41 targets)
# -- UNCHANGED from v3's 2026-09-05 extension. harm_future is not among the
# six intervenable nodes, so this list and its size are unaffected by the
# flip-set change above.
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
            choose(6, 1) + choose(6, 2) + choose(6, 3))  # 6+15+20 = 41

# --- 4: main loop: refit + fit indices + ATEs for every acyclic scenario ----
# (mechanics unchanged from v1/v2/v3 -- only base_edges/flip_candidates
# above changed; the cyclic/fit-failed messaging below is now neutral rather
# than "UNEXPECTED", since a 4th flip edge means we no longer have a strong
# prior that every combination will be acyclic -- explicitly checking, not
# assuming, is the point of this rerun.)

fit_index_names <- c("cfi", "tli", "rmsea", "srmr", "aic", "bic")

full_out <- scenario_list |>
  purrr::imap(function(flip_idx, scen_name) {

    e <- flip_edges(base_edges, flip_idx)
    n_flipped <- length(flip_idx)

    if (!is_acyclic(e)) {
      message(scen_name, " (", scenario_flip_labels[[scen_name]], "): cyclic, skipped.")
      return(list(ate = NULL, fit = tibble::tibble(
        scenario = scen_name, n_flipped = n_flipped,
        flipped_edges = scenario_flip_labels[[scen_name]],
        status = "cyclic_skipped",
        cfi = NA_real_, tli = NA_real_, rmsea = NA_real_, srmr = NA_real_,
        aic = NA_real_, bic = NA_real_
      )))
    }

    model_syntax <- build_lavaan_syntax(e)
    fit <- tryCatch(
      lavaan::sem(model_syntax, data = df_extended, estimator = "MLR", fixed.x = FALSE),
      error = function(err) NULL
    )
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      message(scen_name, " (", scenario_flip_labels[[scen_name]], "): fit failed/non-converged, skipped")
      return(list(ate = NULL, fit = tibble::tibble(
        scenario = scen_name, n_flipped = n_flipped,
        flipped_edges = scenario_flip_labels[[scen_name]],
        status = "fit_failed",
        cfi = NA_real_, tli = NA_real_, rmsea = NA_real_, srmr = NA_real_,
        aic = NA_real_, bic = NA_real_
      )))
    }

    fm <- lavaan::fitmeasures(fit)
    get_fm <- function(nms) {
      hit <- nms[nms %in% names(fm)]
      if (length(hit) == 0) return(NA_real_)
      unname(fm[hit[1]])
    }
    fit_row <- tibble::tibble(
      scenario = scen_name,
      n_flipped = n_flipped,
      flipped_edges = scenario_flip_labels[[scen_name]],
      status = "ok",
      cfi   = get_fm(c("cfi.robust", "cfi")),
      tli   = get_fm(c("tli.robust", "tli")),
      rmsea = get_fm(c("rmsea.robust", "rmsea")),
      srmr  = get_fm(c("srmr")),
      aic   = get_fm(c("aic")),
      bic   = get_fm(c("bic"))
    )

    baseline_sim <- simulate_scm_generic(fit, e, all_nodes, intervene = list(), n = 20000)
    baseline_cb  <- mean(baseline_sim$climate_behavior)

    # Single-node ATEs -- UNCHANGED shape from v3 (same 6 rows per
    # scenario), so 10_figure7_uncertainty_pub.R needs no code changes, only
    # more rows per node (16 scenarios instead of 8).
    ate_rows <- intervene_nodes |>
      purrr::map(function(node) {
        sim <- simulate_scm_generic(fit, e, all_nodes,
                                     intervene = setNames(list(0.5), node), n = 20000)
        tibble::tibble(
          scenario = scen_name,
          node = node,
          ate_climate_behavior = mean(sim$climate_behavior) - baseline_cb
        )
      }) |>
      dplyr::bind_rows()

    # All 41 targets (6 singles + 15 pairs + 20 triples), same baseline_cb
    # and same fitted model as above -- written to a separate csv, unchanged
    # shape from v3's 2026-09-05 extension, just 16 scenarios instead of 8.
    combo_ate_rows <- purrr::imap(intervene_targets, function(combo, label) {
      sim <- simulate_scm_generic(
        fit, e, all_nodes,
        intervene = setNames(as.list(rep(0.5, length(combo))), combo),
        n = 20000
      )
      tibble::tibble(
        scenario = scen_name,
        target_label = label,
        target_size = length(combo),
        ate_climate_behavior = mean(sim$climate_behavior) - baseline_cb
      )
    }) |>
      dplyr::bind_rows()

    list(ate = ate_rows, combo_ate = combo_ate_rows, fit = fit_row)
  })

full_results       <- full_out |> purrr::map("ate")       |> purrr::compact() |> dplyr::bind_rows()
full_combo_results <- full_out |> purrr::map("combo_ate") |> purrr::compact() |> dplyr::bind_rows()
full_fit           <- full_out |> purrr::map("fit")        |> dplyr::bind_rows()

cat("\n--- Scenarios attempted:", length(scenario_list),
    "| acyclic+converged:", sum(full_fit$status == "ok"),
    "| cyclic (excluded):", sum(full_fit$status == "cyclic_skipped"),
    "| fit failed:", sum(full_fit$status == "fit_failed"), "---\n")
cat("(16 attempted. Read the ok/cyclic/failed counts above directly -- v4 does",
    "NOT assume all 16 are valid the way v3 could for its narrower 3-edge set.",
    "If fewer than 16 are 'ok', full_uncertainty_band below is still computed",
    "correctly -- dplyr::bind_rows() on full_results just has fewer scenarios",
    "contributing -- but report the excluded scenario(s) back before trusting",
    "the range as 'the full 16-specification range'.)\n")

print(full_fit)
readr_ok <- requireNamespace("readr", quietly = TRUE)
if (readr_ok) {
  readr::write_csv(full_fit, "tables/orientation_enumeration_fit.csv")
  readr::write_csv(full_results, "tables/orientation_enumeration_ate.csv")
} else {
  write.csv(full_fit, "tables/orientation_enumeration_fit.csv", row.names = FALSE)
  write.csv(full_results, "tables/orientation_enumeration_ate.csv", row.names = FALSE)
}

if (readr_ok) {
  readr::write_csv(full_combo_results, "tables/orientation_enumeration_combo_ate.csv")
} else {
  write.csv(full_combo_results, "tables/orientation_enumeration_combo_ate.csv", row.names = FALSE)
}
message("Wrote tables/orientation_enumeration_combo_ate.csv (",
        nrow(full_combo_results), " rows = ", length(intervene_targets),
        " targets x ", sum(full_fit$status == "ok"), " valid structural specifications).")

# --- 5: headline uncertainty band, now over all VALID (status=="ok") specs -

full_uncertainty_band <- full_results |>
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

cat("\n--- Six single-node intervention ranges across all valid structural",
    "specifications (THE key new result -- compare each node's range here",
    "to its v3/8-specification range to see how much adding",
    "harm_future<->harm_present as a flippable edge widened it) ---\n")
print(full_uncertainty_band)
if (readr_ok) {
  readr::write_csv(full_uncertainty_band, "tables/orientation_uncertainty_band_full.csv")
} else {
  write.csv(full_uncertainty_band, "tables/orientation_uncertainty_band_full.csv", row.names = FALSE)
}

# --- Same treatment for all 41 combo targets (unchanged mechanics from v3) -
full_combo_uncertainty_band <- full_combo_results |>
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

cat("\n--- All 41 intervention targets (baseline effect + range across all",
    "valid structural specifications), sorted by baseline effect size ---\n")
print(as.data.frame(full_combo_uncertainty_band), row.names = FALSE)
if (readr_ok) {
  readr::write_csv(full_combo_uncertainty_band, "tables/orientation_uncertainty_band_combos_full.csv")
} else {
  write.csv(full_combo_uncertainty_band, "tables/orientation_uncertainty_band_combos_full.csv", row.names = FALSE)
}

# --- 6: trust_science check -- unchanged from v3, now over 16 specs --------
# (kept for continuity; trust_science is not the node this v4 rerun is
# actually about, but there is no reason to drop a working diagnostic).
ts_check <- full_results |>
  dplyr::filter(node == "trust_science") |>
  dplyr::inner_join(full_fit, by = "scenario") |>
  dplyr::arrange(ate_climate_behavior)

cat("\n--- trust_science ATE vs. fit, across all valid specifications ---\n")
print(ts_check |> dplyr::select(scenario, flipped_edges, ate_climate_behavior, cfi, tli, rmsea, srmr, bic))

# --- 6b: harm_present check -- NEW in v4. harm_present is one of the six
# intervenable nodes AND one of the two endpoints of the newly-flippable
# edge, so it is the node most directly affected by this rerun (script 27's
# Step 5 substantive-consequence check already found only a +.0046 SD
# difference between the single baseline-vs-reversed comparison; this
# repeats that comparison across the FULL 16-specification space, not just
# the two HF-HP-only endpoints, so cross-effects with the other 3 flippable
# edges are also captured). ---------------------------------------------------
hp_check <- full_results |>
  dplyr::filter(node == "harm_present") |>
  dplyr::inner_join(full_fit, by = "scenario") |>
  dplyr::arrange(ate_climate_behavior)

cat("\n--- harm_present ATE vs. fit, across all valid specifications ---\n")
cat("--- (this is the node most directly affected by adding harm_future<->",
    "harm_present to the flip set -- compare its range here to script 27's",
    "single baseline-vs-reversed diff of +.0046 SD) ---\n")
print(hp_check |> dplyr::select(scenario, flipped_edges, ate_climate_behavior, cfi, tli, rmsea, srmr, bic))

# --- 7: publication-ready Figure 7 ------------------------------------------
if (file.exists("r_patches/10_figure7_uncertainty_pub.R")) {
  source("r_patches/10_figure7_uncertainty_pub.R")
} else {
  message("Figure 7 script not found; full_results is available for plotting.")
}

# --- 8: small supplementary fit table (paste into Overleaf as a table) -----

fit_supplement_table <- full_fit |>
  dplyr::filter(status == "ok") |>
  dplyr::arrange(n_flipped, scenario) |>
  dplyr::select(scenario, flipped_edges, n_flipped, cfi, tli, rmsea, srmr, aic, bic)

print(fit_supplement_table)
if (readr_ok) {
  readr::write_csv(fit_supplement_table, "tables/fit_supplement_table.csv")
} else {
  write.csv(fit_supplement_table, "tables/fit_supplement_table.csv", row.names = FALSE)
}

# =============================================================================
# KEY OUTPUTS FROM THIS RUN (16-edge base, 4-edge/16-combination flip set)
# -----------------
# 1. The console line "Scenarios attempted: 16 | acyclic+converged: N |
#    cyclic (excluded): N | fit failed: N" -- if N is not 16 for
#    acyclic+converged, note which scenario(s) were excluded and why before
#    trusting anything downstream.
# 2. The full_fit console output (up to 16 rows).
# 3. The full_uncertainty_band console output (6 rows) -- THE KEY RESULT.
#    Compare every node's range here to the 8-specification version already
#    in the manuscript, not just harm_present's.
# 4. The ts_check AND hp_check printouts (harm_present is the node most
#    directly touched by the new flippable edge).
# 5. figures/fig7_uncertainty_pub.pdf and tables/fit_supplement_table.csv.
# 6. The "All 41 intervention targets..." console printout,
#    tables/orientation_enumeration_combo_ate.csv (41 targets x up to 16
#    specs), and tables/orientation_uncertainty_band_combos_full.csv (41
#    rows) -- same downstream figures (08b_supp_figure_all_combo_
#    interventions.R) read the latter, no code changes needed there beyond
#    the caption text already updated to say 16.
# =============================================================================

cat("\n--- Full run summary ---\n")
cat("Single-node targets: 6 (unchanged) | Combo targets: ", length(intervene_targets),
    "(6 singles + 15 pairs + 20 triples) | Structural specifications: 16",
    "(up from 8 in v3 -- harm_future<->harm_present added per script 27 + the",
    "newly-adopted cross-alpha sign-consistency criterion)\n")
cat("Total combo simulation cells: ", nrow(full_combo_results),
    "(expect up to", length(intervene_targets) * 16, "if all 16 specs are valid)\n")
