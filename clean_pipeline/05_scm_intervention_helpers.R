# shared scm scaffolding for 06-09 -- base_edges + helper fns used to live
# copy-pasted across 4 scripts (30/31/32/22), pulled into one place here so
# they can't drift apart from each other.


source("clean_pipeline/00_config.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lavaan); library(igraph)
  library(purrr); library(tibble)
})

# ---- The working SCM's 16 edges, in ONE baseline orientation. -------------
# UPDATED 2026-09-11: politics<->belief_concern and policy_support<->
# social_norms now point the way the bootstrap orientation asymmetry
# actually favors (politics asymmetry was -.336, policy_support->social_norms
# was -.325 -- both consistently negative at both alpha levels, not a
# staleness artifact) instead of the old theory-asserted direction. This is
# a real change to the fitted SCM, not just a sensitivity variant -- refit
# and paste the new belief_concern->politics and social_norms->policy_support
# coefficients into r_patches/07_figure5_scm_hierarchical_v3.R's BETA_TR
# before trusting that figure again. See analysis_decisions_log.md.
base_edges <- tibble::tribble(
  ~from,               ~to,
  "belief_concern",    "politics",
  "belief_concern",    "harm_future",
  "belief_concern",    "harm_present",
  "harm_future",       "harm_present",
  "belief_concern",    "trust_science",
  "harm_future",       "trust_science",
  "belief_concern",    "policy_support",
  "trust_science",     "policy_support",
  "politics",          "policy_support",
  "social_norms",      "policy_support",
  "trust_science",     "social_norms",
  "harm_present",      "weather_risk_prep",
  "belief_concern",    "weather_risk_prep",
  "harm_present",      "climate_behavior",
  "social_norms",      "climate_behavior",
  "weather_risk_prep", "climate_behavior"
)
stopifnot(nrow(base_edges) == 16)

# ---- Provenance cross-check: base_edges above should be the SAME 16 edges
# 04_scm_finalize.R wrote to scm_edges_finalized.csv (direction-sensitive,
# order-insensitive). This is the guard against the copy above silently
# going stale relative to whatever 04 last audited -- if 04 has been run at
# least once, this checks agreement; if not, it warns instead of failing so
# 06-09 can still be run standalone (same convention as 05_figures2_3_ggm_
# redesign.R's exists() guards).
scm_csv_path <- file.path(OUTPUT_DIR, "scm_edges_finalized.csv")
if (file.exists(scm_csv_path)) {
  finalized <- utils::read.csv(scm_csv_path, stringsAsFactors = FALSE)
  key <- function(df) sort(paste(df$from, df$to, sep = "->"))
  if (!identical(key(base_edges), key(finalized))) {
    stop(
      "base_edges in 05_scm_intervention_helpers.R does NOT match the 16 edges ",
      "in ", scm_csv_path, " (from 04_scm_finalize.R). One of them is stale -- ",
      "resolve before trusting any intervention ATE."
    )
  }
} else {
  warning(
    scm_csv_path, " not found -- skipping the base_edges provenance check ",
    "(run 04_scm_finalize.R at least once to enable it). Proceeding with ",
    "the hardcoded base_edges as-is."
  )
}

# ---- The 4 directionally-unresolved edges, flippable for the 16-scenario ---
# orientation-sensitivity enumeration. UPDATED 2026-09-11: politics/
# belief_concern and policy_support/social_norms are OUT (their asymmetry is
# large, just pointing the other way -- see base_edges above, they're no
# longer flip candidates, they're just correctly oriented now). IN: the two
# weather_risk_prep edges, whose pooled asymmetry is genuinely near zero
# (+.065 and +.078, both under the .10 band in 00_config.R). social_norms->
# climate_behavior and harm_future->harm_present are unchanged.
flip_candidates <- tibble::tribble(
  ~edge_label,                              ~from,               ~to,
  "belief_concern -> weather_risk_prep",     "belief_concern",    "weather_risk_prep",
  "harm_present -> weather_risk_prep",       "harm_present",      "weather_risk_prep",
  "social_norms -> climate_behavior",        "social_norms",      "climate_behavior",
  "harm_future -> harm_present",             "harm_future",       "harm_present"
)
n_flip <- nrow(flip_candidates)
stopifnot(n_flip == 4)

all_nodes <- unique(c(base_edges$from, base_edges$to))

# ---- Structural helpers -- verbatim from 30/31/32/22 ------------------------
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

# Exact mean propagation for a linear-recursive SCM: E[node] under a given
# intervention assignment is sum(beta_parent * E[parent]) for unintervened
# nodes (residual noise is mean-zero by construction) and the fixed value
# for intervened nodes; exogenous nodes have standardized mean 0. Replaces
# the old Monte Carlo simulate_scm_generic() (n=20000 draws) -- exact rather
# than approximate, and free of the simulation noise that produced a
# spurious small-negative ATE for policy_support in one orientation scenario
# (see 30's header for the full story).
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

target_ate <- function(fit, edges, all_nodes, target_combo) {
  baseline_cb <- scm_mean_propagate(fit, edges, all_nodes, intervene = list())[["climate_behavior"]]
  do_cb <- scm_mean_propagate(
    fit, edges, all_nodes,
    intervene = setNames(as.list(rep(0.5, length(target_combo))), target_combo)
  )[["climate_behavior"]]
  do_cb - baseline_cb
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

# ---- The 16-scenario orientation enumeration -- verbatim from 30/31/32 -----
scenario_list <- setNames(
  lapply(0:(2^n_flip - 1), function(k) which(as.logical(intToBits(k)[1:n_flip]))),
  paste0("combo_", 0:(2^n_flip - 1))
)
scenario_flip_labels <- purrr::imap_chr(scenario_list, function(idx, nm) {
  if (length(idx) == 0) "(none -- baseline)" else paste(flip_candidates$edge_label[idx], collapse = "; ")
})

cat("05_scm_intervention_helpers.R loaded: 16 edges, 4 flip candidates,",
    length(scenario_list), "orientation scenarios.\n")
