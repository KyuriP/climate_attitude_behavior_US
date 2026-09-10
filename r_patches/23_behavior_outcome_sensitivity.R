# =============================================================================
# 23_behavior_outcome_sensitivity.R
#
# Purpose (methodological review, item #4 / Blocker 4): weather_risk_prep
# ("weather risk") is a worry/perceived-risk item, and the six-item
# climate_behavior composite includes two behaviors that are conceptually
# close to weather-related risk (emergency preparedness, relocation). This
# makes the weather_risk_prep -> climate_behavior result potentially partly
# adaptation-specific rather than a general mitigation/engagement finding.
#
# Design, per team decision -- a TARGETED sensitivity, not a second full
# pipeline: rerun only the behavior-facing analyses, on three alternative
# outcome definitions built from the already-existing raw Wave-5 behavior
# items (no new data, no new items):
#   - mitig4   = mean(beh_meat, beh_travel, beh_activ, beh_discuss)  -- the
#                four mitigation/engagement items, excluding both
#                weather-proximal behaviors.
#   - evacuate = beh_evacuate alone (kept separate, NOT averaged into an
#                "adaptation" composite per team decision -- that composite
#                can be added later only if warranted).
#   - move     = beh_move alone, same reasoning.
#
# For each of the three alternative outcomes, this script reruns:
#   (a) the extended GGM (EBICglasso, same spec as the primary analysis),
#       reporting edge weights from every attitude/context node to the
#       alternative outcome;
#   (b) the extended FCI/PC bootstrap existence proportions for the 8
#       attitude-to-behavior adjacencies specifically (same run_one_ext()
#       machinery, same n_boot/alpha settings as the primary extended
#       analysis -- reused as-is, not reimplemented);
#   (c) the behavior equation of the working SCM (climate_behavior ~
#       harm_present + weather_risk_prep + social_norms), refit with the
#       alternative outcome as the dependent variable;
#   (d) the six single-node intervention effects on the alternative outcome,
#       via the same deterministic mean-propagation approach as script 22
#       (full 16-edge working SCM, alternative outcome substituted in place
#       of the original climate_behavior).
#
# The key question this answers: does weather_risk_prep remain an important
# behavioral predictor once the two most weather-proximal behaviors are
# removed from the outcome (mitig4), and does it show up at all for the two
# adaptation items considered on their own?
#
# HOW TO USE: run within the same .qmd session, after Section 4.3 (so
# df_behavior_w5 with beh_meat/beh_travel/beh_activ/beh_discuss/beh_evacuate/
# beh_move exists) AND after Section 7.4's bootstrap setup (so run_one_ext,
# node_order_ext, context_idx, alphas, n_boot are already defined -- this
# script reuses that machinery unmodified rather than reimplementing it).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lavaan)
  library(qgraph)
  library(bootnet)
  library(huge)
  library(purrr)
  library(tibble)
  library(furrr)
  library(future)
})

stopifnot(
  exists("df_main"), exists("df_behavior_w5"),
  exists("node_order_ext"), exists("context_idx"), exists("run_one_ext"),
  exists("main_nodes"), exists("alphas"), exists("n_boot")
)
stopifnot(all(c("beh_meat", "beh_travel", "beh_activ", "beh_discuss",
                "beh_evacuate", "beh_move") %in% names(df_behavior_w5)))

set.seed(2026)

# --- 1. Build the three alternative outcomes (from existing raw items only) -
alt_outcomes <- df_behavior_w5 |>
  dplyr::transmute(
    participant_id,
    mitig4   = rowMeans(dplyr::pick(beh_meat, beh_travel, beh_activ, beh_discuss), na.rm = TRUE),
    evacuate = beh_evacuate,
    move     = beh_move
  )

cat("mitig4 (4-item mitigation/engagement) alpha =",
    round(psych::alpha(df_behavior_w5 |>
      dplyr::select(beh_meat, beh_travel, beh_activ, beh_discuss))$total$raw_alpha, 3), "\n")

# --- Shared helpers (same base_edges / deterministic propagation as script 22)
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
all_nodes <- unique(c(base_edges$from, base_edges$to))
intervene_nodes <- c("belief_concern", "harm_present", "weather_risk_prep",
                      "social_norms", "trust_science", "policy_support")

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

ggm_results   <- list()
fci_adj_results <- list()
scm_eq_results  <- list()
ate_results     <- list()

for (outcome_name in c("mitig4", "evacuate", "move")) {

  cat("\n\n=====================================================\n")
  cat("Alternative outcome:", outcome_name, "\n")
  cat("=====================================================\n")

  df_extended_alt <- df_main |>
    dplyr::inner_join(
      alt_outcomes |> dplyr::select(participant_id, climate_behavior = dplyr::all_of(outcome_name)),
      by = "participant_id"
    )
  cat("N =", nrow(df_extended_alt), "\n")

  # (a) Extended GGM ----------------------------------------------------------
  df_net_ext_alt <- df_extended_alt |>
    dplyr::select(dplyr::all_of(node_order_ext)) |>
    huge::huge.npn(npn.func = "truncation") |>
    as.data.frame()
  colnames(df_net_ext_alt) <- node_order_ext

  net_alt <- bootnet::estimateNetwork(df_net_ext_alt, default = "EBICglasso",
                                      tuning = 0.5, corMethod = "cor")
  wmat_alt <- qgraph::getWmat(net_alt)
  ggm_edges_to_outcome <- tibble::tibble(
    from = main_nodes,
    weight = wmat_alt[main_nodes, "climate_behavior"]
  ) |> dplyr::arrange(dplyr::desc(abs(weight)))
  cat("\n-- (a) EBICglasso edge weights into", outcome_name, "--\n")
  print(as.data.frame(ggm_edges_to_outcome), row.names = FALSE)
  ggm_results[[outcome_name]] <- ggm_edges_to_outcome |> dplyr::mutate(outcome = outcome_name)

  # (b) FCI/PC bootstrap existence proportions, attitude -> behavior edges only
  agg_ext_alt <- as.matrix(df_net_ext_alt[, node_order_ext])
  plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))
  fci_marks_ext <- c("N", "o", ">", "-")
  mark_props_alt <- list()
  for (alph in names(alphas)) {
    cat("Bootstrapping (", outcome_name, ") alpha =", alph, "\n")
    res <- furrr::future_map(
      seq_len(n_boot), run_one_ext,
      data = agg_ext_alt, nms = node_order_ext, alpha = alphas[[alph]],
      ctx_idx = context_idx,
      .options = furrr::furrr_options(seed = TRUE), .progress = TRUE
    )
    p_ext <- length(node_order_ext)
    fci_c <- array(0L, dim = c(p_ext, p_ext, 4L),
                    dimnames = list(node_order_ext, node_order_ext, fci_marks_ext))
    for (r in res) if (!is.null(r)) fci_c <- fci_c + r$fci
    mark_props_alt[[alph]] <- fci_c / n_boot
  }
  plan(sequential)

  pooled_alt <- (mark_props_alt[["0.05"]] + mark_props_alt[["0.01"]]) / 2
  behavior_adj <- tibble::tibble(
    from = main_nodes,
    existence_pct = round(100 * (1 - pooled_alt[main_nodes, "climate_behavior", "N"]), 1)
  ) |> dplyr::arrange(dplyr::desc(existence_pct))
  cat("\n-- (b) Pooled FCI bootstrap existence proportion, attitude ->", outcome_name, "--\n")
  print(as.data.frame(behavior_adj), row.names = FALSE)
  fci_adj_results[[outcome_name]] <- behavior_adj |> dplyr::mutate(outcome = outcome_name)

  # (c) Behavior equation of the working SCM, refit on the alternative outcome
  eq_syntax <- "climate_behavior ~ harm_present + weather_risk_prep + social_norms"
  fit_eq <- lavaan::sem(eq_syntax, data = df_extended_alt, estimator = "MLR", fixed.x = FALSE)
  eq_std <- lavaan::standardizedSolution(fit_eq) |> dplyr::filter(op == "~")
  cat("\n-- (c) Behavior equation, outcome =", outcome_name, "--\n")
  print(as.data.frame(eq_std[, c("lhs", "rhs", "est.std", "se", "pvalue")]), row.names = FALSE)
  scm_eq_results[[outcome_name]] <- eq_std |>
    dplyr::transmute(predictor = rhs, beta = round(est.std, 3), se = round(se, 3),
                      p = round(pvalue, 4), outcome = outcome_name)

  # (d) Six single-node intervention effects, full 16-edge SCM refit
  model_syntax_alt <- build_lavaan_syntax(base_edges)
  fit_full_alt <- tryCatch(
    lavaan::sem(model_syntax_alt, data = df_extended_alt, estimator = "MLR", fixed.x = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit_full_alt) || !lavaan::lavInspect(fit_full_alt, "converged")) {
    cat("\n-- (d) Full 16-edge SCM refit FAILED/non-converged for outcome", outcome_name, "--\n")
    ate_results[[outcome_name]] <- tibble::tibble(node = intervene_nodes, ate = NA_real_,
                                                    outcome = outcome_name)
  } else {
    ates_alt <- single_node_ates(fit_full_alt, base_edges, all_nodes, intervene_nodes)
    cat("\n-- (d) Single-node ATEs on", outcome_name, "(full 16-edge SCM) --\n")
    print(round(ates_alt, 4))
    ate_results[[outcome_name]] <- tibble::tibble(node = names(ates_alt), ate = round(unname(ates_alt), 4),
                                                    outcome = outcome_name)
  }
}

# --- Write outputs -----------------------------------------------------------
dir.create("pipeline_outputs", showWarnings = FALSE)
write.csv(dplyr::bind_rows(ggm_results),     "pipeline_outputs/behavior_sensitivity_ggm.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(fci_adj_results), "pipeline_outputs/behavior_sensitivity_fci_adjacency.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(scm_eq_results),  "pipeline_outputs/behavior_sensitivity_scm_equation.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(ate_results),     "pipeline_outputs/behavior_sensitivity_ate.csv", row.names = FALSE)

cat("\n\nWrote pipeline_outputs/behavior_sensitivity_{ggm,fci_adjacency,scm_equation,ate}.csv\n")
cat("Compare weather_risk_prep's row across mitig4 / evacuate / move against the original\n",
    "six-item climate_behavior result to see whether it survives outcome narrowing.\n", sep = "")
