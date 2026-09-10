# behavior outcome sensitivity, from 23_behavior_outcome_sensitivity.R.
# feeds supp S12 (mitig4 adjacency %s, betas, single-node ates).
#
# weather_risk_prep is a worry item and 2 of the 6 climate_behavior items
# are weather-adjacent (evacuate/move), so this reruns the behavior-facing
# stuff on mitig4 (the other 4 items only) plus evacuate/move alone, to
# check the weather_risk_prep -> behavior result isn't just an adaptation
# artifact.


source("clean_pipeline_alt_will/05_scm_intervention_helpers.R")
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
  library(psych)
})

stopifnot(
  exists("df_main"), exists("df_behavior_w5"),
  exists("node_order_ext"), exists("context_idx"), exists("run_one_ext")
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

intervene_nodes <- c("belief_concern", "harm_present", "weather_risk_prep",
                      "social_norms", "trust_science", "policy_support")

ggm_results     <- list()
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
    from = NODE_ORDER_MAIN,
    weight = wmat_alt[NODE_ORDER_MAIN, "climate_behavior"]
  ) |> dplyr::arrange(dplyr::desc(abs(weight)))
  cat("\n-- (a) EBICglasso edge weights into", outcome_name, "--\n")
  print(as.data.frame(ggm_edges_to_outcome), row.names = FALSE)
  ggm_results[[outcome_name]] <- ggm_edges_to_outcome |> dplyr::mutate(outcome = outcome_name)

  # (b) FCI/PC bootstrap existence proportions, attitude -> behavior edges only
  agg_ext_alt <- as.matrix(df_net_ext_alt[, node_order_ext])
  future::plan(future::multisession, workers = max(1L, parallelly::availableCores() - 1L))
  fci_marks_ext <- c("N", "o", ">", "-")
  mark_props_alt <- list()
  for (alph in names(ALPHAS)) {
    cat("Bootstrapping (", outcome_name, ") alpha =", alph, "\n")
    res <- furrr::future_map(
      seq_len(N_BOOT), run_one_ext,
      data = agg_ext_alt, nms = node_order_ext, alpha = ALPHAS[[alph]],
      ctx_idx = context_idx,
      .options = furrr::furrr_options(seed = TRUE), .progress = TRUE
    )
    p_ext <- length(node_order_ext)
    fci_c <- array(0L, dim = c(p_ext, p_ext, 4L),
                    dimnames = list(node_order_ext, node_order_ext, fci_marks_ext))
    for (r in res) if (!is.null(r)) fci_c <- fci_c + r$fci
    mark_props_alt[[alph]] <- fci_c / N_BOOT
  }
  future::plan(future::sequential)

  pooled_alt <- (mark_props_alt[["0.05"]] + mark_props_alt[["0.01"]]) / 2
  behavior_adj <- tibble::tibble(
    from = NODE_ORDER_MAIN,
    existence_pct = round(100 * (1 - pooled_alt[NODE_ORDER_MAIN, "climate_behavior", "N"]), 1)
  ) |> dplyr::arrange(dplyr::desc(existence_pct))
  cat("\n-- (b) Pooled FCI bootstrap existence proportion, attitude ->", outcome_name, "--\n")
  print(as.data.frame(behavior_adj), row.names = FALSE)
  fci_adj_results[[outcome_name]] <- behavior_adj |> dplyr::mutate(outcome = outcome_name)
  if (outcome_name == "mitig4") {
    cat("Manuscript (S12) reports: present harm 94.3%, social norms 71.5%,",
        "policy support 60.4%, weather risk 38.5%, belief/concern 52.9% -- compare.\n")
  }

  # (c) Behavior equation of the working SCM, refit on the alternative outcome
  eq_syntax <- "climate_behavior ~ harm_present + weather_risk_prep + social_norms"
  fit_eq <- lavaan::sem(eq_syntax, data = df_extended_alt, estimator = "MLR", fixed.x = FALSE)
  eq_std <- lavaan::standardizedSolution(fit_eq) |> dplyr::filter(op == "~")
  cat("\n-- (c) Behavior equation, outcome =", outcome_name, "--\n")
  print(as.data.frame(eq_std[, c("lhs", "rhs", "est.std", "se", "pvalue")]), row.names = FALSE)
  scm_eq_results[[outcome_name]] <- eq_std |>
    dplyr::transmute(predictor = rhs, beta = round(est.std, 3), se = round(se, 3),
                      p = round(pvalue, 4), outcome = outcome_name)
  if (outcome_name == "mitig4") {
    cat("Manuscript (S12) reports weather_risk_prep beta=.135 (SE=.035, p=.0001),",
        "harm_present beta=.348 (SE=.035), social_norms beta=.127 (SE=.034, p=.0002) -- compare.\n")
  }

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
    if (outcome_name == "mitig4") {
      cat("Manuscript (S12) reports weather_risk_prep=.068 SD, social_norms=.064 SD,",
          "present_harm=.204 SD -- compare.\n")
    }
  }
}

# --- Write outputs -----------------------------------------------------------
write.csv(dplyr::bind_rows(ggm_results),     file.path(OUTPUT_DIR, "behavior_sensitivity_ggm.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(fci_adj_results), file.path(OUTPUT_DIR, "behavior_sensitivity_fci_adjacency.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(scm_eq_results),  file.path(OUTPUT_DIR, "behavior_sensitivity_scm_equation.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(ate_results),     file.path(OUTPUT_DIR, "behavior_sensitivity_ate.csv"), row.names = FALSE)

cat("\n\nWrote behavior_sensitivity_{ggm,fci_adjacency,scm_equation,ate}.csv to", OUTPUT_DIR, "\n")
cat("Compare weather_risk_prep's row across mitig4 / evacuate / move against the original\n",
    "six-item climate_behavior result to see whether it survives outcome narrowing.\n", sep = "")
