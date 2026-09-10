# single-node ates across all 16 orientation scenarios -- from
# 30_deterministic_16spec_ates.R, now using the shared helpers/base_edges
# from 05 instead of its own copy.


source("clean_pipeline_alt_will/05_scm_intervention_helpers.R")

stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

intervene_nodes <- c("belief_concern", "harm_present", "weather_risk_prep",
                      "social_norms", "trust_science", "policy_support")

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
      scenario = scen_name, node = intervene_nodes,
      ate_climate_behavior = unname(ates[intervene_nodes])
    )
    list(ate = ate_rows, status = "ok")
  })

det_status  <- purrr::map_chr(det_out, "status")
det_results <- det_out |> purrr::map("ate") |> purrr::compact() |> dplyr::bind_rows()

cat("\n--- Scenarios attempted:", length(scenario_list),
    "| acyclic+converged:", sum(det_status == "ok"),
    "| cyclic (excluded):", sum(det_status == "cyclic_skipped"),
    "| fit failed:", sum(det_status == "fit_failed"), "---\n")

det_uncertainty_band <- det_results |>
  dplyr::group_by(node) |>
  dplyr::summarise(
    ate_baseline = ate_climate_behavior[scenario == "combo_0"],
    ate_min = min(ate_climate_behavior), ate_max = max(ate_climate_behavior),
    range = ate_max - ate_min, n_specs = dplyr::n(), .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(range))

cat("\n--- Six single-node intervention ranges (deterministic, exact mean propagation) ---\n")
print(as.data.frame(det_uncertainty_band), row.names = FALSE)

write.csv(det_uncertainty_band, file.path(TABLES_DIR, "orientation_uncertainty_band_full_deterministic.csv"), row.names = FALSE)
write.csv(det_results, file.path(TABLES_DIR, "orientation_enumeration_ate_deterministic.csv"), row.names = FALSE)
cat("\nWrote", file.path(TABLES_DIR, "orientation_uncertainty_band_full_deterministic.csv"), "and",
    file.path(TABLES_DIR, "orientation_enumeration_ate_deterministic.csv"), "\n")
cat("(Belief/concern and present harm should be the two largest baseline ATEs,",
    "~.206 and ~.195 -- matching the manuscript's Results section, before",
    "trusting anything else here.)\n")
