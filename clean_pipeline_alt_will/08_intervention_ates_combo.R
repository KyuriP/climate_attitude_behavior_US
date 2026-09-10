# combo ates -- all 63 pair/triple targets over the 7 non-political nodes,
# across all 16 scenarios. from 31_deterministic_combo_ates.R.


source("clean_pipeline_alt_will/05_scm_intervention_helpers.R")

stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

intervene_nodes <- c("belief_concern", "harm_present", "harm_future",
                      "weather_risk_prep", "social_norms", "trust_science",
                      "policy_support")

intervene_targets <- build_intervene_targets(intervene_nodes, max_size = 3)
stopifnot(length(intervene_targets) == choose(7, 1) + choose(7, 2) + choose(7, 3))  # 63

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
      do_means <- scm_mean_propagate(fit, e, all_nodes, setNames(as.list(rep(0.5, length(combo))), combo))
      tibble::tibble(scenario = scen_name, target_label = label, target_size = length(combo),
                      ate_climate_behavior = do_means[["climate_behavior"]] - base_cb)
    }) |> dplyr::bind_rows()
    list(ate = combo_rows, status = "ok")
  })

det_combo_status  <- purrr::map_chr(det_combo_out, "status")
det_combo_results <- det_combo_out |> purrr::map("ate") |> purrr::compact() |> dplyr::bind_rows()

cat("\n--- Scenarios attempted:", length(scenario_list),
    "| acyclic+converged:", sum(det_combo_status == "ok"),
    "| cyclic (excluded):", sum(det_combo_status == "cyclic_skipped"),
    "| fit failed:", sum(det_combo_status == "fit_failed"), "---\n")

det_combo_uncertainty_band <- det_combo_results |>
  dplyr::group_by(target_label, target_size) |>
  dplyr::summarise(
    ate_baseline = ate_climate_behavior[scenario == "combo_0"],
    ate_min = min(ate_climate_behavior), ate_max = max(ate_climate_behavior),
    range = ate_max - ate_min, n_specs = dplyr::n(), .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(ate_baseline))

cat("\n--- All 63 intervention targets (deterministic) -- top 5 by baseline: ---\n")
print(as.data.frame(head(det_combo_uncertainty_band, 5)), row.names = FALSE)

write.csv(det_combo_uncertainty_band, file.path(TABLES_DIR, "orientation_uncertainty_band_combos_full_deterministic.csv"), row.names = FALSE)
write.csv(det_combo_results, file.path(TABLES_DIR, "orientation_enumeration_combo_ate_deterministic.csv"), row.names = FALSE)

cat("\n--- Top pair and top triple (for the manuscript's Table 3 / Figure 6): ---\n")
print(as.data.frame(det_combo_uncertainty_band |> dplyr::filter(target_size == 2) |> dplyr::slice(1)), row.names = FALSE)
print(as.data.frame(det_combo_uncertainty_band |> dplyr::filter(target_size == 3) |> dplyr::slice(1)), row.names = FALSE)
cat("\n(Top triple's baseline ATE should read ~.305 -- present harm + weather risk +",
    "social norms -- matching the manuscript's Results section, before trusting",
    "anything else here.)\n")
