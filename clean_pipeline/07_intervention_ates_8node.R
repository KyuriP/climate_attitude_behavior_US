# same as 06 but adds harm_future + politics (fig 7's 8-node version), from
# 32_deterministic_8node_ates.R. politics is a benchmark covariate, not a
# real intervention target -- say so in any caption built off this.


source("clean_pipeline/05_scm_intervention_helpers.R")

stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

node_targets <- c(
  "belief_concern", "harm_present", "weather_risk_prep",
  "social_norms", "trust_science", "policy_support",
  "harm_future", "politics"
)
stopifnot(length(node_targets) == 8)

det_8node_out <- scenario_list |>
  purrr::imap(function(flip_idx, scen_name) {
    e <- flip_edges(base_edges, flip_idx)
    if (!is_acyclic(e)) { message(scen_name, ": cyclic after flip, skipped."); return(NULL) }
    model_syntax <- build_lavaan_syntax(e)
    fit <- tryCatch(
      lavaan::sem(model_syntax, data = df_extended, estimator = "MLR", fixed.x = FALSE),
      error = function(err) NULL
    )
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      message(scen_name, ": fit failed or did not converge, skipped."); return(NULL)
    }
    base_means <- scm_mean_propagate(fit, e, all_nodes, list())
    base_cb    <- base_means[["climate_behavior"]]
    purrr::map_dfr(node_targets, function(node) {
      do_means <- scm_mean_propagate(fit, e, all_nodes, setNames(list(0.5), node))
      tibble::tibble(scenario = scen_name, node = node,
                      ate_climate_behavior = do_means[["climate_behavior"]] - base_cb)
    })
  }) |>
  dplyr::bind_rows()

stopifnot(nrow(det_8node_out) == 16 * 8)

write.csv(det_8node_out, file.path(TABLES_DIR, "orientation_enumeration_ate_deterministic_8node.csv"), row.names = FALSE)
message("Saved ", file.path(TABLES_DIR, "orientation_enumeration_ate_deterministic_8node.csv"),
        " (", nrow(det_8node_out), " rows: 16 specs x 8 nodes).")

cat("\n(Quick check -- baseline (combo_0) values for the two new nodes:)\n")
det_8node_out |> dplyr::filter(scenario == "combo_0", node %in% c("harm_future", "politics")) |> print()
