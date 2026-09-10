# fit indices (cfi/tli/rmsea/srmr/aic/bic) for all 16 scenarios, from
# 02_full_orientation_enumeration_v4.R's fit-index half. that script's ATE
# half used the old monte carlo sim with the known noise bug -- 06/07/08
# replaced that with exact propagation, so this script doesn't touch ates
# at all, just fit.


source("clean_pipeline_alt_will/05_scm_intervention_helpers.R")

stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

full_fit <- scenario_list |>
  purrr::imap(function(flip_idx, scen_name) {
    e <- flip_edges(base_edges, flip_idx)
    n_flipped <- length(flip_idx)

    if (!is_acyclic(e)) {
      message(scen_name, " (", scenario_flip_labels[[scen_name]], "): cyclic, skipped.")
      return(tibble::tibble(
        scenario = scen_name, n_flipped = n_flipped,
        flipped_edges = scenario_flip_labels[[scen_name]], status = "cyclic_skipped",
        cfi = NA_real_, tli = NA_real_, rmsea = NA_real_, srmr = NA_real_,
        aic = NA_real_, bic = NA_real_
      ))
    }

    model_syntax <- build_lavaan_syntax(e)
    fit <- tryCatch(
      lavaan::sem(model_syntax, data = df_extended, estimator = "MLR", fixed.x = FALSE),
      error = function(err) NULL
    )
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      message(scen_name, " (", scenario_flip_labels[[scen_name]], "): fit failed/non-converged, skipped")
      return(tibble::tibble(
        scenario = scen_name, n_flipped = n_flipped,
        flipped_edges = scenario_flip_labels[[scen_name]], status = "fit_failed",
        cfi = NA_real_, tli = NA_real_, rmsea = NA_real_, srmr = NA_real_,
        aic = NA_real_, bic = NA_real_
      ))
    }

    fm <- lavaan::fitmeasures(fit)
    get_fm <- function(nms) {
      hit <- nms[nms %in% names(fm)]
      if (length(hit) == 0) return(NA_real_)
      unname(fm[hit[1]])
    }
    tibble::tibble(
      scenario = scen_name, n_flipped = n_flipped,
      flipped_edges = scenario_flip_labels[[scen_name]], status = "ok",
      cfi   = get_fm(c("cfi.robust", "cfi")),
      tli   = get_fm(c("tli.robust", "tli")),
      rmsea = get_fm(c("rmsea.robust", "rmsea")),
      srmr  = get_fm(c("srmr")),
      aic   = get_fm(c("aic")),
      bic   = get_fm(c("bic"))
    )
  }) |>
  dplyr::bind_rows()

cat("\n--- Scenarios attempted:", length(scenario_list),
    "| acyclic+converged:", sum(full_fit$status == "ok"),
    "| cyclic (excluded):", sum(full_fit$status == "cyclic_skipped"),
    "| fit failed:", sum(full_fit$status == "fit_failed"), "---\n")
print(full_fit)

write.csv(full_fit, file.path(TABLES_DIR, "orientation_enumeration_fit.csv"), row.names = FALSE)

fit_supplement_table <- full_fit |>
  dplyr::filter(status == "ok") |>
  dplyr::arrange(n_flipped, scenario) |>
  dplyr::select(scenario, flipped_edges, n_flipped, cfi, tli, rmsea, srmr, aic, bic)

cat("\n--- Supplementary fit table (Overleaf-ready) ---\n")
print(fit_supplement_table)
write.csv(fit_supplement_table, file.path(TABLES_DIR, "fit_supplement_table.csv"), row.names = FALSE)

cat("\nWrote", file.path(TABLES_DIR, "orientation_enumeration_fit.csv"), "and",
    file.path(TABLES_DIR, "fit_supplement_table.csv"), "\n")
cat("(The baseline/combo_0 row should read CFI=.979, TLI=.962, RMSEA=.086, SRMR=.039",
    "-- matching the manuscript's reported 16-edge model fit -- before trusting",
    "anything else here.)\n")
