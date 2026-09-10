# =============================================================================
# 26_intervention_ci_followup.R
#
# Follow-up to script 22, prompted by the team's precision request on the
# bootstrap-CI wording: (1) get point estimates + 95% CIs for the OTHER two
# pairs the rank-1 frequency table showed tied with belief_concern+
# weather_risk_prep (harm_present+social_norms, harm_present+weather_risk_prep)
# -- script 22 only saved CI summary stats for the single pair that happened
# to win on point estimate; (2) compute the bootstrap CI for the DIFFERENCE
# belief_concern - harm_present directly (paired within each of the 1,000
# resamples), rather than just eyeballing whether the two marginal CIs
# overlap, per the team's explicit request.
#
# No new bootstrap -- reuses pipeline_outputs/intervention_bootstrap_ate_matrix.rds
# (the full 1000 x 41 target matrix already saved by script 22) for the CIs, and
# refits the baseline model once on the real (unresampled) sample for point
# estimates. Requires df_extended in the session (same precondition as 22).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lavaan)
  library(igraph)
  library(purrr)
})

stopifnot(exists("df_extended"))
mat_path <- "pipeline_outputs/intervention_bootstrap_ate_matrix.rds"
if (!file.exists(mat_path)) stop("Run script 22 first -- ", mat_path, " not found.")
boot_ate <- readRDS(mat_path)

# --- Re-derive point estimates on the real sample (same base_edges/helpers as
# script 22) for the two additional pairs. -----------------------------------
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
build_lavaan_syntax <- function(edges) {
  edges |> dplyr::group_by(to) |>
    dplyr::summarise(rhs = paste(from, collapse = " + "), .groups = "drop") |>
    dplyr::mutate(line = paste(to, "~", rhs)) |> dplyr::pull(line) |> paste(collapse = "\n")
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
target_ate <- function(fit, edges, all_nodes, target_combo) {
  base_cb <- scm_mean_propagate(fit, edges, all_nodes, list())[["climate_behavior"]]
  do_cb <- scm_mean_propagate(fit, edges, all_nodes,
             setNames(as.list(rep(0.5, length(target_combo))), target_combo))[["climate_behavior"]]
  do_cb - base_cb
}

fit_pt <- lavaan::sem(build_lavaan_syntax(base_edges), data = df_extended,
                       estimator = "MLR", fixed.x = FALSE)
stopifnot(lavaan::lavInspect(fit_pt, "converged"))

extra_pairs <- list(
  "harm_present+social_norms"      = c("harm_present", "social_norms"),
  "harm_present+weather_risk_prep" = c("harm_present", "weather_risk_prep")
)
point_extra <- purrr::map_dbl(extra_pairs, ~ target_ate(fit_pt, base_edges, all_nodes, .x))

# --- CIs for the three tied pairs, from the already-saved bootstrap matrix --
tied_pairs <- c("belief_concern+weather_risk_prep",
                names(extra_pairs))
stopifnot(all(tied_pairs %in% colnames(boot_ate)))

pair_ci <- purrr::map_dfr(tied_pairs, function(tgt) {
  vals <- boot_ate[, tgt]; vals <- vals[!is.na(vals)]
  pt <- if (tgt %in% names(point_extra)) point_extra[[tgt]] else NA_real_
  tibble::tibble(
    target = tgt,
    point_estimate = round(if (is.na(pt)) mean(vals) else pt, 4),  # belief+weather's point est already in script 22's ci_table
    ci_lower_95 = round(quantile(vals, .025, names = FALSE), 4),
    ci_upper_95 = round(quantile(vals, .975, names = FALSE), 4)
  )
})
cat("=== Three tied pairs: point estimate + 95% CI ===\n")
print(as.data.frame(pair_ci), row.names = FALSE)

# --- Bootstrap CI on the DIFFERENCE belief_concern - harm_present (paired) --
stopifnot(all(c("belief_concern", "harm_present") %in% colnames(boot_ate)))
diff_vals <- boot_ate[, "belief_concern"] - boot_ate[, "harm_present"]
diff_vals <- diff_vals[!is.na(diff_vals)]
diff_ci <- tibble::tibble(
  comparison = "belief_concern - harm_present",
  point_diff = round(0.2059 - 0.1951, 4),  # from script 22's already-reported point estimates
  boot_mean_diff = round(mean(diff_vals), 4),
  ci_lower_95 = round(quantile(diff_vals, .025, names = FALSE), 4),
  ci_upper_95 = round(quantile(diff_vals, .975, names = FALSE), 4),
  ci_excludes_zero = round(quantile(diff_vals, .025, names = FALSE), 4) > 0 ||
                      round(quantile(diff_vals, .975, names = FALSE), 4) < 0
)
cat("\n=== Bootstrap CI on the difference (belief_concern - harm_present) ===\n")
print(as.data.frame(diff_ci), row.names = FALSE)

dir.create("pipeline_outputs", showWarnings = FALSE)
write.csv(pair_ci, "pipeline_outputs/intervention_bootstrap_tied_pairs_ci.csv", row.names = FALSE)
write.csv(diff_ci, "pipeline_outputs/intervention_bootstrap_belief_vs_harm_diff_ci.csv", row.names = FALSE)
cat("\nWrote pipeline_outputs/intervention_bootstrap_tied_pairs_ci.csv and ",
    "intervention_bootstrap_belief_vs_harm_diff_ci.csv\n", sep = "")
