# participant bootstrap cis for the intervention effects, baseline
# orientation only. from 22_intervention_bootstrap_ci.R. this is the
# "1,000 participant bootstrap resamples" in results -- separate bootstrap
# from N_BOOT, uses N_BOOT_INTERVENTION instead so the two never get mixed up.


source("clean_pipeline_alt_will/05_scm_intervention_helpers.R")
suppressPackageStartupMessages({ library(furrr); library(future) })

stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

set.seed(BOOT_SEED_INTERVENTION)

intervene_nodes <- c("belief_concern", "harm_present", "harm_future",
                      "weather_risk_prep", "social_norms", "trust_science",
                      "policy_support")
intervene_targets <- build_intervene_targets(intervene_nodes, max_size = 3)
stopifnot(length(intervene_targets) == choose(7, 1) + choose(7, 2) + choose(7, 3))
# Politics added as its own length-1 target only -- never combined into a
# pair/triple, but still gets the same treatment as every other single.
intervene_targets[["politics"]] <- "politics"
stopifnot(length(intervene_targets) == 64)

# ---- Point estimate on the real (unresampled) sample -----------------------
model_syntax_pt <- build_lavaan_syntax(base_edges)
fit_pt <- lavaan::sem(model_syntax_pt, data = df_extended, estimator = "MLR", fixed.x = FALSE)
stopifnot(lavaan::lavInspect(fit_pt, "converged"))
point_ate <- purrr::imap_dbl(intervene_targets, ~ target_ate(fit_pt, base_edges, all_nodes, .x))
cat("Point-estimate ATEs (deterministic propagation, real sample), all 8 singles:\n")
print(round(point_ate[c(intervene_nodes, "politics")], 4))

# ---- Bootstrap loop (parallelized across resamples) ------------------------
n_row <- nrow(df_extended)

boot_one <- function(b) {
  set.seed(BOOT_SEED_INTERVENTION + b)
  idx  <- sample.int(n_row, n_row, replace = TRUE)
  df_b <- df_extended[idx, ]
  fit_b <- tryCatch(
    lavaan::sem(model_syntax_pt, data = df_b, estimator = "MLR", fixed.x = FALSE),
    error = function(e) NULL
  )
  converged <- !is.null(fit_b) &&
    isTRUE(tryCatch(lavaan::lavInspect(fit_b, "converged"), error = function(e) FALSE))
  if (!converged) return(setNames(rep(NA_real_, length(intervene_targets)), names(intervene_targets)))
  purrr::map_dbl(intervene_targets, ~ target_ate(fit_b, base_edges, all_nodes, .x))
}

plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))
boot_list <- furrr::future_map(
  seq_len(N_BOOT_INTERVENTION), boot_one,
  .options = furrr::furrr_options(seed = TRUE), .progress = TRUE
)
plan(sequential)

boot_ate <- do.call(rbind, boot_list)
rownames(boot_ate) <- NULL
n_converged <- sum(!is.na(boot_ate[, 1]))
cat("\nBootstrap complete:", n_converged, "/", N_BOOT_INTERVENTION, "resamples converged.\n")

# ---- 95% percentile CIs: all 8 singles + dynamically-chosen best pair/triple
pair_labels_pt    <- names(intervene_targets)[purrr::map_int(intervene_targets, length) == 2]
triple_labels_pt  <- names(intervene_targets)[purrr::map_int(intervene_targets, length) == 3]
best_pair_label   <- pair_labels_pt[which.max(point_ate[pair_labels_pt])]
best_triple_label <- triple_labels_pt[which.max(point_ate[triple_labels_pt])]
cat("\nBest pair (real-sample point estimate):", best_pair_label, "(", round(point_ate[[best_pair_label]], 4), ")\n")
cat("Best triple (real-sample point estimate):", best_triple_label, "(", round(point_ate[[best_triple_label]], 4), ")\n")

ci_targets <- c(intervene_nodes, "politics", best_pair_label, best_triple_label)
ci_table <- purrr::map_dfr(ci_targets, function(tgt) {
  vals <- boot_ate[, tgt]; vals <- vals[!is.na(vals)]
  tibble::tibble(target = tgt, point_estimate = round(point_ate[[tgt]], 4),
                 boot_mean = round(mean(vals), 4),
                 ci_lower_95 = round(quantile(vals, .025, names = FALSE), 4),
                 ci_upper_95 = round(quantile(vals, .975, names = FALSE), 4),
                 n_boot = length(vals))
})
cat("\n=== 95% bootstrap CIs (participant resampling, baseline orientation only) ===\n")
print(as.data.frame(ci_table), row.names = FALSE)

rank1_frequency <- function(labels) {
  sub <- boot_ate[, labels, drop = FALSE]
  ok  <- stats::complete.cases(sub)
  sub <- sub[ok, , drop = FALSE]
  winner <- labels[apply(sub, 1, which.max)]
  tibble::tibble(target = labels) |>
    dplyr::left_join(tibble::tibble(target = winner) |> dplyr::count(target, name = "times_ranked_first"), by = "target") |>
    dplyr::mutate(times_ranked_first = tidyr::replace_na(times_ranked_first, 0L),
                   pct_ranked_first = round(100 * times_ranked_first / nrow(sub), 1)) |>
    dplyr::arrange(dplyr::desc(pct_ranked_first))
}
pair_labels   <- names(intervene_targets)[purrr::map_int(intervene_targets, length) == 2]
triple_labels <- names(intervene_targets)[purrr::map_int(intervene_targets, length) == 3]
pair_rank_table   <- rank1_frequency(pair_labels)
triple_rank_table <- rank1_frequency(triple_labels)
cat("\n=== How often each pair ranks first ===\n"); print(as.data.frame(pair_rank_table), row.names = FALSE)
cat("\n=== How often each triple ranks first ===\n"); print(as.data.frame(triple_rank_table), row.names = FALSE)

write.csv(ci_table,          file.path(OUTPUT_DIR, "intervention_bootstrap_ci.csv"), row.names = FALSE)
write.csv(pair_rank_table,   file.path(OUTPUT_DIR, "intervention_bootstrap_pair_rank1_freq.csv"), row.names = FALSE)
write.csv(triple_rank_table, file.path(OUTPUT_DIR, "intervention_bootstrap_triple_rank1_freq.csv"), row.names = FALSE)
saveRDS(boot_ate, file.path(OUTPUT_DIR, "intervention_bootstrap_ate_matrix.rds"))
cat("\nWrote intervention_bootstrap_ci.csv, _pair_rank1_freq.csv, _triple_rank1_freq.csv,",
    "and the full boot_ate matrix to", OUTPUT_DIR, "\n")
cat("(belief_concern CI should read ~.183-.228, present harm ~.161-.229 -- matching",
    "the manuscript's Results section, before trusting anything else here.)\n")
