# ipw sensitivity for wave5 attrition, from 28_ipw_attrition_sensitivity.R.
# feeds supp S11 (weight distribution/ess/truncation, age smd .167->.030,
# largest ate change .0033 for belief_concern).
#
# 12 found age/income/urbanicity over the .10 smd bar even though the 8
# attitude nodes were fine, so this is one targeted ipw check, not a
# replacement for the primary complete-case analysis. not forcing ipw into
# fci/pc-stable either -- just the behavior-facing scm stage.


source("clean_pipeline_alt_will/05_scm_intervention_helpers.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(lavaan)
  library(igraph)
})

stopifnot(exists("df_main"), exists("df_extended"), exists("df_long"))
stopifnot(all(c("dem_age", "dem_educ", "dem_income", "dem_male",
                "dem_race_cat", "dem_urban_cat") %in% names(df_long)))

set.seed(2026)

retained_ids <- df_extended$participant_id
attrited_ids <- setdiff(df_main$participant_id, retained_ids)
stopifnot(length(retained_ids) == 870)
stopifnot(length(retained_ids) + length(attrited_ids) == nrow(df_main))

# --- Step 0: assemble the retention-model dataset (N = 1,987) ---------------
dem_participant <- df_long |>
  dplyr::group_by(participant_id) |>
  dplyr::summarise(
    age       = dplyr::first(na.omit(dem_age)),
    educ      = dplyr::first(na.omit(dem_educ)),
    income    = dplyr::first(na.omit(dem_income)),
    male      = dplyr::first(na.omit(dem_male)),
    race_cat  = dplyr::first(na.omit(as.character(dem_race_cat))),
    urban_cat = dplyr::first(na.omit(as.character(dem_urban_cat))),
    .groups = "drop"
  ) |>
  dplyr::filter(participant_id %in% df_main$participant_id)

race_tab <- table(dem_participant$race_cat)
rare_levels <- names(race_tab)[race_tab < 20]
if (length(rare_levels) > 0) {
  cat("Collapsing rare race_cat levels into 'other' (n < 20 in the N=1,987 sample): ",
      paste(rare_levels, collapse = ", "), "\n", sep = "")
  dem_participant <- dem_participant |>
    dplyr::mutate(race_cat = ifelse(race_cat %in% rare_levels, "other", race_cat))
} else {
  cat("No race_cat levels below the n=20 threshold; using all categories as-is.\n")
}

dem_participant <- dem_participant |>
  dplyr::mutate(
    gender_cat = dplyr::case_when(male == 1 ~ "male", male == 0 ~ "female", TRUE ~ "other"),
    R = as.integer(participant_id %in% retained_ids)
  )

retention_df <- dem_participant |>
  dplyr::inner_join(df_main |> dplyr::select(participant_id, dplyr::all_of(NODE_ORDER_MAIN)),
                     by = "participant_id") |>
  dplyr::mutate(
    educ_f   = factor(educ),
    income_f = factor(income),
    race_f   = factor(race_cat),
    urban_f  = factor(urban_cat),
    gender_f = factor(gender_cat)
  )

cat("\nRetention model sample: N =", nrow(retention_df), "( retained =",
    sum(retention_df$R), ", attrited =", sum(1 - retention_df$R), ")\n")
stopifnot(nrow(retention_df) == nrow(df_main))

# --- Step 1: propensity (retention) model -----------------------------------
node_rhs <- paste(NODE_ORDER_MAIN, collapse = " + ")
retention_formula <- stats::as.formula(
  paste("R ~ age + gender_f + educ_f + income_f + urban_f + race_f +", node_rhs)
)
ret_fit <- stats::glm(retention_formula, data = retention_df, family = stats::binomial())
cat("\n=== Retention (propensity) model: R = 1 if has Wave-5 behavior ===\n")
print(summary(ret_fit)$coefficients)

retention_df$p_hat <- stats::fitted(ret_fit)
cat("\nPredicted P(retained) range: [", round(min(retention_df$p_hat), 4), ", ",
    round(max(retention_df$p_hat), 4), "]\n", sep = "")

# --- Step 2: stabilized IPW weights (retained participants only) -----------
p_marginal <- mean(retention_df$R)
retention_df <- retention_df |>
  dplyr::mutate(sw = ifelse(R == 1, p_marginal / p_hat, NA_real_))

# --- Step 3a: propensity overlap diagnostic ---------------------------------
cat("\n=== (3a) Propensity overlap: P(retained) by actual retention status ===\n")
overlap_summary <- retention_df |>
  dplyr::group_by(R) |>
  dplyr::summarise(
    n      = dplyr::n(),
    min    = round(min(p_hat), 3),
    p25    = round(stats::quantile(p_hat, .25), 3),
    median = round(stats::median(p_hat), 3),
    p75    = round(stats::quantile(p_hat, .75), 3),
    max    = round(max(p_hat), 3),
    .groups = "drop"
  )
print(as.data.frame(overlap_summary), row.names = FALSE)
n_extreme <- sum(retention_df$p_hat < .05 | retention_df$p_hat > .95)
cat("Participants with P(retained) outside [.05, .95]: ", n_extreme,
    " of ", nrow(retention_df), "\n", sep = "")

# --- Step 3b: weight distribution, ESS, truncation --------------------------
ess <- function(w) sum(w)^2 / sum(w^2)

sw_retained <- retention_df$sw[retention_df$R == 1]
cat("\n=== (3b) Stabilized weight distribution (N = 870 retained) ===\n")
print(summary(sw_retained))
cat("Effective sample size (untruncated): ", round(ess(sw_retained), 1),
    " of N = ", length(sw_retained), "\n", sep = "")

trunc_lo <- stats::quantile(sw_retained, .01)
trunc_hi <- stats::quantile(sw_retained, .99)
n_trunc <- sum(sw_retained < trunc_lo | sw_retained > trunc_hi)
sw_trunc <- pmin(pmax(sw_retained, trunc_lo), trunc_hi)
cat("Truncating at [1st, 99th] percentile: [", round(trunc_lo, 3), ", ",
    round(trunc_hi, 3), "] -- ", n_trunc, " weights truncated.\n", sep = "")
cat("Effective sample size (truncated): ", round(ess(sw_trunc), 1), "\n", sep = "")
cat("Manuscript (S11) reports the truncated weights' maximum as 1.97 -- compare.\n")

retention_df$sw_trunc <- NA_real_
retention_df$sw_trunc[retention_df$R == 1] <- sw_trunc

weight_summary <- tibble::tibble(
  stat = c("min", "p25", "median", "mean", "p75", "max", "ess", "n_truncated"),
  untruncated = c(round(min(sw_retained), 3), round(stats::quantile(sw_retained, .25), 3),
                  round(stats::median(sw_retained), 3), round(mean(sw_retained), 3),
                  round(stats::quantile(sw_retained, .75), 3), round(max(sw_retained), 3),
                  round(ess(sw_retained), 1), NA_real_),
  truncated = c(round(min(sw_trunc), 3), round(stats::quantile(sw_trunc, .25), 3),
                round(stats::median(sw_trunc), 3), round(mean(sw_trunc), 3),
                round(stats::quantile(sw_trunc, .75), 3), round(max(sw_trunc), 3),
                round(ess(sw_trunc), 1), n_trunc)
)

# --- Step 3c: balance after weighting ---------------------------------------
wtd_mean <- function(x, w) sum(x * w) / sum(w)
wtd_sd   <- function(x, w) sqrt(sum(w * (x - wtd_mean(x, w))^2) / sum(w))

balance_vars <- c(NODE_ORDER_MAIN, "age")
retained_rows <- retention_df |> dplyr::filter(R == 1)
full_rows     <- retention_df  # N = 1,987

balance_check <- purrr::map_dfr(balance_vars, function(v) {
  x_ret  <- retained_rows[[v]]
  x_full <- full_rows[[v]]
  m1u <- mean(x_ret, na.rm = TRUE);  s1u <- sd(x_ret, na.rm = TRUE)
  m2  <- mean(x_full, na.rm = TRUE); s2  <- sd(x_full, na.rm = TRUE)
  smd_unweighted <- (m1u - m2) / sqrt((s1u^2 + s2^2) / 2)

  m1w <- wtd_mean(x_ret, retained_rows$sw_trunc)
  s1w <- wtd_sd(x_ret, retained_rows$sw_trunc)
  smd_ipw <- (m1w - m2) / sqrt((s1w^2 + s2^2) / 2)

  tibble::tibble(variable = v,
                 smd_unweighted = round(smd_unweighted, 3),
                 smd_ipw_weighted = round(smd_ipw, 3),
                 improved = abs(smd_ipw) < abs(smd_unweighted))
})
cat("\n=== (3c) Balance: retained sample vs. full N=1,987 sample ===\n")
cat("(unweighted vs. IPW-weighted; population being restored is all N=1,987\n",
    "Waves-1-4 participants, not just the attrited group)\n", sep = "")
print(as.data.frame(balance_check), row.names = FALSE)
cat("Manuscript (S11) reports age SMD=.167 before weighting, .030 after -- compare.\n")

# --- Step 4: refit the behavior-facing part of the working SCM -------------
df_extended_w <- df_extended |>
  dplyr::inner_join(
    retention_df |> dplyr::filter(R == 1) |> dplyr::select(participant_id, sw_trunc),
    by = "participant_id"
  )
stopifnot(nrow(df_extended_w) == nrow(df_extended))
stopifnot(!anyNA(df_extended_w$sw_trunc))

extract_eq <- function(fit, label) {
  lavaan::standardizedSolution(fit) |>
    dplyr::filter(op == "~") |>
    dplyr::transmute(predictor = rhs, beta = round(est.std, 3), se = round(se, 3),
                      p = round(pvalue, 4), model = label)
}

# (4a) Behavior equation alone
eq_syntax <- "climate_behavior ~ harm_present + weather_risk_prep + social_norms"
fit_eq_unw <- lavaan::sem(eq_syntax, data = df_extended_w, estimator = "MLR", fixed.x = FALSE)
fit_eq_ipw <- lavaan::sem(eq_syntax, data = df_extended_w, estimator = "MLR", fixed.x = FALSE,
                           sampling.weights = "sw_trunc")
stopifnot(lavaan::lavInspect(fit_eq_unw, "converged"), lavaan::lavInspect(fit_eq_ipw, "converged"))

eq_compare <- dplyr::bind_rows(extract_eq(fit_eq_unw, "unweighted"),
                                extract_eq(fit_eq_ipw, "ipw_weighted"))
cat("\n=== (4a) Behavior equation: unweighted vs. IPW-weighted ===\n")
print(as.data.frame(eq_compare), row.names = FALSE)

# (4b) Full 16-edge SCM, using the SHARED base_edges/helpers from
# 05_scm_intervention_helpers.R (not a local redefinition).
intervene_nodes <- c("belief_concern", "harm_present", "weather_risk_prep",
                      "social_norms", "trust_science", "policy_support")
model_syntax <- build_lavaan_syntax(base_edges)
fit_full_unw <- lavaan::sem(model_syntax, data = df_extended_w, estimator = "MLR", fixed.x = FALSE)
fit_full_ipw <- lavaan::sem(model_syntax, data = df_extended_w, estimator = "MLR", fixed.x = FALSE,
                             sampling.weights = "sw_trunc")
stopifnot(lavaan::lavInspect(fit_full_unw, "converged"), lavaan::lavInspect(fit_full_ipw, "converged"))

ates_unw <- single_node_ates(fit_full_unw, base_edges, all_nodes, intervene_nodes)
ates_ipw <- single_node_ates(fit_full_ipw, base_edges, all_nodes, intervene_nodes)

ate_compare <- tibble::tibble(
  node              = intervene_nodes,
  ate_unweighted    = round(unname(ates_unw[intervene_nodes]), 4),
  ate_ipw_weighted  = round(unname(ates_ipw[intervene_nodes]), 4)
) |>
  dplyr::mutate(diff = round(ate_ipw_weighted - ate_unweighted, 4))

cat("\n=== (4b/5) Six single-node intervention effects: unweighted vs. IPW-weighted ===\n")
print(as.data.frame(ate_compare), row.names = FALSE)
cat("Manuscript (S11) reports the largest absolute change as .0033 SD for\n")
cat("belief/concern, all others <= .0006 SD -- compare the 'diff' column above.\n")

# --- Write outputs -----------------------------------------------------------
write.csv(overlap_summary, file.path(OUTPUT_DIR, "ipw_overlap_summary.csv"), row.names = FALSE)
write.csv(weight_summary,  file.path(OUTPUT_DIR, "ipw_weight_distribution.csv"), row.names = FALSE)
write.csv(balance_check,   file.path(OUTPUT_DIR, "ipw_balance_check.csv"), row.names = FALSE)
write.csv(eq_compare,      file.path(OUTPUT_DIR, "ipw_behavior_equation_compare.csv"), row.names = FALSE)
write.csv(ate_compare,     file.path(OUTPUT_DIR, "ipw_single_node_ate_compare.csv"), row.names = FALSE)

cat("\nWrote ipw_{overlap_summary,weight_distribution,balance_check,",
    "behavior_equation_compare,single_node_ate_compare}.csv to", OUTPUT_DIR, "\n", sep = "")
cat("\nESS (truncated) = ", round(ess(sw_trunc), 1), " of N = ", length(sw_trunc),
    "; ", n_trunc, " weights truncated at [1st,99th] percentile.\n", sep = "")
