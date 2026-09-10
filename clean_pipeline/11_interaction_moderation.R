# politics x belief_concern interaction on policy_support, plus the shift-
# intervention numbers (results "additional check" + supp S9). from
# 11_ida_appendix_and_interaction_test.R, parts B/C/D (part A was already
# dropped from the ms).
#
# builds scm_coefs/resid_sd from base_edges in 05 now instead of its own
# hardcoded copy -- the old copy had silently missed the
# harm_future->trust_science and belief_concern->weather_risk_prep edges
# from the last scm revision, which is what made table S9's climate_behavior
# numbers go stale. building off the shared base_edges means it can't happen
# again next time the scm changes.


source("clean_pipeline/05_scm_intervention_helpers.R")

stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Baseline (no-interaction) 16-edge fit, used both for scm_coefs/
# scm_resid_sd (the simulator's parameters) and as the model this script's
# interaction test is added on top of. Same fitting convention as 06/07/08
# (estimator="MLR", fixed.x=FALSE).
model_syntax_base <- build_lavaan_syntax(base_edges)
fit_base <- lavaan::sem(model_syntax_base, data = df_extended, estimator = "MLR", fixed.x = FALSE)
stopifnot(lavaan::lavInspect(fit_base, "converged"))
std_base <- lavaan::standardizedSolution(fit_base) |> dplyr::filter(op == "~")

# ---- scm_coefs: one named coefficient per edge, using the same short
# abbreviations as climate_analysis_avg_v2_altweather.qmd's own scm_coefs
# list (Section 7.6/7.7) and r_patches/11's p$<name> references inside
# simulate_scm_int()/simulate_scm_shift() -- built programmatically from
# base_edges instead of copy-pasted, so it can't independently go stale.
edge_abbrev <- tibble::tribble(
  ~from,               ~to,                 ~abbrev,
  "politics",          "belief_concern",    "pol_bc",
  "belief_concern",    "harm_future",       "bc_hf",
  "belief_concern",    "harm_present",      "bc_hp",
  "harm_future",       "harm_present",      "hf_hp",
  "belief_concern",    "trust_science",     "bc_ts",
  "harm_future",       "trust_science",     "hf_ts",
  "belief_concern",    "policy_support",    "bc_ps",
  "trust_science",     "policy_support",    "ts_ps",
  "politics",          "policy_support",    "pol_ps",
  "policy_support",    "social_norms",      "ps_sn",
  "trust_science",     "social_norms",      "ts_sn",
  "harm_present",      "weather_risk_prep", "hp_wr",
  "belief_concern",    "weather_risk_prep", "bc_wr",
  "harm_present",      "climate_behavior",  "hp_cb",
  "social_norms",      "climate_behavior",  "sn_cb",
  "weather_risk_prep", "climate_behavior",  "wr_cb"
)
# Guard: edge_abbrev must cover exactly base_edges (direction-sensitive,
# order-insensitive) -- if base_edges ever changes, this stop()s instead of
# silently mapping the wrong coefficient to the wrong name.
key <- function(df) sort(paste(df$from, df$to, sep = "->"))
stopifnot(identical(key(base_edges), key(edge_abbrev[, c("from","to")])))

scm_coefs <- setNames(as.list(rep(NA_real_, nrow(edge_abbrev))), edge_abbrev$abbrev)
for (i in seq_len(nrow(edge_abbrev))) {
  row <- edge_abbrev[i, ]
  val <- std_base$est.std[std_base$lhs == row$to & std_base$rhs == row$from]
  stopifnot(length(val) == 1)
  scm_coefs[[row$abbrev]] <- val
}
cat("\n--- Baseline SCM path coefficients (scm_coefs) ---\n")
purrr::iwalk(scm_coefs, ~cat(sprintf("  %-10s  %.4f\n", .y, .x)))

# ---- scm_resid_sd: politics is exogenous (SD=1 on the standardized scale);
# every endogenous node's residual SD is sqrt(1 - R^2) from fit_base's own
# equation-level R^2 -- same convention as Supplementary S7.
r2_base <- lavaan::lavInspect(fit_base, "rsquare")
endogenous_nodes <- setdiff(all_nodes, "politics")
scm_resid_sd <- c(politics = 1.0, setNames(sqrt(1 - r2_base[endogenous_nodes]), endogenous_nodes))
cat("\n--- Residual SDs (scm_resid_sd) ---\n")
print(round(scm_resid_sd, 4))

# =============================================================================
# PART B -- interaction significance test (mean-centered raw composites,
# group-restandardized) -- verbatim method from r_patches/11 Part B.
# =============================================================================
dat_int <- as.data.frame(df_extended[, all_nodes])
dat_int$bc_c  <- dat_int$belief_concern - mean(dat_int$belief_concern, na.rm = TRUE)
dat_int$pol_c <- dat_int$politics       - mean(dat_int$politics,       na.rm = TRUE)
dat_int$bc_x_politics <- dat_int$bc_c * dat_int$pol_c

model_lines_int <- base_edges |>
  dplyr::group_by(to) |>
  dplyr::summarise(rhs = paste(from, collapse = " + "), .groups = "drop") |>
  dplyr::mutate(rhs = ifelse(to == "policy_support", paste(rhs, "+ bc_x_politics"), rhs),
                line = paste(to, "~", rhs))
model_syntax_int <- paste(model_lines_int$line, collapse = "\n")
cat("\nInteraction-augmented lavaan model:\n", model_syntax_int, "\n")

sem_fit_int <- lavaan::sem(
  model_syntax_int,
  data = as.data.frame(scale(dat_int[, c(all_nodes, "bc_x_politics")])),
  estimator = "MLR"
)
cat("\nFit indices (interaction model):\n")
print(round(lavaan::fitMeasures(sem_fit_int, c("cfi","tli","rmsea","srmr","aic","bic")), 3))

std_int <- lavaan::standardizedSolution(sem_fit_int)
int_row <- std_int[std_int$lhs == "policy_support" & std_int$rhs == "bc_x_politics", ]
cat("\npolitics x belief_concern -> policy_support (Part B, group-restandardized):\n")
print(int_row[, c("lhs","rhs","est.std","se","pvalue")])
cat("(Manuscript reports beta_BP=.080, SE=.022, p<.001 -- compare directly.)\n")

lm_check <- lm(policy_support ~ belief_concern + trust_science + politics + bc_x_politics,
               data = as.data.frame(scale(dat_int[, c(all_nodes, "bc_x_politics")])))
cat("\nPart B lm() cross-check (manuscript reports beta=.0815):\n")
print(summary(lm_check)$coefficients)

# =============================================================================
# PART C -- "clean" single-standardization interaction coefficient -- this is
# the one that actually feeds the simulator (p_int$bc_pol_ps below).
# Manuscript reports beta=.0802 for this version.
# =============================================================================
z_bc  <- as.numeric(scale(df_extended$belief_concern))
z_pol <- as.numeric(scale(df_extended$politics))
bc_x_politics_clean <- z_bc * z_pol

dat_clean <- as.data.frame(scale(df_extended[, all_nodes]))
dat_clean$bc_x_politics <- bc_x_politics_clean

sem_fit_clean <- lavaan::sem(model_syntax_int, data = dat_clean, estimator = "MLR")
std_clean <- lavaan::standardizedSolution(sem_fit_clean)
int_row_clean <- std_clean[std_clean$lhs == "policy_support" & std_clean$rhs == "bc_x_politics", ]
cat("\nPart C clean (single-standardization) interaction coefficient (manuscript: .0802):\n")
print(int_row_clean[, c("lhs","rhs","est.std","se","pvalue")])
stopifnot(nrow(int_row_clean) == 1)

p_int <- scm_coefs
p_int$bc_pol_ps <- int_row_clean$est.std[1]
cat("\nInteraction coefficient feeding the simulator (bc_pol_ps):", round(p_int$bc_pol_ps, 4), "\n")

interaction_summary <- tibble::tibble(
  version = c("part_b_group_restandardized", "part_b_lm_crosscheck", "part_c_clean_single_standardization"),
  beta = c(int_row$est.std[1], unname(coef(lm_check)["bc_x_politics"]), int_row_clean$est.std[1]),
  manuscript_value = c(.080, .0815, .0802)
)
cat("\n--- Interaction coefficient summary (compare 'beta' to 'manuscript_value') ---\n")
print(as.data.frame(interaction_summary), row.names = FALSE)
write.csv(interaction_summary, file.path(TABLES_DIR, "interaction_coefficient.csv"), row.names = FALSE)

# =============================================================================
# PART D -- politics-conditional SHIFT intervention on belief_concern. This
# is the primary comparison the manuscript reports (Results paragraph;
# Supplementary S9, Table~\ref{tab:supp_interaction}): a shift, not a common
# absolute target, because politics -> belief_concern is a real fitted path
# and a common do(belief_concern=0.5) would give the two political groups
# very different-sized treatments. Verbatim method from r_patches/11 Part D
# (already includes the hf_ts/bc_wr patch), parameterized from p_int (built
# above from base_edges) instead of a copy hand-edited in place.
# =============================================================================
simulate_scm_shift <- function(n = 10000, politics_level = 0, belief_shift = 0, seed = 42) {
  set.seed(seed)
  p <- p_int
  pol <- rep(politics_level, n)
  bc  <- p$pol_bc * pol + rnorm(n, 0, scm_resid_sd["belief_concern"]) + belief_shift
  hf  <- p$bc_hf  * bc  + rnorm(n, 0, scm_resid_sd["harm_future"])
  ts  <- p$bc_ts  * bc  + p$hf_ts * hf + rnorm(n, 0, scm_resid_sd["trust_science"])
  hp  <- p$hf_hp * hf + p$bc_hp * bc + rnorm(n, 0, scm_resid_sd["harm_present"])
  ps  <- p$ts_ps * ts + p$bc_ps * bc + p$pol_ps * pol +
           p$bc_pol_ps * bc * pol + rnorm(n, 0, scm_resid_sd["policy_support"])
  wr  <- p$hp_wr * hp + p$bc_wr * bc + rnorm(n, 0, scm_resid_sd["weather_risk_prep"])
  sn  <- p$ts_sn * ts + p$ps_sn * ps + rnorm(n, 0, scm_resid_sd["social_norms"])
  cb  <- p$hp_cb * hp + p$wr_cb * wr + p$sn_cb * sn + rnorm(n, 0, scm_resid_sd["climate_behavior"])
  tibble::tibble(politics = pol, belief_concern = bc, harm_future = hf,
                 trust_science = ts, harm_present = hp, policy_support = ps,
                 weather_risk_prep = wr, social_norms = sn, climate_behavior = cb)
}

N_SIM_SHIFT <- 50000
shift_scenarios <- list(
  low_politics_baseline  = list(politics_level = -1, belief_shift = 0),
  low_politics_shift     = list(politics_level = -1, belief_shift = 0.5),
  high_politics_baseline = list(politics_level =  1, belief_shift = 0),
  high_politics_shift    = list(politics_level =  1, belief_shift = 0.5)
)
shift_results <- purrr::imap_dfr(shift_scenarios, function(s, label) {
  dat <- simulate_scm_shift(n = N_SIM_SHIFT, politics_level = s$politics_level,
                             belief_shift = s$belief_shift, seed = 42)
  tibble::tibble(scenario = label,
                 belief_concern_mean   = mean(dat$belief_concern),
                 policy_support_mean   = mean(dat$policy_support),
                 climate_behavior_mean = mean(dat$climate_behavior))
})
cat("\n=== Monte Carlo shift results (n=50,000 per scenario) ===\n")
print(as.data.frame(shift_results), row.names = FALSE)

shift_delta_ps_low  <- shift_results$policy_support_mean[shift_results$scenario == "low_politics_shift"] -
                       shift_results$policy_support_mean[shift_results$scenario == "low_politics_baseline"]
shift_delta_ps_high <- shift_results$policy_support_mean[shift_results$scenario == "high_politics_shift"] -
                       shift_results$policy_support_mean[shift_results$scenario == "high_politics_baseline"]
shift_delta_cb_low  <- shift_results$climate_behavior_mean[shift_results$scenario == "low_politics_shift"] -
                       shift_results$climate_behavior_mean[shift_results$scenario == "low_politics_baseline"]
shift_delta_cb_high <- shift_results$climate_behavior_mean[shift_results$scenario == "high_politics_shift"] -
                       shift_results$climate_behavior_mean[shift_results$scenario == "high_politics_baseline"]

cat(sprintf("\nSHIFT belief_concern by +0.5 SD -> policy_support: %.3f SD (low politics) vs %.3f SD (high politics)\n",
            shift_delta_ps_low, shift_delta_ps_high))
cat(sprintf("Same, propagated to climate_behavior: %.3f SD (low) vs %.3f SD (high)\n",
            shift_delta_cb_low, shift_delta_cb_high))
cat("Manuscript (Table S9) reports .322 (low) vs .402 (high) for policy_support,\n")
cat("and .202 (low) vs .205 (high) for climate_behavior -- compare directly.\n")

# ---- exact closed-form cross-check, catches any monte carlo noise --------
# politics is held at a fixed constant within each scenario here (unlike the
# free/random politics used in the primary do(X=0.5) ATE table), so
# bc_pol_ps * belief_concern * politics_level is linear in belief_concern
# for fixed politics_level and the whole chain stays linear-recursive: its
# mean is exactly computable by propagation, with no simulation noise or
# seed dependence. This should match the Monte Carlo numbers above to ~3
# decimal places; if it doesn't, STOP and investigate before trusting either.
exact_shift_mean <- function(politics_level, belief_shift) {
  p <- p_int
  e_bc <- p$pol_bc * politics_level + belief_shift
  e_hf <- p$bc_hf * e_bc
  e_ts <- p$bc_ts * e_bc + p$hf_ts * e_hf
  e_hp <- p$hf_hp * e_hf + p$bc_hp * e_bc
  e_ps <- p$ts_ps * e_ts + p$bc_ps * e_bc + p$pol_ps * politics_level +
            p$bc_pol_ps * e_bc * politics_level
  e_wr <- p$hp_wr * e_hp + p$bc_wr * e_bc
  e_sn <- p$ts_sn * e_ts + p$ps_sn * e_ps
  e_cb <- p$hp_cb * e_hp + p$wr_cb * e_wr + p$sn_cb * e_sn
  c(belief_concern = e_bc, policy_support = e_ps, climate_behavior = e_cb)
}

exact_results <- purrr::imap_dfr(shift_scenarios, function(s, label) {
  m <- exact_shift_mean(s$politics_level, s$belief_shift)
  tibble::tibble(scenario = label, belief_concern_mean = m[["belief_concern"]],
                 policy_support_mean = m[["policy_support"]],
                 climate_behavior_mean = m[["climate_behavior"]])
})
cat("\n=== Exact closed-form shift results (cross-check, no simulation) ===\n")
print(as.data.frame(exact_results), row.names = FALSE)

max_diff <- max(abs(shift_results$policy_support_mean - exact_results$policy_support_mean),
                abs(shift_results$climate_behavior_mean - exact_results$climate_behavior_mean),
                abs(shift_results$belief_concern_mean - exact_results$belief_concern_mean))
cat(sprintf("\nMax |Monte Carlo - exact| across all scenarios/outcomes: %.5f\n", max_diff))
if (max_diff > 0.01) {
  stop("Monte Carlo shift simulation and exact closed-form propagation disagree ",
       "by more than 0.01 SD -- investigate before trusting either result.")
} else {
  cat("Monte Carlo and exact closed-form agree (as expected for a linear-recursive\n")
  cat("model with a fixed politics level) -- shift_results above is not an\n")
  cat("artifact of simulation noise.\n")
}

# Closed-form marginal-slope sanity check (same as r_patches/11 Part D).
slope_low  <- p_int$bc_ps + p_int$bc_pol_ps * (-1)
slope_high <- p_int$bc_ps + p_int$bc_pol_ps * ( 1)
cat(sprintf("\nMarginal slope d(policy_support)/d(belief_concern) = bc_ps + bc_pol_ps*politics:\n"))
cat(sprintf("  at politics = -1 (conservative): %.4f\n", slope_low))
cat(sprintf("  at politics = +1 (liberal):      %.4f\n", slope_high))

write.csv(shift_results, file.path(TABLES_DIR, "interaction_shift_results_montecarlo.csv"), row.names = FALSE)
write.csv(exact_results, file.path(TABLES_DIR, "interaction_shift_results_exact.csv"), row.names = FALSE)

cat("\n===== 11_interaction_moderation.R done =====\n")
cat("check shift_results/exact_results above against main9.tex table S9\n")
cat("(.322/.402 policy_support, .205/.207 climate_behavior) and\n")
cat("interaction_coefficient.csv against .080/.0815/.0802.\n")
