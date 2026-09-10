# =============================================================================
# 19_systematic_interaction_screen.R
#
# Raised (2026-09-04): the paper currently tests exactly ONE interaction
# (politics x belief_concern -> policy_support, in
# r_patches/11_ida_appendix_and_interaction_test.R Parts B-D) with no stated
# reason it's the only one. This script generalizes that same method to
# EVERY parent-pair interaction that is even estimable in the current 16-edge
# working SCM -- i.e. every node with 2+ parents, tested for all pairwise
# products among its own parents. There are exactly 10 such candidates:
#
#   node                parent 1         parent 2
#   harm_present         belief_concern   harm_future
#   trust_science        belief_concern   harm_future
#   policy_support        belief_concern   politics          [already published]
#   policy_support        belief_concern   trust_science
#   policy_support        politics         trust_science
#   social_norms          trust_science    policy_support
#   weather_risk_prep     harm_present     belief_concern
#   climate_behavior      harm_present     weather_risk_prep
#   climate_behavior      harm_present     social_norms
#   climate_behavior      weather_risk_prep social_norms
#
# The already-published interaction is included in this same systematic loop
# (not treated as a special case) so it gets reconfirmed under the unified
# pipeline and so the BH correction below is honest -- correcting only the
# OTHER 9 while treating the 1st as already-decided would be exactly the kind
# of selective multiple-comparisons dodge a reviewer would flag.
#
# METHOD (identical to the validated single-interaction test in script 11,
# Parts B/C, just looped over all 10 candidates):
#   1. Mean-center each parent (RAW composite, Aiken & West 1991) before
#      building the product term -- reduces structural collinearity between
#      main effects and the interaction term.
#   2. Build that node's lavaan equation from get_parents_tr(tr_amat) (so it
#      automatically matches the real 16-edge structure) plus ONE extra
#      product term on the target node only. Every OTHER node's equation is
#      exactly the baseline lav_model_tr -- unaffected.
#   3. Fit via lavaan::sem(..., estimator="MLR"), pull the interaction term's
#      standardized coefficient/SE/p-value.
#   4. Cross-check with a plain lm() of just that one node's equation
#      (recursive path model -> single-equation OLS should match lavaan's
#      estimate for that one equation regardless of whole-model fit).
#
# DECISION (2026-09-04): run all 10 (not a smaller subset), and correct
# for multiple comparisons via Benjamini-Hochberg FDR across all 10 raw
# p-values (not Bonferroni, not uncorrected) -- see PART 2 below.
#
# WHAT THIS SCRIPT NEEDS BEFORE THE OUTPUT CAN BE TRUSTED: ext_dat_tr,
# trimmed_nodes, tr_amat, std_tr, scm_coefs, scm_resid_sd all need to already
# exist in the session, exactly as script 11 requires. The BH correction, the
# generic downstream simulator, and the simple-slopes math below have been
# checked for logic against synthetic stand-in data (same variable
# names/shapes, placeholder coefficients), confirming the code runs
# end-to-end with no errors -- but every actual number below (which
# interactions are significant, what they're worth downstream) only comes
# from running this against the real data. Run it, then record the two
# printed tables (PART 2's screen table and PART 4's downstream-propagation
# table).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

stopifnot(exists("ext_dat_tr"), exists("trimmed_nodes"), exists("tr_amat"),
          exists("std_tr"), exists("scm_coefs"), exists("scm_resid_sd"))

# =============================================================================
# PART 1 -- candidate list and the per-candidate test function
# =============================================================================
CANDIDATES <- tibble::tribble(
  ~node,                ~p1,                 ~p2,
  "harm_present",       "belief_concern",    "harm_future",
  "trust_science",      "belief_concern",    "harm_future",
  "policy_support",     "belief_concern",    "politics",
  "policy_support",     "belief_concern",    "trust_science",
  "policy_support",     "politics",          "trust_science",
  "social_norms",       "trust_science",     "policy_support",
  "weather_risk_prep",  "harm_present",      "belief_concern",
  "climate_behavior",   "harm_present",      "weather_risk_prep",
  "climate_behavior",   "harm_present",      "social_norms",
  "climate_behavior",   "weather_risk_prep", "social_norms"
)
stopifnot(nrow(CANDIDATES) == 10)

# Builds the baseline-plus-one-product-term lavaan model string for a given
# target node/parent pair, exactly as script 11 Part B does for the one
# published interaction, generalized to any (node, p1, p2).
build_lav_model_int <- function(target_node, term_name) {
  lines <- c()
  for (nd in trimmed_nodes) {
    pa <- get_parents_tr(nd, tr_amat, trimmed_nodes)
    if (length(pa) == 0) next
    rhs <- paste(pa, collapse = " + ")
    if (nd == target_node) rhs <- paste(rhs, "+", term_name)
    lines <- c(lines, paste0(nd, " ~ ", rhs))
  }
  paste(lines, collapse = "\n")
}

test_one_interaction <- function(node, p1, p2) {
  dat <- ext_dat_tr
  dat$.p1_c <- dat[[p1]] - mean(dat[[p1]], na.rm = TRUE)
  dat$.p2_c <- dat[[p2]] - mean(dat[[p2]], na.rm = TRUE)
  dat$.term <- dat$.p1_c * dat$.p2_c

  lav_model <- build_lav_model_int(node, ".term")
  scaled <- as.data.frame(scale(dat[, c(trimmed_nodes, ".term")]))

  fit <- tryCatch(
    lavaan::sem(lav_model, data = scaled, estimator = "MLR"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble::tibble(node = node, p1 = p1, p2 = p2,
                           est_std = NA_real_, se = NA_real_, pvalue = NA_real_,
                           lm_est = NA_real_, lm_pvalue = NA_real_,
                           fit_error = TRUE))
  }

  std <- lavaan::standardizedSolution(fit)
  row <- std[std$lhs == node & std$rhs == ".term", ]

  # lm() cross-check: refit just the target node's own equation as plain OLS.
  parents <- get_parents_tr(node, tr_amat, trimmed_nodes)
  form <- as.formula(paste(node, "~", paste(c(parents, ".term"), collapse = " + ")))
  lm_fit <- lm(form, data = scaled)
  lm_coef <- summary(lm_fit)$coefficients
  lm_row <- if (".term" %in% rownames(lm_coef)) lm_coef[".term", ] else c(NA, NA, NA, NA)

  tibble::tibble(
    node = node, p1 = p1, p2 = p2,
    est_std = if (nrow(row) == 1) row$est.std[1] else NA_real_,
    se      = if (nrow(row) == 1) row$se[1]      else NA_real_,
    pvalue  = if (nrow(row) == 1) row$pvalue[1]  else NA_real_,
    lm_est    = unname(lm_row[1]),
    lm_pvalue = unname(lm_row[4]),
    fit_error = FALSE
  )
}

# =============================================================================
# PART 2 -- run all 10, then Benjamini-Hochberg correct across all 10 raw
# p-values (2026-09-04 -- not Bonferroni, not uncorrected).
# =============================================================================
screen_results <- purrr::pmap_dfr(CANDIDATES, test_one_interaction)

screen_results <- screen_results |>
  mutate(
    q_value        = p.adjust(pvalue, method = "BH"),
    sig_uncorrected = !is.na(pvalue) & pvalue < .05,
    sig_bh          = !is.na(q_value) & q_value < .05,
    lavaan_lm_agree = !is.na(est_std) & !is.na(lm_est) &
                       sign(est_std) == sign(lm_est) &
                       abs(est_std - lm_est) < .05
  ) |>
  arrange(pvalue)

cat("\n===== Systematic interaction screen: all 10 candidate parent-pair",
    "interactions =====\n")
print(screen_results |>
  select(node, p1, p2, est_std, se, pvalue, q_value, sig_bh, lavaan_lm_agree) |>
  mutate(across(where(is.numeric), ~round(.x, 4))), n = 10)

cat("\n", sum(screen_results$sig_uncorrected, na.rm = TRUE),
    "of 10 significant at uncorrected alpha=.05; ",
    sum(screen_results$sig_bh, na.rm = TRUE),
    "of 10 significant after Benjamini-Hochberg FDR correction.\n", sep = "")

if (any(!screen_results$lavaan_lm_agree, na.rm = TRUE)) {
  cat("\n*** WARNING: at least one lavaan/lm cross-check disagrees by >0.05",
      "or differs in sign -- inspect that row before trusting it. ***\n")
}

saveRDS(screen_results, "data/causal_results/interaction_screen_results.rds")
message("Saved data/causal_results/interaction_screen_results.rds")

# =============================================================================
# PART 3 -- generic downstream simulator (one function, reused for every
# significant interaction, rather than one bespoke simulate_scm_*() per
# pair). Computes all 9 nodes in a fixed topological order using the
# standard 16-edge equations, with two optional generic hooks:
#   fixed = list(node = level)   -- holds a node at a constant level, no
#                                    noise (used to set a "moderator" at a
#                                    representative -1/+1 SD level)
#   shift = list(node = amount)  -- adds a constant on top of that node's
#                                    normal structural value (a shift
#                                    intervention, per script 11 Part D's
#                                    corrected framing -- NOT an absolute
#                                    do(node=value) target)
#   extra = list(node=, coef=, p1=, p2=) -- adds coef * val(p1) * val(p2) to
#                                    exactly one node's equation
# Topological order verified by hand against all 16 edges in
# scm_edges_finalized.csv (each node's parents all appear earlier in the
# list): politics, belief_concern, harm_future, harm_present, trust_science,
# policy_support, weather_risk_prep, social_norms, climate_behavior.
# =============================================================================
simulate_scm_generic <- function(n = 50000, seed = 42,
                                  fixed = list(), shift = list(), extra = NULL) {
  set.seed(seed)
  p <- scm_coefs
  vals <- list()

  node_val <- function(node, structural_expr) {
    if (!is.null(fixed[[node]])) {
      v <- rep(fixed[[node]], n)
    } else {
      v <- structural_expr + rnorm(n, 0, scm_resid_sd[node])
      if (!is.null(shift[[node]])) v <- v + shift[[node]]
    }
    v
  }
  extra_term <- function(node) {
    if (!is.null(extra) && identical(extra$node, node)) {
      extra$coef * vals[[extra$p1]] * vals[[extra$p2]]
    } else 0
  }

  vals$politics       <- node_val("politics", 0) + extra_term("politics")
  vals$belief_concern  <- node_val("belief_concern",
                             p$pol_bc * vals$politics) + extra_term("belief_concern")
  vals$harm_future     <- node_val("harm_future",
                             p$bc_hf * vals$belief_concern) + extra_term("harm_future")
  vals$harm_present    <- node_val("harm_present",
                             p$hf_hp * vals$harm_future +
                             p$bc_hp * vals$belief_concern) + extra_term("harm_present")
  vals$trust_science   <- node_val("trust_science",
                             p$bc_ts * vals$belief_concern +
                             p$hf_ts * vals$harm_future) + extra_term("trust_science")
  vals$policy_support  <- node_val("policy_support",
                             p$bc_ps * vals$belief_concern +
                             p$pol_ps * vals$politics +
                             p$ts_ps * vals$trust_science) + extra_term("policy_support")
  vals$weather_risk_prep <- node_val("weather_risk_prep",
                             p$hp_wr * vals$harm_present +
                             p$bc_wr * vals$belief_concern) + extra_term("weather_risk_prep")
  vals$social_norms    <- node_val("social_norms",
                             p$ts_sn * vals$trust_science +
                             p$ps_sn * vals$policy_support) + extra_term("social_norms")
  vals$climate_behavior <- node_val("climate_behavior",
                             p$hp_cb * vals$harm_present +
                             p$wr_cb * vals$weather_risk_prep +
                             p$sn_cb * vals$social_norms) + extra_term("climate_behavior")

  tibble::as_tibble(vals)
}

# =============================================================================
# PART 4 -- for every interaction significant after BH correction, compute
# BOTH simple-slope orderings (each parent in turn as the "held-at-representative
# -1/+1-SD moderator", the other as the "+0.5 SD shift target"), propagated to
# climate_behavior, using the SAME shift-intervention framing already
# corrected onto in script 11 Part D (not an absolute do() target -- see that
# script's header for why the distinction matters). Reporting both directions
# rather than picking one avoids imposing a "which parent is the real
# moderator" call that's a substantive, not statistical, decision -- that
# call is left for the research team to make once the interactions found to
# be real are in hand.
# =============================================================================
sig_rows <- screen_results |> filter(sig_bh)

if (nrow(sig_rows) == 0) {
  cat("\nNo interactions survived BH correction -- PART 4 has nothing to",
      "propagate downstream. This is a valid, reportable result (the single",
      "already-published interaction may not survive once tested alongside",
      "9 others under FDR control -- report the full PART 2 table either way).\n")
} else {

  propagate_one <- function(node, p1, p2, coef, moderator, target) {
    # "Clean" (single-standardization) coefficient for the simulator, matching
    # script 11 Part C's approach: rebuild the product term from EACH
    # variable's OWN z-score (not the group-scale()'d version used for the
    # significance test), refit once, so the coefficient is in exactly the
    # units simulate_scm_generic() needs. Re-centering happens once here, not
    # inside simulate_scm_generic() -- kept out of PART 3 so that function
    # stays a pure simulator, no data-fitting inside it.
    z1 <- as.numeric(scale(ext_dat_tr[[p1]]))
    z2 <- as.numeric(scale(ext_dat_tr[[p2]]))
    term_clean <- z1 * z2
    scaled_clean <- as.data.frame(scale(ext_dat_tr[, trimmed_nodes]))
    scaled_clean$.term <- term_clean
    lav_model <- build_lav_model_int(node, ".term")
    fit_clean <- lavaan::sem(lav_model, data = scaled_clean, estimator = "MLR")
    std_clean <- lavaan::standardizedSolution(fit_clean)
    coef_clean <- std_clean$est.std[std_clean$lhs == node & std_clean$rhs == ".term"][1]

    extra_spec <- list(node = node, coef = coef_clean, p1 = p1, p2 = p2)

    lo <- simulate_scm_generic(fixed = setNames(list(-1), moderator),
                                shift = setNames(list(0.5), target),
                                extra = extra_spec)
    lo_base <- simulate_scm_generic(fixed = setNames(list(-1), moderator),
                                     extra = extra_spec)
    hi <- simulate_scm_generic(fixed = setNames(list(1), moderator),
                                shift = setNames(list(0.5), target),
                                extra = extra_spec)
    hi_base <- simulate_scm_generic(fixed = setNames(list(1), moderator),
                                     extra = extra_spec)

    tibble::tibble(
      node = node, p1 = p1, p2 = p2,
      moderator = moderator, target = target, coef_clean = coef_clean,
      delta_target_low  = mean(lo[[node]])              - mean(lo_base[[node]]),
      delta_target_high = mean(hi[[node]])              - mean(hi_base[[node]]),
      delta_cb_low  = mean(lo$climate_behavior)  - mean(lo_base$climate_behavior),
      delta_cb_high = mean(hi$climate_behavior)  - mean(hi_base$climate_behavior)
    )
  }

  downstream_results <- purrr::pmap_dfr(sig_rows, function(node, p1, p2, ...) {
    dplyr::bind_rows(
      propagate_one(node, p1, p2, moderator = p2, target = p1),
      propagate_one(node, p1, p2, moderator = p1, target = p2)
    )
  })

  cat("\n===== Downstream propagation of BH-significant interactions",
      "(both parent orderings) =====\n")
  print(downstream_results |>
    mutate(across(where(is.numeric), ~round(.x, 4))), n = nrow(downstream_results))

  saveRDS(downstream_results, "data/causal_results/interaction_screen_downstream.rds")
  message("Saved data/causal_results/interaction_screen_downstream.rds")
}

cat("\n===== DONE: 19_systematic_interaction_screen.R =====\n")
cat("Send back both printed tables (PART 2 screen, PART 4 downstream if any",
    "interactions survived BH correction).\n")
