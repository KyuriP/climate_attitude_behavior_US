# =============================================================================
# 11_ida_appendix_and_interaction_test.R
#
# The working SCM gained 2 edges this round (harm_future -> trust_science,
# belief_concern -> weather_risk_prep -- see 18_finalize_scm_specification_v4.R
# and the qmd's own simulate_scm()). This file has its OWN separate,
# hand-rolled simulation functions -- simulate_scm_int() (Part C) and
# simulate_scm_shift() (Part D, the one actually feeding main.tex's reported
# .080/.402 vs .322/.205 vs .202 numbers) -- which are NOT auto-derived from
# tr_amat and were silently still missing both new edges (same bug class as
# simulate_scm() itself before that got fixed). Patched both functions'
# trust_science and weather_risk_prep equations to include p$hf_ts*hf and
# p$bc_wr*bc respectively, matching the qmd's simulate_scm(). p_int (built as
# p_int <- scm_coefs, then adding bc_pol_ps) already carried hf_ts/bc_wr as
# values -- they just weren't referenced in the equations, so this was a
# silent omission, not a missing-coefficient error. NOT yet re-run against
# live data: Parts B/C's interaction coefficient itself (bc_x_politics/
# bc_pol_ps) should come out fine on a rerun since lav_model_int is built
# from get_parents_tr(tr_amat), which already reflects all 16 edges; it's
# specifically the Part D shift-intervention numbers in main.tex that were
# computed under the stale 14-edge simulation and need a fresh rerun with
# this fix in place.
#
# Two independent, separately-scoped additions:
#
#   PART A -- pcalg::ida()/jointIda() bounds under PC-stable's CPDAG, for the
#             appendix ONLY. This does NOT touch FCI/CCI or the main-text
#             structural-uncertainty story (that stays the orientation
#             enumeration in 02_full_orientation_enumeration.R). ida()/jointIda()
#             only accept a CPDAG/PDAG (type = "cpdag"/"pdag") -- they assume
#             causal sufficiency, the same assumption PC-stable makes and that
#             we've already decided is questionable enough here to demote PC
#             to the appendix. So this block answers "if we trusted PC-stable's
#             assumptions, how far would the equivalence class move the
#             intervention numbers" -- an appendix note *about PC*, not a
#             replacement for anything in the main text. There is no PAG
#             equivalent of ida() available as a package (Malinsky & Spirtes'
#             own implementation, github.com/dmalinsk/lv-ida, is research code,
#             not something to bolt in here) -- that's why this stays PC-only.
#
#   PART B -- adds a politics x belief_concern interaction to the
#             policy_support equation of the working SCM (politics,
#             belief_concern, and trust_science are its three parents per the
#             manuscript's 14-edge model). This is unrelated to the
#             causal-discovery-algorithm question -- it's a downstream
#             refinement of the already-fitted lavaan model -- so it can run
#             regardless of the FCI/CCI/PC framing above.
#
# Run climate_analysis_avg_v2_altweather.qmd first, through Section 7.5/7.6
# (the working SCM + lavaan fit). This script assumes these objects exist:
#   node_order_ext, pc_ext_amat_05, ext_dat_tr, trimmed_nodes, tr_amat,
#   lav_model_tr, sem_fit_tr, std_tr, scm_coefs, scm_resid_sd, simulate_scm()
#
# worth checking -- Part A: console output, especially str(eff)/dim() of what ida()/jointIda()
# actually return in the installed pcalg version -- printed defensively
# (structure before summarizing) since pcalg's exact return shape can differ
# across versions and hasn't been confirmed against a live run here.
# Part B: whether the interaction term's p-value/coefficient printed, and the
# new conditional intervention numbers.
# =============================================================================

suppressPackageStartupMessages({
  library(pcalg)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(lavaan)
})

# =============================================================================
# PART A -- ida() / jointIda() bounds (PC-stable CPDAG, appendix only)
# =============================================================================
# DECISION (documented in the paper-idea assessment notes): this appendix
# computation is DROPPED from the manuscript.
# pc_ext_amat_05 contains a genuine directed cycle (harm_future -> belief_
# concern -> harm_present -> harm_future), so it fails pcalg::isValidGraph()
# and ida()/jointIda() cannot legitimately run on it -- confirmed below, kept
# for the record, but skipped by default so a full-pipeline source() of this
# file doesn't stop on the resulting error. Flip to `if (TRUE)` only if you
# specifically want to re-examine this failure.
if (FALSE) {
stopifnot(exists("pc_ext_amat_05"), exists("node_order_ext"), exists("ext_dat_tr"))

# Rebuild the graph object from pc_ext_amat_05 -- the background-knowledge-
# corrected matrix (climate_behavior forced downstream), the same matrix
# plotPC() already renders -- rather than pc_ext_05@graph directly, which
# predates that correction.
# ---------------------------------------------------------------------------
# DIAGNOSTIC (added after 2nd failed attempt -- type="pdag" didn't fix it,
# same numbers/warning came back, so the problem isn't cpdag-vs-pdag at all.
# Printing the actual matrix + pcalg's own checks instead of guessing again.
# ---------------------------------------------------------------------------
cat("=== DIAGNOSTIC: pc_ext_amat_05 (send this whole block back) ===
")
print(pc_ext_amat_05)

cat("
Direct isValidGraph() calls:
")
cat("  as pdag: ");  print(tryCatch(pcalg::isValidGraph(pc_ext_amat_05, type = "pdag"),
                                     error = function(e) paste("ERROR:", conditionMessage(e))))
cat("  as cpdag:"); print(tryCatch(pcalg::isValidGraph(pc_ext_amat_05, type = "cpdag"),
                                     error = function(e) paste("ERROR:", conditionMessage(e))))

cat("
Manual cycle check on the DIRECTED sub-edges only",
    "(amat[i,j]==1 & amat[j,i]==0 means i -> j):
")
p_ext_local <- nrow(pc_ext_amat_05)
dir_edges <- which(pc_ext_amat_05 == 1 & t(pc_ext_amat_05) == 0, arr.ind = TRUE)
if (nrow(dir_edges) > 0) {
  dir_df <- data.frame(from = node_order_ext[dir_edges[,1]],
                        to   = node_order_ext[dir_edges[,2]])
  print(dir_df)
  if (requireNamespace("igraph", quietly = TRUE)) {
    g_dir <- igraph::graph_from_data_frame(dir_df, vertices = node_order_ext)
    cat("Directed sub-graph acyclic? ", igraph::is_dag(g_dir), "
")
  } else {
    cat("(igraph not installed -- install it for an automatic cycle check,",
        "or eyeball the from/to list above for a cycle.)
")
  }
} else {
  cat("No purely-directed edges found -- check the matrix printed above by eye.
")
}
cat("=== END DIAGNOSTIC ===

")

pc_ext_graph <- as(pc_ext_amat_05, "graphNEL")

# Covariance on the SAME scale simulate_scm()'s SD-based effects are reported
# on: the raw (non-NPN) composite scores the SCM itself is fit on, not the
# NPN-transformed causal-discovery data.
mcov_ext <- cov(ext_dat_tr)

target_y <- which(node_order_ext == "climate_behavior")
stopifnot(length(target_y) == 1)

single_targets <- c("harm_present", "belief_concern", "weather_risk_prep",
                     "social_norms", "trust_science", "policy_support")

cat("=== Part A.1: single-node IDA bounds (PC-stable CPDAG) ===\n")
ida_raw <- purrr::map(single_targets, function(nm) {
  x_pos <- which(node_order_ext == nm)
  eff <- pcalg::ida(x.pos = x_pos, y.pos = target_y, mcov = mcov_ext,
                     graphEst = pc_ext_graph, method = "local", type = "pdag")
  # type="pdag" (not "cpdag"): pc_ext_amat_05 already has background knowledge
  # (climate_behavior forced downstream) applied via addBgKnowledge(), so it's
  # a background-knowledge-augmented PDAG, not guaranteed to pass the stricter
  # CPDAG validity check pcalg::isValidGraph() applies. First run flagged
  # "input graph is not a valid cpdag" on every call and returned a single
  # point estimate per node instead of a real multiset -- this is the fix.
  cat(sprintf("\n--- %s -> climate_behavior ---\n", nm))
  cat("str(eff):\n"); print(str(eff))
  eff
})
names(ida_raw) <- single_targets

ida_bounds <- purrr::imap_dfr(ida_raw, function(eff, nm) {
  tibble::tibble(node = nm, min_effect = min(eff), max_effect = max(eff),
                 n_estimates = length(eff), point_estimate_dag_count = length(unique(eff)))
})

cat("\n=== Part A.1 summary ===\n")
print(ida_bounds)

cat("\n=== Part A.2: joint IDA bounds (2 combined-intervention scenarios) ===\n")
# NOTE: jointIda() returns a matrix with one ROW per intervention variable in
# x.pos and one COLUMN per valid joint-effect vector (see ?jointIda). For a
# genuinely simultaneous +0.5 SD intervention on both variables, the combined
# model-implied effect per valid specification is the column sum (assumes
# additive/linear effects, consistent with simulate_scm()'s own linear SCM --
# if that assumption looks wrong once you see the printed structure below,
# stop and review the output before proceeding rather than trusting the summary).
joint_specs <- list(
  harm_present_plus_norms = c("harm_present", "social_norms"),
  belief_plus_harm        = c("belief_concern", "harm_present")
)

joint_raw <- purrr::map(joint_specs, function(nms) {
  x_pos <- match(nms, node_order_ext)
  eff <- pcalg::jointIda(x.pos = x_pos, y.pos = target_y, mcov = mcov_ext,
                          graphEst = pc_ext_graph, technique = "RRC",
                          type = "pdag")  # see note above single-node ida() call
  cat("\nstr(eff):\n"); print(str(eff)); cat("dim(eff):\n"); print(dim(eff))
  eff
})

joint_bounds <- purrr::imap_dfr(joint_raw, function(eff, label) {
  col_sums <- colSums(eff)
  tibble::tibble(scenario = label, min_effect = min(col_sums),
                 max_effect = max(col_sums), n_estimates = length(col_sums))
})

cat("\n=== Part A.2 summary ===\n")
print(joint_bounds)

dir.create("data/causal_results", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(ida_bounds = ida_bounds, joint_bounds = joint_bounds,
             ida_raw = ida_raw, joint_raw = joint_raw),
        "data/causal_results/ida_bounds_appendix.rds")
message("Saved data/causal_results/ida_bounds_appendix.rds")

cat("\nCompare against the point estimates already in the intervention table:\n")
cat("(these ranges are wider than / different from the working-SCM point\n")
cat(" estimate by construction -- they span the WHOLE PC equivalence class,\n")
cat(" not just the 5 theory-completed edges the orientation enumeration flips.)\n")
}  # end if (FALSE) -- Part A skipped, see decision note above

# =============================================================================
# PART B -- politics x belief_concern interaction on policy_support
# =============================================================================
stopifnot(exists("ext_dat_tr"), exists("trimmed_nodes"), exists("tr_amat"),
          exists("std_tr"))

cat("\n\n=== Part B: testing politics x belief_concern -> policy_support ===\n")

# Mean-center (not full z-score) before the product term, per Aiken & West
# (1991) -- reduces structural collinearity between the main effects and the
# interaction term. ext_dat_tr is the RAW (non-NPN) composite data the
# working SCM is fit on (see Section 7.5/7.6); everything gets scale()'d
# together right before the lavaan fit, exactly as the baseline model does.
ext_dat_int <- ext_dat_tr
ext_dat_int$bc_c  <- ext_dat_int$belief_concern - mean(ext_dat_int$belief_concern, na.rm = TRUE)
ext_dat_int$pol_c <- ext_dat_int$politics       - mean(ext_dat_int$politics,       na.rm = TRUE)
ext_dat_int$bc_x_politics <- ext_dat_int$bc_c * ext_dat_int$pol_c

# Same regression-equation structure as the baseline lav_model_tr
# (get_parents_tr-generated), but with one extra term on policy_support only.
lav_lines_int <- c()
for (nd in trimmed_nodes) {
  pa <- get_parents_tr(nd, tr_amat, trimmed_nodes)
  if (length(pa) == 0) next
  rhs <- paste(pa, collapse = " + ")
  if (nd == "policy_support") rhs <- paste(rhs, "+ bc_x_politics")
  lav_lines_int <- c(lav_lines_int, paste0(nd, " ~ ", rhs))
}
lav_model_int <- paste(lav_lines_int, collapse = "\n")
cat("Interaction-augmented lavaan model:\n", lav_model_int, "\n\n")

sem_fit_int <- lavaan::sem(
  lav_model_int,
  data      = as.data.frame(scale(ext_dat_int[, c(trimmed_nodes, "bc_x_politics")])),
  estimator = "MLR"
)

cat("Fit indices (compare to baseline sem_fit_tr):\n")
print(round(lavaan::fitMeasures(sem_fit_int,
      c("cfi","tli","rmsea","srmr","aic","bic")), 3))

std_int <- lavaan::standardizedSolution(sem_fit_int)
int_row <- std_int[std_int$lhs == "policy_support" & std_int$rhs == "bc_x_politics", ]
cat("\npolitics x belief_concern -> policy_support:\n")
print(int_row[, c("lhs","rhs","est.std","se","pvalue")])

if (nrow(int_row) == 1 && !is.na(int_row$pvalue) && int_row$pvalue < .05) {
  cat("\n*** Interaction is significant at alpha=.05 -- see note below on\n")
  cat("    extending simulate_scm() with this term before re-running the\n")
  cat("    policy_support intervention scenario split by politics level. ***\n")
} else {
  cat("\nInteraction not significant at alpha=.05 -- baseline (no-interaction)\n")
  cat("model for policy_support is not contradicted; you may not need to\n")
  cat("touch simulate_scm() at all. Review this output either way -- it will\n")
  cat("confirm the right next step.\n")
}

# ---------------------------------------------------------------------------
# Cross-check: refit just the policy_support equation as a plain lm(), outside
# lavaan's whole-model-fit framework. In a fully recursive path model (no
# correlated residuals, no feedback loops -- true here), the single-equation
# OLS/lm() estimates should match lavaan's for that one equation regardless of
# what the overall CFI/RMSEA/SRMR say about the full model. This isolates the
# interaction test from the product-term/whole-model-fit artifact noted above.
# ---------------------------------------------------------------------------
scaled_dat_int <- as.data.frame(scale(ext_dat_int[, c(trimmed_nodes, "bc_x_politics")]))
lm_check <- lm(policy_support ~ belief_concern + trust_science + politics + bc_x_politics,
               data = scaled_dat_int)
cat("\n=== Part B cross-check: plain lm() for the policy_support equation ===\n")
print(summary(lm_check)$coefficients)
cat("\nCompare the bc_x_politics row above to lavaan's est.std/se/pvalue",
    "printed earlier -- if they agree, the interaction is not an artifact of",
    "how lavaan handles the product term's own covariance structure in the",
    "whole model, and we're clear to extend simulate_scm() with this term.\n")

# =============================================================================
# PART C -- clean interaction coefficient + simulate_scm() extension +
#           politics-conditional intervention on belief_concern
# =============================================================================
# Part B built bc_x_politics from mean-centered RAW composites, then
# re-standardized it a SECOND time together with everything else in the group
# scale() call. That's fine for a significance test (confirmed twice above:
# lavaan est.std=.08, lm's own beta=.0815, both p<.0001) but awkward to port
# into simulate_scm(), which already works entirely in standardized units.
# Here the product term is built from the two components' OWN z-scores
# directly and is NOT re-scaled again -- so its coefficient is exactly what
# simulate_scm() needs to multiply bc*pol by, no unit-conversion algebra.

z_bc  <- as.numeric(scale(ext_dat_tr$belief_concern))
z_pol <- as.numeric(scale(ext_dat_tr$politics))
bc_x_politics_clean <- z_bc * z_pol

scaled_main_clean <- as.data.frame(scale(ext_dat_tr[, trimmed_nodes]))
scaled_main_clean$bc_x_politics <- bc_x_politics_clean

sem_fit_clean <- lavaan::sem(lav_model_int, data = scaled_main_clean, estimator = "MLR")
std_clean <- lavaan::standardizedSolution(sem_fit_clean)
int_row_clean <- std_clean[std_clean$lhs == "policy_support" & std_clean$rhs == "bc_x_politics", ]

cat("\n=== Part C: clean (single-standardization) interaction coefficient ===\n")
print(int_row_clean[, c("lhs","rhs","est.std","se","pvalue")])
cat("(Should be close to Part B's 0.08/0.022 -- if wildly different, STOP and\n")
cat(" review this output before proceeding before trusting anything below.)\n")

lm_check_clean <- lm(policy_support ~ belief_concern + trust_science + politics + bc_x_politics,
                      data = data.frame(scale(ext_dat_tr[, trimmed_nodes]),
                                        bc_x_politics = bc_x_politics_clean))
cat("\nlm() cross-check (clean version):\n")
print(summary(lm_check_clean)$coefficients)

# ---------------------------------------------------------------------------
# Extend simulate_scm() with this one term. New function so the original
# simulate_scm() (used everywhere else in the paper) is untouched.
# ---------------------------------------------------------------------------
p_int <- scm_coefs
p_int$bc_pol_ps <- int_row_clean$est.std[1]
cat("\nInteraction coefficient going into the simulator (bc_pol_ps):",
    round(p_int$bc_pol_ps, 4), "\n")

simulate_scm_int <- function(n = 10000, intervene = list(), seed = 42) {
  set.seed(seed)
  p <- p_int
  noise <- function(node, n) {
    if (node %in% names(intervene)) rep(0, n) else rnorm(n, 0, scm_resid_sd[node])
  }
  pin <- function(node, expr) {
    if (node %in% names(intervene)) rep(intervene[[node]], n) else expr
  }
  pol  <- pin("politics", rnorm(n, 0, scm_resid_sd["politics"]))
  bc   <- pin("belief_concern", p$pol_bc * pol + noise("belief_concern", n))
  hf   <- pin("harm_future",    p$bc_hf  * bc  + noise("harm_future", n))
  ts   <- pin("trust_science",  p$bc_ts  * bc  + p$hf_ts * hf + noise("trust_science", n))
  hp   <- pin("harm_present",   p$hf_hp * hf + p$bc_hp * bc + noise("harm_present", n))
  # ---- only line that differs from the baseline simulate_scm() ----
  ps   <- pin("policy_support",
              p$ts_ps * ts + p$bc_ps * bc + p$pol_ps * pol +
                p$bc_pol_ps * bc * pol + noise("policy_support", n))
  wr   <- pin("weather_risk_prep", p$hp_wr * hp + p$bc_wr * bc + noise("weather_risk_prep", n))
  sn   <- pin("social_norms", p$ts_sn * ts + p$ps_sn * ps + noise("social_norms", n))
  cb   <- pin("climate_behavior",
              p$hp_cb * hp + p$wr_cb * wr + p$sn_cb * sn + noise("climate_behavior", n))
  tibble::tibble(politics = pol, belief_concern = bc, harm_future = hf,
                 trust_science = ts, harm_present = hp, policy_support = ps,
                 weather_risk_prep = wr, social_norms = sn, climate_behavior = cb)
}

obs_baseline_int <- simulate_scm_int(n = 50000)
cat("\nSimulated marginal SDs WITH interaction term (compare to ~1.0 -- a big\n")
cat("departure would mean the product term is destabilizing implied variances):\n")
print(round(apply(obs_baseline_int, 2, sd), 3))

# ---------------------------------------------------------------------------
# Politics-conditional scenario. Important framing note: intervening
# (pinning) policy_support directly makes the interaction moot -- a
# do(policy_support = v) replaces its whole equation, interaction term
# included. The substantively meaningful version of "does politics moderate
# this path" is: does intervening on belief_concern (upstream of
# policy_support) move policy_support -- and downstream climate_behavior --
# differently depending on where someone sits on politics? That's a
# moderated-mediation-style question (belief_concern -> policy_support,
# moderated by politics -> social_norms -> climate_behavior), and it's what's
# tested below: fix politics at two representative levels (-1 SD / +1 SD,
# i.e., simulating two subgroups) and compare do(belief_concern=+0.5) against
# each subgroup's own baseline.
# ---------------------------------------------------------------------------
N_SIM_INT <- 50000
cond_scenarios <- list(
  low_politics_baseline   = list(politics = -1),
  low_politics_intervene  = list(politics = -1, belief_concern = 0.5),
  high_politics_baseline  = list(politics =  1),
  high_politics_intervene = list(politics =  1, belief_concern = 0.5)
)

cond_results <- purrr::imap_dfr(cond_scenarios, function(intv, label) {
  dat <- simulate_scm_int(n = N_SIM_INT, intervene = intv, seed = 42)
  tibble::tibble(scenario = label,
                 policy_support_mean   = mean(dat$policy_support),
                 climate_behavior_mean = mean(dat$climate_behavior))
})

cat("\n=== Politics-conditional do(belief_concern = +0.5 SD) ===\n")
print(cond_results)

delta_low  <- cond_results$policy_support_mean[cond_results$scenario == "low_politics_intervene"] -
              cond_results$policy_support_mean[cond_results$scenario == "low_politics_baseline"]
delta_high <- cond_results$policy_support_mean[cond_results$scenario == "high_politics_intervene"] -
              cond_results$policy_support_mean[cond_results$scenario == "high_politics_baseline"]
delta_cb_low  <- cond_results$climate_behavior_mean[cond_results$scenario == "low_politics_intervene"] -
                 cond_results$climate_behavior_mean[cond_results$scenario == "low_politics_baseline"]
delta_cb_high <- cond_results$climate_behavior_mean[cond_results$scenario == "high_politics_intervene"] -
                 cond_results$climate_behavior_mean[cond_results$scenario == "high_politics_baseline"]

cat(sprintf("\ndo(belief_concern=+0.5) effect on policy_support: %.3f SD at low politics (-1SD), %.3f SD at high politics (+1SD)\n",
            delta_low, delta_high))
cat(sprintf("Same, propagated to climate_behavior: %.3f SD (low) vs %.3f SD (high)\n",
            delta_cb_low, delta_cb_high))
cat("If these two numbers are noticeably different, that's the moderated\n")
cat("intervention story for the paper -- the belief-to-policy-support pathway\n")
cat("(and its downstream reach to behavior) is politics-dependent, not fixed.\n")

saveRDS(list(p_int = p_int, cond_results = cond_results), "data/causal_results/interaction_simulation.rds")
message("Saved data/causal_results/interaction_simulation.rds")

cat("\n===== DONE =====\n")
cat("New objects: pc_ext_graph, mcov_ext, ida_raw, ida_bounds, joint_raw,\n")
cat("joint_bounds, ext_dat_int, sem_fit_int, std_int, int_row.\n")
cat("Record the full console output above (esp. the str()/dim() prints\n")
cat("in Part A, and the fit/interaction printout in Part B) before drafting\n")
cat("the next step (main.tex text for Part A, or the simulate_scm()\n")
cat("extension + conditional-effect re-simulation for Part B).\n")


# =============================================================================
# PART D -- CORRECTED framing: shift intervention, not a common absolute target
# =============================================================================
# DECISION (documented in the paper-idea assessment notes): the
# cond_results simulation directly above (Part B/C's conditional scenario)
# answers a different question than the one the paper wants to make. It
# sets belief_concern to the SAME
# ABSOLUTE value (+0.5 SD) for both political subgroups via do(belief_concern
# = 0.5). Because politics -> belief_concern is a real, fitted path
# (p_int$pol_bc), the two subgroups start from very different baseline
# belief_concern levels, so "set everyone to +0.5" is a much bigger jump for
# conservatives than for liberals. The resulting asymmetry (+.660 vs -.027 in
# cond_results above) is therefore a starting-point/target-level heterogeneity
# result, NOT evidence that the causal effect of belief_concern on downstream
# outcomes is stronger for conservatives -- and it is the opposite of what the
# positive interaction coefficient (bc_pol_ps > 0) implies for the marginal
# slope d(policy_support)/d(belief_concern) = bc_ps + bc_pol_ps * politics,
# which is LARGER (more positive) at higher (more liberal) politics.
#
# This block implements the corrected comparison: a SHIFT intervention,
# belief_concern <- belief_concern + 0.5 applied on top of each unit's own
# structural value (preserving that unit's idiosyncratic deviation and the
# politics -> belief_concern path), evaluated at the same two representative
# politics levels (-1 SD / +1 SD) used above. This answers the actual
# moderation question the paper wants: does the SAME increase in belief/
# concern propagate differently downstream depending on political
# orientation? cond_results above is NOT deleted or overwritten -- it may
# still be worth reporting, correctly re-framed (see cat() note below), as a
# secondary "same absolute target" comparison, but shift_results is the
# corrected primary comparison to write into the Results section.
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

cat("\n=== CORRECTED: belief_concern <- belief_concern + 0.5 (shift, not absolute target) ===\n")
print(shift_results)

shift_delta_ps_low  <- shift_results$policy_support_mean[shift_results$scenario == "low_politics_shift"] -
                       shift_results$policy_support_mean[shift_results$scenario == "low_politics_baseline"]
shift_delta_ps_high <- shift_results$policy_support_mean[shift_results$scenario == "high_politics_shift"] -
                       shift_results$policy_support_mean[shift_results$scenario == "high_politics_baseline"]
shift_delta_cb_low  <- shift_results$climate_behavior_mean[shift_results$scenario == "low_politics_shift"] -
                       shift_results$climate_behavior_mean[shift_results$scenario == "low_politics_baseline"]
shift_delta_cb_high <- shift_results$climate_behavior_mean[shift_results$scenario == "high_politics_shift"] -
                       shift_results$climate_behavior_mean[shift_results$scenario == "high_politics_baseline"]

cat(sprintf("\nSHIFT belief_concern by +0.5 SD, effect on policy_support: %.3f SD at low politics (-1SD), %.3f SD at high politics (+1SD)\n",
            shift_delta_ps_low, shift_delta_ps_high))
cat(sprintf("Same, propagated to climate_behavior: %.3f SD (low) vs %.3f SD (high)\n",
            shift_delta_cb_low, shift_delta_cb_high))
cat("Expected direction, given bc_pol_ps > 0: the direct belief->policy_support\n")
cat("slope (bc_ps + bc_pol_ps*politics) should be LARGER at high (liberal)\n")
cat("politics. Downstream climate_behavior can differ from this because of\n")
cat("the rest of the SCM (e.g., weather_risk_prep and social_norms pathways).\n")

# Closed-form marginal-slope sanity check, independent of simulation noise:
slope_low  <- p_int$bc_ps + p_int$bc_pol_ps * (-1)
slope_high <- p_int$bc_ps + p_int$bc_pol_ps * ( 1)
cat(sprintf("\nClosed-form marginal slope d(policy_support)/d(belief_concern) = bc_ps + bc_pol_ps*politics:\n"))
cat(sprintf("  at politics = -1 (conservative): %.4f\n", slope_low))
cat(sprintf("  at politics = +1 (liberal):      %.4f\n", slope_high))
cat("If shift_delta_ps_low/high above don't track this closed-form ordering,\n")
cat("STOP and review the output before proceeding -- something would be inconsistent.\n")

saveRDS(list(p_int = p_int, shift_results = shift_results,
             slope_low = slope_low, slope_high = slope_high),
        "data/causal_results/interaction_shift_simulation.rds")
message("Saved data/causal_results/interaction_shift_simulation.rds")

cat("\n===== PART D DONE =====\n")
cat("New objects: simulate_scm_shift, shift_scenarios, shift_results,\n")
cat("slope_low, slope_high.\n")
cat("Send back: the printed shift_results table, the two SHIFT-effect lines,\n")
cat("and the two closed-form slope lines above.\n")
