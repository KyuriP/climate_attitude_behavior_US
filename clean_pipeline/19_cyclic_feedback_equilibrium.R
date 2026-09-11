# 19_cyclic_feedback_equilibrium.R -- feedback/equilibrium extension for the
# 4 structurally-cyclic orientation combinations (combo_1/5/9/13) that
# 07_intervention_ates_8node.R and 08/09/10 skip. New 2026-09-11, per Kyuri's
# plan (analysis_decisions_log.md Section 39) -- proposed as Supplementary
# Section S15 ("Feedback extension for cyclic orientation specifications",
# \label{supp:feedback}), NOT yet added to run_all.R's scripts_in_order
# (optional/separate, same convention as 15/16/17/18).
#
# ---- The single question this script exists to answer ---------------------
# If the four cyclic orientation combinations are interpreted as equilibrium
# systems (x = Bx + eps) rather than discarded, do they materially change the
# intervention conclusions relative to the 12 recursive (acyclic) combos?
#
# All four share the same core 3-node loop (same mechanism documented in
# 07's header): flipping belief_concern -> weather_risk_prep alone (leaving
# harm_present -> weather_risk_prep and the fixed belief_concern ->
# harm_present edge untouched) always closes
#   weather_risk_prep -> belief_concern -> harm_present -> weather_risk_prep
# regardless of the other two flip candidates (social_norms<->climate_
# behavior, harm_future<->harm_present), which is why exactly 4 of the 16
# combinations are cyclic and why all 4 are "the same loop" in different
# surrounding contexts.
#
# ---- Two-tier estimation strategy per combo --------------------------------
# PRIMARY: fit the whole 16-edge system for that combo as ONE simultaneous
# lavaan model (build_lavaan_syntax() from 05 already just emits "to ~ from1
# + from2 + ..." per node with no acyclicity assumption baked in, so handing
# it a cyclic edge list produces exactly the reciprocal-paths non-recursive
# SEM syntax this needs -- no new syntax builder required). Checked for:
# convergence, finite/identified standard errors, no negative residual
# variances (Heywood cases), det(I-B) != 0, condition number of (I-B), and
# spectral radius rho(B). If rho(B) < 1, solve the full-system equilibrium
# for do(node=0.5) by removing that node's own equation and solving the
# reduced linear system, general and exact for whatever that combo's actual
# climate_behavior parent set is (handles combo_5/13 losing social_norms as
# a climate_behavior parent automatically, no special-casing needed).
#
# FALLBACK (only if the joint fit fails to converge/identify, or rho(B)>=1):
# per Kyuri's plan item 8 -- do NOT force identification. Instead fit every
# node's regression equation SEPARATELY via lm() on the combo's actual
# parent structure (always well-defined for observed, non-latent variables,
# and correct for every non-loop edge and for the two non-swept loop edges
# belief_concern->harm_present-or-harm_future and harm_present->
# weather_risk_prep, which don't depend on the contested feedback
# coefficient). The one coefficient NOT separately identified this way --
# weather_risk_prep -> belief_concern, the flipped edge that closes the loop
# -- is swept across a plausible range instead of estimated, and the
# equilibrium is solved at each swept value (a feedback-strength sensitivity
# analysis rather than a point estimate). This is explicitly an
# approximation (single-equation OLS on a jointly-endogenous loop is not the
# same as a properly re-identified simultaneous-equations estimate -- no
# instrument for the loop is currently available, same caveat as this
# session's earlier by-hand version of this calculation) and is flagged as
# such in the output.

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
n_all <- length(all_nodes)

# ---- Identify the 4 cyclic combos and each one's "nearest acyclic" analog --
# Nearest acyclic = same flip set with edge #1 (belief_concern ->
# weather_risk_prep, the edge whose flip is what closes the loop -- see
# 07's header) un-flipped, everything else identical. This isolates "what
# changes when only the loop-closing flip is added," which is the
# comparison the amplification ratio below is meant to capture.
is_combo_acyclic <- vapply(scenario_list, function(idx) is_acyclic(flip_edges(base_edges, idx)), logical(1))
cyclic_names  <- names(scenario_list)[!is_combo_acyclic]
stopifnot(length(cyclic_names) == 4)

nearest_acyclic_for <- function(nm) {
  this_idx <- scenario_list[[nm]]
  target_idx <- sort(setdiff(this_idx, 1L))
  hit <- vapply(scenario_list, function(s) identical(sort(s), target_idx), logical(1))
  stopifnot(sum(hit) == 1)
  names(scenario_list)[hit]
}
nearest_acyclic_map <- setNames(vapply(cyclic_names, nearest_acyclic_for, character(1)), cyclic_names)

cat("Cyclic combos and their nearest-acyclic analogs (differ only in the\n",
    "belief_concern<->weather_risk_prep flip, everything else held fixed):\n", sep = "")
print(data.frame(cyclic = names(nearest_acyclic_map), nearest_acyclic = unname(nearest_acyclic_map)),
      row.names = FALSE)

# ---- Read 07's already-computed 12-acyclic-combo ATEs for the comparison ---
acyclic_ate_path <- file.path(TABLES_DIR, "orientation_enumeration_ate_deterministic_8node.csv")
if (!file.exists(acyclic_ate_path)) {
  stop("can't find ", acyclic_ate_path, " -- run 07_intervention_ates_8node.R first, ",
       "this script needs it for the nearest-acyclic comparison.")
}
acyclic_ate <- utils::read.csv(acyclic_ate_path, stringsAsFactors = FALSE)

get_nearest_acyclic_ate <- function(nearest_combo, node) {
  row <- acyclic_ate[acyclic_ate$scenario == nearest_combo & acyclic_ate$node == node, ]
  if (nrow(row) != 1) return(NA_real_)
  row$ate_climate_behavior[1]
}

# ---- Build a full B matrix (B[to, from] = standardized beta) from a set of
# edges + a lookup of coefficients (named "from->to" -> beta). Any edge in
# `edges` missing from `coef_lookup` is left at 0 (shouldn't happen if
# coef_lookup is built correctly for that same edge list). -----------------
build_B <- function(edges, coef_lookup) {
  B <- matrix(0, n_all, n_all, dimnames = list(all_nodes, all_nodes))
  for (k in seq_len(nrow(edges))) {
    key <- paste(edges$from[k], edges$to[k], sep = "->")
    b <- coef_lookup[[key]]
    if (is.null(b) || is.na(b)) {
      warning("no coefficient found for edge ", key, " -- left at 0 in B, check coef_lookup.")
      b <- 0
    }
    B[edges$to[k], edges$from[k]] <- b
  }
  B
}

# Equilibrium single-node ATE from a full B matrix: remove node j's own
# equation, fix x_j = 0.5, solve the reduced linear system for everything
# else, read off climate_behavior. Baseline (no intervention) mean is 0 for
# every node by construction (standardized, mean-zero exogenous residuals),
# so the equilibrium value IS the ATE, no separate baseline subtraction
# needed -- same convention scm_mean_propagate() uses elsewhere in 05.
equilibrium_ate <- function(B, node) {
  j <- which(all_nodes == node)
  other <- setdiff(seq_len(n_all), j)
  B_oo <- B[other, other, drop = FALSE]
  B_oj <- B[other, j]
  IminusB_oo <- diag(length(other)) - B_oo
  x_other <- tryCatch(solve(IminusB_oo, B_oj * 0.5), error = function(e) NULL)
  if (is.null(x_other)) return(NA_real_)
  names(x_other) <- all_nodes[other]
  if (node == "climate_behavior") return(0.5)  # do(climate_behavior) is not a target anywhere else either
  unname(x_other["climate_behavior"])
}

spectral_radius <- function(B) max(Mod(eigen(B, only.values = TRUE)$values))

diagnostics <- list()
ate_rows <- list()
sweep_rows <- list()

for (nm in cyclic_names) {

  flip_idx <- scenario_list[[nm]]
  e <- flip_edges(base_edges, flip_idx)
  nearest_nm <- nearest_acyclic_map[[nm]]

  cat("\n=== ", nm, " (", scenario_flip_labels[[nm]], ") ===\n", sep = "")

  # ---- PRIMARY: joint non-recursive lavaan fit -----------------------------
  model_syntax <- build_lavaan_syntax(e)
  fit <- tryCatch(
    lavaan::sem(model_syntax, data = df_extended, estimator = "MLR", fixed.x = FALSE),
    error = function(err) err
  )

  fit_error   <- inherits(fit, "error") || inherits(fit, "simpleError")
  converged   <- !fit_error && isTRUE(tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE))
  se_finite   <- FALSE
  no_heywood  <- FALSE
  std_sol     <- NULL

  if (converged) {
    std_sol <- lavaan::standardizedSolution(fit) |> dplyr::filter(op == "~")
    se_finite <- all(is.finite(std_sol$se)) && all(std_sol$se < 5)  # >=5 SD SEs treated as practically unidentified
    # Heywood check done directly off parameterEstimates' own variance rows
    # (lhs==rhs, op=="~~") rather than a specific internal matrix slot, since
    # which slot (theta vs psi) holds observed-variable residual variances
    # depends on whether lavaan treats them as indicator-less "y" variables
    # or not -- parameterEstimates() is unambiguous regardless of that.
    pe <- lavaan::parameterEstimates(fit)
    resid_var <- pe[pe$op == "~~" & pe$lhs == pe$rhs, ]
    no_heywood <- nrow(resid_var) > 0 && all(resid_var$est > 0)
  }

  identified <- converged && se_finite && no_heywood

  method_used <- NA_character_
  rho_B <- NA_real_; detIB <- NA_real_; condIB <- NA_real_; stable <- NA

  if (identified) {
    coef_lookup <- setNames(std_sol$est.std, paste(std_sol$rhs, std_sol$lhs, sep = "->"))
    B <- build_B(e, coef_lookup)
    detIB <- det(diag(n_all) - B)
    condIB <- tryCatch(kappa(diag(n_all) - B, exact = TRUE), error = function(e) NA_real_)
    rho_B <- spectral_radius(B)
    stable <- isTRUE(rho_B < 1)

    if (stable) {
      method_used <- "joint_lavaan"
      cat("  joint lavaan fit converged and identified (finite SEs, no Heywood cases).\n",
          "  det(I-B) =", signif(detIB, 4), " | cond(I-B) =", signif(condIB, 4),
          " | rho(B) =", signif(rho_B, 4), "-> stable equilibrium.\n")
      for (node in node_targets) {
        dy <- equilibrium_ate(B, node)
        near_ate <- get_nearest_acyclic_ate(nearest_nm, node)
        ate_rows[[length(ate_rows) + 1]] <- tibble::tibble(
          combo = nm, node = node, method = method_used,
          ate_equilibrium = dy, nearest_acyclic_combo = nearest_nm,
          nearest_acyclic_ate = near_ate,
          amplification_ratio = ifelse(abs(near_ate) > 1e-6, dy / near_ate, NA_real_),
          amplification_diff = dy - near_ate
        )
      }
    } else {
      cat("  joint lavaan fit converged and identified, but rho(B) =", signif(rho_B, 4),
          ">= 1 -- NOT a stable equilibrium under this fit. Falling back to the",
          "feedback-strength sweep instead of reporting an unstable point estimate.\n")
    }
  } else {
    reason <- if (fit_error) "lavaan::sem() errored"
      else if (!converged) "did not converge"
      else if (!se_finite) "standard errors not finite / implausibly large (>=5 SD) -- likely underidentified"
      else "negative residual variance (Heywood case)"
    cat("  joint lavaan fit not usable (", reason, "). Falling back per Kyuri's plan item 8:\n",
        "  single-equation OLS for every node's regression on this combo's actual parent\n",
        "  structure, feedback-strength sweep for the one coefficient that closes the loop\n",
        "  (weather_risk_prep -> belief_concern) instead of trying to force identification.\n",
        sep = "")
  }

  if (!identified || !isTRUE(stable)) {

    method_used <- "fallback_sweep"

    # single-equation OLS per node, on THIS combo's actual parent set --
    # correct for every non-loop edge and for the loop's two non-contested
    # edges; the loop-closing edge itself is swept, not estimated here.
    fit_one_equation <- function(node) {
      parents <- e$from[e$to == node]
      if (length(parents) == 0) return(numeric(0))
      dat <- as.data.frame(scale(df_extended[, c(node, parents), drop = FALSE]))
      form <- stats::as.formula(paste(node, "~ 0 +", paste(parents, collapse = " + ")))
      m <- stats::lm(form, data = dat)
      stats::setNames(stats::coef(m), parents)
    }

    coef_lookup_fallback <- list()
    for (node in all_nodes) {
      cf <- fit_one_equation(node)
      for (p in names(cf)) coef_lookup_fallback[[paste(p, node, sep = "->")]] <- unname(cf[[p]])
    }

    # the flipped edge that closes the loop, in THIS combo's orientation
    loop_from <- "weather_risk_prep"; loop_to <- "belief_concern"
    loop_key <- paste(loop_from, loop_to, sep = "->")
    point_estimate_c <- stats::cor(df_extended[[loop_from]], df_extended[[loop_to]])  # bivariate corr = single-predictor std beta, same fact used earlier this session
    sweep_vals <- sort(unique(c(seq(-0.6, 0.6, by = 0.02), round(point_estimate_c, 4))))

    for (c_val in sweep_vals) {
      coef_lookup_this <- coef_lookup_fallback
      coef_lookup_this[[loop_key]] <- c_val
      B_sweep <- build_B(e, coef_lookup_this)
      rho_sweep <- spectral_radius(B_sweep)
      stable_sweep <- isTRUE(rho_sweep < 1)
      for (node in node_targets) {
        dy <- if (stable_sweep) equilibrium_ate(B_sweep, node) else NA_real_
        sweep_rows[[length(sweep_rows) + 1]] <- tibble::tibble(
          combo = nm, feedback_coef_c = c_val, is_point_estimate = isTRUE(all.equal(c_val, round(point_estimate_c, 4))),
          rho_B = rho_sweep, stable = stable_sweep, node = node, ate_equilibrium = dy
        )
      }
    }

    # also file the point-estimate row (c = raw empirical correlation) into
    # the main comparison table, clearly flagged as the fallback/approximate
    # method rather than the joint-fit method
    coef_lookup_point <- coef_lookup_fallback
    coef_lookup_point[[loop_key]] <- point_estimate_c
    B_point <- build_B(e, coef_lookup_point)
    rho_B <- spectral_radius(B_point)
    detIB <- det(diag(n_all) - B_point)
    condIB <- tryCatch(kappa(diag(n_all) - B_point, exact = TRUE), error = function(err) NA_real_)
    stable <- isTRUE(rho_B < 1)
    cat("  fallback point estimate (feedback coef = raw correlation =", signif(point_estimate_c, 4),
        "): rho(B) =", signif(rho_B, 4), "-> stable:", stable, "\n")
    if (stable) {
      for (node in node_targets) {
        dy <- equilibrium_ate(B_point, node)
        near_ate <- get_nearest_acyclic_ate(nearest_nm, node)
        ate_rows[[length(ate_rows) + 1]] <- tibble::tibble(
          combo = nm, node = node, method = method_used,
          ate_equilibrium = dy, nearest_acyclic_combo = nearest_nm,
          nearest_acyclic_ate = near_ate,
          amplification_ratio = ifelse(abs(near_ate) > 1e-6, dy / near_ate, NA_real_),
          amplification_diff = dy - near_ate
        )
      }
    }
  }

  diagnostics[[length(diagnostics) + 1]] <- tibble::tibble(
    combo = nm, nearest_acyclic_combo = nearest_nm, method = method_used,
    joint_fit_converged = converged, joint_fit_se_finite = se_finite,
    joint_fit_no_heywood = no_heywood, joint_fit_identified = identified,
    det_I_minus_B = detIB, cond_I_minus_B = condIB, rho_B = rho_B, stable = stable
  )
}

diagnostics_out <- dplyr::bind_rows(diagnostics)
ate_out <- dplyr::bind_rows(ate_rows)
sweep_out <- if (length(sweep_rows) > 0) dplyr::bind_rows(sweep_rows) else NULL

cat("\n=== per-combo diagnostics ===\n")
print(as.data.frame(diagnostics_out), row.names = FALSE)

cat("\n=== equilibrium single-node ATEs vs nearest-acyclic (amplification) ===\n")
print(as.data.frame(ate_out), row.names = FALSE)

write.csv(diagnostics_out, file.path(TABLES_DIR, "cyclic_feedback_diagnostics.csv"), row.names = FALSE)
write.csv(ate_out, file.path(TABLES_DIR, "cyclic_feedback_ate_equilibrium.csv"), row.names = FALSE)
if (!is.null(sweep_out)) {
  write.csv(sweep_out, file.path(TABLES_DIR, "cyclic_feedback_strength_sweep.csv"), row.names = FALSE)
  cat("\nWrote", file.path(TABLES_DIR, "cyclic_feedback_strength_sweep.csv"),
      "(", nrow(sweep_out), "rows -- one or more combos needed the feedback-strength",
      "sensitivity sweep instead of a single joint-fit point estimate; see diagnostics",
      "above for which and why).\n")
}

# ---- Unified 16-combo table for the proposed Supp. figure (Kyuri's item 6:
# diamonds = baseline, light points = 11 other acyclic, different symbol =
# 4 cyclic equilibrium) -- combines 07's 12 acyclic rows with this script's
# 4 equilibrium rows into one long table, tagged by combo_type. -------------
acyclic_tagged <- acyclic_ate |>
  dplyr::mutate(combo_type = ifelse(scenario == "combo_0", "baseline", "acyclic")) |>
  dplyr::rename(ate = ate_climate_behavior) |>
  dplyr::select(scenario, node, ate, combo_type)

cyclic_tagged <- if (nrow(ate_out) > 0) {
  ate_out |>
    dplyr::group_by(combo, node) |>
    dplyr::slice(1) |>  # one row per (combo, node): the joint-fit or fallback-point-estimate row
    dplyr::ungroup() |>
    dplyr::transmute(scenario = combo, node = node, ate = ate_equilibrium, combo_type = "cyclic_equilibrium")
} else {
  warning("No cyclic combo produced a stable equilibrium (neither joint-fit nor fallback ",
          "point estimate) -- see diagnostics above. The all-16 table below has only the ",
          "12 acyclic rows; nothing to add from this run.")
  acyclic_tagged[0, ]
}

all16_out <- dplyr::bind_rows(acyclic_tagged, cyclic_tagged)
write.csv(all16_out, file.path(TABLES_DIR, "orientation_enumeration_ate_all16_with_equilibrium.csv"), row.names = FALSE)
cat("\nWrote", file.path(TABLES_DIR, "orientation_enumeration_ate_all16_with_equilibrium.csv"),
    "(", nrow(all16_out), "rows: 12 acyclic combos x 8 nodes + 4 cyclic combos x 8 nodes --",
    "input for the proposed Supp. figure comparing recursive vs equilibrium intervention effects.)\n")

cat("\n===== 19_cyclic_feedback_equilibrium.R done =====\n")
