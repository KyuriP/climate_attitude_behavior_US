# =============================================================================
# 15_sensitivity_mixed_ci_fci.R
#
# SENSITIVITY CHECK #2 of 2, companion
# to 14_sensitivity_rank_copula_fci.R. Where script 14 keeps every node
# continuous (just changes how the correlation matrix is estimated), this
# script instead treats the four genuinely single-item nodes as what they
# actually are -- ordered categorical variables -- and uses a conditional-
# Gaussian likelihood-ratio CI test built for exactly that mix, rather than
# forcing everything through a continuous-margin correlation matrix at all.
#
# Single-item (ordinal) nodes, confirmed directly against
# climate_analysis_avg_v2_altweather.qmd's own prose (Section 4, composite-
# validation summary, 2026-09-03): "trust_science, social_norms, and
# politics are retained as substantively specific single-item nodes," plus
# weather_risk_prep per the earlier-established ew5-alone fix (also already
# baked into this .qmd's df_main/df_main_pre5, per r_patches/10's header):
#   trust_science      (= cvcc9_cc)
#   social_norms       (= cvcc4_should)
#   politics           (= pol_ideology)
#   weather_risk_prep  (= ew5)
# Everything else (belief_concern, harm_present, harm_future, policy_support)
# is a genuine multi-item continuous composite and is left as-is.
#
# Originally SCOPED TO THE MAIN NETWORK ONLY (same update as script 14 --
# see its header for the full account), because Section 7.4's
# extended-network (+ climate_behavior) analysis applies a
# JCI/background-knowledge constraint (contextVars/jci="1" for FCI;
# addBgKnowledge() for PC), which was not yet accounted for. Confirmed
# directly from climate_analysis_avg_v2_altweather.qmd -- it's native to
# pcalg's fci()/pc()+addBgKnowledge(), independent of which indepTest is
# passed (gaussCItest or micd::mixCItest alike), so it carries over here the
# same way it did to script 14's rank-copula correlation. Added the
# extended-network section below (same motivation as script 14) to speak to
# the social_norms->climate_behavior edge (flagged elsewhere as not
# clean-cut), which the main-network-only scope couldn't reach.
#
# r_patches/10-13 already exist for unrelated work (CCI as a third
# causal-discovery algorithm, an IDA appendix, an all-wave sensitivity
# check, an uncertain-cell diagnostic) -- this script is numbered 15 (after
# script 14) to not collide with any of that.
#
# CI TEST: micd::mixCItest() (CRAN package `micd`, actively maintained) -- a
# conditional-Gaussian likelihood-ratio test for conditional independence
# between mixed continuous/discrete variables, documented by its own authors
# as designed "to be used within pcalg::skeleton, pcalg::pc or pcalg::fci" as
# a drop-in indepTest, exactly like gaussCItest is used in the primary
# pipeline. Passing suffStat as a plain data.frame (continuous columns
# numeric, discrete columns as ordered factors) is its documented interface
# -- no correlation-matrix construction step at all for this one.
#
# ============================ TESTING STATUS ================================
# micd could not be installed in the environment this script was drafted in
# (no apt package, and CRAN/Bioconductor/GitHub/r-universe/conda-forge were
# all unreachable there), and neither could pcalg (see script 14's header
# for that same limitation). That means nothing involving an actual
# mixCItest() or fci() call has been run against the real packages yet --
# unlike the other sensitivity/analysis scripts here. What WAS tested, with
# synthetic data built to mimic the real skew of these items (9-pt
# trust_science, 5-pt social_norms with a rare bottom category, 7-pt
# politics, 5-pt weather_risk_prep with a rare bottom category, N=870):
#   - the factor-coding step itself (numeric -> ordered factor) runs cleanly
#     and produces the right column classes;
#   - 1000 bootstrap resamples of that synthetic data, refactored back onto
#     the ORIGINAL factor levels each time (see refactor step below), NEVER
#     dropped an observed category entirely, though the smallest single-
#     category cell count seen across all reps/variables was as low as 3 --
#     confirming the sparse-cell risk micd's own documentation warns about
#     ("conditional on each combination of values of the discrete
#     variables") is real enough here to guard against, not just a
#     theoretical caveat.
#   - additionally, tiny STUB "pcalg" and "micd" packages were built (real
#     fci()/pc()/mixCItest() function signatures and the real fciAlgo/
#     pcAlgo S4 classes incl. an `amat` slot, but random internals) and ran
#     this ENTIRE script (the 9-node extended-network version that predated
#     the scope cut above) against them end to end with synthetic data
#     matching the real skew described above. It ran clean with zero R
#     errors: single-run fits, the bootstrap loop (tryCatch skip-counting,
#     the refactor-onto-original-levels step), endpoint-mark aggregation,
#     CSV writing, and the comparison-to-primary join all produced
#     correctly-shaped output, joined correctly against a copy of the real
#     bootstrap_fci_stability.csv. This is the same stub approach used for
#     script 14 (see its header), and caught one real bug there (a
#     dimnames-assignment error) -- nothing analogous surfaced here.
#   - the EXTENDED-NETWORK section below was tested the same way, after the
#     scope note above, against the same updated stub
#     "pcalg" (accepts contextVars/jci/addBgKnowledge without erroring, but
#     does NOT faithfully model JCI orientation) and a stub "micd" -- ran
#     clean end to end (single-run FCI-JCI + PC-with-background-knowledge
#     fits, the extended bootstrap, CSV output, comparison against a
#     synthetic stand-in for fci_props_ext). As with script 14, this confirms
#     the R plumbing only, not that real micd::mixCItest() + real JCI
#     orientation behaves correctly on the dataset.
# What this does NOT establish: that real micd::mixCItest() (as opposed to
# the stub's random placeholder) will actually run on the dataset, whether the
# category counts above are sparse enough to make it fail on many bootstrap
# replicates, or that the statistical results are correct -- only that the R
# code AROUND those calls is sound. BEFORE trusting this script's output:
#   1. In a fresh R session, confirm `install.packages("micd")` works and
#      `micd::mixCItest(1, 2, integer(0), suffStat = some_toy_df)` runs on a
#      trivial 3-column toy data.frame (2 continuous, 1 factor).
#   2. Run the SINGLE-run fci() block below first (not the bootstrap) and
#      sanity-check the printed PAG against the primary NPN-based fci_05
#      before letting the 1000-rep bootstrap loop go.
#   3. Expect some bootstrap replicates to fail/warn on sparse discrete-cell
#      combinations (micd's own documented caveat) -- these are caught and
#      skipped below, with a final count printed; if a large fraction fail,
#      that's itself informative (it says the single-item categories are too
#      sparse for this particular check to be reliable) and worth reporting
#      rather than working around.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

if (!requireNamespace("pcalg", quietly = TRUE)) {
  stop("pcalg is required (as the fci()/pc() driver) and isn't installed -- ",
       "install it first.")
}
if (!requireNamespace("micd", quietly = TRUE)) {
  stop("micd is required for mixCItest() and isn't installed -- ",
       "install.packages(\"micd\") first.")
}

if (exists("node_order_cd")) {
  MAIN_NODES <- node_order_cd
} else {
  MAIN_NODES <- c(
    "belief_concern", "harm_present", "harm_future", "policy_support",
    "trust_science", "social_norms", "politics", "weather_risk_prep"
  )
  message("node_order_cd not found in the session -- falling back to a ",
          "hardcoded 8-node list that matched it when last checked ",
          "(2026-09-03). Confirm this still matches Section 7.1's ",
          "node_order_cd before trusting anything below.")
}
DISCRETE_NODES <- c("trust_science", "social_norms", "politics", "weather_risk_prep")

stopifnot(exists("df_main"))

make_mixed_df <- function(raw, nodes, discrete_nodes) {
  d <- raw |> dplyr::select(dplyr::all_of(nodes)) |> as.data.frame()
  for (v in discrete_nodes) d[[v]] <- factor(d[[v]], ordered = TRUE)
  d
}

mixed_main <- make_mixed_df(df_main, MAIN_NODES, DISCRETE_NODES)

cat("Column classes going into the mixed suffStat:\n")
print(sapply(mixed_main, class))
cat("\nLevel counts for the four discrete nodes -- eyeball this for\n")
cat("categories with very few observations before trusting the bootstrap below:\n")
for (v in DISCRETE_NODES) { cat(" ", v, ": "); print(table(mixed_main[[v]])) }

ALPHAS <- c(a05 = .05, a01 = .01)

fit_one_mixed <- function(dat, nodes, alpha, method = c("fci", "pc")) {
  method <- match.arg(method)
  f <- if (method == "fci") pcalg::fci else pcalg::pc
  f(dat, indepTest = micd::mixCItest, labels = nodes, alpha = alpha, verbose = FALSE)
}

# ---- Single-run fits (sanity-check these first; see header) ----------------
fci_mixed_a05 <- fit_one_mixed(mixed_main, MAIN_NODES, ALPHAS[["a05"]], "fci")
fci_mixed_a01 <- fit_one_mixed(mixed_main, MAIN_NODES, ALPHAS[["a01"]], "fci")
pc_mixed_a05  <- fit_one_mixed(mixed_main, MAIN_NODES, ALPHAS[["a05"]], "pc")
pc_mixed_a01  <- fit_one_mixed(mixed_main, MAIN_NODES, ALPHAS[["a01"]], "pc")

message("Single-run mixed-CI fits done. Sanity-check e.g.:")
message('  pcalg::plot(fci_mixed_a05, main = "mixed-CI FCI, alpha=.05")')
message("against the primary NPN-based fci_05 before running the bootstrap below.")

# ---- Endpoint-mark extraction -- identical to script 14, repeated here so
# this script stands alone; convention verified against pcalg's own docs
# ("the edgemark-code refers to the column index" for amat.pag; example
# "amat[a,b]=2, amat[b,a]=3 implies a --> b") against a hand-built toy case.
extract_pair_marks <- function(amat, nodeA, nodeB) {
  mark_at_B <- amat[nodeA, nodeB]
  mark_at_A <- amat[nodeB, nodeA]
  c(adjacent = as.numeric(mark_at_A != 0 || mark_at_B != 0),
    arrow_at_A  = as.numeric(mark_at_A == 2),
    tail_at_A   = as.numeric(mark_at_A == 3),
    circle_at_A = as.numeric(mark_at_A == 1),
    arrow_at_B  = as.numeric(mark_at_B == 2),
    tail_at_B   = as.numeric(mark_at_B == 3),
    circle_at_B = as.numeric(mark_at_B == 1))
}
all_pairs <- function(nodes) combn(nodes, 2, simplify = FALSE)

# Converts a the standard pooled mark-proportion array format used elsewhere in this project (dimnames = list(nodes,
# nodes, c("N","o",">","-")), e.g. fci_props_combined / fci_props_ext) to the
# same long A,B,p_adjacent,... shape this script's own bootstrap produces --
# identical to script 14's helper of the same name, repeated here so this
# script stands alone.
array_to_edge_table <- function(props_array) {
  nodes <- dimnames(props_array)[[1]]
  pairs <- all_pairs(nodes)
  out <- lapply(pairs, function(pr) {
    A <- pr[1]; B <- pr[2]
    mark_at_B <- props_array[A, B, ]
    mark_at_A <- props_array[B, A, ]
    data.frame(A = A, B = B,
               p_adjacent = 1 - mark_at_A[["N"]],
               p_arrow_at_A = mark_at_A[[">"]], p_tail_at_A = mark_at_A[["-"]], p_circle_at_A = mark_at_A[["o"]],
               p_arrow_at_B = mark_at_B[[">"]], p_tail_at_B = mark_at_B[["-"]], p_circle_at_B = mark_at_B[["o"]])
  })
  dplyr::bind_rows(out) |> dplyr::arrange(dplyr::desc(p_adjacent))
}

# BOOT_ALPHA: match to whatever alpha the primary bootstrap_fci_stability.csv
# used (assumed .05 here -- see the same note in script 14).
#
# N_BOOT: starting at 10, not 1000 -- neither micd nor this fci()/mixCItest()
# combination is tested against the real packages (see header), and this CI
# test is heavier per-call than gaussCItest and more exposed to the sparse-
# discrete-cell failure mode. Run 10 first, look at the skip-count message
# and the resulting CSV, then bump to 1000 (or whatever matches the primary
# bootstrap) once that looks sane.
BOOT_ALPHA <- 0.05
N_BOOT <- 10
set.seed(20260903)

run_boot_network_mixed <- function(mixed_dat, nodes, discrete_nodes, n_boot, alpha) {
  n <- nrow(mixed_dat)
  pairs <- all_pairs(nodes)
  orig_levels <- lapply(discrete_nodes, function(v) levels(mixed_dat[[v]]))
  names(orig_levels) <- discrete_nodes

  fit_rep <- function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    rs <- mixed_dat[idx, , drop = FALSE]
    # Refactor onto the ORIGINAL level set every replicate -- if a resample
    # happens to drop the rarest category of a sparse item, this keeps the
    # factor's level set (and therefore mixCItest's expected cell structure)
    # consistent with the full-sample fit rather than silently shrinking it.
    for (v in discrete_nodes) rs[[v]] <- factor(rs[[v]], levels = orig_levels[[v]], ordered = TRUE)
    fit <- tryCatch(
      pcalg::fci(rs, indepTest = micd::mixCItest, labels = nodes, alpha = alpha, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    fit@amat
  }

  amats <- if (requireNamespace("furrr", quietly = TRUE) && requireNamespace("future", quietly = TRUE)) {
    furrr::future_map(seq_len(n_boot), fit_rep, .options = furrr::furrr_options(seed = TRUE))
  } else {
    message("furrr/future not available -- running bootstrap sequentially (slower, and ",
            "this CI test is itself heavier per-call than gaussCItest, so budget more time). ",
            "Install furrr and call future::plan(future::multisession) first to parallelize.")
    lapply(seq_len(n_boot), fit_rep)
  }

  n_failed <- sum(vapply(amats, is.null, logical(1)))
  if (n_failed > 0) {
    # NOTE (fixed 2026-09-03): R does not auto-concatenate adjacent quoted
    # string literals -- the split literals below were each being consumed as
    # a SUBSTITUTION VALUE for "%d" rather than as format-string continuation,
    # which crashes as soon as n_failed > 0 actually occurs on real data
    # ("invalid format '%d'; use format %s for character objects"). Found via
    # a deliberate audit after the same bug turned up in script 18; this
    # script had never been tested against a real replicate failure before.
    # paste0() now builds one complete format string before sprintf() sees it.
    message(sprintf(
      paste0("  %d / %d bootstrap replicates failed (mixCItest/fci error, most likely a sparse ",
             "discrete-cell combination per micd's own documented caveat) and were skipped. ",
             "A large fraction here is itself a finding worth reporting, not just a nuisance."),
      n_failed, n_boot))
  }
  amats <- amats[!vapply(amats, is.null, logical(1))]
  if (length(amats) == 0) {
    warning("Every bootstrap replicate failed -- returning NULL. Check the single-run ",
            "fit above works at all before debugging the bootstrap.")
    return(NULL)
  }

  out <- lapply(pairs, function(pr) {
    marks <- t(vapply(amats, function(am) extract_pair_marks(am, pr[1], pr[2]), numeric(7)))
    props <- colMeans(marks)
    data.frame(A = pr[1], B = pr[2],
               p_adjacent = props[["adjacent"]],
               p_arrow_at_A = props[["arrow_at_A"]], p_tail_at_A = props[["tail_at_A"]], p_circle_at_A = props[["circle_at_A"]],
               p_arrow_at_B = props[["arrow_at_B"]], p_tail_at_B = props[["tail_at_B"]], p_circle_at_B = props[["circle_at_B"]])
  })
  dplyr::bind_rows(out) |> dplyr::arrange(dplyr::desc(p_adjacent))
}

message(sprintf("Running %d-replicate mixed-CI bootstrap (alpha=%.2f) ...", N_BOOT, BOOT_ALPHA))
mixed_boot_main <- run_boot_network_mixed(mixed_main, MAIN_NODES, DISCRETE_NODES, N_BOOT, BOOT_ALPHA)
if (!is.null(mixed_boot_main)) write.csv(mixed_boot_main, "pipeline_outputs/bootstrap_fci_stability_mixedCI_main.csv", row.names = FALSE)

# ---- Compare MAIN network against the primary (NPN-based) bootstrap -------
# Prefers the real in-session object (fci_props_combined) over the CSV
# fallback, same as script 14 -- see its header for why (no file literally
# named bootstrap_fci_stability.csv actually exists in the project).
primary_main_edges <- if (exists("fci_props_combined")) {
  array_to_edge_table(fci_props_combined)
} else if (file.exists("bootstrap_fci_stability.csv")) {
  read.csv("bootstrap_fci_stability.csv", stringsAsFactors = FALSE)
} else {
  NULL
}

if (!is.null(mixed_boot_main) && !is.null(primary_main_edges)) {
  cmp <- primary_main_edges |>
    dplyr::inner_join(mixed_boot_main, by = c("A", "B"), suffix = c("_npn", "_mixedCI")) |>
    dplyr::mutate(abs_diff_p_adjacent = abs(p_adjacent_npn - p_adjacent_mixedCI)) |>
    dplyr::arrange(dplyr::desc(abs_diff_p_adjacent))
  write.csv(cmp, "pipeline_outputs/sensitivity_mixedCI_vs_primary_main.csv", row.names = FALSE)
  cat("\n===== MAIN network: edges where mixed-CI vs. primary NPN bootstrap disagree most =====\n")
  print(utils::head(cmp[, c("A", "B", "p_adjacent_npn", "p_adjacent_mixedCI", "abs_diff_p_adjacent")], 10))
} else if (is.null(mixed_boot_main)) {
  message("Bootstrap produced no usable results -- see warning above.")
} else {
  message("Neither fci_props_combined (in session) nor bootstrap_fci_stability.csv ",
          "(on disk) found -- skipping the main-network comparison step.")
}

# =============================================================================
# EXTENDED NETWORK (+ climate_behavior) -- added after reading Section 7.4's
# real JCI/background-knowledge code directly. See the header note above for
# why this exists and what it does and doesn't establish yet. climate_behavior
# itself is continuous (not one of the four
# single-item discrete nodes), so DISCRETE_NODES is unchanged -- only the
# node list and the JCI plumbing differ from the main-network section above.
# =============================================================================

if (exists("df_extended")) {
  EXTENDED_NODES <- c(MAIN_NODES, "climate_behavior")
  mixed_ext <- make_mixed_df(df_extended, EXTENDED_NODES, DISCRETE_NODES)

  # Attitude nodes are context (upstream) variables for FCI-JCI -- identical
  # construction to Section 7.4's `context_idx <- which(node_order_ext %in%
  # node_order_cd)`, and to script 14's context_idx_ext.
  context_idx_ext <- which(EXTENDED_NODES %in% MAIN_NODES)

  # ---- Single-run FCI-JCI + PC-with-background-knowledge fits --------------
  # Mirrors Section 7.4's ext-fci/ext-pc chunks exactly, just with
  # micd::mixCItest swapped in for gaussCItest and mixed_ext (data.frame,
  # discrete columns as ordered factors) swapped in for suffStat_ext.
  fci_mixed_ext_a05 <- pcalg::fci(
    mixed_ext, indepTest = micd::mixCItest, labels = EXTENDED_NODES,
    alpha = ALPHAS[["a05"]], contextVars = context_idx_ext, jci = "1",
    selectionBias = FALSE, verbose = FALSE
  )
  fci_mixed_ext_a01 <- pcalg::fci(
    mixed_ext, indepTest = micd::mixCItest, labels = EXTENDED_NODES,
    alpha = ALPHAS[["a01"]], contextVars = context_idx_ext, jci = "1",
    selectionBias = FALSE, verbose = FALSE
  )

  .pc_with_bg_mixed <- function(alpha) {
    pc_fit <- pcalg::pc(mixed_ext, indepTest = micd::mixCItest, alpha = alpha,
                         labels = EXTENDED_NODES, skel.method = "stable", verbose = FALSE)
    amat <- as(pc_fit@graph, "matrix")
    bg_pairs <- Filter(function(a)
      amat[a, "climate_behavior"] == 1 && amat["climate_behavior", a] == 1,
      MAIN_NODES)
    if (length(bg_pairs) > 0) {
      bg_result <- pcalg::addBgKnowledge(
        gInput = pc_fit@graph, x = bg_pairs, y = rep("climate_behavior", length(bg_pairs))
      )
      amat <- as(bg_result, "matrix")
    }
    amat
  }
  pc_mixed_ext_amat_a05 <- .pc_with_bg_mixed(ALPHAS[["a05"]])
  pc_mixed_ext_amat_a01 <- .pc_with_bg_mixed(ALPHAS[["a01"]])

  message("Extended-network single-run fits done (mixed-CI, FCI-JCI + PC-with-background-knowledge).")
  message('Sanity-check e.g.: pcalg::plot(fci_mixed_ext_a05, main = "mixed-CI FCI-JCI extended, alpha=.05")')
  message("against the primary fci_ext_05 before running the extended bootstrap below.")

  # ---- Extended-network bootstrap (FCI-JCI only, matching script 14's and
  # the primary bootstrap's shape -- PC-with-background-knowledge is
  # single-run-only here, same as PC is single-run-only for the main network
  # above) ---------------------------------------------------------------
  run_boot_network_mixed_ext <- function(mixed_dat, nodes, discrete_nodes, n_boot, alpha,
                                          contextVars, jci) {
    n <- nrow(mixed_dat)
    pairs <- all_pairs(nodes)
    orig_levels <- lapply(discrete_nodes, function(v) levels(mixed_dat[[v]]))
    names(orig_levels) <- discrete_nodes

    fit_rep <- function(b) {
      idx <- sample.int(n, n, replace = TRUE)
      rs <- mixed_dat[idx, , drop = FALSE]
      for (v in discrete_nodes) rs[[v]] <- factor(rs[[v]], levels = orig_levels[[v]], ordered = TRUE)
      fit <- tryCatch(
        pcalg::fci(rs, indepTest = micd::mixCItest, labels = nodes, alpha = alpha, verbose = FALSE,
                   contextVars = contextVars, jci = jci, selectionBias = FALSE),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      fit@amat
    }

    amats <- if (requireNamespace("furrr", quietly = TRUE) && requireNamespace("future", quietly = TRUE)) {
      furrr::future_map(seq_len(n_boot), fit_rep, .options = furrr::furrr_options(seed = TRUE))
    } else {
      message("furrr/future not available -- running bootstrap sequentially (slower). ",
              "Install furrr and call future::plan(future::multisession) first to parallelize.")
      lapply(seq_len(n_boot), fit_rep)
    }

    n_failed <- sum(vapply(amats, is.null, logical(1)))
    if (n_failed > 0) {
      # Same sprintf()/paste0() fix as the main-network block above -- see
      # that comment for the full explanation.
      message(sprintf(
        paste0("  %d / %d EXTENDED-network bootstrap replicates failed (mixCItest/fci error, ",
               "most likely a sparse discrete-cell combination) and were skipped."),
        n_failed, n_boot))
    }
    amats <- amats[!vapply(amats, is.null, logical(1))]
    if (length(amats) == 0) {
      warning("Every EXTENDED-network bootstrap replicate failed -- returning NULL.")
      return(NULL)
    }

    out <- lapply(pairs, function(pr) {
      marks <- t(vapply(amats, function(am) extract_pair_marks(am, pr[1], pr[2]), numeric(7)))
      props <- colMeans(marks)
      data.frame(A = pr[1], B = pr[2],
                 p_adjacent = props[["adjacent"]],
                 p_arrow_at_A = props[["arrow_at_A"]], p_tail_at_A = props[["tail_at_A"]], p_circle_at_A = props[["circle_at_A"]],
                 p_arrow_at_B = props[["arrow_at_B"]], p_tail_at_B = props[["tail_at_B"]], p_circle_at_B = props[["circle_at_B"]])
    })
    dplyr::bind_rows(out) |> dplyr::arrange(dplyr::desc(p_adjacent))
  }

  message(sprintf("Running %d-replicate mixed-CI bootstrap, EXTENDED network (alpha=%.2f) ...", N_BOOT, BOOT_ALPHA))
  mixed_boot_ext <- run_boot_network_mixed_ext(mixed_ext, EXTENDED_NODES, DISCRETE_NODES, N_BOOT, BOOT_ALPHA,
                                                contextVars = context_idx_ext, jci = "1")
  if (!is.null(mixed_boot_ext)) write.csv(mixed_boot_ext, "pipeline_outputs/bootstrap_fci_stability_mixedCI_ext.csv", row.names = FALSE)

  # ---- Compare EXTENDED network against the primary (NPN-based) bootstrap --
  primary_ext_edges <- if (exists("fci_props_ext")) {
    array_to_edge_table(fci_props_ext)
  } else if (file.exists("pipeline_outputs/fci_props_ext.rds")) {
    array_to_edge_table(readRDS("pipeline_outputs/fci_props_ext.rds"))
  } else {
    NULL
  }

  if (!is.null(mixed_boot_ext) && !is.null(primary_ext_edges)) {
    cmp_ext <- primary_ext_edges |>
      dplyr::inner_join(mixed_boot_ext, by = c("A", "B"), suffix = c("_npn", "_mixedCI")) |>
      dplyr::mutate(abs_diff_p_adjacent = abs(p_adjacent_npn - p_adjacent_mixedCI)) |>
      dplyr::arrange(dplyr::desc(abs_diff_p_adjacent))
    write.csv(cmp_ext, "pipeline_outputs/sensitivity_mixedCI_vs_primary_ext.csv", row.names = FALSE)
    cat("\n===== EXTENDED network: edges where mixed-CI vs. primary NPN bootstrap disagree most =====\n")
    print(utils::head(cmp_ext[, c("A", "B", "p_adjacent_npn", "p_adjacent_mixedCI", "abs_diff_p_adjacent")], 10))
    sn_cb <- cmp_ext[(cmp_ext$A == "social_norms" & cmp_ext$B == "climate_behavior") |
                        (cmp_ext$A == "climate_behavior" & cmp_ext$B == "social_norms"), ]
    if (nrow(sn_cb) > 0) {
      cat("\n----- social_norms - climate_behavior specifically -----\n")
      print(sn_cb)
    }
  } else if (is.null(mixed_boot_ext)) {
    message("Extended bootstrap produced no usable results -- see warning above.")
  } else {
    message("Neither fci_props_ext (in session) nor fci_props_ext.rds (on disk) found -- ",
            "skipping the extended-network comparison step.")
  }
} else {
  message("df_extended not found in the session -- skipping the extended-network section entirely.")
}

message("Done.")
