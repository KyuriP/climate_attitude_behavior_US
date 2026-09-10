# =============================================================================
# 16_sensitivity_rcot_fci.R
#
# SENSITIVITY CHECK #3, companion to 14_sensitivity_rank_copula_fci.R and
# 15_sensitivity_mixed_ci_fci.R. Where script 14 still assumes an underlying
# Gaussian copula (semiparametric) and script 15 turned out to be unusable on
# this data (micd::mixCItest()'s discrete G^2 sub-test hits a hard
# sample-size floor -- see script 15's header and the project doc for the
# diagnosis), this script uses a GENUINELY nonparametric CI test: RCoT
# (Randomized Conditional independence Test), a kernel/random-Fourier-
# feature-based test (Strobl, Zhang & Visweswaran 2019, JMLR) that makes no
# distributional assumption on any variable -- continuous, ordinal, or
# binary alike are just embedded via the same kernel. This directly
# addresses whether RCoT/RCIT could resolve the social_norms ->
# climate_behavior edge, which came out messier than the rest under the
# primary test.
#
# WHY THIS EXISTS AS ITS OWN SCRIPT RATHER THAN A FIX TO SCRIPT 15: RCoT
# doesn't need a continuous/discrete split at all (no DISCRETE_NODES, no
# ordered-factor coding, no sparse-cell failure mode) -- every node is just a
# numeric column. That's a different enough data-prep path that mirroring
# script 15's structure directly would mean a lot of dead code, so this is a
# clean companion script instead, following the same conventions (node lists,
# ALPHAS, endpoint-mark extraction, bootstrap/comparison shape) as 14/15 so
# it's easy to read side by side with them.
#
# CI TEST: RCoT() from the `RCIT` package (github.com/ericstrobl/RCIT --
# NOT on CRAN). Depends on `MASS` (base R) and `momentchi2` (small CRAN
# package providing the sw/hbe/lpb4 null-distribution approximations RCoT
# uses). License: CC BY-NC 4.0 (non-commercial) -- flag to the paper's coauthors
# before using this beyond an internal sensitivity check, separately from the
# statistical question.
#
# INSTALLATION: RCIT is GitHub-only, so the normal `install.packages()` path
# doesn't apply.
#   install.packages("momentchi2")   # CRAN, should just work
#   install.packages("remotes")      # if not already installed
#   remotes::install_github("ericstrobl/RCIT")
# If install_github is blocked, RCIT is small and pure-R (no compiled code)
# -- source()-ing its 8 .R files directly (RCoT.R, RCIT.R, RIT.R, matrix2.R,
# random_fourier_features.R, normalize.R, Sta_perm.R, repmat.R, from the
# GitHub repo's R/ folder) works identically.
#
# ============================ TESTING STATUS ================================
# RCIT/RCoT has been verified on synthetic data (not yet on the live
# df_main/df_extended):
#   - RCoT() correctly returns a large p-value for true independence (both
#     unconditional and conditional) and a near-zero p-value for true
#     dependence (linear, nonlinear, and conditional), across several
#     hand-built synthetic cases;
#   - tested specifically on data shaped like the single-item nodes -- a
#     5-level ordinal variable and a binary outcome with a continuous
#     confounder -- runs without error and gives sensible answers on both
#     the unconditional and conditional tests;
#   - the rcotCItest() wrapper below (matching pcalg's documented
#     function(x, y, S, suffStat) indepTest contract exactly) plugs cleanly
#     into a stub pcalg::fci() call -- i.e. the same "swap the indepTest
#     argument" pattern already used for rank-copula (script 14) and
#     mixed-CI (script 15) works for RCoT too;
#   - timed at ~6ms per CI-test call on n=300 (a plausible bootstrap-resample
#     size) -- a 1000-rep bootstrap across several edges should be
#     comfortably feasible time-wise, unlike script 15's sample-size wall.
# What this does NOT establish: that RCoT behaves well on the actual dataset
# (real skew/sparsity of the single-item nodes, real n), or that pcalg's
# real fci()/pc() drives it correctly end to end. BEFORE trusting this
# script's output:
#   1. In a fresh R session, confirm `RCIT::RCoT(rnorm(100), rnorm(100))`
#      runs and returns a p-value close to what's expected (large, since
#      those are independent) -- a trivial smoke test that the package
#      installed correctly.
#   2. Run the SINGLE-run fci() block below first (not the bootstrap) and
#      sanity-check the printed PAG against the primary NPN-based fci_05
#      before letting the full bootstrap loop go.
#   3. RCoT's num_f/num_f2 (random-feature counts) and the "lpd4" null
#      approximation are the package's own recommended defaults (used
#      as-is below) -- if you see many warnings/NaNs from lpb4(), the
#      package's own fallback to hbe() should catch it silently, but worth
#      watching the console the first time this runs on real data.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

if (!requireNamespace("pcalg", quietly = TRUE)) {
  stop("pcalg is required (as the fci()/pc() driver) and isn't installed -- ",
       "install it first.")
}
if (!requireNamespace("RCIT", quietly = TRUE)) {
  stop("RCIT is required for RCoT() and isn't installed -- RCIT is GitHub-only, ",
       "not on CRAN. Try:\n",
       "  install.packages(\"momentchi2\")\n",
       "  install.packages(\"remotes\")\n",
       "  remotes::install_github(\"ericstrobl/RCIT\")\n",
       "See header for a source()-based fallback if install_github is blocked.")
}
# IMPORTANT: library(RCIT), not just requireNamespace() above -- RCIT's own
# DESCRIPTION lists momentchi2 as a Depends, not an Imports, and its NAMESPACE
# doesn't formally import it either, so RCoT()'s internal (unqualified) calls
# to hbe()/sw()/lpb4() only resolve once momentchi2 has actually been
# ATTACHED to the search path. Only library()/require() do that --
# requireNamespace() and RCIT::RCoT() (namespace-qualified, without a prior
# library(RCIT)) do NOT, and fail with "could not find function \"hbe\"" as
# soon as RCoT() hits its default lpd4/hbe null-approximation branch.
# Confirmed as the exact cause: building RCIT and momentchi2 as two separate
# installed packages (matching a normal
# `remotes::install_github("ericstrobl/RCIT")` install) reproduces this
# precise error with a namespace-qualified call; library(RCIT) fixes it.
suppressPackageStartupMessages(library(RCIT))

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

stopifnot(exists("df_main"))

# No continuous/discrete split needed -- RCoT treats every column as a
# numeric variable regardless of how many distinct values it takes, so the
# single-item ordinal/binary nodes go in as their raw numeric-coded values,
# unlike script 15's ordered-factor coding.
make_rcot_df <- function(raw, nodes) {
  raw |> dplyr::select(dplyr::all_of(nodes)) |> as.data.frame() |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
}

rcot_main <- make_rcot_df(df_main, MAIN_NODES)

cat("Column classes going into the RCoT suffStat (all numeric by design):\n")
print(sapply(rcot_main, class))
cat("\nDistinct-value counts -- eyeball this to see which columns are the\n")
cat("coarse single-item nodes vs. the continuous composites:\n")
print(sapply(rcot_main, function(x) length(unique(x))))

ALPHAS <- c(a05 = .05, a01 = .01)

# indepTest contract: function(x, y, S, suffStat) -> p-value, matching
# pcalg's documented interface exactly (same contract gaussCItest,
# micd::mixCItest, and script 14's rank-copula wrapper all satisfy). suffStat
# here is a plain list holding the data.frame, x/y/S are column indices into
# it -- pcalg's own convention (see e.g. ?gaussCItest, where suffStat$C /
# suffStat$n play the analogous role).
rcotCItest <- function(x, y, S, suffStat, approx = "lpd4", num_f = 100, num_f2 = 5, seed = NULL) {
  dat <- suffStat$data
  xv <- dat[[x]]
  yv <- dat[[y]]
  if (length(S) == 0) {
    out <- RCIT::RCoT(xv, yv, approx = approx, num_f2 = num_f2, seed = seed)
  } else {
    zv <- as.matrix(dat[, S, drop = FALSE])
    out <- RCIT::RCoT(xv, yv, zv, approx = approx, num_f = num_f, num_f2 = num_f2, seed = seed)
  }
  out$p
}

suffStat_rcot_main <- list(data = rcot_main)

fit_one_rcot <- function(suffStat, nodes, alpha, method = c("fci", "pc")) {
  method <- match.arg(method)
  f <- if (method == "fci") pcalg::fci else pcalg::pc
  f(suffStat, indepTest = rcotCItest, labels = nodes, alpha = alpha, verbose = FALSE)
}

# ---- Single-run fits (sanity-check these first; see header) ----------------
fci_rcot_a05 <- fit_one_rcot(suffStat_rcot_main, MAIN_NODES, ALPHAS[["a05"]], "fci")
fci_rcot_a01 <- fit_one_rcot(suffStat_rcot_main, MAIN_NODES, ALPHAS[["a01"]], "fci")
pc_rcot_a05  <- fit_one_rcot(suffStat_rcot_main, MAIN_NODES, ALPHAS[["a05"]], "pc")
pc_rcot_a01  <- fit_one_rcot(suffStat_rcot_main, MAIN_NODES, ALPHAS[["a01"]], "pc")

message("Single-run RCoT fits done. Sanity-check e.g.:")
message('  pcalg::plot(fci_rcot_a05, main = "RCoT FCI, alpha=.05")')
message("against the primary NPN-based fci_05 before running the bootstrap below.")

# ---- Endpoint-mark extraction -- identical to scripts 14/15, repeated here
# so this script stands alone; convention verified against pcalg's own docs
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
# identical to scripts 14/15's helper of the same name, repeated here so this
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
# used (assumed .05 here -- see the same note in scripts 14/15).
#
# N_BOOT: starting at 10, not 1000 -- same "start small" pattern as 14/15,
# since this hasn't been run against the live data or real pcalg yet. RCoT
# itself is fast (~6ms/call in testing on n=300), so once the single-
# run sanity check looks right, bumping N_BOOT up should be cheap compared to
# script 15's mixed-CI bootstrap.
BOOT_ALPHA <- 0.05
N_BOOT <- 10
set.seed(20260903)

run_boot_network_rcot <- function(dat, nodes, n_boot, alpha) {
  n <- nrow(dat)
  pairs <- all_pairs(nodes)

  fit_rep <- function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    rs <- dat[idx, , drop = FALSE]
    suffStat_rs <- list(data = rs)
    fit <- tryCatch(
      pcalg::fci(suffStat_rs, indepTest = rcotCItest, labels = nodes, alpha = alpha, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    fit@amat
  }

  amats <- if (requireNamespace("furrr", quietly = TRUE) && requireNamespace("future", quietly = TRUE)) {
    furrr::future_map(seq_len(n_boot), fit_rep, .options = furrr::furrr_options(seed = TRUE))
  } else {
    message("furrr/future not available -- running bootstrap sequentially (RCoT itself is ",
            "fast, so this should still be manageable at N_BOOT=10; parallelize before ",
            "bumping to 1000 -- install furrr and call future::plan(future::multisession) first).")
    lapply(seq_len(n_boot), fit_rep)
  }

  n_failed <- sum(vapply(amats, is.null, logical(1)))
  if (n_failed > 0) {
    message(sprintf("  %d / %d bootstrap replicates failed (fci/RCoT error) and were skipped.",
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

message(sprintf("Running %d-replicate RCoT bootstrap (alpha=%.2f) ...", N_BOOT, BOOT_ALPHA))
rcot_boot_main <- run_boot_network_rcot(rcot_main, MAIN_NODES, N_BOOT, BOOT_ALPHA)
if (!is.null(rcot_boot_main)) write.csv(rcot_boot_main, "pipeline_outputs/bootstrap_fci_stability_RCoT_main.csv", row.names = FALSE)

# ---- Compare MAIN network against the primary (NPN-based) bootstrap -------
# Prefers the real in-session object (fci_props_combined) over the CSV
# fallback, same as scripts 14/15 -- see their headers for why (no file
# literally named bootstrap_fci_stability.csv actually exists in the project).
primary_main_edges <- if (exists("fci_props_combined")) {
  array_to_edge_table(fci_props_combined)
} else if (file.exists("bootstrap_fci_stability.csv")) {
  read.csv("bootstrap_fci_stability.csv", stringsAsFactors = FALSE)
} else {
  NULL
}

if (!is.null(rcot_boot_main) && !is.null(primary_main_edges)) {
  cmp <- primary_main_edges |>
    dplyr::inner_join(rcot_boot_main, by = c("A", "B"), suffix = c("_npn", "_RCoT")) |>
    dplyr::mutate(abs_diff_p_adjacent = abs(p_adjacent_npn - p_adjacent_RCoT)) |>
    dplyr::arrange(dplyr::desc(abs_diff_p_adjacent))
  write.csv(cmp, "pipeline_outputs/sensitivity_RCoT_vs_primary_main.csv", row.names = FALSE)
  cat("\n===== MAIN network: edges where RCoT vs. primary NPN bootstrap disagree most =====\n")
  print(utils::head(cmp[, c("A", "B", "p_adjacent_npn", "p_adjacent_RCoT", "abs_diff_p_adjacent")], 10))
} else if (is.null(rcot_boot_main)) {
  message("Bootstrap produced no usable results -- see warning above.")
} else {
  message("Neither fci_props_combined (in session) nor bootstrap_fci_stability.csv ",
          "(on disk) found -- skipping the main-network comparison step.")
}

# =============================================================================
# EXTENDED NETWORK (+ climate_behavior) -- same JCI/background-knowledge
# mechanism as scripts 14/15's extended sections, read directly from
# climate_analysis_avg_v2_altweather.qmd Section 7.4: it's
# native to pcalg's fci()/pc()+addBgKnowledge(), independent of which
# indepTest/suffStat is used, so it carries over to RCoT the same way it did
# to rank-copula (script 14) and mixed-CI (script 15).
# =============================================================================

if (exists("df_extended")) {
  EXTENDED_NODES <- c(MAIN_NODES, "climate_behavior")
  rcot_ext <- make_rcot_df(df_extended, EXTENDED_NODES)
  suffStat_rcot_ext <- list(data = rcot_ext)

  # Attitude nodes are context (upstream) variables for FCI-JCI -- identical
  # construction to Section 7.4's `context_idx <- which(node_order_ext %in%
  # node_order_cd)`, and to scripts 14/15's context_idx_ext.
  context_idx_ext <- which(EXTENDED_NODES %in% MAIN_NODES)

  # ---- Single-run FCI-JCI + PC-with-background-knowledge fits --------------
  # Mirrors Section 7.4's ext-fci/ext-pc chunks exactly, just with rcotCItest
  # swapped in for gaussCItest and suffStat_rcot_ext (list(data=...)) swapped
  # in for suffStat_ext.
  fci_rcot_ext_a05 <- pcalg::fci(
    suffStat_rcot_ext, indepTest = rcotCItest, labels = EXTENDED_NODES,
    alpha = ALPHAS[["a05"]], contextVars = context_idx_ext, jci = "1",
    selectionBias = FALSE, verbose = FALSE
  )
  fci_rcot_ext_a01 <- pcalg::fci(
    suffStat_rcot_ext, indepTest = rcotCItest, labels = EXTENDED_NODES,
    alpha = ALPHAS[["a01"]], contextVars = context_idx_ext, jci = "1",
    selectionBias = FALSE, verbose = FALSE
  )

  .pc_with_bg_rcot <- function(alpha) {
    pc_fit <- pcalg::pc(suffStat_rcot_ext, indepTest = rcotCItest, alpha = alpha,
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
  pc_rcot_ext_amat_a05 <- .pc_with_bg_rcot(ALPHAS[["a05"]])
  pc_rcot_ext_amat_a01 <- .pc_with_bg_rcot(ALPHAS[["a01"]])

  message("Extended-network single-run fits done (RCoT, FCI-JCI + PC-with-background-knowledge).")
  message('Sanity-check e.g.: pcalg::plot(fci_rcot_ext_a05, main = "RCoT FCI-JCI extended, alpha=.05")')
  message("against the primary fci_ext_05 before running the extended bootstrap below.")

  # ---- Extended-network bootstrap (FCI-JCI only, matching scripts 14/15's
  # and the primary bootstrap's shape -- PC-with-background-knowledge is
  # single-run-only here, same as PC is single-run-only for the main network
  # above) ---------------------------------------------------------------
  run_boot_network_rcot_ext <- function(dat, nodes, n_boot, alpha, contextVars, jci) {
    n <- nrow(dat)
    pairs <- all_pairs(nodes)

    fit_rep <- function(b) {
      idx <- sample.int(n, n, replace = TRUE)
      rs <- dat[idx, , drop = FALSE]
      suffStat_rs <- list(data = rs)
      fit <- tryCatch(
        pcalg::fci(suffStat_rs, indepTest = rcotCItest, labels = nodes, alpha = alpha, verbose = FALSE,
                   contextVars = contextVars, jci = jci, selectionBias = FALSE),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      fit@amat
    }

    amats <- if (requireNamespace("furrr", quietly = TRUE) && requireNamespace("future", quietly = TRUE)) {
      furrr::future_map(seq_len(n_boot), fit_rep, .options = furrr::furrr_options(seed = TRUE))
    } else {
      message("furrr/future not available -- running bootstrap sequentially. ",
              "Install furrr and call future::plan(future::multisession) first to parallelize.")
      lapply(seq_len(n_boot), fit_rep)
    }

    n_failed <- sum(vapply(amats, is.null, logical(1)))
    if (n_failed > 0) {
      # Uses sprintf() wrapping a paste0()-joined string, same as the other
      # scripts in this batch (14, 15) -- R doesn't auto-concatenate adjacent
      # quoted string literals, so a split literal here crashes with
      # "invalid format '%d'; use format %s for character objects" as soon as
      # n_failed > 0. The main-network instance of this message (above, in
      # run_boot_network_rcot()) already uses a single unsplit literal and is
      # fine as-is.
      message(sprintf(
        paste0("  %d / %d EXTENDED-network bootstrap replicates failed (fci/RCoT error) ",
               "and were skipped."),
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

  message(sprintf("Running %d-replicate RCoT bootstrap, EXTENDED network (alpha=%.2f) ...", N_BOOT, BOOT_ALPHA))
  rcot_boot_ext <- run_boot_network_rcot_ext(rcot_ext, EXTENDED_NODES, N_BOOT, BOOT_ALPHA,
                                              contextVars = context_idx_ext, jci = "1")
  if (!is.null(rcot_boot_ext)) write.csv(rcot_boot_ext, "pipeline_outputs/bootstrap_fci_stability_RCoT_ext.csv", row.names = FALSE)

  # ---- Compare EXTENDED network against the primary (NPN-based) bootstrap --
  primary_ext_edges <- if (exists("fci_props_ext")) {
    array_to_edge_table(fci_props_ext)
  } else if (file.exists("pipeline_outputs/fci_props_ext.rds")) {
    array_to_edge_table(readRDS("pipeline_outputs/fci_props_ext.rds"))
  } else {
    NULL
  }

  if (!is.null(rcot_boot_ext) && !is.null(primary_ext_edges)) {
    cmp_ext <- primary_ext_edges |>
      dplyr::inner_join(rcot_boot_ext, by = c("A", "B"), suffix = c("_npn", "_RCoT")) |>
      dplyr::mutate(abs_diff_p_adjacent = abs(p_adjacent_npn - p_adjacent_RCoT)) |>
      dplyr::arrange(dplyr::desc(abs_diff_p_adjacent))
    write.csv(cmp_ext, "pipeline_outputs/sensitivity_RCoT_vs_primary_ext.csv", row.names = FALSE)
    cat("\n===== EXTENDED network: edges where RCoT vs. primary NPN bootstrap disagree most =====\n")
    print(utils::head(cmp_ext[, c("A", "B", "p_adjacent_npn", "p_adjacent_RCoT", "abs_diff_p_adjacent")], 10))
    sn_cb <- cmp_ext[(cmp_ext$A == "social_norms" & cmp_ext$B == "climate_behavior") |
                        (cmp_ext$A == "climate_behavior" & cmp_ext$B == "social_norms"), ]
    if (nrow(sn_cb) > 0) {
      cat("\n----- social_norms - climate_behavior specifically -----\n")
      print(sn_cb)
    }
  } else if (is.null(rcot_boot_ext)) {
    message("Extended bootstrap produced no usable results -- see warning above.")
  } else {
    message("Neither fci_props_ext (in session) nor fci_props_ext.rds (on disk) found -- ",
            "skipping the extended-network comparison step.")
  }
} else {
  message("df_extended not found in the session -- skipping the extended-network section entirely.")
}

message("Done.")
