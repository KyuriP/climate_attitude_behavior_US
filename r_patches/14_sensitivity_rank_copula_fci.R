# =============================================================================
# 14_sensitivity_rank_copula_fci.R
#
# SENSITIVITY CHECK #1 of 2, responding
# to: several nodes going into FCI/PC (trust_science, social_norms, politics,
# weather_risk_prep) are single ordinal Likert items, not continuous
# composites, but the primary pipeline's huge::huge.npn() nonparanormal
# transform formally assumes continuous margins (ties occur with probability
# zero -- Harris & Drton 2013, JMLR 14, sec. 2). This script is an
# ALTERNATIVE correlation estimator for the exact same gaussCItest-based
# FCI/PC machinery already used in Section 7.1-7.4: instead of
# npn-transforming then taking Pearson correlations, it takes Spearman's
# rho directly on the observed (untransformed) data and maps it to a latent
# Gaussian-copula correlation via rho_latent = 2*sin(pi/6 * rho_spearman)
# (Harris & Drton 2013's OTHER worked example, RPC -- reported there as
# computationally about as cheap as ordinary Pearson PC). This sidesteps the
# continuous-margins assumption entirely rather than trying to satisfy it via
# a monotonic transform.
#
# Originally this script covered the MAIN network only -- Section 7.4's
# extended-network (+ climate_behavior) analysis encodes measurement order
# as background knowledge (a JCI formulation for FCI via
# `contextVars`/`jci="1"`; `addBgKnowledge()` for PC) so climate_behavior
# can't be inferred as a cause of the earlier nodes -- that mechanism is
# confirmed to be a native pcalg feature, independent of which
# indepTest/suffStat is passed -- so it carries over to the rank-copula
# correlation just as easily as to gaussCItest+Pearson. Added the
# extended-network section below, reproducing that exact mechanism, to
# speak to the social_norms->climate_behavior edge (flagged elsewhere as
# not clean-cut) which the main-network-only scope couldn't reach.
#
# r_patches/10-13 already exist in this folder for unrelated work (CCI as a
# third causal-discovery algorithm, an IDA appendix, an all-wave
# sensitivity check, an uncertain-cell diagnostic) -- none of that existed
# yet when these two sensitivity scripts were started, so this one is
# numbered 14 to avoid colliding with it.
#
# WHAT CHANGES vs. the primary pipeline: only the correlation matrix fed into
# suffStat$C. Everything downstream -- gaussCItest, fci()/pc(), the alpha
# grid, the JCI background-knowledge mechanism, the bootstrap-resampling idea
# -- is unchanged, so a big divergence between this script's results and the
# primary NPN-based results tells you the graph is sensitive to that specific
# modeling choice; close agreement tells you it isn't.
#
# ============================ TESTING STATUS ================================
# The MAIN-network portion of this script (everything through the first
# "Compare against the primary" block) has been run against real data with
# pcalg installed (N_BOOT=10 smoke test, 2026-09-03) and produced correctly-
# shaped, internally-consistent output (arrow+tail+circle proportions summed
# to p_adjacent on both ends of every pair) -- so that part is confirmed
# working against real pcalg, not just drafted against a stub. The
# EXTENDED-NETWORK section below is NOT yet confirmed against real pcalg --
# it was drafted and tested against a stub package that accepts
# `contextVars`/`jci`/`addBgKnowledge()` calls without erroring, but the
# stub's JCI orientation logic is NOT a faithful reimplementation (it
# ignores contextVars functionally), so this only confirms the surrounding
# R plumbing (argument-passing, the bootstrap loop, CSV output, the
# array-based comparison against fci_props_ext) is sound -- not that the
# real JCI mechanism behaves correctly here. Also verified before writing
# the original (main-network) version:
#   (a) the rho_latent = 2*sin(pi/6*rho_spearman) transform, applied to data
#       with 4 coarsened/discretized columns (mimicking the real single-item
#       nodes), reliably returns a valid correlation matrix (symmetric, unit
#       diagonal, positive semi-definite) -- confirmed over 200 bootstrap
#       resamples with no failures -- AND lands closer to the true
#       (pre-coarsening) latent correlation matrix than naive Pearson-on-
#       coarsened-data does (Frobenius distance 0.144 vs. 0.196 in the test
#       run; max per-cell error 0.037 vs. 0.051).
#   (b) the endpoint-mark extraction logic below (extract_pair_marks()) was
#       checked against pcalg's own DOCUMENTED amat.pag convention (fetched
#       directly from pcalg's help pages, not assumed): "the edgemark-code
#       refers to the COLUMN index" for amat.pag, with the worked example
#       "amat[a,b] = 2 and amat[b,a] = 3 implies a --> b" -- i.e. amat[i,j]
#       is the mark AT NODE j (the column), codes 0=none/1=circle/2=arrow/
#       3=tail. Verified against a hand-built 3-node toy amat.pag; also
#       reconfirmed indirectly by matching the live fci_props_ext array
#       (same N/o/>/- convention, dimname-indexed, not positional).
# RUN the extended-network single-run block first (not its bootstrap)
# and sanity-check fci_rank_ext_a05 against the primary fci_ext_05 before
# trusting the extended bootstrap.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

if (!requireNamespace("pcalg", quietly = TRUE)) {
  stop(
    "pcalg is required for this script and isn't installed -- ",
    "install it (install.packages(\"pcalg\")) before running. The fci()/pc() ",
    "calls below have not been run without pcalg available -- verify output ",
    "before trusting it."
  )
}

# Prefer the live node_order_cd from Section 7.1 if it's in the session --
# falls back to a hardcoded copy that matched it exactly when checked
# directly against climate_analysis_avg_v2_altweather.qmd on 2026-09-03.
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
RAW_MAIN <- df_main |> dplyr::select(dplyr::all_of(MAIN_NODES)) |> as.data.frame()

# ---- Harris & Drton (2013) rank-based latent correlation -------------------
# Deliberately applied to the RAW composite scores (df_main/df_extended), NOT
# to df_net_main/df_net_ext/agg_data/agg_ext (those are already
# huge::huge.npn()-transformed for the *primary* pipeline) -- the whole point
# of this check is to bypass the npn continuous-margin assumption, not stack
# another transform on top of it.
rank_latent_cor <- function(mat) {
  rho_s <- stats::cor(mat, method = "spearman", use = "pairwise.complete.obs")
  rho_l <- 2 * sin((pi / 6) * rho_s)
  diag(rho_l) <- 1
  rho_l[rho_l >  1] <-  1   # numerical guard; shouldn't trigger for |rho_s|<=1
  rho_l[rho_l < -1] <- -1
  # rho_l already carries the right p x p dimnames from cor()/arithmetic on
  # rho_s (colnames(mat) on both dimensions) -- do NOT reassign dimnames(mat)
  # here, that's the n x p data frame's OWN dimnames (row *labels* + column
  # names), not a p x p pair, and will error ("length of 'dimnames' [1] not
  # equal to array extent"). Caught via the stub-package plumbing test before
  # this ever reached pcalg.
  rho_l
}

C_main <- rank_latent_cor(RAW_MAIN)
n_main <- nrow(RAW_MAIN)

ALPHAS <- c(a05 = .05, a01 = .01)

fit_one <- function(C, n, nodes, alpha, method = c("fci", "pc")) {
  method <- match.arg(method)
  suffStat <- list(C = C, n = n)
  f <- if (method == "fci") pcalg::fci else pcalg::pc
  f(suffStat, indepTest = pcalg::gaussCItest, labels = nodes, alpha = alpha, verbose = FALSE)
}

# ---- Single-run fits (sanity-check these first; see header) ----------------
fci_rank_a05 <- fit_one(C_main, n_main, MAIN_NODES, ALPHAS[["a05"]], "fci")
fci_rank_a01 <- fit_one(C_main, n_main, MAIN_NODES, ALPHAS[["a01"]], "fci")
pc_rank_a05  <- fit_one(C_main, n_main, MAIN_NODES, ALPHAS[["a05"]], "pc")
pc_rank_a01  <- fit_one(C_main, n_main, MAIN_NODES, ALPHAS[["a01"]], "pc")

message("Single-run rank-copula fits done. Sanity-check e.g.:")
message('  pcalg::plot(fci_rank_a05, main = "rank-copula FCI, alpha=.05")')
message("against the primary NPN-based fci_05 before running the bootstrap below.")

# ---- Endpoint-mark extraction (amat.pag convention verified -- see header) --
extract_pair_marks <- function(amat, nodeA, nodeB) {
  mark_at_B <- amat[nodeA, nodeB]  # amat[i,j]: edgemark code refers to column j
  mark_at_A <- amat[nodeB, nodeA]  # amat[j,i]: edgemark code refers to column i
  c(adjacent = as.numeric(mark_at_A != 0 || mark_at_B != 0),
    arrow_at_A  = as.numeric(mark_at_A == 2),
    tail_at_A   = as.numeric(mark_at_A == 3),
    circle_at_A = as.numeric(mark_at_A == 1),
    arrow_at_B  = as.numeric(mark_at_B == 2),
    tail_at_B   = as.numeric(mark_at_B == 3),
    circle_at_B = as.numeric(mark_at_B == 1))
}

all_pairs <- function(nodes) combn(nodes, 2, simplify = FALSE)

# Converts the standard pooled mark-proportion array format used elsewhere in this project (dimnames = list(nodes,
# nodes, c("N","o",">","-")), e.g. fci_props_combined / fci_props_ext) to
# the same long A,B,p_adjacent,... shape this script's own bootstrap produces
# -- so the two are always directly comparable regardless of which primary
# object is in the session. Indexed BY NAME (not position), since that's
# what's actually safe here -- confirmed the real objects carry these exact
# dimnames, confirmed by inspecting fci_props_ext.rds directly (2026-09-03).
array_to_edge_table <- function(props_array) {
  nodes <- dimnames(props_array)[[1]]
  pairs <- all_pairs(nodes)
  out <- lapply(pairs, function(pr) {
    A <- pr[1]; B <- pr[2]
    mark_at_B <- props_array[A, B, ]  # column B -> mark at B
    mark_at_A <- props_array[B, A, ]  # column A -> mark at A
    data.frame(A = A, B = B,
               p_adjacent = 1 - mark_at_A[["N"]],
               p_arrow_at_A = mark_at_A[[">"]], p_tail_at_A = mark_at_A[["-"]], p_circle_at_A = mark_at_A[["o"]],
               p_arrow_at_B = mark_at_B[[">"]], p_tail_at_B = mark_at_B[["-"]], p_circle_at_B = mark_at_B[["o"]])
  })
  dplyr::bind_rows(out) |> dplyr::arrange(dplyr::desc(p_adjacent))
}

# ---- Bootstrap: resample rows, recompute rank-latent correlation, refit FCI,
# aggregate endpoint-mark proportions -- same output shape as the primary
# bootstrap results (A,B,p_adjacent,p_arrow_at_A,p_tail_at_A,p_circle_at_A,
# p_arrow_at_B,p_tail_at_B,p_circle_at_B), so the two are directly comparable
# row-by-row.
#
# BOOT_ALPHA: set to whichever single alpha the primary bootstrap used
# (this script assumes .05 -- confirm/adjust if yours was run at a different
# alpha).
#
# N_BOOT: starting at 10, not 1000 -- confirmed clean against real pcalg at
# N_BOOT=10 for the main network (2026-09-03); the extended-network section
# below hasn't had that same real-data confirmation yet, so keep starting
# small there too. Once a
# 10-rep run completes cleanly and looks sane, bump this to 1000 (or
# whatever matches the primary bootstrap) for the real run.
BOOT_ALPHA <- 0.05
N_BOOT <- 10
set.seed(20260903)

run_boot_network <- function(raw_data, nodes, n_boot, alpha,
                              contextVars = NULL, jci = "0") {
  n <- nrow(raw_data)
  pairs <- all_pairs(nodes)
  fit_rep <- function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    C_b <- tryCatch(rank_latent_cor(raw_data[idx, , drop = FALSE]), error = function(e) NULL)
    if (is.null(C_b)) return(NULL)
    fit <- tryCatch(
      pcalg::fci(list(C = C_b, n = n), indepTest = pcalg::gaussCItest,
                 labels = nodes, alpha = alpha, verbose = FALSE,
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
    # NOTE (fixed 2026-09-03): R does not auto-concatenate adjacent quoted
    # string literals the way C/Python do -- writing them as separate
    # sprintf() arguments treats the second literal as a SUBSTITUTION VALUE
    # for the first "%d", not as format-string continuation, and crashes
    # ("invalid format '%d'; use format %s for character objects") the first
    # time n_failed > 0 actually happens. Wrapping in paste0() first builds
    # one complete format string before it reaches sprintf(). Found via a
    # deliberate audit after this exact bug turned up in script 18; confirmed
    # this script had never been tested against a real replicate failure.
    message(sprintf(
      paste0("  %d / %d bootstrap replicates failed (fci() error or degenerate ",
             "resample) and were skipped."),
      n_failed, n_boot))
  }
  amats <- amats[!vapply(amats, is.null, logical(1))]

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

message(sprintf("Running %d-replicate rank-copula bootstrap, MAIN network (alpha=%.2f) ...", N_BOOT, BOOT_ALPHA))
rank_boot_main <- run_boot_network(RAW_MAIN, MAIN_NODES, N_BOOT, BOOT_ALPHA)
write.csv(rank_boot_main, "pipeline_outputs/bootstrap_fci_stability_rankcopula_main.csv", row.names = FALSE)

# ---- Compare MAIN network against the primary (NPN-based) bootstrap -------
# Prefers the real in-session object (fci_props_combined, from Section 7.3 /
# r_patches/cache_fig4_fig9_objects.R -- confirmed its real shape by
# inspecting fig4_fig9_cache.rds directly, 2026-09-03: an 8x8x4 array,
# dimnames = list(node_order_cd, node_order_cd, c("N","o",">","-"))) over the
# CSV fallback, since no file literally named bootstrap_fci_stability.csv
# actually exists anywhere in the project as far as can be determined.
primary_main_edges <- if (exists("fci_props_combined")) {
  array_to_edge_table(fci_props_combined)
} else if (file.exists("bootstrap_fci_stability.csv")) {
  read.csv("bootstrap_fci_stability.csv", stringsAsFactors = FALSE)
} else {
  NULL
}

if (!is.null(primary_main_edges)) {
  cmp_main <- primary_main_edges |>
    dplyr::inner_join(rank_boot_main, by = c("A", "B"), suffix = c("_npn", "_rankcopula")) |>
    dplyr::mutate(abs_diff_p_adjacent = abs(p_adjacent_npn - p_adjacent_rankcopula)) |>
    dplyr::arrange(dplyr::desc(abs_diff_p_adjacent))
  write.csv(cmp_main, "sensitivity_rankcopula_vs_primary_main.csv", row.names = FALSE)
  cat("\n===== MAIN network: edges where rank-copula vs. primary NPN bootstrap disagree most =====\n")
  print(utils::head(cmp_main[, c("A", "B", "p_adjacent_npn", "p_adjacent_rankcopula", "abs_diff_p_adjacent")], 10))
} else {
  message("Neither fci_props_combined (in session) nor bootstrap_fci_stability.csv ",
          "(on disk) found -- skipping the main-network comparison step.")
}

# =============================================================================
# EXTENDED NETWORK (+ climate_behavior) -- added after reading Section 7.4's
# real JCI/background-knowledge code directly. See the header note above for
# why this exists and what it does and doesn't establish yet.
# =============================================================================

if (exists("df_extended")) {
  EXTENDED_NODES <- c(MAIN_NODES, "climate_behavior")
  RAW_EXT <- df_extended |> dplyr::select(dplyr::all_of(EXTENDED_NODES)) |> as.data.frame()
  C_ext <- rank_latent_cor(RAW_EXT)
  n_ext <- nrow(RAW_EXT)

  # Attitude nodes are context (upstream) variables for FCI-JCI -- identical
  # construction to Section 7.4's `context_idx <- which(node_order_ext %in%
  # node_order_cd)`.
  context_idx_ext <- which(EXTENDED_NODES %in% MAIN_NODES)

  # ---- Single-run FCI-JCI + PC-with-background-knowledge fits --------------
  # Mirrors Section 7.4's `ext-fci-05`/`ext-fci-01`/`ext-pc-05`/`ext-pc-01`
  # chunks exactly, just with C_ext (rank-copula) swapped in for
  # suffStat_ext$C (npn+Pearson).
  fci_rank_ext_a05 <- pcalg::fci(
    suffStat = list(C = C_ext, n = n_ext), indepTest = pcalg::gaussCItest,
    alpha = ALPHAS[["a05"]], labels = EXTENDED_NODES,
    contextVars = context_idx_ext, jci = "1", selectionBias = FALSE, verbose = FALSE
  )
  fci_rank_ext_a01 <- pcalg::fci(
    suffStat = list(C = C_ext, n = n_ext), indepTest = pcalg::gaussCItest,
    alpha = ALPHAS[["a01"]], labels = EXTENDED_NODES,
    contextVars = context_idx_ext, jci = "1", selectionBias = FALSE, verbose = FALSE
  )

  .pc_with_bg <- function(alpha) {
    pc_fit <- pcalg::pc(suffStat = list(C = C_ext, n = n_ext), indepTest = pcalg::gaussCItest,
                         alpha = alpha, labels = EXTENDED_NODES, skel.method = "stable", verbose = FALSE)
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
  pc_rank_ext_amat_a05 <- .pc_with_bg(ALPHAS[["a05"]])
  pc_rank_ext_amat_a01 <- .pc_with_bg(ALPHAS[["a01"]])

  message("Extended-network single-run fits done (rank-copula, FCI-JCI + PC-with-background-knowledge).")
  message('Sanity-check e.g.: pcalg::plot(fci_rank_ext_a05, main = "rank-copula FCI-JCI extended, alpha=.05")')
  message("against the primary fci_ext_05 before running the extended bootstrap below.")

  # ---- Extended-network bootstrap (FCI-JCI only, matching the shape of the
  # main-network bootstrap above and of the primary bootstrap_fci_stability
  # output -- PC-with-background-knowledge is single-run-only here, same as
  # PC is single-run-only for the main network above) ------------------------
  message(sprintf("Running %d-replicate rank-copula bootstrap, EXTENDED network (alpha=%.2f) ...", N_BOOT, BOOT_ALPHA))
  rank_boot_ext <- run_boot_network(RAW_EXT, EXTENDED_NODES, N_BOOT, BOOT_ALPHA,
                                     contextVars = context_idx_ext, jci = "1")
  write.csv(rank_boot_ext, "bootstrap_fci_stability_rankcopula_ext.csv", row.names = FALSE)

  # ---- Compare EXTENDED network against the primary (NPN-based) bootstrap --
  # Prefers fci_props_ext (confirmed real shape by inspecting
  # fci_props_ext.rds directly, 2026-09-03: a 9x9x4 array, dimnames =
  # list(node_order_ext, node_order_ext, c("N","o",">","-"))).
  primary_ext_edges <- if (exists("fci_props_ext")) {
    array_to_edge_table(fci_props_ext)
  } else if (file.exists("pipeline_outputs/fci_props_ext.rds")) {
    array_to_edge_table(readRDS("pipeline_outputs/fci_props_ext.rds"))
  } else {
    NULL
  }

  if (!is.null(primary_ext_edges)) {
    cmp_ext <- primary_ext_edges |>
      dplyr::inner_join(rank_boot_ext, by = c("A", "B"), suffix = c("_npn", "_rankcopula")) |>
      dplyr::mutate(abs_diff_p_adjacent = abs(p_adjacent_npn - p_adjacent_rankcopula)) |>
      dplyr::arrange(dplyr::desc(abs_diff_p_adjacent))
    write.csv(cmp_ext, "sensitivity_rankcopula_vs_primary_ext.csv", row.names = FALSE)
    cat("\n===== EXTENDED network: edges where rank-copula vs. primary NPN bootstrap disagree most =====\n")
    print(utils::head(cmp_ext[, c("A", "B", "p_adjacent_npn", "p_adjacent_rankcopula", "abs_diff_p_adjacent")], 10))
    sn_cb <- cmp_ext[(cmp_ext$A == "social_norms" & cmp_ext$B == "climate_behavior") |
                        (cmp_ext$A == "climate_behavior" & cmp_ext$B == "social_norms"), ]
    if (nrow(sn_cb) > 0) {
      cat("\n----- social_norms - climate_behavior specifically -----\n")
      print(sn_cb)
    }
  } else {
    message("Neither fci_props_ext (in session) nor fci_props_ext.rds (on disk) found -- ",
            "skipping the extended-network comparison step.")
  }
} else {
  message("df_extended not found in the session -- skipping the extended-network section entirely.")
}

message("Done.")
