# =============================================================================
# 10_cci_addition_and_figure4_pc_appendix.R
#
# UPDATED per the CCI orientation-instability finding (r_patches/11c and
# 11d): a synthetic-data permutation check showed CCI.KP::cci()'s ORIENTATION
# marks (which end gets tail/arrowhead/circle) are order-dependent -- the
# same true causal structure, same data, produced different (and in one
# case reversed) edge orientations depending purely on variable/column
# order. Skeleton recovery (which pairs are connected at all) was stable
# across that same test, so CCI is kept here ONLY as a skeleton-level
# sensitivity check (Part 3), not as a main-text structure-discovery method.
# Figure 4 is FCI-only in main text; PC-stable (cycle-caveated, see
# r_patches/11) and CCI (orientation-unreliable) are BOTH appendix-level
# checks now, for different, independently-discovered reasons.
#
# Originally: added CCI as a third causal-discovery algorithm alongside FCI,
# restructuring Figure 4 so FCI + CCI were the main-text comparison and
# PC-stable moved to its own supplementary figure -- prompted by a reviewer
# question about why PC was used instead of CCI, and a follow-up about
# moving PC to the appendix. That FCI+CCI design is superseded by the note
# above.
#
# WHERE THIS FITS IN THE PIPELINE
# --------------------------------
# Run climate_analysis_avg_v2_altweather.qmd first (interactively or knit),
# at least through Section 7.2's single-run fits. This script assumes the
# following objects already exist in the session:
#   node_order_cd, abbr, p, agg_data, suffStat_gauss
#   fci_05, pc_amat_05                        (Section 7.2 single-run fits)
# weather_risk_prep = ew5 alone is already baked into df_main / df_main_pre5
# in the _altweather copy, so agg_data/suffStat_gauss already reflect the
# changed variable set -- this script is purely the CCI addition on top of
# that, nothing here touches the composite definitions again.
#
# UPDATED 2026-09-11: the old Part 4 (the FCI-detail + PC-stable stability-
# matrix figures, figS_stability_matrix_detailed.pdf / figS_stability_pc.pdf)
# moved to clean_pipeline/20_stability_matrix_figures.R, which is fully
# qmd-independent (it reads node_order_cd/fci_props_main/pc_props_main/
# fci_marks/pc_marks from clean_pipeline/03_bootstrap_causal_discovery.R
# instead). Regenerate those two figures from there, not from here. The CCI
# logic below (Parts 1-3) stays qmd-dependent by design -- see
# claude/pipeline-refactor-plan.md's "Triaged and explicitly NOT migrated"
# list, which names the CCI diagnostic chain specifically -- so this script
# still needs the qmd's plotAG()/fci_props_combined/pc_props_combined/
# fci_marks/pc_marks for that reason, not out of oversight.
#
# WHAT THIS SCRIPT DOES
# ----------------------
# 1. Single-run CCI graphs at alpha = .05 / .01 (mirrors the FCI/PC
#    single-run panels in Section 7.2).
# 2. Extends the bootstrap logic to also fit CCI at both alphas x 1000
#    resamples, giving cci_props_combined in the same [p x p x 4] N/o/>/-
#    format as fci_props_combined. Cached to .rds so you don't have to
#    rerun the bootstrap just to re-plot.
# 3. Three-way skeleton comparison (FCI vs PC vs CCI), extending 7.4.
#
# WHAT TO CHECK
# --------------
# The console output (esp. the three-way skeleton comparison in Part 3).
# The Part 4 figures (figS_stability_matrix_detailed.pdf, figS_stability_pc.pdf)
# now come from clean_pipeline/20_stability_matrix_figures.R -- see there.
# =============================================================================

suppressPackageStartupMessages({
  library(pcalg)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(patchwork)
  library(furrr)
  library(future)
  if (requireNamespace("CCI.KP", quietly = TRUE)) {
    library(CCI.KP)
  } else if (requireNamespace("CCI", quietly = TRUE)) {
    library(CCI)
    message("CCI.KP not found -- using CCI. Results may differ slightly ",
            "(CCI.KP is the maintained fork with a corrected $maag extraction).")
  } else {
    stop("Neither CCI.KP nor CCI is installed. Install one before running ",
         "this script, e.g. install.packages('CCI.KP').")
  }
})

stopifnot(
  exists("node_order_cd"), exists("agg_data"), exists("suffStat_gauss"),
  exists("abbr"), exists("plotAG"), exists("p"),
  exists("fci_05"), exists("pc_amat_05")
  # fci_props_combined/pc_props_combined/fci_marks/pc_marks are no longer
  # required here -- they were only needed by the old Part 4 matrix figures,
  # which moved to clean_pipeline/20_stability_matrix_figures.R.
)

source("r_patches/03_figure_style.R")

# ---- helper: run CCI once, return an amat-style 0/1/2/3 matrix -------------
fit_cci <- function(suffStat, alpha, labels, p) {
  if (requireNamespace("CCI.KP", quietly = TRUE)) {
    fit <- CCI.KP::cci(
      suffStat = suffStat, indepTest = pcalg::gaussCItest,
      labels = labels, alpha = alpha, p = p,
      skel.method = "stable.fast", numCores = 1L, verbose = FALSE)
    matrix(as.integer(fit$maag), p, p, dimnames = list(labels, labels))
  } else {
    fit <- CCI::cci(suffStat = suffStat, indepTest = pcalg::gaussCItest,
                     alpha = alpha, p = p)
    matrix(as.integer(as(fit, "matrix")), p, p, dimnames = list(labels, labels))
  }
}

# =============================================================================
# 1. Single-run CCI graphs (mirrors Section 7.2's fci_05/fci_01/pc_05/pc_01)
# =============================================================================
cci_amat_05 <- fit_cci(suffStat_gauss, alpha = 0.05, labels = node_order_cd, p = p)
cci_amat_01 <- fit_cci(suffStat_gauss, alpha = 0.01, labels = node_order_cd, p = p)

dir.create("figures", showWarnings = FALSE)

# NOTE: no plotAG()-based mark figure is generated from cci_amat_05/01 here.
# 11c/11d showed CCI's specific tail/arrowhead/circle marks are order-
# dependent (not invariant to a causally-meaningless column permutation),
# so a figure decorating edges with those marks would misrepresent how
# reliable they are. cci_amat_05/01 are still used below for the
# skeleton-only comparison (Part 3), which IS stable under that same test.

# =============================================================================
# 2. Bootstrap: fit CCI at both alphas x 1000 resamples, same scheme as
#    Section 7.3's run_one() for FCI/PC.
#
#    NOTE: given the orientation-instability finding above, treat
#    cci_props_combined's mark-specific proportions (o/>/-) as unreliable.
#    Only the aggregate "is there an edge at all between these two nodes"
#    rate (i.e. 1 - the "N" proportion) reflects something checked as
#    stable. This output is no longer fed into Figure 4 (Part 4) -- it's
#    kept here in case you want an adjacency-only bootstrap stability
#    statistic for CCI in the supplement. Skip this whole Part 2 block if
#    it would be preferable to not spend the compute on a rerun; nothing downstream in
#    this script depends on it once Part 4's CCI panel is removed.
# =============================================================================
n_boot    <- 1000
alphas    <- c("0.05" = 0.05, "0.01" = 0.01)
cci_marks <- c("N", "o", ">", "-")   # same 4-slot encoding as fci_marks

run_one_cci <- function(b, data, nms, alpha) {
  p_local <- length(nms)
  idx     <- sample(nrow(data), replace = TRUE)
  suff_b  <- list(C = cor(as.matrix(data[idx, ])), n = length(idx))
  cci_arr <- array(0L, dim = c(p_local, p_local, 4L))
  tryCatch({
    am <- fit_cci(suff_b, alpha = alpha, labels = nms, p = p_local)
    for (i in seq_len(p_local)) for (j in seq_len(p_local)) {
      if (i == j) next
      slot <- am[i, j] + 1L
      cci_arr[i, j, slot] <- cci_arr[i, j, slot] + 1L
    }
  }, error = function(e) NULL)
  cci_arr
}

n_cores <- max(1L, parallelly::availableCores() - 1L)
plan(multisession, workers = n_cores)

set.seed(42)  # same seed as Section 7.3's FCI/PC bootstrap
cci_results <- list()
for (alph in names(alphas)) {
  cat("CCI bootstrap alpha =", alph, "\n")
  cci_results[[alph]] <- furrr::future_map(
    seq_len(n_boot), run_one_cci,
    data = agg_data, nms = node_order_cd, alpha = alphas[[alph]],
    .options  = furrr::furrr_options(seed = TRUE),
    .progress = TRUE)
}
plan(sequential)

cci_mark_props <- list()
for (alph in names(alphas)) {
  cci_c <- array(0L, dim = c(p, p, 4L),
                 dimnames = list(node_order_cd, node_order_cd, cci_marks))
  for (res in cci_results[[alph]]) {
    if (is.null(res)) next
    cci_c <- cci_c + res
  }
  cci_sums <- apply(cci_c, c(1, 2), sum); diag(cci_sums) <- n_boot
  if (any(abs(cci_sums - n_boot) > 1))
    warning(sprintf(paste("CCI alpha %s: cell sums deviate from %d --",
                           "some bootstrap fits likely failed silently;",
                           "check cci_results[['%s']] for NULLs."),
                     alph, n_boot, alph))
  cci_mark_props[[alph]] <- cci_c / n_boot
}

cci_props_combined <- (cci_mark_props[["0.05"]] * n_boot +
                        cci_mark_props[["0.01"]] * n_boot) / (2 * n_boot)
dimnames(cci_props_combined) <- list(node_order_cd, node_order_cd, cci_marks)

dir.create("data/causal_results", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(cci_props_combined = cci_props_combined,
             cci_amat_05 = cci_amat_05, cci_amat_01 = cci_amat_01),
        "data/causal_results/cci_main_pipeline_rerun.rds")
message("Saved data/causal_results/cci_main_pipeline_rerun.rds -- reload ",
        "with readRDS() to re-plot without rerunning the bootstrap.")

# =============================================================================
# 3. Three-way skeleton comparison (extends Section 7.4)
#    This is the recommended, sole use of CCI's output in the paper: which
#    pairs are connected at all, not how they're oriented (see header note).
# =============================================================================
pc_skel  <- (pc_amat_05 + t(pc_amat_05)) > 0
fci_skel <- (fci_05@amat + t(fci_05@amat)) > 0
cci_skel <- (cci_amat_05 + t(cci_amat_05)) > 0

edge_pairs <- function(skel, nms) {
  which(skel, arr.ind = TRUE) |>
    as.data.frame() |>
    dplyr::filter(row < col) |>
    dplyr::mutate(from = nms[row], to = nms[col]) |>
    dplyr::select(from, to)
}

cat("=== FCI skeleton edges (alpha=.05) ===\n");  print(edge_pairs(fci_skel, node_order_cd))
cat("\n=== PC skeleton edges (alpha=.05) ===\n");  print(edge_pairs(pc_skel,  node_order_cd))
cat("\n=== CCI skeleton edges (alpha=.05) ===\n"); print(edge_pairs(cci_skel, node_order_cd))
cat("\n=== In all three ===\n")
print(edge_pairs(fci_skel & pc_skel & cci_skel, node_order_cd))
cat("\n=== FCI & CCI but not PC (candidate latent-confounding edges) ===\n")
print(edge_pairs(fci_skel & cci_skel & !pc_skel, node_order_cd))

# Part 4 (the FCI-detail + PC-stable stability-matrix figures) moved to
# clean_pipeline/20_stability_matrix_figures.R on 2026-09-11 -- see the
# header note above. Run that script (after clean_pipeline/00-03) to
# regenerate figS_stability_matrix_detailed.pdf / figS_stability_pc.pdf;
# it no longer depends on this file or the qmd.

# =============================================================================
# 5. NOT generating a CCI mark figure here (see header note): 11c/11d showed
#    CCI's tail/arrowhead/circle marks are order-dependent, so a figure
#    decorating edges with them would misrepresent their reliability. Kept
#    as a numbered, empty section so the script's section numbering still
#    matches its history/console messages below.
# =============================================================================

cat("\n===== DONE =====\n")
cat("New objects in the session: cci_amat_05, cci_amat_01, cci_props_combined,\n")
cat("cci_skel, cci_results, cci_mark_props.\n")
cat("New/changed files:\n")
cat("  (none from this script -- the stability-matrix figures moved to\n")
cat("   clean_pipeline/20_stability_matrix_figures.R, run that separately)\n")
cat("No CCI mark/orientation figures produced -- see header note on the\n")
cat("11c/11d order-instability finding. CCI's skeleton comparison (Part 3\n")
cat("console output above) is its sole recommended use in the paper.\n")
cat("  data/causal_results/cci_main_pipeline_rerun.rds (cached bootstrap output)\n")
