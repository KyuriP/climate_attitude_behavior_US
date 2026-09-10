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
# at least through Section 7.3's bootstrap. This script assumes the
# following objects already exist in the session:
#   node_order_cd, abbr, node_labels, p, agg_data, suffStat_gauss
#   fci_05, pc_amat_05                        (Section 7.2 single-run fits)
#   fci_props_combined, pc_props_combined     (Section 7.3 bootstrap output)
#   plotAG(), plotPC()                        (Section 7.1 plot helpers)
# weather_risk_prep = ew5 alone is already baked into df_main / df_main_pre5
# in the _altweather copy, so agg_data/suffStat_gauss already reflect the
# changed variable set -- this script is purely the CCI addition on top of
# that, nothing here touches the composite definitions again.
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
# 4. Rebuilds Figure 4 as FCI-ONLY (main text); PC-stable AND CCI are both
#    appendix-level checks now (see header note above). PC-stable's own
#    supplementary figure (figS_stability_pc.pdf) is unchanged.
# 5. CCI single-run mark figures are NOT generated (see header note) --
#    this step is now a no-op, left in place so the file numbering/
#    section structure below still matches the script's history.
#
# WHAT TO CHECK
# --------------
# The console output (esp. the three-way skeleton comparison in Part 3), plus:
#   figures/figS_stability_matrix_detailed.pdf (FCI per-end-mark detail, Supplement)
#   figures/figS_stability_pc.pdf              (PC-stable, appendix)
# No CCI mark/orientation figures are produced (see header note).
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
  exists("fci_05"), exists("pc_amat_05"),
  exists("fci_props_combined"), exists("pc_props_combined"),
  exists("fci_marks"), exists("pc_marks")   # needed by the Part 4 matrix redesign
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

# =============================================================================
# 4. Supplement-only edge-matrix figures (FCI detail + PC-stable appendix).
#    The main-text Figure 3/fig:causal now comes from 06_figure4_stability_
#    redesign.R's continuous dot-matrix (figures/fig4_stability_dotmatrix.pdf)
#    instead of the categorical design built here -- this section's output
#    (figS_stability_matrix_detailed.pdf, figS_stability_pc.pdf) is wired
#    into the Supplement only. PC-stable and CCI remain appendix-level checks
#    -- PC for the directed-cycle finding in r_patches/11, CCI for the
#    orientation-instability finding in 11c/11d.
#
#    Restyled (v4) to actually match the rest of the paper's visual system.
#    v3 got the STRUCTURE right -- one lower-triangle cell per pair, explicit
#    direction glyph, printed confidence number -- verified correct against
#    the live PDF. But v3 was never restyled to use 03_figure_style.R's
#    theme_pub/COL/FIG_FONT: it had its own ad hoc palette
#    (#2a78d6/#eb6834/#1baf7a/#B0B0B0), plain theme_minimal, and
#    alpha-blended the confidence proportion into fill darkness, which
#    produces a lot of muddy, hard-to-name in-between shades. Sitting next
#    to Figure 1's clean tint-and-colorbar look, that reads as off-brand and
#    unfinished. This version fixes the STYLE, not the logic.
#
#    What changed from v3:
#      - Fill is now a flat, light tint per category -- COL$blue_light for
#        directed, a matching light tint of COL$rust for bidirected (mixed
#        at the same ~88% white ratio that makes blue_light out of
#        blue_dark, so the two categories read as a matched pair), COL$panel
#        for absent -- instead of an alpha-blended saturated color. The
#        confidence number is still printed in every cell, so the fill no
#        longer needs to double as a magnitude encoding; that double-duty is
#        what produced the muddy gradient.
#      - Glyph + confidence number are now colored ink (COL$blue_dark /
#        COL$rust) on the light tint, the same "tint fill + saturated
#        stroke/text" grammar Figure 1 already uses for its diagonal density
#        insets -- instead of white text on a solid block, which is what
#        made v3 read like conditional-formatting-in-a-spreadsheet rather
#        than a designed figure.
#      - Real ggplot legend (scale_fill_manual + guide_legend) at the
#        bottom, matching Figure 1's bottom-legend convention, instead of a
#        separately-drawn swatch panel glued on underneath with patchwork.
#      - Row axis now shows full node names (node_labels_oneline) instead of
#        the same 2-3 letter code used on both axes -- matches Figure 1's
#        own convention (abbreviated codes on one axis, full names on the
#        other) instead of forcing the reader to hold 8 codes in memory.
#      - theme_pub / FIG_FONT / PT() throughout, so the type matches every
#        other figure in the paper instead of ggplot's default look.
#
#    Category logic (which endpoint-mark combination maps to which category)
#    is UNCHANGED from v3 -- only the rendering changed. The "uncertain"/
#    undirected bucket is still available for the PC-stable panel (FCI has
#    zero cells in it, confirmed by 13_uncertain_cell_diagnostic.R).
# =============================================================================

# Light tint of COL$rust, mixed at the same ~88% white / 12% color ratio
# that produces COL$blue_light from COL$blue_dark, so "directed" (blue) and
# "bidirected" (rust) read as a matched pair rather than two unrelated
# palettes bolted together.
rust_light <- "#F4EBE8"

edge_cat_colors <- c(absent = COL$panel, directed = COL$blue_light,
                      bidirected = rust_light, uncertain = "#E4EFE9")
edge_cat_text   <- c(absent = COL$ink_mid, directed = COL$blue_dark,
                      bidirected = COL$rust, uncertain = "#2F6B4F")
edge_cat_labels <- c(absent = "absent", directed = "directed (→ / ←)",
                      bidirected = "bidirected (↔)", uncertain = "undirected (○)")

plot_edge_matrix <- function(props, marks, title, nodes = node_order_cd,
                              include_undirected = TRUE) {
  p_local <- length(nodes)
  cells <- expand.grid(i = seq_len(p_local), j = seq_len(p_local),
                        stringsAsFactors = FALSE) |>
    dplyr::filter(i > j) |>   # ONE cell per pair (lower triangle only)
    dplyr::mutate(row_node = nodes[i], col_node = nodes[j]) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      mark_at_col = marks[which.max(props[i, j, ])],
      mark_at_row = marks[which.max(props[j, i, ])],
      prop_at_col = max(props[i, j, ]),
      prop_at_row = max(props[j, i, ]),
      edge_cat = dplyr::case_when(
        mark_at_row == "N" & mark_at_col == "N" ~ "absent",
        mark_at_row == "-" & mark_at_col == ">" ~ "directed",   # row -> col
        mark_at_row == ">" & mark_at_col == "-" ~ "directed",   # col -> row
        mark_at_row == ">" & mark_at_col == ">" ~ "bidirected",
        TRUE ~ "uncertain"
      ),
      glyph = dplyr::case_when(
        edge_cat == "absent" ~ "",
        edge_cat == "bidirected" ~ "↔",
        edge_cat == "uncertain" ~ "○",
        mark_at_row == "-" ~ "→",   # row -> col
        TRUE ~ "←"                  # col -> row
      ),
      # Confidence is now shown for EVERY cell, including "absent" (dropping
      # it for absent cells would imply those pairs are 100% certainly
      # unconnected, which isn't right -- absence has its own dominant-mark
      # proportion just like directed/bidirected does).
      conf  = pmin(prop_at_row, prop_at_col),
      label = dplyr::if_else(edge_cat == "absent",
                              sprintf("%.2f", conf),
                              paste0(glyph, "\n", sprintf("%.2f", conf))),
      edge_cat = factor(edge_cat, levels = c("absent", "directed", "bidirected", "uncertain"))
    ) |> dplyr::ungroup()

  legend_cats <- if (include_undirected) names(edge_cat_colors) else
    setdiff(names(edge_cat_colors), "uncertain")

  # No diagonal identity labels needed: with lower-triangle-only cells, the
  # row and column axis tick labels below already name every node exactly
  # once each (same convention as Figure 1's correlation matrix).
  ggplot() +
    geom_tile(data = cells, aes(x = col_node, y = row_node, fill = edge_cat),
              colour = "white", linewidth = 1.0, width = .98, height = .98) +
    geom_text(data = cells,
              aes(x = col_node, y = row_node, label = label, colour = edge_cat),
              family = FIG_FONT, size = PT(8), lineheight = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = edge_cat_colors, breaks = legend_cats,
                       labels = edge_cat_labels[legend_cats], name = NULL, drop = FALSE) +
    scale_colour_manual(values = edge_cat_text, breaks = legend_cats, drop = FALSE) +
    scale_x_discrete(limits = nodes[-p_local], labels = abbr[nodes[-p_local]],
                      position = "top") +
    scale_y_discrete(limits = rev(nodes[-1]), labels = node_labels_oneline[rev(nodes[-1])]) +
    coord_fixed(clip = "off") +
    labs(title = title) +
    theme_pub +
    theme(
      axis.title  = element_blank(),
      axis.line   = element_blank(),
      axis.ticks  = element_blank(),
      axis.text.x = element_text(face = "bold", size = 8.6),
      axis.text.y = element_text(size = 8.8),
      panel.grid  = element_blank(),
      plot.title  = element_text(face = "bold", size = 10.5, hjust = .5),
      legend.position       = "bottom",
      legend.justification  = "center",
      plot.margin = margin(4, 8, 2, 4)
    )
}

# Repurposed: the MAIN-text Figure 3 is the continuous dot-matrix (adjacency
# = size, orientation asymmetry = color) in 06_figure4_stability_redesign.R,
# not this per-end-mark-proportion matrix --
# the mode-based directed/bidirected/absent bucketing here can hide real
# asymmetry (e.g. future harm-trust science reads as bidirected here despite
# clearing the .01 asymmetry threshold comfortably). This matrix is still
# genuinely useful, though: it's the most precise, fully transparent view of
# "what proportion of bootstrap runs put which exact mark at which exact
# end," which is exactly what a careful reader would want after seeing the
# main figure's continuous summary. So it moves to the Supplement instead of
# being discarded, alongside the PC-stable matrix below.
fig_detail <- plot_edge_matrix(fci_props_combined, fci_marks, "FCI", include_undirected = FALSE)
save_ms_figure(fig_detail, "figures/figS_stability_matrix_detailed.pdf", FIG_DIMS_MM$fig4)
message("Saved figures/figS_stability_matrix_detailed.pdf (FCI per-end-mark detail, Supplement)")

figS_pc <- plot_edge_matrix(pc_props_combined, pc_marks, "PC-stable", include_undirected = TRUE)
save_ms_figure(figS_pc, "figures/figS_stability_pc.pdf", FIG_DIMS_MM$fig4)
message("Saved figures/figS_stability_pc.pdf (PC-stable, restyled v4, appendix)")

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
cat("  figures/figS_stability_matrix_detailed.pdf   (FCI per-end-mark detail, Supplement)\n")
cat("  figures/figS_stability_pc.pdf       (PC-stable, appendix)\n")
cat("No CCI mark/orientation figures produced -- see header note on the\n")
cat("11c/11d order-instability finding. CCI's skeleton comparison (Part 3\n")
cat("console output above) is its sole recommended use in the paper.\n")
cat("  data/causal_results/cci_main_pipeline_rerun.rds (cached bootstrap output)\n")
