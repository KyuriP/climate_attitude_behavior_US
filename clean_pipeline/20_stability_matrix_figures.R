# =============================================================================
# 20_stability_matrix_figures.R
#
# Split out of r_patches/10_cci_addition_and_figure4_pc_appendix.R's old
# "Part 4" (2026-09-11), specifically so these two Supplement figures no
# longer depend on climate_analysis_avg_v2_altweather.qmd at all. Everything
# below is a pure relocation, not a re-derivation: plot_edge_matrix() and its
# two call sites are unchanged from r_patches/10 (including the earlier
# lazy-evaluation fix -- force(props)/force(marks)/force(nodes), the
# unconditional mark_at_col/mark_at_row/prop_at_col/prop_at_row computation,
# dplyr::if_else() for conf, full lower-triangle-including-diagonal via
# i >= j, and full 8-node axes). Only the inputs changed name:
#   fci_props_combined -> fci_props_main   (from clean_pipeline/03's bootstrap)
#   pc_props_combined  -> pc_props_main    (from clean_pipeline/03's bootstrap)
# clean_pipeline/03_bootstrap_causal_discovery.R's bootstrap-combination
# formula ((mark_props[["0.05"]]$fci_counts + mark_props[["0.01"]]$fci_counts)
# / (2 * N_BOOT), same dimnames convention) was confirmed line-for-line
# identical in structure to the qmd's fci_props_combined/pc_props_combined --
# this is a rename, not new math.
#
# WHAT THIS SCRIPT PRODUCES
# --------------------------
#   figures/figS_stability_matrix_detailed.pdf (FCI per-end-mark detail, Supplement)
#   figures/figS_stability_pc.pdf              (PC-stable, appendix)
#
# WHERE THIS FITS IN THE PIPELINE (qmd-independent)
# ---------------------------------------------------
# Run, in ONE R session, in order:
#   source("clean_pipeline/00_config.R")
#   source("clean_pipeline/01_data_prep.R")
#   source("clean_pipeline/02_ggm.R")
#   source("clean_pipeline/03_bootstrap_causal_discovery.R")
#   source("clean_pipeline/20_stability_matrix_figures.R")
# (01/02/03 are the same prerequisite chain run_all.R already uses -- this
# script just needs 03's bootstrap output: node_order_cd, fci_props_main,
# pc_props_main, fci_marks, pc_marks.) If those objects are already in the
# session from an earlier 03 run today, this script reuses them as-is
# rather than re-running the 1000x2-alpha bootstrap; otherwise it falls
# back to loading the most recently saved fci_props_main_nboot*.rds /
# pc_props_main_nboot*.rds from pipeline_outputs/ (03 saves both there via
# its tag() helper on every run) so re-plotting never requires re-running
# the bootstrap.
#
# The CCI logic (Parts 1-3 of the old r_patches/10) is NOT part of this
# migration and stays qmd-dependent by design -- see
# claude/pipeline-refactor-plan.md's "Triaged and explicitly NOT migrated"
# list, which names the CCI diagnostic chain specifically. r_patches/10
# still exists for that purpose; its own stopifnot() and header have been
# trimmed to note that Part 4 moved here.
# =============================================================================

source("clean_pipeline/00_config.R")
source("r_patches/03_figure_style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---- load node_order_cd / fci_props_main / pc_props_main / fci_marks /
#      pc_marks from the current session if 03 already ran, otherwise from
#      the most recently saved .rds (avoids re-running the bootstrap just
#      to re-plot) --------------------------------------------------------
.load_from_latest_rds <- function(name) {
  pattern <- sprintf("^%s_nboot.*\\.rds$", name)
  cands <- list.files(OUTPUT_DIR, pattern = pattern, full.names = TRUE)
  if (length(cands) == 0) {
    stop("`", name, "` is not in this session and no ", pattern, " file was ",
         "found in ", OUTPUT_DIR, ". Run clean_pipeline/03_bootstrap_causal_",
         "discovery.R first (in this session, or on a prior run whose .rds ",
         "output is still in ", OUTPUT_DIR, ").")
  }
  latest <- cands[order(file.mtime(cands), decreasing = TRUE)][1]
  message("Loading ", name, " from ", latest, " (not found in current session).")
  readRDS(latest)
}

if (!exists("node_order_cd") || !exists("fci_marks") || !exists("pc_marks")) {
  stop("node_order_cd / fci_marks / pc_marks not in session -- these are ",
       "cheap to define (no bootstrap needed), so run ",
       "clean_pipeline/03_bootstrap_causal_discovery.R's setup lines (or the ",
       "whole script) first rather than reloading them from disk.")
}
if (!exists("fci_props_main")) fci_props_main <- .load_from_latest_rds("fci_props_main")
if (!exists("pc_props_main"))  pc_props_main  <- .load_from_latest_rds("pc_props_main")

# abbr: verbatim from the qmd's causal-setup chunk (climate_analysis_avg_v2_
# altweather.qmd, line ~1619) -- duplicated here rather than sourced, since
# duplicating one 8-line named vector is a smaller, more legible dependency
# than sourcing the whole qmd for it.
abbr <- c(
  belief_concern       = "BC",
  harm_present         = "HP",
  harm_future          = "HF",
  policy_support       = "PS",
  trust_science        = "TS",
  social_norms         = "SN",
  politics             = "POL",
  weather_risk_prep    = "WR"
)

dir.create("figures", showWarnings = FALSE)

# =============================================================================
# Supplement-only edge-matrix figures (FCI detail + PC-stable appendix).
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
                      bidirected = rust_light, uncertain = "#E4EFE9",
                      self = COL$panel)
edge_cat_text   <- c(absent = COL$ink_mid, directed = COL$blue_dark,
                      bidirected = COL$rust, uncertain = "#2F6B4F",
                      self = COL$panel)
edge_cat_labels <- c(absent = "absent", directed = "directed (→ / ←)",
                      bidirected = "bidirected (↔)", uncertain = "undirected (○)",
                      self = "")

plot_edge_matrix <- function(props, marks, title, nodes = node_order_cd,
                              include_undirected = TRUE) {
  force(props); force(marks); force(nodes)   # resolve promises up front, before rowwise()
  p_local <- length(nodes)
  cells <- expand.grid(i = seq_len(p_local), j = seq_len(p_local),
                        stringsAsFactors = FALSE) |>
    dplyr::filter(i >= j) |>   # full lower triangle, including the diagonal
    dplyr::mutate(row_node = nodes[i], col_node = nodes[j]) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      # Computed unconditionally for every row -- including the diagonal --
      # rather than guarded by a scalar if(i == j): the bootstrap loop leaves
      # props[i, i, ] at all-zeros, so which.max()/max() on it resolve
      # harmlessly (first slot, 0), and doing it this way avoids a base-R
      # if/else inside rowwise() forcing the lazy `props` promise on
      # whichever row happens to hit it first -- that pattern was losing the
      # environment chain back to the caller's fci_props_main /
      # pc_props_main and erroring "object not found". case_when()/
      # if_else() below are real vectorized dplyr verbs and don't have that
      # problem, so the diagonal is blanked out downstream instead.
      mark_at_col = marks[which.max(props[i, j, ])],
      mark_at_row = marks[which.max(props[j, i, ])],
      prop_at_col = max(props[i, j, ]),
      prop_at_row = max(props[j, i, ]),
      edge_cat = dplyr::case_when(
        i == j ~ "self",
        mark_at_row == "N" & mark_at_col == "N" ~ "absent",
        mark_at_row == "-" & mark_at_col == ">" ~ "directed",   # row -> col
        mark_at_row == ">" & mark_at_col == "-" ~ "directed",   # col -> row
        mark_at_row == ">" & mark_at_col == ">" ~ "bidirected",
        TRUE ~ "uncertain"
      ),
      glyph = dplyr::case_when(
        edge_cat %in% c("absent", "self") ~ "",
        edge_cat == "bidirected" ~ "↔",
        edge_cat == "uncertain" ~ "○",
        mark_at_row == "-" ~ "→",   # row -> col
        TRUE ~ "←"                  # col -> row
      ),
      # Confidence is shown for every real pair, including "absent" (dropping
      # it for absent cells would imply those pairs are 100% certainly
      # unconnected, which isn't right -- absence has its own dominant-mark
      # proportion just like directed/bidirected does). Diagonal ("self")
      # cells carry no bootstrap information and are always blanked here.
      conf  = dplyr::if_else(i == j, NA_real_, pmin(prop_at_row, prop_at_col)),
      label = dplyr::case_when(
        edge_cat == "self" ~ "",
        edge_cat == "absent" ~ sprintf("%.2f", conf),
        TRUE ~ paste0(glyph, "\n", sprintf("%.2f", conf))
      ),
      edge_cat = factor(edge_cat, levels = c("absent", "directed", "bidirected", "uncertain", "self"))
    ) |> dplyr::ungroup()

  legend_cats <- if (include_undirected) names(edge_cat_colors) else
    setdiff(names(edge_cat_colors), "uncertain")
  legend_cats <- setdiff(legend_cats, "self")   # diagonal placeholder never appears in the legend

  # Diagonal cells are included (as blank "self" tiles) purely so both axes
  # list every node and line up with their own row/column -- same convention
  # as Figure 1's correlation matrix, but shown rather than implied.
  ggplot() +
    geom_tile(data = cells, aes(x = col_node, y = row_node, fill = edge_cat),
              colour = "white", linewidth = 1.0, width = .98, height = .98) +
    geom_text(data = cells,
              aes(x = col_node, y = row_node, label = label, colour = edge_cat),
              family = FIG_FONT, size = PT(8), lineheight = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = edge_cat_colors, breaks = legend_cats,
                       labels = edge_cat_labels[legend_cats], name = NULL, drop = FALSE) +
    scale_colour_manual(values = edge_cat_text, breaks = legend_cats, drop = FALSE) +
    scale_x_discrete(limits = nodes, labels = abbr[nodes],
                      position = "top") +
    scale_y_discrete(limits = rev(nodes), labels = node_labels_oneline[rev(nodes)]) +
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
fig_detail <- plot_edge_matrix(fci_props_main, fci_marks, "FCI", include_undirected = FALSE)
save_ms_figure(fig_detail, "figures/figS_stability_matrix_detailed.pdf", FIG_DIMS_MM$fig_stability)
message("Saved figures/figS_stability_matrix_detailed.pdf (FCI per-end-mark detail, Supplement)")

figS_pc <- plot_edge_matrix(pc_props_main, pc_marks, "PC-stable", include_undirected = TRUE)
save_ms_figure(figS_pc, "figures/figS_stability_pc.pdf", FIG_DIMS_MM$fig_stability)
message("Saved figures/figS_stability_pc.pdf (PC-stable, restyled v4, appendix)")

cat("\n===== DONE (20_stability_matrix_figures.R) =====\n")
cat("New/changed files:\n")
cat("  figures/figS_stability_matrix_detailed.pdf   (FCI per-end-mark detail, Supplement)\n")
cat("  figures/figS_stability_pc.pdf       (PC-stable, appendix)\n")
cat("No qmd dependency -- see header for the clean_pipeline/00-03 prerequisite chain.\n")
