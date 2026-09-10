# =============================================================================
# 05_figures2_3_ggm_redesign.R  -- uses qgraph's own native plot()/qgraph()
# rendering, the look from the original analysis .qmd, rather than a custom
# ggplot2 geom_segment/geom_label reimplementation. Network estimation
# (bootnet EBICglasso, tuning=.5) is unaffected by this; only the rendering
# engine and its house-style parameters are involved.
#
# We call qgraph::qgraph() directly on network_main$graph / network_ext$graph
# (the estimated partial-correlation weight matrices), rather than
# plot(network_main, ...) on the bootnet estimateNetwork object. Visually
# this is the same thing -- bootnet's own plot method for these objects is
# itself a thin wrapper that calls qgraph::qgraph() on $graph -- but calling
# qgraph() directly means every argument below is tested explicitly (against
# a synthetic stand-in network, not real numbers) rather than depended on
# undocumented S3-forwarding behavior.
#
# Preserved from the previous ggplot2 version (both are still true here, just
# achieved through qgraph's own arguments instead of manual rescaling):
#   - identical node positions for shared nodes across panels A and B
#     (fixed_layout_matrix() from 03_figure_style.R, not qgraph::averageLayout(),
#      so the extended network's one extra node -- climate_behavior -- doesn't
#      perturb the eight shared nodes' positions between panels)
#   - one common edge-width scale across both panels (qgraph's own `maximum`
#     argument, set to the same COMMON_MAX for both plots, so a given line
#     width means the same partial-correlation magnitude in both panels;
#     `cut = 0` keeps this continuous rather than qgraph's default two-tier
#     thick/thin dichotomization)
#   - solid/dashed encoding of edge sign, not color alone, matching this
#     figure's existing caption ("solid blue positive, dashed vermillion
#     negative"). Implemented as a MATRIX (not a flat vector) so qgraph reads
#     the style for a given (row, col) node pair directly off that pair,
#     rather than depending on flat-vector edge ordering -- an easy place to
#     silently mislabel which edge is dashed.
#   - the manuscript's own blue/rust palette (COL$blue_dark / COL$rust) for
#     edge color, layered on top of qgraph's theme="colorblind" via explicit
#     posCol/negCol -- qgraph's colorblind theme still governs everything
#     else (node border, background, legend). To use qgraph's own default
#     colorblind-theme edge hues instead of the manuscript's blue/rust, just
#     delete the posCol/negCol lines below.
#
# Node fill colors and node labels come from the SAME node_family_colors /
# node_labels objects used everywhere else in the manuscript (03_figure_style.R)
# -- these already give weather_risk_prep its own distinct amber/orange fill,
# same grouping logic as the original .qmd's node_colors (just different
# specific hex values).
#
# Tested against a synthetic stand-in network built with qgraph::EBICglasso()
# directly on random correlated data (not bootnet, and not the live
# df_net_main/df_net_ext) -- confirms the script runs end to end and that
# shared layout + shared edge scale + dashed-negative-as-matrix all render
# correctly. Needs network_main / network_ext (or df_net_main/df_net_ext to
# re-estimate them) already in scope.
# =============================================================================

suppressPackageStartupMessages({
  library(bootnet)
  library(qgraph)
  library(huge)
  library(dplyr)
})
source("r_patches/03_figure_style.R")

# Nonparanormal transform, matching the original "network-setup" chunk
# exactly (huge::huge.npn(npn.func = "truncation"), applied separately to the
# main-node and extended-node subsets before estimateNetwork()) -- this is
# also what the Methods text (sec:ggm) already claims happens. Built here
# directly from df_main/df_extended (the two objects the other figure
# scripts already depend on existing) rather than assuming
# df_net_main/df_net_ext already exist, pre-transformed, from somewhere
# upstream. Guarded by exists() so it's a no-op if df_net_main/df_net_ext
# are already in scope from elsewhere.
MAIN_NODES <- c(
  "belief_concern", "harm_present", "harm_future", "policy_support",
  "trust_science", "social_norms", "politics", "weather_risk_prep"
)
EXTENDED_NODES <- c(MAIN_NODES, "climate_behavior")

if (!exists("df_net_main")) {
  stopifnot(exists("df_main"))
  df_net_main <- df_main |>
    dplyr::select(dplyr::all_of(MAIN_NODES)) |>
    huge::huge.npn(npn.func = "truncation") |>
    as.data.frame()
  colnames(df_net_main) <- MAIN_NODES
}
if (!exists("df_net_ext")) {
  stopifnot(exists("df_extended"))
  df_net_ext <- df_extended |>
    dplyr::select(dplyr::all_of(EXTENDED_NODES)) |>
    huge::huge.npn(npn.func = "truncation") |>
    as.data.frame()
  colnames(df_net_ext) <- EXTENDED_NODES
}

if (!exists("network_main")) {
  network_main <- bootnet::estimateNetwork(
    df_net_main, default = "EBICglasso", tuning = .5, corMethod = "cor"
  )
}
if (!exists("network_ext")) {
  network_ext <- bootnet::estimateNetwork(
    df_net_ext, default = "EBICglasso", tuning = .5, corMethod = "cor"
  )
}

W_main <- network_main$graph
W_ext  <- network_ext$graph
COMMON_MAX <- max(abs(c(W_main, W_ext)), na.rm = TRUE)

# Layout: qgraph::averageLayout() computed independently for each network --
# this is what the original exploratory .qmd used (Section 6.3/6.5:
# `layout_main <- qgraph::averageLayout(getWmat(network_main))`, `layout_ext
# <- qgraph::averageLayout(getWmat(network_ext))`), rather than the
# fixed-grid table below, which would put shared nodes at IDENTICAL
# positions across both panels. That cross-panel guarantee doesn't hold
# here: because each network's layout is solved independently (and the
# extended network has a ninth node, climate_behavior, that the main
# network doesn't), a shared node's position can shift slightly between
# panels A and B, same as in the original .qmd. The manuscript caption
# (fig:ggm) is worded to match -- see main.tex. To restore identical
# cross-panel positions instead, swap to:
#   layout_for <- function(W) fixed_layout_matrix(rownames(W))
layout_for <- function(W) qgraph::averageLayout(W)

# Per-edge line type as a MATRIX aligned to the input adjacency matrix's own
# dimnames -- see header comment above.
lty_matrix_for <- function(W) {
  matrix(ifelse(W < 0, 2, 1), nrow(W), ncol(W), dimnames = dimnames(W))
}

plot_network_qgraph <- function(W, title) {
  nodes <- rownames(W)
  qgraph::qgraph(
    W,
    layout    = layout_for(W),
    labels    = node_labels[nodes],
    color     = node_fill_for(nodes),
    theme     = "colorblind",
    posCol    = COL$blue_dark,
    negCol    = COL$rust,
    lty       = lty_matrix_for(W),
    maximum   = COMMON_MAX,
    cut       = 0,
    vsize     = 14,
    label.cex = 1,
    title     = title,
    title.cex = 1.2,
    mar       = c(3, 3, 6, 3)
  )
}

save_ms_basegraphics(
  function() plot_network_qgraph(W_main, "Main belief network"),
  "figures/fig2_ggm_main_fixed.pdf", FIG_DIMS_MM$fig2
)
save_ms_basegraphics(
  function() plot_network_qgraph(W_ext, "Extended network"),
  "figures/fig3_ggm_extended_fixed.pdf", FIG_DIMS_MM$fig3
)
message("Saved Figures 2 and 3 (qgraph-native rendering)")
