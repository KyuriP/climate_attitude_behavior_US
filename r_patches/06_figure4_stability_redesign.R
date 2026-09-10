# =============================================================================
# 06_figure4_stability_redesign.R
#
# Figure 3 / fig:causal, FCI-only. Replaces the earlier tile-matrix design
# (which lived in 10_cci_addition_and_figure4_pc_appendix.R) after
# identifying a real accuracy problem with it: that design bucketed every
# pair into "directed"/"bidirected"/"absent" using which single endpoint
# mark was most common at each end (a mode-based, categorical call). That
# bucketing can call a pair "bidirected" even when the continuous asymmetry
# statistic -- the actual quantity Table 2 uses to decide which edges are
# bootstrap-resolved vs. flip-tested -- clearly favors one direction (e.g.
# future harm-trust science shows as bidirected in the tile matrix but has
# asymmetry +.29, well past the .01 resolved threshold in
# 18_finalize_scm_specification_v4.R). A categorical bucket can't show that;
# a continuous color can.
#
# This dot-matrix design drops the PC-stable panel entirely -- FCI is
# main-text only, PC-stable is Supplement-only (a locked methodology
# decision, sec 4 of the confirmed study design) -- and keeps ONLY the FCI
# panel:
#   Dot SIZE  = adjacency stability (1 - avg "N" proportion across both ends).
#   Dot FILL  = orientation asymmetry (P(arrowhead at destination) -
#               P(arrowhead at source)), the SAME statistic and SAME .01
#               threshold used in Table 2 -- shown only when adjacency
#               stability >= .50 (below that, an asymmetry number isn't
#               meaningful -- the pair might not even be a real edge).
# This is rendered as its own PDF and combined with a single-run FCI
# example panel via LaTeX minipages in main2.tex (the same two-panel
# technique already used for Figure 2's GGM main/extended panels), rather
# than composing them on the R side.
#
# Uses the 9-node EXTENDED set (fci_props_ext/node_order_ext/abbr_ext),
# which includes climate_behavior, rather than the 8-node core set
# (fci_props_combined/node_order_cd) -- Figure 3 needs to include behavior.
# Both scale_x_discrete() and scale_y_discrete() set drop = FALSE: this is a
# lower-triangle-only combn() layout, so the very first node
# (belief_concern) never appears as a row value and the very last node
# never appears as a col value -- ggplot's default drop = TRUE would
# silently drop those unused factor levels from the axis labels, losing the
# diagonal. drop = FALSE forces both axes to show the full node set
# regardless. Output dimensions are 118x140mm (up from 105x124mm) to
# accommodate the larger 9x9 grid (estimate -- nudge further if labels
# crowd).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})
source("r_patches/03_figure_style.R")

stopifnot(exists("fci_props_ext"), exists("node_order_ext"), exists("abbr_ext"))

make_stability_df <- function(props, nodes) {
  pairs <- utils::combn(nodes, 2, simplify = FALSE)
  bind_rows(lapply(pairs, function(z) {
    a <- z[1]; b <- z[2]
    p_none <- mean(c(props[a, b, "N"], props[b, a, "N"]), na.rm = TRUE)
    tibble(
      row = b,
      col = a,
      adjacency = 1 - p_none,
      # Same statistic, same direction convention as Table 2 / 18_finalize_
      # scm_specification_v4.R's get_pair_stats(): P(arrowhead at b) -
      # P(arrowhead at a), i.e. destination-minus-source for the a -> b
      # reading of this cell.
      asymmetry = props[a, b, ">"] - props[b, a, ">"]
    )
  }))
}

stab <- make_stability_df(fci_props_ext, node_order_ext) |>
  mutate(
    orient_fill = if_else(adjacency >= .50, asymmetry, NA_real_),
    row = factor(row, levels = rev(node_order_ext)),
    col = factor(col, levels = node_order_ext)
  )

fig4_dots <- ggplot(stab, aes(x = col, y = row)) +
  # Faint cell scaffold gives the eye a matrix without creating a blocky heatmap.
  geom_tile(width = .92, height = .92, fill = "white", colour = "#EFEFEF", linewidth = .35) +
  geom_point(
    aes(size = adjacency, fill = orient_fill),
    shape = 21, stroke = .45, colour = "#7A7A7A"
  ) +
  scale_size_continuous(
    limits = c(0, 1), range = c(.8, 7.2), breaks = c(0, .25, .50, .60, .75, 1.00),
    name = "Adjacency stability"
  ) +
  scale_fill_gradient2(
    low = COL$blue_dark,
    mid = "#F5F5F5",
    high = COL$rust,
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "#E7E7E7",
    breaks = c(-1, -.5, 0, .5, 1),
    name = "Orientation asymmetry\n(destination − source arrowhead prob.)"
  ) +
  scale_x_discrete(labels = abbr_ext[node_order_ext], expand = c(0, 0), position = "top", drop = FALSE) +
  scale_y_discrete(labels = node_labels_oneline[rev(node_order_ext)], expand = c(0, 0), drop = FALSE) +
  coord_fixed() +
  labs(title = "A   FCI bootstrap stability", x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    axis.text.x = element_text(face = "bold", size = 8.4),
    axis.text.y = element_text(size = 8.6),
    plot.title  = element_text(size = 10.5, face = "bold", margin = margin(b = 4)),
    legend.position = "bottom",
    legend.box      = "vertical",
    legend.title    = element_text(size = 7.6),
    legend.text     = element_text(size = 7.2),
    plot.margin = margin(4, 6, 2, 4)
  ) +
  guides(
    size = guide_legend(order = 1, nrow = 1, override.aes = list(fill = COL$ink_mid)),
    fill = guide_colorbar(order = 2, barwidth = grid::unit(34, "mm"), barheight = grid::unit(3, "mm"))
  )

save_ms_figure(fig4_dots, "figures/fig4_stability_dotmatrix.pdf",
               c(width = 122, height = 140))
message("Saved figures/fig4_stability_dotmatrix.pdf (FCI-only dot matrix, Panel A of Figure 3)")
message("v3: now 9-node extended set incl. climate_behavior; drop=FALSE so both axes")
message("show the full diagonal starting at belief/concern.")
message("NOTE: this is Panel A only. Panel B (single-run FCI example) and the")
message("combined two-panel Figure 3 are handled separately -- see the project plan.")

