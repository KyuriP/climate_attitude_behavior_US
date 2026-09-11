# =============================================================================
# 04_figure1_heatmap.R
#
# Adds each variable's marginal distribution on the diagonal, turning the
# lower-triangle-only correlation matrix into a standard lower-triangle +
# diagonal correlogram.
#
# - The plotting grid uses numeric axes (1..8 per node) with custom
#   breaks/labels rather than discrete factor axes, so that a density curve
#   can be drawn as a small inset inside each diagonal cell -- discrete
#   ggplot axes can't host a freely-positioned continuous curve the way
#   numeric axes can. Visual node order and reading direction match the
#   off-diagonal layout.
# - Off-diagonal (lower-triangle) cells: tile fill on the sequential 0..1
#   Pearson-r scale, numeric label in each cell.
# - Diagonal cells show that node's own univariate density (each panel
#   independently rescaled to fill its cell -- the standard convention for a
#   correlogram diagonal, e.g. GGally::ggpairs), on a neutral gray fill that
#   deliberately does NOT participate in the Pearson-r colour scale/legend
#   (a density isn't a correlation, so it shouldn't borrow that scale).
# - Because the diagonal is populated, both axes show all 8 node names
#   (Belief through Social norms) instead of 7 each -- this is expected and
#   correct; there's no "Belief" row or "Social norms" column in a
#   lower-triangle-only layout because neither ever appears in such a cell.
#
# Needs df_main in scope. Tested against a structurally identical stand-in
# dataset (same 8 columns, same N) rather than the live df_main, which isn't
# available outside a live analysis session -- confirms the geometry, math,
# and axis logic are correct; only the real numbers will differ from that
# test render, which is expected and fine.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(scales)
})
source("r_patches/03_figure_style.R")

stopifnot(exists("df_main"))

NODES <- c(
  "belief_concern", "harm_present", "harm_future",
  "trust_science", "policy_support", "politics",
  "weather_risk_prep", "social_norms"
)
stopifnot(all(NODES %in% names(df_main)))
N <- length(NODES)

cor_mat <- df_main |>
  dplyr::select(dplyr::all_of(NODES)) |>
  stats::cor(use = "pairwise.complete.obs")

# --- Off-diagonal (lower triangle): unchanged from before ------------------
cor_df <- as.data.frame(cor_mat) |>
  tibble::rownames_to_column("row") |>
  tidyr::pivot_longer(-row, names_to = "col", values_to = "r") |>
  dplyr::mutate(
    row_i = match(row, NODES),
    col_i = match(col, NODES)
  ) |>
  dplyr::filter(row_i > col_i) |>
  dplyr::mutate(
    x  = col_i,
    y  = N + 1 - row_i,
    text_col = if_else(r >= .55, "white", COL$ink)
  )

# --- Diagonal cell backgrounds (neutral, not on the r scale) ---------------
diag_bg <- tibble::tibble(
  node = NODES,
  i    = seq_len(N)
) |>
  dplyr::mutate(x = i, y = N + 1 - i)

# --- Diagonal marginal densities, each independently rescaled to fill its
#     own cell (standard correlogram-diagonal convention) -------------------
CELL_HALF <- 0.44
diag_density <- purrr::map_dfr(seq_len(N), function(i) {
  v <- df_main[[NODES[i]]]
  v <- v[is.finite(v)]
  d <- stats::density(v, n = 128, from = min(v), to = max(v))
  x0 <- i
  y0 <- N + 1 - i
  tibble::tibble(
    node = NODES[i],
    x    = scales::rescale(d$x, to = c(x0 - CELL_HALF, x0 + CELL_HALF)),
    ymin = y0 - CELL_HALF,
    ymax = scales::rescale(d$y, to = c(y0 - CELL_HALF, y0 + CELL_HALF), from = c(0, max(d$y)))
  )
})

fig1 <- ggplot() +
  # off-diagonal tiles (Pearson r)
  geom_tile(
    data = cor_df, aes(x = x, y = y, fill = r),
    colour = "white", linewidth = 1.0, width = .98, height = .98
  ) +
  geom_text(
    data = cor_df,
    aes(x = x, y = y, label = sprintf("%.2f", r), colour = text_col),
    family = FIG_FONT, size = PT(8.4), show.legend = FALSE
  ) +
  # diagonal cell background (neutral gray, outside the r scale)
  geom_tile(
    data = diag_bg, aes(x = x, y = y),
    fill = COL$panel, colour = "white", linewidth = 1.0, width = .98, height = .98
  ) +
  # diagonal marginal-distribution insets
  geom_ribbon(
    data = diag_density,
    aes(x = x, ymin = ymin, ymax = ymax, group = node),
    fill = COL$blue_light, colour = COL$blue_dark, linewidth = .32
  ) +
  scale_colour_identity() +
  scale_fill_gradient(
    low = "#F2F2F2",
    high = COL$blue_dark,
    limits = c(0, 1),
    breaks = c(0, .25, .5, .75, 1),
    name = "Pearson r",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = .5,
      barwidth = grid::unit(40, "mm"),
      barheight = grid::unit(3.6, "mm")
    )
  ) +
  scale_x_continuous(
    breaks = seq_len(N), labels = node_labels_short[NODES],
    expand = c(0, 0), limits = c(0.5, N + 0.5)
  ) +
  scale_y_continuous(
    breaks = seq_len(N), labels = node_labels_oneline[rev(NODES)],
    expand = c(0, 0), limits = c(0.5, N + 0.5)
  ) +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_pub +
  theme(
    text            = element_text(size = 10.8),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, size = 9.2),
    axis.text.y = element_text(size = 9.4),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 8.8),
    legend.title = element_text(size = 9.4),
    plot.margin = margin(4, 8, 2, 4)
  )

save_ms_figure(fig1, "figures/fig1_heatmap_construct.pdf",
               c(width = 205, height = 138))
message("Saved Figure 1: figures/fig1_heatmap_construct.pdf (now with diagonal marginal distributions)")
