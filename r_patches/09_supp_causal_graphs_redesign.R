# =============================================================================
# 09_supp_causal_graphs_redesign.R
# Supplementary Figure 8: four single-run causal-discovery graphs with exact
# endpoint semantics and one shared layout.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(pcalg)
})
source("r_patches/03_figure_style.R")

stopifnot(exists("suffStat_gauss"), exists("node_order_cd"))

if (!exists("fci_05")) {
  fci_05 <- pcalg::fci(
    suffStat = suffStat_gauss, indepTest = pcalg::gaussCItest,
    alpha = .05, labels = node_order_cd, verbose = FALSE
  )
}
if (!exists("fci_01")) {
  fci_01 <- pcalg::fci(
    suffStat = suffStat_gauss, indepTest = pcalg::gaussCItest,
    alpha = .01, labels = node_order_cd, verbose = FALSE
  )
}
if (!exists("pc_05")) {
  pc_05 <- pcalg::pc(
    suffStat = suffStat_gauss, indepTest = pcalg::gaussCItest,
    alpha = .05, labels = node_order_cd, skel.method = "stable", verbose = FALSE
  )
}
if (!exists("pc_01")) {
  pc_01 <- pcalg::pc(
    suffStat = suffStat_gauss, indepTest = pcalg::gaussCItest,
    alpha = .01, labels = node_order_cd, skel.method = "stable", verbose = FALSE
  )
}

fci_amat_05 <- fci_05@amat
fci_amat_01 <- fci_01@amat
pc_amat_05 <- as(pc_05@graph, "matrix")
pc_amat_01 <- as(pc_01@graph, "matrix")

stopifnot(
  all(fci_amat_05 %in% 0:3), all(fci_amat_01 %in% 0:3),
  all(pc_amat_05 %in% 0:1), all(pc_amat_01 %in% 0:1)
)

node_layout <- fixed_layout_coords |>
  filter(name %in% node_order_cd) |>
  mutate(label = node_labels[name])
node_layout <- node_layout[match(node_order_cd, node_layout$name), ]
stopifnot(!anyNA(node_layout$x), !anyNA(node_layout$y))

# FCI pcalg endpoint coding: 1 circle, 2 arrowhead, 3 tail.
decode_fci_mark <- function(x) {
  case_when(x == 1 ~ "circle", x == 2 ~ "arrow", x == 3 ~ "tail", TRUE ~ "none")
}

fci_to_edges <- function(amat) {
  nodes <- rownames(amat); out <- list(); k <- 1L
  for (i in seq_len(length(nodes) - 1L)) {
    for (j in (i + 1L):length(nodes)) {
      if (amat[i, j] == 0 && amat[j, i] == 0) next
      out[[k]] <- tibble(
        from = nodes[i], to = nodes[j],
        mark_from = decode_fci_mark(amat[j, i]),
        mark_to   = decode_fci_mark(amat[i, j])
      )
      k <- k + 1L
    }
  }
  bind_rows(out)
}

pc_to_edges <- function(amat) {
  nodes <- rownames(amat); out <- list(); k <- 1L
  for (i in seq_len(length(nodes) - 1L)) {
    for (j in (i + 1L):length(nodes)) {
      ij <- amat[i, j]; ji <- amat[j, i]
      if (ij == 0 && ji == 0) next
      if (ij == 1 && ji == 1) {
        mf <- "tail"; mt <- "tail"
      } else if (ij == 1 && ji == 0) {
        mf <- "tail"; mt <- "arrow"
      } else if (ij == 0 && ji == 1) {
        mf <- "arrow"; mt <- "tail"
      } else stop("Unexpected PC matrix coding")
      out[[k]] <- tibble(from = nodes[i], to = nodes[j], mark_from = mf, mark_to = mt)
      k <- k + 1L
    }
  }
  bind_rows(out)
}

prepare_geometry <- function(edge_df, trim = .18, arrow_fraction = .42) {
  edge_df |>
    left_join(node_layout |> select(from = name, x0 = x, y0 = y), by = "from") |>
    left_join(node_layout |> select(to = name, x1c = x, y1c = y), by = "to") |>
    rowwise() |>
    mutate(
      dx = x1c - x0, dy = y1c - y0,
      d = sqrt(dx^2 + dy^2), ux = dx/d, uy = dy/d,
      x1 = x0 + trim * ux, y1 = y0 + trim * uy,
      x2 = x1c - trim * ux, y2 = y1c - trim * uy,
      tx = x2 - arrow_fraction * (x2 - x1),
      ty = y2 - arrow_fraction * (y2 - y1),
      fx = x1 + arrow_fraction * (x2 - x1),
      fy = y1 + arrow_fraction * (y2 - y1)
    ) |>
    ungroup()
}

plot_single_run <- function(edge_df, title) {
  ee <- prepare_geometry(edge_df)

  p <- ggplot() +
    geom_segment(
      data = ee,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      linewidth = .62, colour = "#555555", lineend = "round"
    )

  if (any(ee$mark_to == "arrow")) {
    p <- p + geom_segment(
      data = filter(ee, mark_to == "arrow"),
      aes(x = tx, y = ty, xend = x2, yend = y2),
      linewidth = 1.15, colour = COL$ink, lineend = "round", linejoin = "mitre",
      arrow = grid::arrow(type = "closed", length = grid::unit(3.6, "mm"), angle = 22)
    )
  }
  if (any(ee$mark_from == "arrow")) {
    p <- p + geom_segment(
      data = filter(ee, mark_from == "arrow"),
      aes(x = fx, y = fy, xend = x1, yend = y1),
      linewidth = 1.15, colour = COL$ink, lineend = "round", linejoin = "mitre",
      arrow = grid::arrow(type = "closed", length = grid::unit(3.6, "mm"), angle = 22)
    )
  }
  if (any(ee$mark_from == "circle")) {
    p <- p + geom_point(
      data = filter(ee, mark_from == "circle"),
      aes(x = x1, y = y1), shape = 21, size = 2.8, stroke = .85,
      fill = "white", colour = COL$ink
    )
  }
  if (any(ee$mark_to == "circle")) {
    p <- p + geom_point(
      data = filter(ee, mark_to == "circle"),
      aes(x = x2, y = y2), shape = 21, size = 2.8, stroke = .85,
      fill = "white", colour = COL$ink
    )
  }

  p +
    geom_point(
      data = node_layout, aes(x, y), shape = 21, size = 15.5, stroke = .8,
      fill = "#F2F2F2", colour = COL$ink_mid
    ) +
    geom_text(
      data = node_layout, aes(x, y, label = label),
      family = FIG_FONT, size = PT(6.3), lineheight = .86, colour = COL$ink
    ) +
    labs(title = title) +
    coord_equal(
      xlim = range(fixed_layout_coords$x) + c(-.55, .55),
      ylim = c(-2.75, 2.05),
      clip = "off"
    ) +
    theme_void(base_family = FIG_FONT) +
    theme(
      plot.title = element_text(size = 8.8, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.margin = margin(2, 2, 2, 2)
    )
}

fig_s1 <- plot_single_run(fci_to_edges(fci_amat_05), "A  FCI, α = .05")
fig_s2 <- plot_single_run(pc_to_edges(pc_amat_05),  "B  PC-stable, α = .05")
fig_s3 <- plot_single_run(fci_to_edges(fci_amat_01), "C  FCI, α = .01")
fig_s4 <- plot_single_run(pc_to_edges(pc_amat_01),  "D  PC-stable, α = .01")

fig_supp_causal <- (fig_s1 | fig_s2) / (fig_s3 | fig_s4)

save_ms_figure(
  fig_supp_causal,
  "figures/fig_supp_causal_redesign.pdf",
  FIG_DIMS_MM$supp
)
message("Saved Supplementary Figure 8: figures/fig_supp_causal_redesign.pdf")
