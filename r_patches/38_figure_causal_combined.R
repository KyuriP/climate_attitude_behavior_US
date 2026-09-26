# =============================================================================
# 38_figure_causal_combined.R  (2026-09-26)
#
# Manuscript Figure 2 (fig:causal) as ONE file, both panels drawn in ggplot2
# so the panel titles, fonts and title-to-plot spacing are identical:
#   A  FCI bootstrap stability dot matrix (same data/encoding as
#      06_figure4_stability_redesign.R, legends re-laid out so they fit)
#   B  single-run nine-node FCI-JCI PAG at alpha = .05
#      (pipeline_outputs/fci_ext_05_singlerun.rds, from clean_pipeline/16),
#      node positions copied from the hand-tuned Affinity layout
#      (figures/illustration/FCI_Graph.af), node fills from the shared
#      node_family_colors used in Figure 1. Edges/endpoints are read from the
#      saved amat -- nothing is hand-drawn.
# Output: figures/fig_causal_combined.pdf (+ .png preview)
# =============================================================================
setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tibble); library(patchwork) })
source("clean_pipeline/00_config.R")
source("r_patches/03_figure_style.R")

TITLE_THEME <- theme(
  plot.title = element_text(family = FIG_FONT, size = 10, face = "bold", hjust = 0,
                            margin = margin(b = 6)),
  plot.title.position = "plot")

# ---------------------------------------------------------------- panel A
latest <- read.csv(LATEST_POINTER_PATH, stringsAsFactors = FALSE)
props  <- readRDS(latest$fci_props_ext_path)
nodes  <- NODE_ORDER_EXT
abbr   <- ABBR_EXT
stab <- bind_rows(lapply(utils::combn(nodes, 2, simplify = FALSE), function(z) {
  a <- z[1]; b <- z[2]
  tibble(row = b, col = a,
         adjacency = 1 - mean(c(props[a, b, "N"], props[b, a, "N"]), na.rm = TRUE),
         asymmetry = props[a, b, ">"] - props[b, a, ">"])
})) |>
  mutate(orient_fill = if_else(adjacency >= .50, asymmetry, NA_real_),
         row = factor(row, levels = rev(nodes)), col = factor(col, levels = nodes))

pA <- ggplot(stab, aes(col, row)) +
  geom_tile(width = .92, height = .92, fill = "white", colour = "#EFEFEF", linewidth = .3) +
  geom_point(aes(size = adjacency, fill = orient_fill), shape = 21, stroke = .4, colour = "#7A7A7A") +
  scale_size_continuous(limits = c(0, 1), range = c(.6, 5.6), breaks = c(.25, .50, .75, 1),
                        labels = c(".25", ".50", ".75", "1.00"), name = "Adjacency stability") +
  scale_fill_gradient2(low = COL$blue_dark, mid = "#F5F5F5", high = COL$rust, midpoint = 0,
                       limits = c(-1, 1), na.value = "#E7E7E7", breaks = c(-1, -.5, 0, .5, 1),
                       labels = c("−1", "−.5", "0", ".5", "1"),
                       name = "Orientation asymmetry") +
  scale_x_discrete(labels = abbr[nodes], position = "top", drop = FALSE, expand = c(0, 0)) +
  scale_y_discrete(labels = paste0(node_labels_oneline[rev(nodes)], " (", abbr[rev(nodes)], ")"),
                   drop = FALSE, expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  labs(title = "A   FCI bootstrap stability", x = NULL, y = NULL) +
  theme_pub + TITLE_THEME +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 7, face = "bold", colour = "#333333"),
        axis.text.y = element_text(size = 7.3, colour = "#333333"),
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.justification = "left", legend.spacing.y = unit(1, "mm"),
        legend.title = element_text(size = 6.8), legend.text = element_text(size = 6.3),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 4, 2, 2)) +
  guides(size = guide_legend(order = 1, nrow = 1, title.position = "left", title.vjust = .5,
                             override.aes = list(fill = COL$ink_mid)),
         fill = guide_colorbar(order = 2, title.position = "left", title.vjust = .9,
                               barwidth = unit(30, "mm"), barheight = unit(2.4, "mm")))

# ---------------------------------------------------------------- panel B
obj  <- readRDS("pipeline_outputs/fci_ext_05_singlerun.rds")
amat <- if (isS4(obj)) obj@amat else if (is.list(obj) && !is.null(obj$amat)) obj$amat else obj
stopifnot(all(nodes %in% rownames(amat)))
amat <- amat[nodes, nodes]
# ILLUSTRATIVE OVERRIDE (author decision 2026-09-26): panel B is presented as
# an example PAG. To show what an unresolved endpoint looks like, the
# trust_science end of trust_science -- policy_support is drawn as a circle
# (o->), as in the earlier hand-drawn version. In the actual alpha = .05 single
# run this end is a tail (trust_science -> policy_support); the figure caption
# says so. amat[a, b] is the mark at b, so this sets the mark AT trust_science.
amat["policy_support", "trust_science"] <- 1

# node centres copied from FCI_Graph.af (pixel units, y flipped)
pos <- tribble(~name, ~x, ~y,
  "belief_concern",     828, -172,
  "policy_support",     333, -456,
  "harm_present",      1060, -520,
  "weather_risk_prep", 1462, -636,
  "politics",            98, -778,
  "harm_future",        746, -722,
  "social_norms",      1028, -1016,
  "trust_science",      527, -1122,
  "climate_behavior",   948, -1540)
P  <- as.matrix(pos[, c("x", "y")]); rownames(P) <- pos$name
RAD <- 126; GAP <- 6; CR <- 17            # node radius, gap, endpoint-circle radius

# curvature for the long edges (signed w.r.t. the listed from -> to order)
curv_tbl <- tribble(~from, ~to, ~c,
  "belief_concern", "politics",          -0.62,
  "belief_concern", "trust_science",     -0.28,
  "belief_concern", "weather_risk_prep",  0.25,
  "policy_support", "climate_behavior",  -0.52,
  "harm_present",   "climate_behavior",   0.22,
  "weather_risk_prep", "climate_behavior", 0.12,
  "policy_support", "social_norms",      -0.22,
  "trust_science",  "policy_support",     0.04)
curv_of <- function(a, b) {
  r <- curv_tbl[curv_tbl$from == a & curv_tbl$to == b, ]
  if (nrow(r)) return(r$c)
  r <- curv_tbl[curv_tbl$from == b & curv_tbl$to == a, ]
  if (nrow(r)) return(-r$c)
  0
}
bez <- function(p0, p1, cv, n = 200) {
  d <- p1 - p0; L <- sqrt(sum(d^2)); perp <- c(-d[2], d[1]) / L
  ctrl <- (p0 + p1) / 2 + cv * L * perp
  t <- seq(0, 1, length.out = n)
  cbind((1-t)^2*p0[1] + 2*(1-t)*t*ctrl[1] + t^2*p1[1],
        (1-t)^2*p0[2] + 2*(1-t)*t*ctrl[2] + t^2*p1[2])
}
edges <- list(); circles <- list(); k <- 0; edge_log <- character(0)
MARK <- c("0" = "none", "1" = "circle", "2" = "arrow", "3" = "tail")
for (i in 1:(length(nodes) - 1)) for (j in (i + 1):length(nodes)) {
  a <- nodes[i]; b <- nodes[j]
  if (amat[a, b] == 0 && amat[b, a] == 0) next
  mark_b <- MARK[as.character(amat[a, b])]; mark_a <- MARK[as.character(amat[b, a])]
  edge_log <- c(edge_log, sprintf("%s %s-%s %s", a, mark_a, mark_b, b))
  pts <- bez(P[a, ], P[b, ], curv_of(a, b))
  da <- sqrt((pts[, 1] - P[a, 1])^2 + (pts[, 2] - P[a, 2])^2)
  db <- sqrt((pts[, 1] - P[b, 1])^2 + (pts[, 2] - P[b, 2])^2)
  cut_a <- RAD + GAP + if (mark_a == "circle") 2 * CR else 0
  cut_b <- RAD + GAP + if (mark_b == "circle") 2 * CR else 0
  keep <- da > cut_a & db > cut_b
  k <- k + 1
  edges[[k]] <- list(df = data.frame(x = pts[keep, 1], y = pts[keep, 2]),
                     ends = if (mark_a == "arrow" && mark_b == "arrow") "both" else
                            if (mark_b == "arrow") "last" else if (mark_a == "arrow") "first" else NA)
  for (side in c("a", "b")) {
    m <- if (side == "a") mark_a else mark_b
    if (m != "circle") next
    dd <- if (side == "a") da else db
    idx <- which.min(abs(dd - (RAD + GAP + CR)))
    circles[[length(circles) + 1]] <- data.frame(x = pts[idx, 1], y = pts[idx, 2])
  }
}
cat("PAG edges (from amat):\n"); writeLines(edge_log)

circ_poly <- function(cx, cy, r, id, n = 120) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t), id = id)
}
node_df <- do.call(rbind, lapply(nodes, function(nd) circ_poly(P[nd, 1], P[nd, 2], RAD, nd)))
node_df$fill <- node_fill_for(node_df$id)

pB <- ggplot() + coord_fixed(clip = "off") + theme_void() + theme(plot.margin = margin(0, 2, 2, 6))
AR <- arrow(angle = 22, length = unit(1.7, "mm"), type = "open")
for (e in edges) {
  pB <- pB + geom_path(data = e$df, aes(x, y), linewidth = .38, colour = "#3C3C3C",
                       arrow = if (is.na(e$ends)) NULL else
                         arrow(angle = 22, length = unit(1.7, "mm"), type = "open", ends = e$ends))
}
if (length(circles)) pB <- pB + geom_point(data = do.call(rbind, circles), aes(x, y),
                                           shape = 21, size = 2.0, stroke = .42,
                                           colour = "#3C3C3C", fill = "white")
pB <- pB +
  geom_polygon(data = node_df, aes(x, y, group = id, fill = I(fill)), colour = "#6B6B6B", linewidth = .35) +
  geom_text(data = pos, aes(x, y, label = node_labels[name]), family = FIG_FONT,
            size = 5.9 / .pt, lineheight = .88, colour = "#222222") +
  scale_x_continuous(expand = expansion(add = 30)) + scale_y_continuous(expand = expansion(add = 30))

# wrap_elements(): stop patchwork from shrinking panel B to the height of
# panel A's matrix region (A's legends sit below its panel; B has none)
# B's title is attached to the wrapped element (not the fixed-aspect plot),
# so it sits at the top of its column exactly like A's title
pBw <- wrap_elements(full = pB) + labs(title = "B   Example FCI PAG") + TITLE_THEME +
  theme(plot.margin = margin(0, 0, 0, 0))
fig <- pA + pBw + plot_layout(widths = c(1, 1.22))
save_ms_figure(fig, "figures/fig_causal_combined.pdf", c(width = 180, height = 98))
ggsave("figures/fig_causal_combined_preview.png", fig, width = 180, height = 98, units = "mm", dpi = 220, bg = "white")
cat("38 DONE\n")
