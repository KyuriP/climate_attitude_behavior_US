# =============================================================================
# 51_figure_confound_sensitivity.R  (2026-09-27, NCC Figure 5 replacement)
#
# Main-text Figure 5 for main_ncc.tex, drawn as the companion to Figure 4
# (r_patches/34_figure_bootstrap_lvida_ridge.R, panel A). Figure 4 varies the
# DIRECTION of the four weakly oriented relationships; this figure asks what
# happens when relationships that FCI most often returned as bidirected at
# both alpha levels (the eight-relationship rule, clean_pipeline/15) are no
# longer treated as directed pathways.
#
#   diamond  = working SCM (spec "baseline")
#   circles  = the eight single replacements (specs "confound_1of8_*"), each
#              directed path replaced by a residual covariance, one at a time
#   square   = all eight replaced together (spec "confound_all8")
#
# Input : pipeline_outputs/tables/confound_sensitivity_singlenode_ate.csv
#         (written by clean_pipeline/15_confound_sensitivity_diagnostic.R)
# Output: figures/fig_confound_sensitivity.pdf (140 x 110 mm, same as Fig. 4)
# No new estimation: this only re-displays existing point estimates.
# The separate present harm -> climate behavior check (r_patches/37) is not
# part of the eight-relationship rule and is deliberately not drawn here.
# =============================================================================

REPO <- "/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US"
if (dir.exists(REPO)) setwd(REPO)  # otherwise run from the repository root

suppressPackageStartupMessages(library(ggplot2))
source("r_patches/03_figure_style.R")

X_LAB <- "Effect on climate behavior (ΔY, SD)"
NODES <- c("belief_concern", "harm_present", "harm_future", "weather_risk_prep",
           "social_norms", "trust_science", "policy_support", "politics")
NOT_A_LEVER <- "politics"

d <- utils::read.csv("pipeline_outputs/tables/confound_sensitivity_singlenode_ate.csv",
                     stringsAsFactors = FALSE)
d <- d[d$node %in% NODES, ]
single_specs <- grep("^confound_1of8_", unique(d$spec), value = TRUE)
stopifnot(length(single_specs) == 8,
          all(c("baseline", "confound_all8") %in% d$spec))

base <- d[d$spec == "baseline", ]
node_order <- base$node[order(-base$ate_climate_behavior, match(base$node, NODES))]
lab_levels <- rev(node_labels_oneline[node_order])
ypos <- setNames(seq_along(lab_levels), lab_levels)

KIND <- c(base = "Working SCM",
          one  = "One of the eight relationships replaced",
          all8 = "All eight replaced")

d$kind <- ifelse(d$spec == "baseline", KIND[["base"]],
          ifelse(d$spec == "confound_all8", KIND[["all8"]],
          ifelse(d$spec %in% single_specs, KIND[["one"]], NA)))
d <- d[!is.na(d$kind), ]
d$label <- node_labels_oneline[d$node]
d$y <- unname(ypos[d$label])

# spread the eight single replacements evenly around the row (deterministic,
# so identical effects stay visible as eight separate dots)
one <- d[d$kind == KIND[["one"]], ]
one$y <- one$y + (match(one$spec, sort(single_specs)) - 4.5) * 0.045
pol <- function(x) x$node == NOT_A_LEVER

p <- ggplot() +
  geom_vline(xintercept = 0, linewidth = .45, linetype = "dashed", colour = COL$ink_light) +
  geom_point(data = one[!pol(one), ], aes(x = ate_climate_behavior, y = y),
             shape = 21, size = 1.7, stroke = 0, fill = COL$blue_dark, alpha = .25) +
  geom_point(data = one[pol(one), ], aes(x = ate_climate_behavior, y = y),
             shape = 21, size = 1.7, stroke = 0, fill = COL$ink_light, alpha = .35) +
  geom_point(data = d[d$kind == KIND[["all8"]], ],
             aes(x = ate_climate_behavior, y = y),
             shape = 22, size = 2.7, stroke = .35, fill = COL$blue_light, colour = COL$blue_dark, alpha = .7) +
  geom_point(data = d[d$kind == KIND[["base"]] & !pol(d), ],
             aes(x = ate_climate_behavior, y = y),
             shape = 23, size = 2.6, stroke = .3, fill = COL$rust, colour = COL$ink, alpha = .65) +
  geom_point(data = d[d$kind == KIND[["base"]] & pol(d), ],
             aes(x = ate_climate_behavior, y = y),
             shape = 23, size = 2.6, stroke = .3, fill = COL$ink_mid, colour = COL$ink, alpha = .65) +
  # legend keys (drawn off-panel, NA coordinates are dropped)
  geom_point(data = data.frame(kind = factor(unname(KIND), levels = unname(KIND)), x = NA_real_, y = NA_real_),
             aes(x = x, y = y, shape = kind, fill = kind, colour = kind), na.rm = TRUE) +
  scale_shape_manual(name = NULL, values = setNames(c(23, 21, 22), unname(KIND)), drop = FALSE) +
  scale_fill_manual(name = NULL, values = setNames(c(COL$rust, COL$blue_dark, COL$blue_light), unname(KIND)), drop = FALSE) +
  scale_colour_manual(name = NULL, values = setNames(c(COL$ink, COL$blue_dark, COL$blue_dark), unname(KIND)), drop = FALSE) +
  guides(shape = guide_legend(override.aes = list(size = c(3, 2.2, 3), stroke = c(.3, 0, .35), alpha = c(.65, .35, .7)))) +
  scale_y_continuous(breaks = ypos, labels = names(ypos), limits = c(.5, length(ypos) + .5),
                     expand = expansion(0)) +
  scale_x_continuous(limits = c(-.03, .25), breaks = c(0, .1, .2),
                     expand = expansion(mult = c(0, .01))) +
  labs(x = X_LAB, y = NULL) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 9.2),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = COL$grid, linewidth = .35),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7.6),
    plot.title = element_text(size = 10, hjust = 0),
    plot.subtitle = element_text(size = 8.2, colour = COL$ink_mid, hjust = 0),
    plot.margin = margin(6, 6, 4, 4)
  )

dir.create("figures", showWarnings = FALSE)
save_ms_figure(p, "figures/fig_confound_sensitivity.pdf", c(width = 140, height = 110))
ggsave("figures/fig_confound_sensitivity_preview.png", p, width = 140, height = 110,
       units = "mm", dpi = 200, bg = "white")
message("Saved: figures/fig_confound_sensitivity.pdf")

tab <- reshape(d[, c("spec", "node", "ate_climate_behavior")], idvar = "node",
               timevar = "spec", direction = "wide")
print(round(tab[match(node_order, tab$node), -1], 3))
