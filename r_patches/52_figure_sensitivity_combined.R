# =============================================================================
# 52_figure_sensitivity_combined.R  (2026-09-27, NCC Figure 4, panels A and B)
#
# Combines the two structural-sensitivity figures of main_ncc.tex into one
# two-panel figure that shares its rows and its x-axis:
#
#   A  Uncertain directions   (was Figure 4, r_patches/34 panel A)
#        diamond = working SCM (scenario "combo_0")
#        circles = the 11 other acyclic direction choices
#        line    = 95% participant-bootstrap interval for the working SCM
#   B  Possible shared causes (was Figure 5, r_patches/51)
#        diamond = working SCM (spec "baseline")
#        circles = one of the eight relationships replaced by a covariance
#        square  = all eight replaced together
#
# Inputs : pipeline_outputs/tables/orientation_enumeration_ate_deterministic_8node.csv
#          pipeline_outputs/intervention_bootstrap_ci.csv
#          pipeline_outputs/tables/confound_sensitivity_singlenode_ate.csv
# Output : figures/fig_sensitivity_combined.pdf (180 x 100 mm) and a preview png
# Legend : one shared legend. Panel A circles are blue, panel B circles and
#          squares are dark green, so each marker type can be read from the legend.
# No new estimation: this only re-displays existing point estimates and
# intervals. Scripts 34 and 51 are left unchanged (the longer manuscript
# still uses their separate outputs).
# =============================================================================

REPO <- "/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US"
if (dir.exists(REPO)) setwd(REPO)  # otherwise run from the repository root

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("r_patches/03_figure_style.R")

X_LAB <- "Effect on climate behavior (ΔY, SD)"
NODES <- c("belief_concern", "harm_present", "harm_future", "weather_risk_prep",
           "social_norms", "trust_science", "policy_support", "politics")
NOT_A_LEVER <- "politics"
X_LIM <- c(-.03, .25)
X_BREAKS <- c(0, .1, .2)

# ---- data -------------------------------------------------------------------
ori <- utils::read.csv("pipeline_outputs/tables/orientation_enumeration_ate_deterministic_8node.csv",
                       stringsAsFactors = FALSE)
ori <- ori[ori$node %in% NODES, ]
stopifnot("combo_0" %in% ori$scenario, length(unique(ori$scenario)) == 12)

ci <- utils::read.csv("pipeline_outputs/intervention_bootstrap_ci.csv", stringsAsFactors = FALSE)
ci <- ci[ci$target %in% NODES, ]
names(ci)[names(ci) == "target"] <- "node"

con <- utils::read.csv("pipeline_outputs/tables/confound_sensitivity_singlenode_ate.csv",
                       stringsAsFactors = FALSE)
con <- con[con$node %in% NODES, ]
single_specs <- sort(grep("^confound_1of8_", unique(con$spec), value = TRUE))
stopifnot(length(single_specs) == 8, all(c("baseline", "confound_all8") %in% con$spec))

# Both analyses start from the same working SCM, so their baselines must agree.
b_ori <- ori[ori$scenario == "combo_0", ]
b_con <- con[con$spec == "baseline", ]
stopifnot(max(abs(b_ori$ate_climate_behavior -
                  b_con$ate_climate_behavior[match(b_ori$node, b_con$node)])) < 1e-6)

# Shared row order: working-SCM effect, largest at the top.
node_order <- b_con$node[order(-b_con$ate_climate_behavior, match(b_con$node, NODES))]
lab_levels <- rev(node_labels_oneline[node_order])
ypos <- setNames(seq_along(lab_levels), lab_levels)
y_of <- function(node) unname(ypos[node_labels_oneline[node]])

stopifnot(all(c(ori$ate_climate_behavior, con$ate_climate_behavior,
                ci$ci_lower_95, ci$ci_upper_95) >= X_LIM[1]),
          all(c(ori$ate_climate_behavior, con$ate_climate_behavior,
                ci$ci_lower_95, ci$ci_upper_95) <= X_LIM[2]))

ori$y <- y_of(ori$node)
con$y <- y_of(con$node)
ci$y  <- y_of(ci$node)

# Deterministic vertical spread of the circles (no random jitter).
alt_specs <- sort(setdiff(unique(ori$scenario), "combo_0"))
alt <- ori[ori$scenario != "combo_0", ]
alt$y <- alt$y + (match(alt$scenario, alt_specs) - (length(alt_specs) + 1) / 2) * 0.03
base_a <- ori[ori$scenario == "combo_0", ]

one <- con[con$spec %in% single_specs, ]
one$y <- one$y + (match(one$spec, single_specs) - (length(single_specs) + 1) / 2) * 0.045
all8 <- con[con$spec == "confound_all8", ]
base_b <- con[con$spec == "baseline", ]

pol <- function(x) x$node == NOT_A_LEVER

# ---- shared styling ---------------------------------------------------------
Y_SCALE <- scale_y_continuous(breaks = ypos, labels = names(ypos),
                              limits = c(.5, length(ypos) + .5), expand = expansion(0))
X_SCALE <- scale_x_continuous(limits = X_LIM, breaks = X_BREAKS,
                              expand = expansion(mult = c(0, .01)))
# Panel B gets its own hue so one shared legend can tell the two kinds of
# circle apart. This dark green was checked against rust and blue_dark with the
# dataviz palette validator: normal-vision dE >= 16 for every pair, colour-blind
# dE 6.1 against the rust diamonds, which is acceptable here because the two
# also differ in shape (diamond versus circle/square).
GREEN       <- "#12805F"
GREEN_DARK  <- "#0B5E45"
GREEN_LIGHT <- "#9FD0B5"

panel_theme <- theme_pub +
  theme(
    axis.text.y = element_text(size = 9.2),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = COL$grid, linewidth = .35),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7.6),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    plot.title.position = "plot",
    plot.margin = margin(4, 6, 2, 4)
  )

draw_circles <- function(df, fill, alpha = .3) list(
  geom_point(data = df[!pol(df), ], aes(x = ate_climate_behavior, y = y),
             shape = 21, size = 1.7, stroke = 0, fill = fill, alpha = alpha),
  geom_point(data = df[pol(df), ], aes(x = ate_climate_behavior, y = y),
             shape = 21, size = 1.7, stroke = 0, fill = COL$ink_light, alpha = .35)
)
draw_diamonds <- function(df) list(
  geom_point(data = df[!pol(df), ], aes(x = ate_climate_behavior, y = y),
             shape = 23, size = 2.6, stroke = .3, fill = COL$rust, colour = COL$ink, alpha = .65),
  geom_point(data = df[pol(df), ], aes(x = ate_climate_behavior, y = y),
             shape = 23, size = 2.6, stroke = .3, fill = COL$ink_mid, colour = COL$ink, alpha = .65)
)
vline0 <- geom_vline(xintercept = 0, linewidth = .45, linetype = "dashed", colour = COL$ink_light)

# One shared legend. Both panels carry the same (NA-coordinate) legend keys and
# identical scales, so patchwork's guides = "collect" merges them into one.
KIND  <- c("Working model", "Alternative directions (A)",
           "One relationship replaced (B)", "All eight replaced (B)")
K_SHP <- c(23, 21, 21, 22)
K_FIL <- c(COL$rust, COL$blue_dark, GREEN, GREEN_LIGHT)
K_COL <- c(COL$ink, COL$blue_dark, GREEN, GREEN_DARK)
K_OVR <- list(size = c(3, 2.2, 2.2, 3), stroke = c(.3, 0, 0, .35), alpha = c(.65, .45, .65, .85))
CI_LAB <- "95% bootstrap interval (A)"
shared_legend <- list(
  geom_point(data = data.frame(kind = factor(KIND, levels = KIND), x = NA_real_, y = NA_real_),
             aes(x = x, y = y, shape = kind, fill = kind, colour = kind), na.rm = TRUE),
  geom_segment(data = data.frame(x = NA_real_, y = NA_real_),
               aes(x = x, xend = x, y = y, yend = y, linetype = CI_LAB), na.rm = TRUE),
  scale_shape_manual(name = NULL, values = setNames(K_SHP, KIND), drop = FALSE),
  scale_fill_manual(name = NULL, values = setNames(K_FIL, KIND), drop = FALSE),
  scale_colour_manual(name = NULL, values = setNames(K_COL, KIND), drop = FALSE),
  scale_linetype_manual(name = NULL, values = setNames("solid", CI_LAB)),
  guides(shape    = guide_legend(order = 1, nrow = 2, byrow = TRUE, override.aes = K_OVR),
         fill     = guide_legend(order = 1, nrow = 2, byrow = TRUE, override.aes = K_OVR),
         colour   = guide_legend(order = 1, nrow = 2, byrow = TRUE, override.aes = K_OVR),
         linetype = guide_legend(order = 2, override.aes = list(colour = COL$ink_light, linewidth = .8)))
)

# ---- panel A: uncertain directions -------------------------------------------
pA <- ggplot() + vline0 +
  geom_segment(data = ci[ci$ci_upper_95 > ci$ci_lower_95, ],
               aes(x = ci_lower_95, xend = ci_upper_95, y = y, yend = y),
               linewidth = .8, colour = COL$ink_light) +
  draw_circles(alt, COL$blue_dark) + draw_diamonds(base_a) +
  shared_legend + Y_SCALE + X_SCALE +
  labs(x = X_LAB, y = NULL, title = "A   Uncertain directions") +
  panel_theme

# ---- panel B: possible shared causes -----------------------------------------
pB <- ggplot() + vline0 +
  draw_circles(one, GREEN, alpha = .5) +
  geom_point(data = all8[!pol(all8), ], aes(x = ate_climate_behavior, y = y),
             shape = 22, size = 2.7, stroke = .35, fill = GREEN_LIGHT, colour = GREEN_DARK, alpha = .85) +
  geom_point(data = all8[pol(all8), ], aes(x = ate_climate_behavior, y = y),
             shape = 22, size = 2.7, stroke = .35, fill = "#EDEDED", colour = COL$ink_light, alpha = .8) +
  draw_diamonds(base_b) +
  shared_legend + Y_SCALE + X_SCALE +
  labs(x = X_LAB, y = NULL, title = "B   Possible shared causes") +
  panel_theme + theme(axis.text.y = element_blank())

fig <- (pA + pB) +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.spacing.x = unit(1.5, "mm"), legend.box.spacing = unit(2.5, "mm"),
        legend.key.height = unit(3.8, "mm"))

dir.create("figures", showWarnings = FALSE)
save_ms_figure(fig, "figures/fig_sensitivity_combined.pdf", c(width = 180, height = 100))
ggsave("figures/fig_sensitivity_combined_preview.png", fig, width = 180, height = 100,
       units = "mm", dpi = 200, bg = "white")
message("Saved: figures/fig_sensitivity_combined.pdf")
