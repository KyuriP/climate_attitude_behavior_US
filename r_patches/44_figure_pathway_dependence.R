# =============================================================================
# 44_figure_pathway_dependence.R  (2026-09-26, NCC framing revision)
#
# Main-text figure for main_ncc.tex: how the predicted effect on climate
# behavior changes when one of two specific working-SCM paths is treated as
# non-transmitting (the directed path is replaced by a residual covariance,
# so the association is kept but an intervention cannot pass through it).
#
#   Panel A: belief/concern -> future harm replaced
#            (pipeline_outputs/tables/confound_sensitivity_singlenode_ate.csv,
#             spec "confound_1of8_belief_concern_harm_future")
#   Panel B: present harm -> climate behavior replaced
#            (pipeline_outputs/tables/presentharm_cb_sensitivity_singlenode.csv,
#             column "hp_cb_cov"; written by 37_presentharm_behavior_sensitivity_and_S8.R)
#
# No new estimation: this only re-displays existing point estimates.
# Output: figures/fig_pathway_dependence.pdf (180 x 74 mm)
# Labels (2026-09-26, readability pass): plain-language legend/subtitles; the
# alternative model still replaces the directed path with a residual covariance.
# =============================================================================

REPO <- "/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US"
if (dir.exists(REPO)) setwd(REPO)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("r_patches/03_figure_style.R")

X_LAB   <- "Effect on climate behavior (ΔY, SD)"
TARGETS <- c("belief_concern", "harm_present", "harm_future",
             "weather_risk_prep", "social_norms", "trust_science")
# policy support and political orientation have no directed path to behavior
# in any of these specifications (effect 0 throughout), so they are omitted.

conf <- read.csv("pipeline_outputs/tables/confound_sensitivity_singlenode_ate.csv")
hpcb <- read.csv("pipeline_outputs/tables/presentharm_cb_sensitivity_singlenode.csv")

get_conf <- function(spec) {
  x <- conf[conf$spec == spec, ]
  setNames(x$ate_climate_behavior, x$node)[TARGETS]
}

panel_data <- function(base, alt, focal) {
  data.frame(
    node     = TARGETS,
    label    = factor(node_labels_oneline[TARGETS], levels = rev(node_labels_oneline[TARGETS])),
    baseline = unname(base),
    alt      = unname(alt),
    focal    = TARGETS == focal
  )
}

dA <- panel_data(get_conf("baseline"),
                 get_conf("confound_1of8_belief_concern_harm_future"),
                 focal = "belief_concern")
dB <- panel_data(setNames(hpcb$baseline, hpcb$node)[TARGETS],
                 setNames(hpcb$hp_cb_cov, hpcb$node)[TARGETS],
                 focal = "harm_present")

# sanity checks against the values reported in the manuscript
stopifnot(abs(dA$baseline[1] - .201) < .0005, abs(dA$alt[1] - .091) < .0005,
          abs(dB$baseline[2] - .187) < .0005, abs(dB$alt[2] - .055) < .0005)

KIND_BASE <- "Working model"
KIND_ALT  <- "Link attributed to a shared cause"

make_panel <- function(d, title, show_y = TRUE) {
  long <- rbind(
    data.frame(label = d$label, x = d$baseline, kind = KIND_BASE),
    data.frame(label = d$label, x = d$alt,      kind = KIND_ALT)
  )
  long$kind <- factor(long$kind, levels = c(KIND_BASE, KIND_ALT))
  moved <- d[abs(d$alt - d$baseline) > .004, ]
  # stop the arrow just short of both markers so the arrowhead stays visible
  GAP <- .0062
  moved$x0 <- moved$baseline + sign(moved$alt - moved$baseline) * GAP
  moved$x1 <- moved$alt      - sign(moved$alt - moved$baseline) * GAP
  focal <- d[d$focal, ]
  focal$pct <- sprintf("%+.0f%%", 100 * (focal$alt / focal$baseline - 1))
  focal$pct <- sub("-", "−", focal$pct)
  fy <- as.numeric(focal$label)

  ggplot() +
    geom_blank(data = long, aes(x = x, y = label)) +  # establishes the discrete y scale
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = fy - .5, ymax = fy + .5,
             fill = COL$blue_light, alpha = .75) +
    geom_vline(xintercept = 0, linewidth = .45, linetype = "dashed", colour = COL$ink_light) +
    geom_segment(
      data = moved,
      aes(x = x0, xend = x1, y = label, yend = label),
      colour = COL$ink_mid, linewidth = .45,
      arrow = grid::arrow(length = unit(1.5, "mm"), type = "closed")
    ) +
    # alternative (circle) drawn first and slightly larger, baseline (diamond) on
    # top, so an unchanged effect reads as a diamond sitting inside a circle
    geom_point(
      data = long[long$kind == KIND_ALT, ],
      aes(x = x, y = label, shape = kind, fill = kind, colour = kind),
      size = 2.9, stroke = .45
    ) +
    geom_point(
      data = long[long$kind == KIND_BASE, ],
      aes(x = x, y = label, shape = kind, fill = kind, colour = kind),
      size = 2.4, stroke = .45
    ) +
    geom_text(
      data = focal,
      aes(x = (baseline + alt) / 2, y = label, label = pct),
      nudge_y = .32, size = PT(7.8), colour = COL$ink, fontface = "bold", family = FIG_FONT
    ) +
    scale_shape_manual(name = NULL, values = setNames(c(23, 21), c(KIND_BASE, KIND_ALT)), breaks = c(KIND_BASE, KIND_ALT)) +
    scale_fill_manual(name = NULL, values = setNames(c(COL$rust, COL$blue_dark), c(KIND_BASE, KIND_ALT)), breaks = c(KIND_BASE, KIND_ALT)) +
    scale_colour_manual(name = NULL, values = setNames(c(COL$ink, COL$blue_dark), c(KIND_BASE, KIND_ALT)), breaks = c(KIND_BASE, KIND_ALT)) +
    scale_x_continuous(limits = c(-.012, .235), breaks = c(0, .1, .2),
                       expand = expansion(mult = c(0, .01))) +
    labs(x = X_LAB, y = NULL, title = title,
         subtitle = "link attributed to a shared unmeasured cause") +
    theme_pub +
    theme(
      axis.text.y  = if (show_y) element_text(size = 9, colour = COL$ink) else element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      panel.grid.major.x = element_line(colour = COL$grid, linewidth = .35),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1.5, "mm"),
      plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8.2, colour = COL$ink_mid, hjust = 0),
      plot.margin = margin(6, 6, 4, 4)
    )
}

pA <- make_panel(dA, "A   Belief and concern → Future harm", show_y = TRUE)
pB <- make_panel(dB, "B   Present harm → Climate behavior", show_y = FALSE)

fig <- (pA | pB) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box.spacing = unit(0, "mm"))

dir.create("figures", showWarnings = FALSE)
save_ms_figure(fig, "figures/fig_pathway_dependence.pdf", c(width = 180, height = 74))
ggsave("figures/fig_pathway_dependence_preview.png", fig, width = 180, height = 74,
       units = "mm", dpi = 200, bg = "white")
message("Saved: figures/fig_pathway_dependence.pdf")
print(dA); print(dB)
