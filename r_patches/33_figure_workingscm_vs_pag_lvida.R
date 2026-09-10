# =============================================================================
# 33_figure_workingscm_vs_pag_lvida.R
#
# Main-text intervention/uncertainty figure -- replaces figures/fig7_uncertainty_pub.pdf
# in the main text. r_patches/10_figure7_uncertainty_pub.R itself is
# unchanged and untouched by this file.
#
# Panel A keeps the old fig7 design (baseline diamond + jittered cloud of the
# 15 alternative orientation specifications) rather than a plain forest plot
# of point + bootstrap CI only, since the jittered cloud carries information
# a plain interval doesn't.
#
# Panel B plots one point per individual compatible MAG, rather than one
# point per distinct identified value with an "N MAGs" text count, so
# multiplicity is visible directly instead of via a label. Points sit at a
# deterministic vertical stack within each node/alpha's row-band (dynamically
# spaced so any bin from 1 to 30 MAGs fits without bleeding into the next
# band or row), plus a small fixed horizontal jitter so repeated values read
# as a little cloud rather than a single stick (see the data-prep section
# below for the exact width and why it's safe). No density/violin smoothing
# is used here: the MAGs in a bin are compatible structures in an equivalence
# class, not draws from a distribution, so a smoothed curve would visually
# imply a probability that isn't there. "Unidentified" MAGs render in a fully
# separate panel (panelB_side) with no numeric x-axis at all -- a categorical
# column of x-marks -- so they can never be misread as points on the effect
# scale.
#
# ------------------------------------------------------------------ Panel A
# "A   Working SCM" / "Baseline and 15 alternative orientations" -- same
# design as r_patches/10_figure7_uncertainty_pub.R:
#   - diamond = baseline working-SCM effect (combo_0, the 16-edge model)
#   - small semi-transparent points = the other 15 acyclic orientation
#     specifications (the 4 genuinely uncertain edges' alternate directions)
#   - thin horizontal line = the participant-bootstrap 95% CI for the
#     baseline (09_intervention_bootstrap_ci.R) -- explained in the caption
#     (see bottom of file) rather than added to the legend
#   - politics kept visually distinct (own muted colour/shape), no separate
#     legend entries -- same convention as the old script
# Data: pipeline_outputs/tables/orientation_enumeration_ate_deterministic_8node.csv
# (16 specs, r_patches/32_deterministic_8node_ates.R) +
# pipeline_outputs/intervention_bootstrap_ci.csv (CI line only).
#
# ------------------------------------------------------------------ Panel B
# "B   PAG-compatible effects" / "LV-IDA across compatible MAGs" -- one
# point per INDIVIDUAL compatible mag, x = its exact effect value (only a
# small fixed jitter around it, never enough to blur two nodes' values
# together), y = a deterministic vertical stack within that node/alpha's
# row-band. Unidentified mags render in a fully separate adjoining panel
# (panelB_side), one x-mark per unidentified mag, same vertical-stack logic,
# no numeric x-axis.
# Data: pipeline_outputs/tables/extended_pag_lvida_full_vectors.csv
# (16_extended_pag_lvida.R output, with NA preserved for unidentified MAGs).
#
# Node ordering: all three panels share one order, ranked by Panel A's
# baseline (combo_0) working-SCM effect, descending, computed from the data
# (set NODE_ORDER_OVERRIDE below to force a different order instead).
#
# Not yet rendered -- exact spacing/margins between the three combined
# panels are the most likely thing left to need tuning once viewed.
#
# Output: figures/fig_workingscm_vs_pag_lvida.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
})
source("r_patches/03_figure_style.R")

NODE_ORDER_OVERRIDE <- NULL

TARGET_ORDER <- c(
  "belief_concern", "harm_present", "harm_future", "weather_risk_prep",
  "social_norms", "trust_science", "policy_support", "politics"
)
NOT_A_LEVER <- "politics"

# =============================================================================
# Panel A data: 16 orientation specifications (same source as the old fig7
# script) + the bootstrap CI for the baseline only.
# =============================================================================
full_results <- utils::read.csv(
  "pipeline_outputs/tables/orientation_enumeration_ate_deterministic_8node.csv",
  stringsAsFactors = FALSE
)
stopifnot(all(TARGET_ORDER %in% full_results$node))

boot_ci <- utils::read.csv("pipeline_outputs/intervention_bootstrap_ci.csv",
                            stringsAsFactors = FALSE)
stopifnot(all(TARGET_ORDER %in% boot_ci$target))

panelA_all <- full_results |>
  filter(node %in% TARGET_ORDER) |>
  mutate(is_baseline = scenario == "combo_0", is_politics = node == NOT_A_LEVER)

baseline_df_tmp <- filter(panelA_all, is_baseline)

# node order: baseline (combo_0) working-SCM effect, descending -- computed
# fresh from the data, not hand-typed. Strict magnitude order puts
# belief_concern (.2059) above harm_present (.1951). Set NODE_ORDER_OVERRIDE
# above to a specific 8-name vector for a different top-to-bottom order.
node_order <- if (!is.null(NODE_ORDER_OVERRIDE)) {
  stopifnot(setequal(NODE_ORDER_OVERRIDE, TARGET_ORDER))
  NODE_ORDER_OVERRIDE
} else {
  baseline_df_tmp$node[order(-baseline_df_tmp$ate_climate_behavior)]
}
node_label_levels <- rev(node_labels_oneline[node_order])
node_positions <- setNames(seq_along(node_label_levels), node_label_levels)

panelA_all <- panelA_all |>
  mutate(
    node_label = factor(node_labels_oneline[node], levels = node_label_levels),
    kind = factor(
      case_when(
        is_politics &  is_baseline ~ "Politics (baseline)",
        is_politics & !is_baseline ~ "Politics (alternative)",
        is_baseline                ~ "Baseline specification",
        TRUE                       ~ "Alternative orientation"
      ),
      levels = c("Baseline specification", "Alternative orientation",
                 "Politics (baseline)", "Politics (alternative)")
    )
  )
baseline_df    <- filter(panelA_all, is_baseline)
alternative_df <- filter(panelA_all, !is_baseline)

ci_df <- boot_ci |>
  filter(target %in% TARGET_ORDER) |>
  rename(node = target) |>
  mutate(node_label = factor(node_labels_oneline[node], levels = node_label_levels))

# =============================================================================
# Panel B data: every individual compatible mag as its own point, x anchored
# at its exact effect value plus a small fixed jitter, y a deterministic
# vertical stack.
# =============================================================================
lvida_long <- utils::read.csv(
  "pipeline_outputs/tables/extended_pag_lvida_full_vectors.csv",
  stringsAsFactors = FALSE
)
stopifnot(all(TARGET_ORDER %in% lvida_long$node))
stopifnot(all(c("zero", "identified_nonzero", "unidentified") %in% unique(lvida_long$status)))

# alpha reads back numeric (0.05/0.01) even though the csv quotes it as a
# string -- read.csv's type.convert coerces a quoted-but-numeric field.
# match on the numeric value (with a tolerance), not a hardcoded string.
alpha_label_of <- function(a) ifelse(abs(a - 0.05) < 1e-9, "alpha = .05",
                               ifelse(abs(a - 0.01) < 1e-9, "alpha = .01",
                                      paste0("alpha = ", a)))
ALPHA_LEVELS <- c("alpha = .05", "alpha = .01")

lvida_long <- lvida_long |>
  filter(node %in% TARGET_ORDER) |>
  mutate(
    node_label = factor(node_labels_oneline[node], levels = node_label_levels),
    alpha_label = factor(alpha_label_of(alpha), levels = ALPHA_LEVELS),
    # separate row-bands for the two alphas within a node's row -- same
    # fixed-offset idea as Panel A's baseline/alternative separation.
    y_row_nudge = ifelse(alpha_label == "alpha = .05", 0.20, -0.20)
  )

# A purely deterministic vertical-only stack (no horizontal jitter at all)
# made dense bins read as thin vertical sticks, losing the sense that a bin
# holds multiple compatible mags. A small horizontal jitter (+/-.008,
# tight enough that a cluster still visually sits at its true value -- .2321
# and .2545 are ~.022 apart, so +/-.008 leaves a clean gap between the two
# clusters even at their closest jittered points) plus a small vertical
# jitter turns each value into a visible little cloud instead of a line.
# Reproducible: one fixed seed, one ordered sequence of draws, not re-seeded
# per group. A half-width up to about .009 stays safely inside that same
# margin (.022 - 2*.009 = .004 remaining gap) if more spread is ever wanted.
JITTER_X_HALFWIDTH <- 0.008
JITTER_Y_HALFWIDTH <- 0.11    # stays inside the +/-0.20 row-band nudge

set.seed(20260908)
identified_df <- lvida_long |>
  filter(status %in% c("zero", "identified_nonzero")) |>
  mutate(
    y_num  = node_positions[as.character(node_label)] + y_row_nudge +
             stats::runif(n(), -JITTER_Y_HALFWIDTH, JITTER_Y_HALFWIDTH),
    x_plot = effect_ate_scale + stats::runif(n(), -JITTER_X_HALFWIDTH, JITTER_X_HALFWIDTH)
  )

unidentified_df <- lvida_long |>
  filter(status == "unidentified") |>
  mutate(
    y_num  = node_positions[as.character(node_label)] + y_row_nudge +
             stats::runif(n(), -JITTER_Y_HALFWIDTH, JITTER_Y_HALFWIDTH),
    # panelB_side's x has no numeric meaning at all (see below) -- centred
    # on 1 with a small jitter purely so the x-marks read as a cloud too,
    # not a stack. no "preserve the true value" concern here since there is
    # no true x value to preserve.
    x_plot = 1 + stats::runif(n(), -0.18, 0.18)
  )

# =============================================================================
# x-axis range for the two quantitative panels (A and B-main) -- shared, so
# the effect scale is directly comparable across them. The "Unidentified"
# column (panelB_side) is a separate categorical panel with no numeric axis
# at all, so it needs no x-range here.
# =============================================================================
all_effect_x <- c(baseline_df$ate_climate_behavior, alternative_df$ate_climate_behavior,
                   ci_df$ci_lower_95, ci_df$ci_upper_95,
                   identified_df$effect_ate_scale)
x_lo <- min(-0.02, min(all_effect_x, na.rm = TRUE) - 0.02)
x_hi <- max(all_effect_x, na.rm = TRUE) + 0.02
x_breaks <- scales::pretty_breaks(n = 5)(c(x_lo, x_hi))
X_LAB <- "Effect on climate behavior (ΔY, SD)"

Y_SCALE_B <- scale_y_continuous(
  breaks = node_positions, labels = names(node_positions),
  limits = c(0.35, length(node_order) + 0.65), expand = c(0, 0)
)

# =============================================================================
# Panel A: working-SCM baseline + 15 alternative orientation specifications,
# same visual design as r_patches/10_figure7_uncertainty_pub.R, plus the
# bootstrap CI line for the baseline drawn behind the jitter cloud (see
# caption for what the CI line is).
# =============================================================================
panelA <- ggplot() +
  geom_vline(xintercept = 0, linewidth = .45, linetype = "dashed", colour = "#9A9A9A") +
  geom_segment(
    data = ci_df,
    aes(x = ci_lower_95, xend = ci_upper_95, y = node_label, yend = node_label),
    linewidth = .8, colour = COL$ink_light
  ) +
  geom_point(
    data = alternative_df,
    aes(x = ate_climate_behavior, y = node_label, shape = kind, fill = kind, colour = kind),
    size = 2.4, stroke = 0, alpha = .3,
    # horizontal jitter width for the alternative-orientation cloud
    position = position_jitter(height = .16, width = 0.014, seed = 20260830)
  ) +
  geom_point(
    data = baseline_df,
    aes(x = ate_climate_behavior, y = node_label, shape = kind, fill = kind, colour = kind),
    size = 2.6, stroke = .5, alpha = .9
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Baseline specification" = 23, "Alternative orientation" = 21,
               "Politics (baseline)" = 23, "Politics (alternative)" = 21),
    breaks = c("Baseline specification", "Alternative orientation")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Baseline specification" = COL$rust, "Alternative orientation" = COL$blue_dark,
               "Politics (baseline)" = COL$ink_mid, "Politics (alternative)" = COL$ink_light),
    breaks = c("Baseline specification", "Alternative orientation")
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Baseline specification" = COL$ink, "Alternative orientation" = COL$blue_dark,
               "Politics (baseline)" = COL$ink, "Politics (alternative)" = COL$ink_light),
    breaks = c("Baseline specification", "Alternative orientation")
  ) +
  guides(
    shape  = guide_legend(override.aes = list(size = 3, stroke = .4, alpha = .85)),
    fill   = guide_legend(override.aes = list(size = 3, stroke = .4, alpha = .85)),
    colour = guide_legend(override.aes = list(size = 3, stroke = .4, alpha = .85))
  ) +
  scale_x_continuous(limits = c(x_lo, x_hi), breaks = x_breaks,
                      expand = expansion(mult = c(0, .01))) +
  labs(x = X_LAB, y = NULL, title = "A   Working SCM",
       subtitle = "Baseline and 15 alternative orientations") +
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
    plot.margin = margin(10, 4, 4, 4)
  )

# =============================================================================
# Panel B-main: every compatible mag as its own point, small +/-.008
# horizontal jitter around its true value (tight enough that a cluster still
# visually sits at its true value; see data-prep section above) plus a small
# vertical jitter, so multiplicity reads as a little cloud rather than a
# stick.
# =============================================================================
panelB_main <- ggplot() +
  geom_vline(xintercept = 0, linewidth = .45, linetype = "dashed", colour = "#9A9A9A") +
  geom_point(
    data = identified_df,
    aes(x = x_plot, y = y_num, colour = alpha_label, fill = alpha_label),
    shape = 21, size = 1.7, stroke = .25, alpha = .4
  ) +
  scale_colour_manual(name = NULL, values = c("alpha = .05" = COL$blue_dark, "alpha = .01" = COL$rust)) +
  scale_fill_manual(name = NULL, values = c("alpha = .05" = COL$blue_dark, "alpha = .01" = COL$rust)) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
  scale_x_continuous(limits = c(x_lo, x_hi), breaks = x_breaks,
                      expand = expansion(mult = c(0, .01))) +
  Y_SCALE_B +
  labs(x = X_LAB, y = NULL, title = "B   PAG-compatible effects",
       subtitle = "LV-IDA across compatible MAGs") +
  theme_pub +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = COL$grid, linewidth = .35),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7.6),
    plot.title = element_text(size = 10, hjust = 0),
    plot.subtitle = element_text(size = 8.2, colour = COL$ink_mid, hjust = 0),
    plot.margin = margin(10, 2, 4, 2)
  )

# =============================================================================
# Panel B-side: "Unidentified" -- a fully separate categorical strip, no
# numeric x-axis at all, so it can never be misread as a point on the effect
# scale. One x-mark per unidentified mag, same vertical-stack logic as
# panelB_main, all at a single fixed (meaningless) x position per alpha.
# =============================================================================
panelB_side <- ggplot(unidentified_df, aes(x = x_plot, y = y_num, colour = alpha_label)) +
  geom_point(shape = 4, size = 1.9, stroke = .7, alpha = .7) +
  scale_colour_manual(name = NULL, values = c("alpha = .05" = COL$blue_dark, "alpha = .01" = COL$rust),
                       guide = "none") +
  scale_x_continuous(limits = c(0.4, 1.6), breaks = NULL) +
  Y_SCALE_B +
  labs(x = NULL, y = NULL, title = "Unidentified") +
  theme_pub +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 8.4, hjust = .5, colour = COL$ink_mid, face = "italic"),
    plot.margin = margin(10, 4, 4, 0)
  )

# =============================================================================
# combine + save -- NOT using patchwork's automatic tag_levels, since the
# "A"/"B" letters are already baked into panelA/panelB_main's own titles and
# panelB_side is a sub-column of B, not its own lettered panel.
# =============================================================================
fig_combined <- panelA + panelB_main + panelB_side +
  patchwork::plot_layout(widths = c(1, 0.95, 0.22))

save_ms_figure(fig_combined, "figures/fig_workingscm_vs_pag_lvida.pdf",
               c(width = 210, height = 122))
message("Saved: figures/fig_workingscm_vs_pag_lvida.pdf")
message("Node order (top to bottom): ", paste(rev(node_order), collapse = " > "))

# =============================================================================
# Caption -- NOT inserted into main9.tex; kept here for reference only.
# main9.tex is not touched by this project until explicitly decided.
#
# \caption{Single-node intervention effects under the working SCM and under
# PAG-level structural uncertainty. \textbf{(A)} Effects of setting each
# attitude/context variable to $+0.5$ SD in the baseline working SCM
# (diamonds), with effects from the 15 alternative directional
# specifications shown as points. Horizontal intervals show 95\%
# participant-bootstrap confidence intervals for the baseline model.
# Political orientation is included as a structural comparison rather than
# as a realistic intervention target. \textbf{(B)} LV-IDA effects across
# MAGs compatible with the FCI PAGs at $\alpha=.05$ and $\alpha=.01$, shown
# on the same $+0.5$-SD intervention scale. Each point represents one
# compatible MAG. Repeated points at the same value indicate that multiple
# MAGs imply the same identified effect and should not be interpreted as
# probabilities. MAGs for which the total effect was not identified are
# shown separately in the ``Unidentified'' column.}
# =============================================================================
