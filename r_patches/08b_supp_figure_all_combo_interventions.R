# =============================================================================
# 08b_supp_figure_all_combo_interventions.R
#
# The full set of combined-intervention results is not too much for a
# chart -- just too much for the compact main-text figure it would have sat
# alongside. This is the full picture, for the appendix.
#
# Restricted to PAIRS + TRIPLES only (35 rows), not all 41: the 6 single-node
# targets live in the merged main-text figure (10_figure7_uncertainty_pub.R,
# which absorbed the old bar chart's value labels), so repeating them here
# would just duplicate that figure. This appendix figure is exactly "the
# other 35" the main text doesn't show. tables/supp_table_all_combo_
# interventions.csv (the CSV table) still lists all 41 including singles,
# for anyone cross-referencing by rank -- only this FIGURE is narrowed.
#
# Each shown with its baseline (point-estimate SCM) effect AND its range
# across the same 16 structural-uncertainty specifications used in Figure 7
# (harm_future<->harm_present is the 4th flippable edge, see
# r_patches/02_full_orientation_enumeration_v4.R) -- same underlying data as
# the old Figure 6 barplot's "best pair"/"best triple" bars, just every
# combo instead of the top 2. 35 rows still needs real vertical space
# (230mm here) to stay legible; text sizes below are absolute pt, not
# scaled by figure size, so the extra height buys row-spacing, not
# distorted labels.
#
# Data source: tables/orientation_uncertainty_band_combos_full_deterministic.csv,
# written by r_patches/31_deterministic_combo_ates.R (16 structural
# specifications, same set used everywhere else in the intervention
# analysis). This uses exact mean propagation rather than the
# 20,000-draw Monte Carlo version (tables/orientation_uncertainty_band_
# combos_full.csv, from 02_full_orientation_enumeration_v4.R) -- same
# reasoning as script 30's single-node fix, just extended to all 41
# targets: the SCM is linear-recursive, so mean propagation gives the exact
# answer with no simulation noise at all. Every target's baseline/min/max
# shifts by at most ~.01 SD versus the Monte Carlo version (checked in
# tables/orientation_combo_ate_deterministic_vs_montecarlo.csv), and the top
# pair and top triple by baseline effect are unchanged. Run
# r_patches/31_deterministic_combo_ates.R first if this file doesn't exist
# yet or looks stale -- no other script needs to be re-run for this one.
#
# Design: one row per target, sorted by baseline effect (matches
# tables/supp_table_all_combo_interventions.csv's ranking so the reader can
# cross-reference the two). Filled diamond = baseline/point-estimate SCM
# effect (same structure used everywhere else in the paper). Horizontal
# segment = range across all 16 structural specifications (min to max) --
# same range-across-specifications idea as Figure 7, just without the
# individual jittered per-specification dots Figure 7 shows for its 6 rows;
# at 41 rows that would be 328 dots and too dense to read. Colour encodes
# target size (single/pair/triple), matching Figure 6's colour scheme
# exactly so the two figures read as one system.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
source("r_patches/03_figure_style.R")

COMBO_BAND_CSV <- "pipeline_outputs/tables/orientation_uncertainty_band_combos_full_deterministic.csv"  # this writes under pipeline_outputs/tables now, not bare tables/
if (!file.exists(COMBO_BAND_CSV)) {
  stop("'", COMBO_BAND_CSV, "' not found -- run ",
       "08_intervention_ates_combo.R first.")
}
combo_band <- utils::read.csv(COMBO_BAND_CSV, stringsAsFactors = FALSE)
stopifnot(nrow(combo_band) == 63)  # 7 singles + 21 pairs + 35 triples (was 41 over 6 nodes)

# Drop the 6 singles for THIS FIGURE ONLY -- they're already shown, with
# their own value labels, in the merged 10_figure7_uncertainty_pub.R.
combo_band <- combo_band |> dplyr::filter(target_size %in% c(2, 3))
stopifnot(nrow(combo_band) == 56)  # 21 pairs + 35 triples (was 35 over 6 nodes)

friendly_combo_label <- function(target_label) {
  nodes <- strsplit(target_label, "\\+", fixed = FALSE)[[1]]
  paste(node_labels_oneline[nodes], collapse = " + ")
}

plot_df <- combo_band |>
  mutate(
    label = vapply(target_label, friendly_combo_label, character(1)),
    size_group = factor(
      c("1" = "Single node", "2" = "Pair", "3" = "Triple")[as.character(target_size)],
      levels = c("Single node", "Pair", "Triple")
    )
  ) |>
  arrange(ate_baseline) |>
  mutate(label = factor(label, levels = label))  # sorted top-to-bottom, largest at top

# Same 3 colours as Figure 6 (Single node / Best pair / Best triple), so the
# two figures share one visual vocabulary.
GROUP_COLOR <- c(
  "Single node" = COL$blue_dark,
  "Pair"        = COL$rust,
  "Triple"      = "#5B7553"
)

fig_supp_combo <- ggplot(plot_df, aes(y = label)) +
  geom_vline(xintercept = 0, colour = COL$ink_light, linewidth = .35) +
  geom_segment(
    aes(x = ate_min, xend = ate_max, yend = label, colour = size_group),
    linewidth = 1.1, lineend = "round", alpha = .75
  ) +
  geom_point(
    aes(x = ate_baseline, colour = size_group),
    size = 2.6, shape = 18
  ) +
  scale_colour_manual(values = GROUP_COLOR, name = NULL) +
  scale_x_continuous(
    name = "Model-implied change in climate behavior (SD units)",
    expand = expansion(mult = c(.02, .05))
  ) +
  scale_y_discrete(name = NULL) +
  # The explanatory clause used to be crammed into the axis title itself as
  # one long string -- rendering confirmed it ran right off the page edge
  # ("range across a" then cut off), since axis titles don't wrap. Moved to
  # a caption instead, which DOES wrap.
  labs(caption = paste0(
    "Diamond = baseline structure. Segment = range across all 16 structural\n",
    "specifications (see Figure 4)."
  )) +
  theme_pub +
  theme(
    text            = element_text(size = 9.5),
    axis.text.y     = element_text(size = 8.0, colour = COL$ink),
    axis.text.x     = element_text(size = 8.4),
    axis.title.x    = element_text(size = 8.8),
    axis.line.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    panel.grid.major.x = element_line(colour = COL$grid, linewidth = .3),
    panel.grid.major.y = element_line(colour = COL$grid, linewidth = .15),
    legend.position = "top",
    legend.text     = element_text(size = 8.6),
    plot.caption    = element_text(size = 7.6, colour = COL$ink_mid, hjust = 0, margin = margin(t = 6)),
    plot.margin     = margin(6, 14, 4, 6)
  )

# Height is 275mm (up from 230mm) -- 56 rows (up from 35) need more
# vertical space to stay legible at the same per-row spacing.
save_ms_figure(fig_supp_combo, "figures/figS_all_combo_interventions.pdf",
               c(width = 172, height = 275))
message("Saved supplementary figure: figures/figS_all_combo_interventions.pdf (",
        nrow(plot_df), " rows [pairs+triples only], sorted by baseline effect size)")
