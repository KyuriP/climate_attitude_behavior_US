# =============================================================================
# 34_figure_bootstrap_lvida_ridge.R
#
# EXPERIMENTAL SIBLING of r_patches/33_figure_workingscm_vs_pag_lvida.R --
# does NOT overwrite that file or its output. Panel A is copied verbatim
# (unchanged). Panel B is replaced entirely: instead of one point per
# compatible MAG from a SINGLE full-sample PAG (33's design), this draws a
# weighted point cloud per node built from the bootstrap-resampled
# PAG/MAG/LV-IDA results (clean_pipeline/18_bootstrap_lvida.R), so Panel B now
# reflects resampling/structural uncertainty in addition to the
# within-PAG completion uncertainty that 33's Panel B already showed.
#
# The MAG/bootstrap material (Panel B) is too technical for the main text and
# is moving to Supplementary Materials, while Panel A (Working SCM) stays in
# the main text as Figure 4 -- so the two are no longer one combined figure.
# Each panel is saved as its own standalone file (see the save section near
# the bottom): Panel A -> figures/fig4_workingscm.pdf (main text), Panel B ->
# figures/figS_bootstrap_lvida_pointcloud.pdf (Supplementary). Both panels'
# in-plot titles dropped their "A"/"B" prefixes accordingly, since neither is
# a lettered sub-panel of a shared figure anymore.
#
# ------------------------------------------------------- PANEL B DESIGN -----
# The goal is to show, per node, the distribution of LV-IDA effects across
# MAGs compatible with PAGs refit on 1,000 bootstrap resamples -- i.e.
# resampling/structural uncertainty layered on top of the within-PAG
# completion uncertainty that 33's Panel B already shows for a single PAG.
#
# Raw MAG counts are not probabilities: a resample with many compatible MAGs
# would otherwise inject that many points, over-weighting large equivalence
# classes relative to resamples with few MAGs. Per-resample quantities ARE
# legitimate frequentist quantities (each resample is one exchangeable draw),
# so every resample must count equally overall regardless of its MAG count.
# Collapsing each resample to one summary value (e.g. median identified
# effect) would fix that weighting problem but throws away the
# within-resample completion-uncertainty spread that is the whole point of
# this bootstrap. The design below instead keeps every individual MAG's raw
# effect value and controls for MAG count through a weight of
# 1 / n_mags_in_resample, so every resample's total contribution sums to 1
# regardless of how many MAGs it has.
#
# A weighted kernel density using that per-MAG weighting was tried first and
# rejected: most nodes are overwhelmingly zero-mass (many "all-agree-zero"
# MAGs, plus the unidentified-to-zero recoding described below), and a
# density curve normalizes its own peak height to 1, so the zero-spike
# becomes the peak and any real non-zero signal gets squashed toward the
# baseline and effectively disappears.
#
# A weighted point cloud (opacity instead of density mass) was tried next,
# using the same per-MAG values and weighting: this kept non-zero points
# individually visible, but hit a different failure mode from the same root
# cause (zero/recoded-zero MAGs vastly outnumbering non-zero ones) --
# thousands of near-zero points stacked inside the narrow jitter band around
# x=0, and even faint (low-alpha) points saturate to fully opaque once enough
# of them overlap in a small area. The result was a solid opaque block at
# x=0 in every row: real data, not a rendering bug, but useless as a visual,
# since once saturated the block can't show whether a node is 70% zero-weight
# or 99% zero-weight, and it dominates the panel's visual space at the
# expense of the actually-informative non-zero MAGs.
#
# Splitting the zero/unidentified share into its own separate encoding (a
# weighted-percentage bar beside the point cloud) solved the saturation
# problem, but added a second panel and a side-by-side layout that didn't
# suit the final page layout, so that path was dropped in favor of folding
# zero/unidentified points back into the same cloud at a much lower, flat
# opacity, small size, and thin stroke.
#
# CURRENT VERSION: one flat, low alpha/small size/thin stroke for every
# point. Identified_positive/identified_negative MAGs are plotted in full.
# Zero/unidentified MAGs are drawn from a WEIGHTED RANDOM SUBSAMPLE
# (dplyr::slice_sample(prop = ZERO_SUBSAMPLE_FRAC, weight_by = weight))
# instead of every row -- because each resample's own weights already sum to
# 1, sampling with probability proportional to `weight` gives every resample
# the same expected number of surviving dots regardless of its MAG count, so
# the "one resample, one equal draw" principle that used to live in the alpha
# channel now lives in which dots get drawn. Layout is Panel A stacked above
# Panel B (patchwork's `/`), not side by side.
#
# UNIDENTIFIED MAGS: MAGs for which LV-IDA could not identify an effect
# (status == "unidentified", i.e. lv.ida()/lvida_from_mags() returned NA for
# that MAG) are recoded to effect_ate_scale = 0 for the purposes of this
# plot -- so "zero" and "unidentified" are visually indistinguishable
# throughout (both thinned and plotted together in panelB at the same flat,
# low opacity). This is a real, flagged simplification, not a neutral
# default: it merges "identified as exactly no effect" with "could not be
# determined from this MAG's structure." See the caption note at the bottom
# -- this recoding needs to be stated in the manuscript methods/caption text,
# not left implicit.
#
# Not yet rendered against actual data -- this is a draft to look at, not a
# finished figure. ZERO_SUBSAMPLE_FRAC (see data-prep section below) is a
# first guess: even a fairly aggressive subsample may still leave visible
# mass at zero for an "all-agree-zero" node -- that's expected and even
# somewhat desirable (it should still read as "there's a lot of zero here"),
# the open question is just whether the current value lands that at a
# reasonable, non-block-like density once rendered, or whether it cuts so far
# that a genuinely zero-heavy node starts looking sparse instead. This is the
# one constant most likely needing adjustment.
#
# Data: pipeline_outputs/tables/bootstrap_lvida_full_vectors.csv
# (18_bootstrap_lvida.R output, 1000-resample manuscript run, NMAGS_CAP =
# 2000L, PER_BOOT_TIMEOUT_SECS = 300L).
# Output: figures/fig4_workingscm.pdf (main text) and
# figures/figS_bootstrap_lvida_pointcloud.pdf (Supplementary), saved separately.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})
source("r_patches/03_figure_style.R")

NODE_ORDER_OVERRIDE <- NULL

TARGET_ORDER <- c(
  "belief_concern", "harm_present", "harm_future", "weather_risk_prep",
  "social_norms", "trust_science", "policy_support", "politics"
)
NOT_A_LEVER <- "politics"

# =============================================================================
# Panel A data -- copied verbatim from r_patches/33_figure_workingscm_vs_pag_lvida.R.
# Unchanged: same source files, same node-ordering logic, same specification.
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
# Panel B data: bootstrap-resampled PAG/MAG/LV-IDA results
# (clean_pipeline/18_bootstrap_lvida.R output), one row per (resample, MAG,
# node) -- NOT a single fixed PAG's MAGs like 33's panel B.
# =============================================================================
boot_full <- utils::read.csv(
  "pipeline_outputs/tables/bootstrap_lvida_full_vectors.csv",
  stringsAsFactors = FALSE
)
stopifnot(all(TARGET_ORDER %in% union(boot_full$node, NOT_A_LEVER)))
stopifnot(all(c("zero", "identified_positive", "identified_negative", "unidentified") %in%
              unique(boot_full$status)))

# Same defensive numeric-with-tolerance alpha matching as 33 (read.csv's
# type.convert can coerce a quoted-but-numeric alpha column either way
# depending on how write.csv typed it upstream -- don't assume either form).
alpha_label_of <- function(a) {
  a_num <- suppressWarnings(as.numeric(a))
  ifelse(!is.na(a_num) & abs(a_num - 0.05) < 1e-9, "alpha = .05",
  ifelse(!is.na(a_num) & abs(a_num - 0.01) < 1e-9, "alpha = .01",
         paste0("alpha = ", a)))
}
ALPHA_LEVELS <- c("alpha = .05", "alpha = .01")

boot_full <- boot_full |>
  filter(node %in% TARGET_ORDER) |>
  mutate(
    node_label = factor(node_labels_oneline[node], levels = node_label_levels),
    alpha_label = factor(alpha_label_of(alpha), levels = ALPHA_LEVELS),
    # Fold "could not be identified" into "zero effect", per the header
    # comment above -- this is a real methodological choice, not a neutral
    # cleanup step.
    effect_for_density = dplyr::if_else(status == "unidentified", 0, effect_ate_scale)
  )
stopifnot(!anyNA(boot_full$effect_for_density))

# Weight scheme: every bootstrap resample contributes total weight 1, split
# evenly across however many MAGs THAT resample has (1 / n_mags_b each). A
# resample with 372 MAGs still shows all 372 of its raw values (full
# within-resample completion-uncertainty spread preserved), but they only sum
# to the same total mass as a resample with 1 MAG -- every resample counts
# equally overall. n_mags_b is recomputed directly here (count of rows per
# alpha/node/boot_id) rather than joined from
# bootstrap_lvida_resample_summary.csv, so this script has no dependency on
# that file being in sync.
boot_full <- boot_full |>
  dplyr::group_by(alpha_label, node_label, boot_id) |>
  dplyr::mutate(n_mags_b = dplyr::n(), weight = 1 / n_mags_b) |>
  dplyr::ungroup()

# =============================================================================
# x-axis range for the two quantitative panels (A and B) -- shared, so the
# effect scale is directly comparable across them, same convention as 33.
# =============================================================================
all_effect_x <- c(baseline_df$ate_climate_behavior, alternative_df$ate_climate_behavior,
                   ci_df$ci_lower_95, ci_df$ci_upper_95,
                   boot_full$effect_for_density)
x_lo <- min(-0.02, min(all_effect_x, na.rm = TRUE) - 0.02)
x_hi <- max(all_effect_x, na.rm = TRUE) + 0.02
x_breaks <- scales::pretty_breaks(n = 5)(c(x_lo, x_hi))
X_LAB <- "Effect on climate behavior (ΔY, SD)"

Y_SCALE_B <- scale_y_continuous(
  breaks = node_positions, labels = names(node_positions),
  limits = c(0.35, length(node_order) + 0.65), expand = c(0, 0)
)

# =============================================================================
# Weighted point-cloud data prep. Same design language as 33's original panel
# B (one point per individual compatible MAG, small jitter so multiplicity
# reads as a cloud rather than a stick -- JITTER_X_HALFWIDTH/JITTER_Y_HALFWIDTH
# values copied from 33 rather than re-derived, since that spacing was
# already tuned there).
#
# Panel B's aesthetics use one flat low alpha/small size/thin stroke for
# every point (no per-point alpha driven by `weight`). Flattening alpha to a
# constant drops the "alpha = weight" mechanism that would otherwise keep a
# resample's visual contribution equal regardless of its MAG count -- that
# mechanism is reinstated here differently for the zero/unidentified pile:
# rather than varying opacity, a WEIGHTED RANDOM SUBSAMPLE of the
# zero/unidentified rows is drawn (dplyr::slice_sample(..., weight_by =
# weight)), keeping only ZERO_SUBSAMPLE_FRAC of them. Because each resample's
# own weights already sum to exactly 1 (weight = 1/n_mags_b, summed over that
# resample's n_mags_b rows), sampling with selection probability proportional
# to `weight` gives every resample the SAME expected number of surviving dots
# regardless of its MAG count (expected survivors from resample r is
# proportional to sum(weight_r) = 1, the same for every resample) -- so the
# same "one resample, one equal draw" principle is preserved, just through
# which dots get drawn rather than through their opacity.
# Identified_positive/identified_negative points are NOT thinned -- there are
# comparatively few of them per row (typically far short of the thousands
# seen for zero/unidentified), so the flat alpha used for them doesn't
# reproduce the same over-weighting risk in practice.
#
# ZERO_SUBSAMPLE_FRAC = 0.02 -- the one constant to move if the pile still
# looks too dense or too sparse once actually rendered.
# =============================================================================
JITTER_X_HALFWIDTH <- 0.006
JITTER_Y_HALFWIDTH <- 0.11
ZERO_SUBSAMPLE_FRAC <- 0.02
POINT_ALPHA <- 0.2
POINT_SIZE  <- 1
POINT_STROKE <- 0.1

set.seed(20260910)
point_df_all <- boot_full |>
  dplyr::mutate(
    y_row_nudge = ifelse(alpha_label == "alpha = .05", 0.20, -0.20),
    y_num = node_positions[as.character(node_label)] + y_row_nudge +
            stats::runif(dplyr::n(), -JITTER_Y_HALFWIDTH, JITTER_Y_HALFWIDTH),
    x_plot = effect_for_density +
             stats::runif(dplyr::n(), -JITTER_X_HALFWIDTH, JITTER_X_HALFWIDTH)
  )

identified_rows <- point_df_all |>
  dplyr::filter(status %in% c("identified_positive", "identified_negative"))

zero_rows_thinned <- point_df_all |>
  dplyr::filter(status %in% c("zero", "unidentified")) |>
  dplyr::group_by(alpha_label, node_label) |>
  dplyr::slice_sample(prop = ZERO_SUBSAMPLE_FRAC, weight_by = weight) |>
  dplyr::ungroup()

point_df <- dplyr::bind_rows(identified_rows, zero_rows_thinned)

# Console-only diagnostic (not plotted): how much of each node/alpha's total
# resample-weight is zero/unidentified -- kept as a printed sanity check
# further down even though it no longer has its own panel, so you can still
# see the number behind however solid or faint that part of the cloud ends
# up looking.
zero_summary <- boot_full |>
  dplyr::group_by(alpha_label, node_label) |>
  dplyr::summarise(
    pct_zero_weight = 100 * sum(weight[status %in% c("zero", "unidentified")]) / sum(weight),
    .groups = "drop"
  )

# =============================================================================
# Panel A: unchanged, copied verbatim from 33.
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
  labs(x = X_LAB, y = NULL, title = "Working SCM",
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
# Panel B: bootstrap point cloud -- identified_positive/identified_negative
# MAGs shown in full, zero/unidentified MAGs shown as a weighted subsample
# (ZERO_SUBSAMPLE_FRAC of them, see data-prep section above for why this
# preserves equal per-resample representation). One row per node, two alphas
# per row (offset +/-0.20, same convention as 33). Flat, low alpha/small
# size/thin stroke for every point -- no per-point alpha column, so no
# scale_alpha_identity() needed here.
#
# Stacked vertically with panelA above, not side by side -- so panelB gets
# its own y-axis labels back (no longer redundant with an adjacent panelA
# sharing the same rows).
# =============================================================================
panelB <- ggplot() +
  geom_vline(xintercept = 0, linewidth = .45, linetype = "dashed", colour = "#9A9A9A") +
  geom_point(
    data = point_df,
    aes(x = x_plot, y = y_num, colour = alpha_label, fill = alpha_label),
    shape = 21, size = POINT_SIZE, stroke = POINT_STROKE, alpha = POINT_ALPHA
  ) +
  scale_colour_manual(name = NULL, values = c("alpha = .05" = COL$blue_dark, "alpha = .01" = COL$rust)) +
  scale_fill_manual(name = NULL, values = c("alpha = .05" = COL$blue_dark, "alpha = .01" = COL$rust)) +
  guides(
    colour = guide_legend(override.aes = list(size = 2.6, alpha = 1)),
    # fill kept suppressed here -- otherwise it draws a second, redundant
    # fill legend beside the colour one for the same alpha_label variable.
    fill = "none"
  ) +
  scale_x_continuous(limits = c(x_lo, x_hi), breaks = x_breaks,
                      expand = expansion(mult = c(0, .01))) +
  Y_SCALE_B +
  labs(x = X_LAB, y = NULL, title = "PAG-compatible effects (bootstrap)",
       subtitle = "Point cloud (zero/unidentified MAGs subsampled)") +
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
# save -- separate files, not a combined figure. Panel B (the bootstrap point
# cloud) moves to Supplementary Materials, so it no longer shares a figure
# with Panel A (Working SCM, staying in the main text as Figure 4) -- there
# is no longer a single "Figure 4" that contains both. Each panel now stands
# alone as its own figure, sized and titled accordingly (title text above no
# longer has an "A"/"B" prefix, since neither is a lettered sub-panel of a
# shared figure anymore).
#
# Dimensions are first guesses, not measured against a rendered output --
# panelA_fig is sized for a main-text single-panel figure, panelB_fig for a
# supplementary one; both are plausible starting points for a portrait,
# 8-row forest-plot-style figure but may need adjusting once actually
# rendered.
# =============================================================================
save_ms_figure(panelA, "figures/fig4_workingscm.pdf",
               c(width = 140, height = 110))
message("Saved: figures/fig4_workingscm.pdf (main text, Figure 4)")

save_ms_figure(panelB, "figures/figS_bootstrap_lvida_pointcloud.pdf",
               c(width = 150, height = 130))
message("Saved: figures/figS_bootstrap_lvida_pointcloud.pdf (Supplementary)")

message("Node order (top to bottom): ", paste(rev(node_order), collapse = " > "))

# Quick sanity printout -- how many resamples actually contributed at least
# one identified-non-zero MAG to each node/alpha (note: filters point_df
# down to the identified statuses first -- point_df itself now holds ALL
# statuses since zero/unidentified were folded back in above, so grouping it
# directly without this filter would just count every resample). Compare
# against zero_summary's pct_zero_weight and against
# bootstrap_lvida_node_summary.csv's own pct_resamples_any_positive /
# pct_resamples_any_negative, which are resample-level "any," not
# MAG-weighted like these.
weight_check <- point_df |>
  dplyr::filter(status %in% c("identified_positive", "identified_negative")) |>
  dplyr::group_by(alpha_label, node_label) |>
  dplyr::summarise(n_resamples_with_identified = dplyr::n_distinct(boot_id), .groups = "drop") |>
  dplyr::full_join(zero_summary, by = c("alpha_label", "node_label")) |>
  dplyr::mutate(n_resamples_with_identified = ifelse(is.na(n_resamples_with_identified), 0L,
                                                       n_resamples_with_identified)) |>
  dplyr::select(alpha_label, node_label, n_resamples_with_identified, pct_zero_weight) |>
  dplyr::arrange(node_label, alpha_label)
cat("\n--- Panel B (bootstrap point cloud) sanity check ---\n")
print(as.data.frame(weight_check), row.names = FALSE)

# =============================================================================
# Caption (DRAFT -- not inserted into main9.tex; kept here for reference
# only, main9.tex is not touched by this project until explicitly decided).
# Adapted from 33's Panel B caption language to describe the bootstrap
# version. The unidentified-to-zero recoding is called out explicitly per
# the discussion above -- this needs to survive into whatever caption
# actually ships, not just live in this comment.
#
# \caption{Single-node intervention effects under the working SCM and under
# bootstrap-resampled structural uncertainty. \textbf{(A)} Effects of setting
# each attitude/context variable to $+0.5$ SD in the baseline working SCM
# (diamonds), with effects from the 15 alternative directional
# specifications shown as points. Horizontal intervals show 95\%
# participant-bootstrap confidence intervals for the baseline model.
# Political orientation is included as a structural comparison rather than
# as a realistic intervention target. \textbf{(B)} LV-IDA effects across
# MAGs compatible with FCI PAGs refit on 1,000 participant-resampled
# datasets at $\alpha=.05$ and $\alpha=.01$, shown on the same $+0.5$-SD
# intervention scale as (A). Each point represents one MAG for which LV-IDA
# identified a non-zero effect. MAGs for which the total effect was
# identified as exactly zero, or could not be identified at all from the
# MAG's structure, are shown at zero as a random subsample (weighted so that
# resamples with many compatible MAGs are not over-represented relative to
# resamples with few), so the figure conveys how common this category is
# without every such MAG being individually plotted.}
# =============================================================================
