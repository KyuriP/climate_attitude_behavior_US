# =============================================================================
# 10_figure7_uncertainty_pub.R
#
# Merged main-text figure: for each candidate intervention node, shows (1)
# the model-implied effect on climate_behavior under the baseline working
# SCM (combo_0), labeled directly with its value, and (2) how that effect
# moves across the other 15 acyclic, converged alternative orientations of
# the four genuinely uncertain edges (politics->belief_concern,
# policy_support->social_norms, social_norms->climate_behavior,
# harm_future->harm_present). This absorbs both the old intervention barplot
# and the old separate uncertainty-range chart into one figure. (The "7" in
# the filename is a leftover from when only 3 edges were flippable, 2^3=8
# specs, 7 alternatives; the 4th flip edge (harm_future<->harm_present) was
# added later in r_patches/02_full_orientation_enumeration_v4.R, making it
# 16 specs / 15 alternatives -- the plotted data was always correct either
# way, since it just filters on scenario != "combo_0" rather than
# hardcoding a count.)
#
# Uses a plain vertical jitter of the 15 alternative-orientation points per
# node rather than a smoothed density -- several other encodings were tried
# (beeswarm, KDE hump, binned histogram, capped-bandwidth half-density hump)
# but on reflection the jitter is the more honest picture: 15 discrete
# points is not enough to support any smoothed density without implying a
# continuous distribution that isn't there. 10_figure7_uncertainty_pub_jitter.R
# carries the same reasoning for the side-by-side comparison it was built
# for; that file's plotting code is now this script's canonical code
# (output filename unchanged, so no manuscript-side references need to
# change).
#
# Reads tables/orientation_enumeration_ate_deterministic_8node.csv (see Data
# source below), not the older tables/orientation_enumeration_ate.csv
# 20,000-draw Monte Carlo table -- that older file predates script 30's
# exact mean-propagation values and would leave this figure's baseline
# labels quietly out of step with the manuscript text and Table 3.
#
# Covers all 8 non-outcome nodes in the working SCM, using
# tables/orientation_enumeration_ate_deterministic_8node.csv from
# r_patches/32_deterministic_8node_ates.R (see that script's header for why
# harm_future and politics never had single-node ATEs computed before now).
# harm_future is added as an ordinary row -- it's a belief construct like
# present harm, "do(harm_future=0.5)" reads the same way as "do(harm_present=
# 0.5)". politics is NOT an ordinary row: it isn't an actionable intervention
# (nobody's political identity is a messaging target), so it's kept visually
# distinct -- own colour, own shape, its own legend entry and a caption
# sentence saying explicitly that it's a descriptive benchmark, not a
# candidate lever, so the figure can't be skimmed as implying otherwise.
#
# Data source: tables/orientation_enumeration_ate_deterministic_8node.csv,
# produced by r_patches/32_deterministic_8node_ates.R (16 base edges, 16
# structural specifications, exact mean propagation -- no simulation noise;
# superset of tables/orientation_enumeration_ate_deterministic.csv's 6 nodes).
# Falls back to reading that CSV directly if full_results isn't already in
# the session, so this script also runs standalone.
#
# Keeping all 8 rows in the main figure rather than reverting to the
# original 6 (moving future harm/politics to Supplement-only, which was
# considered): Table S7 reports the same 8 nodes for cross-reference, and
# this figure and that table should not disagree in scope.
#
# Design notes -- alternative-orientation encoding
# -------------------------------------------------
# Each node has exactly 15 alternative-orientation values (the non-baseline
# structural specifications), a complete enumeration rather than a sample.
# Encoding: a plain vertical jitter of the 15 real points (blue, semi-
# transparent, no shape beyond the point itself -- no smoothing, no implied
# density). The baseline specification is marked separately with a rust
# diamond and its own "+X.XXX SD" label -- this is the number reported in
# text/Table 3 for the 6 policy-relevant nodes. politics uses a third,
# visually distinct colour/shape pair for both its baseline and alternative
# points (see "REVISION 2026-09-07 (b)" above).
#
# Output: figures/fig7_uncertainty_pub.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})
source("r_patches/03_figure_style.R")

# This used to be `if (!exists("full_results")) {...}`, a convenience so the
# script could reuse an already-loaded object instead of hitting disk. That
# silently backfired the first time this script grew a new data source (the
# 8-node file): a `full_results` object left over in the R session from an
# earlier script (e.g. 02/30, which only ever build the 6-node version)
# satisfied `exists()`, so this script quietly reused the stale 6-node
# object instead of loading the new file -- no error, no warning, just 6
# rows on the plot instead of 8. Always reading fresh from disk now; the
# disk read is cheap (a few KB) and "definitely current data" is worth more
# here than "skip a fast read sometimes."
full_results <- utils::read.csv("pipeline_outputs/tables/orientation_enumeration_ate_deterministic_8node.csv",  # this writes under pipeline_outputs/tables now, not bare tables/
                                 stringsAsFactors = FALSE)
stopifnot(all(c("harm_future", "politics") %in% full_results$node))

TARGET_ORDER <- c(
  "harm_present","belief_concern","harm_future","weather_risk_prep",
  "social_norms","trust_science","policy_support","politics"
)
NOT_A_LEVER <- "politics"  # descriptive benchmark only -- see header note.

KIND_LEVELS <- c(
  "Baseline specification", "Alternative orientation",
  "Politics (baseline)", "Politics (alternative)"
)

plot_df <- full_results |>
  filter(node %in% TARGET_ORDER) |>
  mutate(
    is_baseline = scenario == "combo_0",
    is_politics = node == NOT_A_LEVER,
    node_label = factor(
      node_labels_oneline[node],
      levels = rev(node_labels_oneline[TARGET_ORDER])
    ),
    # Baseline/alternative still drives shape (diamond vs circle) for every
    # node, including politics -- politics just gets its own two colours
    # within that same shape scheme, so it reads as "same kind of point,
    # different (muted) group" rather than a wholly separate visual system.
    kind = factor(
      dplyr::case_when(
        is_politics &  is_baseline ~ "Politics (baseline)",
        is_politics & !is_baseline ~ "Politics (alternative)",
        is_baseline                ~ "Baseline specification",
        TRUE                       ~ "Alternative orientation"
      ),
      levels = KIND_LEVELS
    )
  )

baseline_df    <- filter(plot_df, is_baseline)
alternative_df <- filter(plot_df, !is_baseline)

# Data-driven axis range: computed from the data itself (rather than a
# hardcoded literal) so a future change in the value range can't silently
# clip a point off the edge of the plot.
data_min <- min(plot_df$ate_climate_behavior, na.rm = TRUE)
data_max <- max(plot_df$ate_climate_behavior, na.rm = TRUE)
data_span <- data_max - data_min
x_lo <- min(-.015, data_min - .05 * data_span)
x_hi <- data_max + .12 * data_span

fig7 <- ggplot() +
  geom_vline(xintercept=0, linewidth=.45, linetype="dashed", colour="#9A9A9A") +
  geom_point(
    data=alternative_df,
    aes(x=ate_climate_behavior,y=node_label,shape=kind,fill=kind,colour=kind),
    size=2.6, stroke=0, alpha=.3,
    # height bumped .09 -> .16: the 15 alternative points per row were
    # sitting too close together at .09 to read as separated jitter rather
    # than a smear; .16 still leaves clear separation from the adjacent
    # node's row (rows are 1 unit apart, so .16 is under a fifth of that gap
    # on either side).
    position=position_jitter(height=.16,width=0.006,seed=20260830)
  ) +
  geom_point(
    data=baseline_df,
    aes(x=ate_climate_behavior,y=node_label,shape=kind,fill=kind,colour=kind),
    size=2.6, stroke=.5, alpha=.72
  ) +
  geom_text(
    data=baseline_df,
    aes(x=ate_climate_behavior,y=node_label,
        label=sprintf("+%.3f SD", ate_climate_behavior)),
    hjust=-.35, vjust=-1.1, size=PT(8.2), colour=COL$ink, family=FIG_FONT
  ) +
  # politics keeps its own muted-grey colour/shape values
  # (still visually distinct from the other 7 rows) but is no longer given
  # its own legend entries -- `breaks=` below limits the drawn legend keys to
  # just the 2 real categories, even though `values=` still defines all 4 so
  # politics' points render correctly. The politics caveat itself moved out
  # of the figure entirely (was a plot.caption here, now just described in
  # the LaTeX figure caption in the manuscript) -- per-viewer request that
  # the figure stay visually simple and put the explanation in prose instead.
  scale_shape_manual(
    name=NULL,
    values=c("Baseline specification"=23, "Alternative orientation"=21,
             "Politics (baseline)"=23,
             "Politics (alternative)"=21),
    breaks=c("Baseline specification", "Alternative orientation")
  ) +
  scale_fill_manual(
    name=NULL,
    values=c("Baseline specification"=COL$rust, "Alternative orientation"=COL$blue_dark,
             "Politics (baseline)"=COL$ink_mid,
             "Politics (alternative)"=COL$ink_light),
    breaks=c("Baseline specification", "Alternative orientation")
  ) +
  scale_colour_manual(
    name=NULL,
    values=c("Baseline specification"=COL$ink, "Alternative orientation"=COL$blue_dark,
             "Politics (baseline)"=COL$ink, "Politics (alternative)"=COL$ink_light),
    breaks=c("Baseline specification", "Alternative orientation")
  ) +
  # Back to a single legend row now that only 2 breaks are drawn (politics'
  # 2 entries dropped above). Scalar override.aes kept as-is: a scalar here
  # (rather than a per-level vector) is what survives across ggplot2
  # versions' differing guide-merge behavior.
  guides(
    shape=guide_legend(override.aes=list(size=3.2, stroke=0.4, alpha=0.85)),
    fill=guide_legend(override.aes=list(size=3.2, stroke=0.4, alpha=0.85)),
    colour=guide_legend(override.aes=list(size=3.2, stroke=0.4, alpha=0.85))
  ) +
  scale_x_continuous(
    limits=c(x_lo, x_hi), breaks=scales::pretty_breaks(n = 5)(c(x_lo, x_hi)),
    expand=expansion(mult=c(0,.01))
  ) +
  labs(x="Model-implied intervention effect on climate behavior (SD)", y=NULL) +
  theme_pub +
  theme(
    text            = element_text(size = 10.5),
    axis.text.y=element_text(size=10.2),
    axis.text.x=element_text(size=9.2),
    axis.title.x=element_text(size=9.8),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    panel.grid.major.x=element_line(colour=COL$grid,linewidth=.35),
    panel.grid.major.y=element_blank(),
    legend.position="bottom",
    legend.text=element_text(size=8.4),
    plot.margin=margin(6,16,4,4)
  )

# Height: 140mm from the 8-row version, trimmed slightly to 132mm now that
# the caption line and second legend row are both gone (one row of 2 legend
# items + no caption text needs less vertical footprint than 2 rows + a
# 2-line caption did).
save_ms_figure(fig7, "figures/fig7_uncertainty_pub.pdf",
               c(width = 172, height = 132))
message("Saved Figure 7 (merged single-node intervention + uncertainty figure): ",
        "figures/fig7_uncertainty_pub.pdf (x-axis: ", round(x_lo,3), " to ", round(x_hi,3), ")")

