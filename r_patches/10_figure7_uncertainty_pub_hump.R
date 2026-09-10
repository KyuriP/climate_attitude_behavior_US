# =============================================================================
# 10_figure7_uncertainty_pub_hump.R  -- OPTION C for the alternative-
# orientation encoding, a third comparison point alongside
# 10_figure7_uncertainty_pub_violin.R and 10_figure7_uncertainty_pub_jitter.R
# (see those files' headers, and 10_figure7_uncertainty_pub.R's header, for
# the full back-and-forth that led here).
#
# This brings back the earlier half-density hump per node (not a full
# symmetric violin), with the specific problem fixed this time: earlier
# attempts looked lopsided -- a slanted, triangle-like taper (clearest on
# present_harm, whose 7 points have a big gap between a cluster near
# baseline and one outlier near .24).
# The previous fix only clamped the curve's ENDPOINTS to the real data range
# (from/to = range(x)); it didn't address the cause of the wedge shape
# itself, which is bandwidth: stats::bw.nrd0() picks ONE global bandwidth
# per node, and when 6 points cluster tightly and 1 sits far away, that
# bandwidth ends up wide enough to smear the whole thing into one long
# diagonal ramp connecting the cluster to the outlier, instead of showing
# a distinct bump where the 6 actually are. Fixed by capping the bandwidth
# to a fraction of the node's own range (diff(range(x))/6, floored so it
# never collapses to ~0 for a tight cluster) -- this keeps the curve
# responsive to real local clustering (visible bumps near where points
# actually sit) rather than one smoothed-over wedge across the whole span.
# Endpoints are still clamped to range(x), so no extrapolation past the data.
#
# Uses the same y_num numeric-position workaround as
# 10_figure7_uncertainty_pub.R (see that file's note on the "discrete values
# supplied to continuous scale" error) since this layer, like the
# histogram/hump attempts there, needs raw numeric y offsets rather than the
# node_label factor directly.
#
# Output: figures/fig7_uncertainty_pub_hump.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})
source("r_patches/03_figure_style.R")

if (!exists("full_results")) {
  full_results <- utils::read.csv("pipeline_outputs/tables/orientation_enumeration_ate_deterministic_8node.csv",  # was reading the old monte carlo file (script 02 v4) by mistake -- pointing at the same deterministic file the main figure uses now
                                   stringsAsFactors = FALSE)
  message("full_results not found in session -- loaded from ",
          "pipeline_outputs/tables/orientation_enumeration_ate_deterministic_8node.csv.")
}

TARGET_ORDER <- c(
  "harm_present","belief_concern","weather_risk_prep",
  "social_norms","trust_science","policy_support"
)

plot_df <- full_results |>
  filter(node %in% TARGET_ORDER) |>
  mutate(
    is_baseline = scenario == "combo_0",
    node_label = factor(
      node_labels_oneline[node],
      levels = rev(node_labels_oneline[TARGET_ORDER])
    )
  ) |>
  mutate(
    kind = factor(
      if_else(is_baseline, "Baseline specification", "Alternative orientation"),
      levels = c("Baseline specification", "Alternative orientation")
    ),
    # Numeric row position, used by every layer below -- see header note.
    y_num = as.numeric(node_label)
  )

NODE_LEVELS <- levels(plot_df$node_label)

baseline_df    <- filter(plot_df, is_baseline)
alternative_df <- filter(plot_df, !is_baseline)

# Half-density hump per node, bandwidth capped so it can't smear a tight
# cluster + a lone outlier into one diagonal wedge -- see header note.
HUMP_HEIGHT <- 0.34
dens_hump <- alternative_df |>
  dplyr::group_by(node_label) |>
  dplyr::group_modify(function(d, key) {
    x <- d$ate_climate_behavior
    row_y <- as.numeric(key$node_label)
    rng <- range(x)
    if (diff(rng) < 1e-9) {
      return(tibble::tibble(x = numeric(0), y = numeric(0)))
    }
    bw_default <- stats::bw.nrd0(x)
    bw_cap     <- diff(rng) / 6
    bw_floor   <- diff(rng) / 20
    bw <- max(min(bw_default, bw_cap), bw_floor)
    dd <- stats::density(x, bw = bw, n = 128, from = rng[1], to = rng[2])
    hump_y <- row_y + (dd$y / max(dd$y)) * HUMP_HEIGHT
    tibble::tibble(
      x = c(dd$x[1], dd$x, dd$x[length(dd$x)]),
      y = c(row_y,   hump_y, row_y)
    )
  }, .keep = TRUE) |>
  dplyr::ungroup()

# One short tick per actual alternative-orientation value, below the row
# line, so the real discrete values stay visible under the hump.
tick_df <- alternative_df |>
  dplyr::mutate(row_y = y_num)

# Data-driven axis range (see 10_figure7_uncertainty_pub.R's note on this).
data_min <- min(plot_df$ate_climate_behavior, na.rm = TRUE)
data_max <- max(plot_df$ate_climate_behavior, na.rm = TRUE)
data_span <- data_max - data_min
x_lo <- min(-.015, data_min - .05 * data_span)
x_hi <- data_max + .12 * data_span

fig7_hump <- ggplot() +
  geom_vline(xintercept=0, linewidth=.45, linetype="dashed", colour="#9A9A9A") +
  geom_polygon(
    data=dens_hump,
    aes(x=x,y=y,group=node_label),
    fill=COL$blue_dark, colour=NA, alpha=.30
  ) +
  geom_segment(
    data=tick_df,
    aes(x=ate_climate_behavior,xend=ate_climate_behavior,
        y=row_y-0.16,yend=row_y-0.04),
    colour=COL$blue_dark, alpha=.55, linewidth=.6, lineend="round"
  ) +
  geom_point(
    data=alternative_df,
    aes(x=ate_climate_behavior,y=y_num,shape=kind,fill=kind,colour=kind),
    size=0, stroke=0, alpha=0, show.legend=TRUE
  ) +
  geom_point(
    data=baseline_df,
    aes(x=ate_climate_behavior,y=y_num,shape=kind,fill=kind,colour=kind),
    size=3, stroke=.8, alpha=.72
  ) +
  geom_text(
    data=baseline_df,
    aes(x=ate_climate_behavior,y=y_num,
        label=sprintf("+%.3f SD", ate_climate_behavior)),
    hjust=-.35, vjust=-1.1, size=PT(8.2), colour=COL$ink, family=FIG_FONT
  ) +
  scale_shape_manual(
    name=NULL,
    values=c("Baseline specification"=23, "Alternative orientation"=124)
  ) +
  scale_fill_manual(
    name=NULL,
    values=c("Baseline specification"=COL$blue_dark, "Alternative orientation"=COL$blue_dark)
  ) +
  scale_colour_manual(
    name=NULL,
    values=c("Baseline specification"=COL$ink, "Alternative orientation"=COL$blue_dark)
  ) +
  guides(
    # "Alternative orientation" is drawn on the chart only as a short tick
    # (geom_segment) under the hump -- shape 124 ("|") in the legend key
    # matches that, replacing the earlier open-circle key (shape 21) that
    # implied a dot mark which is never actually plotted.
    shape=guide_legend(override.aes=list(size=c(3,4.3), stroke=c(1.0,1.0), alpha=c(1,.65))),
    fill=guide_legend(override.aes=list(size=c(3,4.3), stroke=c(1.0,1.0), alpha=c(1,.65))),
    colour=guide_legend(override.aes=list(size=c(3,4.3), stroke=c(1.0,1.0), alpha=c(1,.65)))
  ) +
  scale_x_continuous(
    limits=c(x_lo, x_hi), breaks=scales::pretty_breaks(n = 5)(c(x_lo, x_hi)),
    expand=expansion(mult=c(0,.01))
  ) +
  scale_y_continuous(
    breaks=seq_along(NODE_LEVELS), labels=NODE_LEVELS,
    expand=expansion(add=0.6)
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
    plot.margin=margin(6,16,4,4)
  )

save_ms_figure(fig7_hump, "figures/fig7_uncertainty_pub_hump.pdf",
               c(width = 172, height = 108))
message("Saved Figure 7 [OPTION C: half-density hump, bandwidth-capped]: ",
        "figures/fig7_uncertainty_pub_hump.pdf")
