# =============================================================================
# 10_figure7_uncertainty_pub_violin.R  -- OPTION A for the alternative-
# orientation encoding, produced alongside 10_figure7_uncertainty_pub_jitter.R
# after several rounds of iteration on this one design choice didn't land
# (beeswarm -> KDE hump -> 6-bin histogram -> clamped KDE hump -> 14-bin
# histogram -- see 10_figure7_uncertainty_pub.R's header for the full history).
# A full-distribution violin and a plain jitter were requested side by side
# for comparison, rather than iterating blind again. This is the "full
# distribution" option: a standard symmetric ggplot2 violin per node.
#
# Uses geom_violin()'s own `trim = TRUE` default (not manual density code
# this time) -- ggplot2 already clips the estimated density to each group's
# actual data range out of the box, which is exactly the "don't extend past
# where the real points are" property we were hand-rolling before. With only
# n=7 per node this is still a rough density estimate (7 points, not a
# sample), so the shape should be read as "roughly where the 7 structural
# alternatives cluster," not as a precise smooth population -- same caveat
# as every version of this figure has carried.
#
# Everything else (baseline diamond + value label, axis range, legend,
# theme) is identical to 10_figure7_uncertainty_pub.R. Output:
# figures/fig7_uncertainty_pub_violin.pdf
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
    ),
    kind = factor(
      if_else(is_baseline, "Baseline specification", "Alternative orientation"),
      levels = c("Baseline specification", "Alternative orientation")
    )
  )

baseline_df    <- filter(plot_df, is_baseline)
alternative_df <- filter(plot_df, !is_baseline)

# Data-driven axis range (see 10_figure7_uncertainty_pub.R's note on this).
data_min <- min(plot_df$ate_climate_behavior, na.rm = TRUE)
data_max <- max(plot_df$ate_climate_behavior, na.rm = TRUE)
data_span <- data_max - data_min
x_lo <- min(-.015, data_min - .05 * data_span)
x_hi <- data_max + .12 * data_span

fig7_violin <- ggplot() +
  geom_vline(xintercept=0, linewidth=.45, linetype="dashed", colour="#9A9A9A") +
  geom_violin(
    data=alternative_df,
    aes(x=ate_climate_behavior, y=node_label, group=node_label),
    orientation="y", trim=TRUE, scale="width", width=0.7,
    fill=COL$blue_dark, colour=NA, alpha=.30
  ) +
  geom_point(
    data=alternative_df,
    aes(x=ate_climate_behavior,y=node_label,shape=kind,fill=kind,colour=kind),
    size=0, stroke=0, alpha=0, show.legend=TRUE
  ) +
  geom_point(
    data=baseline_df,
    aes(x=ate_climate_behavior,y=node_label,shape=kind,fill=kind,colour=kind),
    size=3.6, stroke=.8, alpha=.72
  ) +
  geom_text(
    data=baseline_df,
    aes(x=ate_climate_behavior,y=node_label,
        label=sprintf("+%.3f SD", ate_climate_behavior)),
    hjust=-.35, vjust=-1.1, size=PT(8.2), colour=COL$ink, family=FIG_FONT
  ) +
  scale_shape_manual(
    name=NULL,
    values=c("Baseline specification"=23, "Alternative orientation"=21)
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
    shape=guide_legend(override.aes=list(size=c(4.6,3.4), stroke=c(1.0,0), alpha=c(1,.65))),
    fill=guide_legend(override.aes=list(size=c(4.6,3.4), stroke=c(1.0,0), alpha=c(1,.65))),
    colour=guide_legend(override.aes=list(size=c(4.6,3.4), stroke=c(1.0,0), alpha=c(1,.65)))
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
    plot.margin=margin(6,16,4,4)
  )

save_ms_figure(fig7_violin, "figures/fig7_uncertainty_pub_violin.pdf",
               c(width = 172, height = 108))
message("Saved Figure 7 [OPTION A: violin]: figures/fig7_uncertainty_pub_violin.pdf")
