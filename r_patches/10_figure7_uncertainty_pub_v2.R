# figure5_uncertainty (single-node intervention effects + spread across the
# 16 structural specs). v2 - just swapping the data source.
#
# v1 read from tables/orientation_enumeration_ate.csv, the Monte Carlo ATEs
# (script 02 v4, n=20,000 draws). During the manuscript pass I realized the
# figure still had those old MC numbers baked in while the text (S10, Table
# 3, section 3.5) had already moved to the deterministic ATEs from script 30
# (exact mean propagation - no more spurious negative floor on policy_support,
# and no spurious positive floor on trust_science either). That mismatch was
# the biggest thing off between the figure and the text, so this is just the
# fix: point it at tables/orientation_enumeration_ate_deterministic.csv
# instead.
#
# everything else below is untouched from v1 - bandwidth calc, silhouette
# geometry, ticks, baseline diamond + label, axis scaling, all identical.
# added two stopifnot checks after loading so it fails loudly instead of
# quietly plotting whatever's in the csv if the row counts are ever off.
#
# output path is unchanged (figures/fig7_uncertainty_pub.pdf) so nothing
# needs to change on the main2.tex side.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})
source("r_patches/03_figure_style.R")

DETERMINISTIC_ATE_PATH <- "pipeline_outputs/tables/orientation_enumeration_ate_deterministic.csv"  # this writes under pipeline_outputs/tables now, not bare tables/
stopifnot(file.exists(DETERMINISTIC_ATE_PATH))
full_results <- utils::read.csv(DETERMINISTIC_ATE_PATH, stringsAsFactors = FALSE)

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
    # keep everything on this numeric y so the whole plot is on one
    # continuous scale - mixing the factor and raw numeric offsets across
    # layers throws a "discrete values supplied to continuous scale" error
    y_num = as.numeric(node_label)
  )

NODE_LEVELS <- levels(plot_df$node_label)

baseline_df    <- filter(plot_df, is_baseline)
alternative_df <- filter(plot_df, !is_baseline)

stopifnot(nrow(baseline_df) == length(TARGET_ORDER))
stopifnot(nrow(alternative_df) == length(TARGET_ORDER) * 15)  # 16 specs - 1 baseline

# half-density silhouette per node
HUMP_HEIGHT <- 0.34
dens_hump <- alternative_df |>
  dplyr::group_by(node_label) |>
  dplyr::group_modify(function(d, key) {
    x <- d$ate_climate_behavior
    row_y <- as.numeric(key$node_label)
    rng <- range(x)
    if (diff(rng) < 1e-9) {
      # all 15 alt values identical, nothing to draw - tick mark below still
      # shows it
      return(tibble::tibble(x = numeric(0), y = numeric(0)))
    }
    bw_default <- stats::bw.nrd0(x)
    bw_cap     <- diff(rng) / 6
    bw_floor   <- diff(rng) / 20
    bw <- max(min(bw_default, bw_cap), bw_floor)
    # pin from/to to the actual data range so the curve can't spill past
    # where the 15 points actually are
    dd <- stats::density(x, bw = bw, n = 128, from = rng[1], to = rng[2])
    hump_y <- row_y + (dd$y / max(dd$y)) * HUMP_HEIGHT
    tibble::tibble(
      x = c(dd$x[1], dd$x, dd$x[length(dd$x)]),
      y = c(row_y,   hump_y, row_y)
    )
  }, .keep = TRUE) |>
  dplyr::ungroup()

# one tick per alternative-orientation value, just under the row line
tick_df <- alternative_df |>
  dplyr::mutate(row_y = y_num)

# axis range from the data itself, not hardcoded, so it can't silently clip
data_min <- min(plot_df$ate_climate_behavior, na.rm = TRUE)
data_max <- max(plot_df$ate_climate_behavior, na.rm = TRUE)
data_span <- data_max - data_min
x_lo <- min(-.015, data_min - .05 * data_span)
x_hi <- data_max + .12 * data_span

fig7 <- ggplot() +
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
    # "Alternative orientation" only ever shows up as a tick (geom_segment),
    # so give it the "|" key (124) in the legend instead of a dot that's
    # never actually plotted
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

save_ms_figure(fig7, "figures/fig7_uncertainty_pub.pdf",
               c(width = 172, height = 108))
message("saved: figures/fig7_uncertainty_pub.pdf (x-axis ", round(x_lo,3), " to ", round(x_hi,3), ")")
