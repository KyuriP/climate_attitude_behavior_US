# 35_figure_cyclic_feedback_equilibrium.R
#
# Supplementary figure for \label{supp:feedback} (Section titled "Feedback
# extension for cyclic orientation specifications" -- deliberately not
# hard-coded as "S15": Supplementary Section S9 is being deleted, which
# shifts every later section's number, so this section is referenced by
# label only until the final numbering is settled -- see
# analysis_decisions_log.md Section 43).
#
# Per Kyuri's spec (2026-09-11): one row per intervention node, diamonds =
# baseline working SCM, light circles = the other 11 acyclic directional
# specifications, triangles = the 4 cyclic specifications analyzed as
# stable linear equilibrium systems (clean_pipeline/19_cyclic_feedback_
# equilibrium.R). Reads that script's already-written unified table
# (pipeline_outputs/tables/orientation_enumeration_ate_all16_with_
# equilibrium.csv -- 07's 12 acyclic rows + 19's 4 equilibrium rows, one
# row per node x specification, tagged combo_type). Deliberately simple
# (3 shape categories, no politics-distinctness styling like fig7 uses) --
# this is a compact supplementary figure, not the main-text one.
#
# Kept as a plain jitter for the acyclic points, same reasoning as
# 10_figure7_uncertainty_pub.R's header: a handful of discrete points isn't
# enough to support a smoothed density without implying a continuous
# distribution that isn't there.
#
# Output: figures/fig_supp_cyclic_feedback.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})
source("r_patches/03_figure_style.R")

ate_path <- file.path("pipeline_outputs", "tables", "orientation_enumeration_ate_all16_with_equilibrium.csv")
if (!file.exists(ate_path)) {
  stop("can't find ", ate_path, " -- run clean_pipeline/07_intervention_ates_8node.R ",
       "and clean_pipeline/19_cyclic_feedback_equilibrium.R first (in that order).")
}
all16 <- utils::read.csv(ate_path, stringsAsFactors = FALSE)
stopifnot(all(c("scenario", "node", "ate", "combo_type") %in% names(all16)))
stopifnot(setequal(unique(all16$combo_type), c("baseline", "acyclic", "cyclic_equilibrium")))

# Same row order as 10_figure7_uncertainty_pub.R's TARGET_ORDER (that
# script defines it locally rather than in 03_figure_style.R, so it's
# repeated here rather than depended on cross-script) and the same
# one-line label set (node_labels_oneline, from 03_figure_style.R) --
# this is a companion figure to fig7, should read consistently with it.
node_targets <- c(
  "harm_present", "belief_concern", "harm_future", "weather_risk_prep",
  "social_norms", "trust_science", "policy_support", "politics"
)
stopifnot(setequal(unique(all16$node), node_targets))

plot_df <- all16 |>
  dplyr::mutate(
    node_label = factor(node_labels_oneline[node], levels = rev(node_labels_oneline[node_targets])),
    kind = factor(
      dplyr::recode(combo_type,
        baseline = "Baseline specification",
        acyclic = "Other acyclic specification",
        cyclic_equilibrium = "Cyclic (equilibrium) specification"
      ),
      levels = c("Baseline specification", "Other acyclic specification",
                 "Cyclic (equilibrium) specification")
    )
  )

baseline_df <- plot_df |> dplyr::filter(kind == "Baseline specification")
acyclic_df  <- plot_df |> dplyr::filter(kind == "Other acyclic specification")
cyclic_df   <- plot_df |> dplyr::filter(kind == "Cyclic (equilibrium) specification")
stopifnot(nrow(baseline_df) == length(node_targets))
stopifnot(nrow(acyclic_df) == 11 * length(node_targets))
stopifnot(nrow(cyclic_df) <= 4 * length(node_targets))  # <= in case a combo had no stable equilibrium at all -- see 19's header

data_min <- min(plot_df$ate, na.rm = TRUE)
data_max <- max(plot_df$ate, na.rm = TRUE)
data_span <- data_max - data_min
x_lo <- min(-.015, data_min - .05 * data_span)
x_hi <- data_max + .12 * data_span

fig_cyclic <- ggplot() +
  geom_vline(xintercept = 0, linewidth = .45, linetype = "dashed", colour = "#9A9A9A") +
  geom_point(
    data = acyclic_df,
    aes(x = ate, y = node_label, shape = kind, fill = kind, colour = kind),
    size = 2.4, stroke = 0, alpha = .28,
    position = position_jitter(height = .16, width = 0.006, seed = 20260911)
  ) +
  geom_point(
    data = cyclic_df,
    aes(x = ate, y = node_label, shape = kind, fill = kind, colour = kind),
    size = 2.6, stroke = .4, alpha = .8,
    position = position_jitter(height = .10, width = 0.004, seed = 20260911)
  ) +
  geom_point(
    data = baseline_df,
    aes(x = ate, y = node_label, shape = kind, fill = kind, colour = kind),
    size = 2.6, stroke = .5, alpha = .8
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Baseline specification" = 23,
               "Other acyclic specification" = 21,
               "Cyclic (equilibrium) specification" = 24)
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Baseline specification" = COL$rust,
               "Other acyclic specification" = COL$blue_dark,
               "Cyclic (equilibrium) specification" = COL$ink_mid)
  ) +
  scale_colour_manual(
    name = NULL,
    values = c("Baseline specification" = COL$ink,
               "Other acyclic specification" = COL$blue_dark,
               "Cyclic (equilibrium) specification" = COL$ink)
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 3.0, stroke = 0.4, alpha = 0.85)),
    fill = guide_legend(override.aes = list(size = 3.0, stroke = 0.4, alpha = 0.85)),
    colour = guide_legend(override.aes = list(size = 3.0, stroke = 0.4, alpha = 0.85))
  ) +
  scale_x_continuous(
    limits = c(x_lo, x_hi), breaks = scales::pretty_breaks(n = 5)(c(x_lo, x_hi)),
    expand = expansion(mult = c(0, .01))
  ) +
  labs(x = "Model-implied intervention effect on climate behavior (SD)", y = NULL) +
  theme_pub +
  theme(
    text = element_text(size = 10.5),
    axis.text.y = element_text(size = 10.2),
    axis.text.x = element_text(size = 9.2),
    axis.title.x = element_text(size = 9.8),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = COL$grid, linewidth = .35),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8.2),
    plot.margin = margin(6, 16, 4, 4)
  )

save_ms_figure(fig_cyclic, "figures/fig_supp_cyclic_feedback.pdf", c(width = 172, height = 132))
message("Saved supplementary feedback-extension figure: figures/fig_supp_cyclic_feedback.pdf ",
        "(x-axis: ", round(x_lo, 3), " to ", round(x_hi, 3), ")")
