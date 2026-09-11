# =============================================================================
# 03_figure_style.R
#
# Shared visual system for every main-text and supplementary figure. Source
# this file at the top of every other 03x_figure_*.R script.
#
# Adopts an explicit global theme (theme_pub), explicit mm-based final
# dimensions per figure, and a tight font-size spec (9pt axis/node text,
# 10.5pt panel titles, 8pt legends). Everything below matches that spec
# exactly.
#
# PACKAGES: this file only needs ggplot2. Individual figure scripts load
# whatever else they need (ggraph, tidygraph, qgraph, patchwork, ggrepel).
#
# COL (the shared palette object), PT() (a pt-to-mm size helper for
# geom_text), and FIG_FONT below are what figures 1, 2/3, 4, 6, 7, and the
# supplementary causal-discovery panels depend on. The palette was
# reconstructed by sampling the actual rendered colors out of the existing
# fig1/fig6/fig7 PDFs (pixel-exact for blue_dark, blue_light, panel, ink,
# grid, rust), plus the hardcoded fallback hex values already used elsewhere
# in the pipeline (ink_mid/ink_light, which match #4A4A4A/#9A9A9A already
# hardcoded in 10_figure7_uncertainty_pub.R). PT() uses ggplot2's own .pt
# constant, the standard definition. Sanity-check the reconstructed swatches
# against the originally intended values before relying on them for a final
# render.
# ==============================================================================

suppressPackageStartupMessages(library(ggplot2))

# -----------------------------------------------------------------------------
# 0. Shared palette (COL) and small helpers used across figure scripts.
#    Reconstructed from the rendered PDFs -- see note above.
# -----------------------------------------------------------------------------
COL <- list(
  blue_dark  = "#244E68",  # sampled from fig1/fig7 (density stroke, diamond fill, gradient high end)
  blue_light = "#DDEAF2",  # sampled from fig1 (density fill)
  rust       = "#A75A43",  # sampled from fig6 ("Combined" bars)
  ink        = "#202020",  # sampled from fig7 (axis/category text)
  ink_mid    = "#4A4A4A",  # matches the hardcoded fallback already used in 10_figure7_uncertainty_pub.R
  ink_light  = "#9A9A9A",  # matches the hardcoded fallback already used in 10_figure7_uncertainty_pub.R
  grid       = "#E8E8E8",  # sampled from fig7 (vertical gridlines)
  panel      = "#F6F6F6"   # sampled from fig1 (diagonal panel background)
)

# Font family used throughout (matches theme_pub's text family below).
FIG_FONT <- "Arial"

# ggplot2's geom_text()/geom_label() "size" aesthetic is in mm, not points --
# PT(x) converts a point size (e.g. "8.4pt") to the mm value ggplot expects,
# using ggplot2's own .pt constant (72.27/25.4). This is the standard
# conversion, not a stylistic choice.
PT <- function(x) x / .pt

# -----------------------------------------------------------------------------
# 1. Global theme (project visual spec, verbatim)
# -----------------------------------------------------------------------------
# NOTE on "Arial": base R / ggplot2 don't ship Arial. If it's not installed
# as a system font, this silently falls back to the device default (usually
# Helvetica/Nimbus Sans on Linux, which is metrically near-identical) --
# fine for a first pass; for true Arial, load it via the `extrafont` or
# `systemfonts` package before running these scripts. "Arial" is kept here
# per the project's explicit font choice, and most journal production
# systems substitute an equivalent sans anyway.
theme_pub <- theme_classic(base_size = 9) +
  theme(
    text            = element_text(family = "Arial"),
    axis.title      = element_text(size = 9),
    axis.text       = element_text(size = 8.5),
    plot.title      = element_text(size = 10.5, face = "bold"),
    plot.subtitle   = element_text(size = 8.5),
    legend.title    = element_text(size = 8.5),
    legend.text     = element_text(size = 8),
    legend.position = "right",
    plot.margin     = margin(5, 7, 5, 5)
  )

# Kept as an alias so any earlier-written code calling theme_ms() still
# works without edits -- theme_ms is now just theme_pub.
theme_ms <- function(base_size = 9) theme_pub

# -----------------------------------------------------------------------------
# 2. Node-family color system (used in Figures 2, 3, 5/SCM, 6, and 7)
#
#    These 5 hexes are copied verbatim from the .qmd's node_colors
#    (Section 1): every figure that colors nodes by family -- Figure 2's
#    GGM, Figure 3 Panel B's FCI-JCI PAG (plotAG, via node_colors directly),
#    and Figure 5's SCM (via node_fill_for() below) -- draws from this one
#    identical palette, so a node reads the same color everywhere it appears.
# -----------------------------------------------------------------------------
# Lightened further toward white (more pastel, matching the GGM look) --
# same hues throughout, blended ~35% toward white so fills read as soft
# pastel; kept byte-for-byte identical to node_colors in the .qmd (see that
# file) since the two must match across Fig 2 / SCM / FCI panels.
node_family_colors <- c(
  "core_attitude" = "#D5ECF9",  # belief_concern, harm_present, harm_future,
                                 # trust_science, policy_support
  "politics"      = "#D5F1E2",
  "weather_risk"  = "#FCEBD1",
  "social_norms"  = "#EDDEF2",
  "behavior"      = "#F8D7D2"
)

node_family <- c(
  belief_concern     = "core_attitude",
  harm_present       = "core_attitude",
  harm_future        = "core_attitude",
  trust_science      = "core_attitude",
  policy_support     = "core_attitude",
  politics           = "politics",
  weather_risk_prep  = "weather_risk",
  social_norms       = "social_norms",
  climate_behavior   = "behavior"
)

node_labels <- c(
  belief_concern     = "Belief/\nConcern",
  harm_present       = "Present\nHarm",
  harm_future        = "Future\nHarm",
  trust_science      = "Trust\nScience",
  policy_support     = "Policy\nSupport",
  politics           = "Politics",
  weather_risk_prep  = "Weather\nrisk",
  social_norms       = "Social\nNorms",
  climate_behavior   = "Climate\nBehavior"
)

node_labels_oneline <- c(
  belief_concern     = "Belief/concern",
  harm_present       = "Present harm",
  harm_future        = "Future harm",
  trust_science      = "Trust in science",
  policy_support     = "Policy support",
  politics           = "Political orientation",
  weather_risk_prep  = "Weather risk",
  social_norms       = "Social norms",
  climate_behavior   = "Climate behavior"
)

# Short/abbreviated labels for tight axes (e.g. Figure 1's angled x-axis,
# Figure 4's stability-matrix axes). Reconstructed by reading the exact
# abbreviations off the already-rendered fig1_heatmap_construct.pdf x-axis,
# so these match the approved rendering rather than a fresh guess.
node_labels_short <- c(
  belief_concern     = "Belief",
  harm_present       = "Present harm",
  harm_future        = "Future harm",
  trust_science      = "Trust",
  policy_support     = "Policy",
  politics           = "Politics",
  weather_risk_prep  = "Weather risk",
  social_norms       = "Norms",
  climate_behavior   = "Behavior"
)

node_fill_for <- function(node_names) unname(node_family_colors[node_family[node_names]])

# -----------------------------------------------------------------------------
# 3. FIXED node layout shared across Figures 2, 3, and the supplementary
#    causal-discovery graphs (requirement: these should read as the same
#    figure before and after adding behavior, with identical coordinates
#    across all four supplementary panels). One set of coordinates, reused
#    everywhere a node appears, so comparing panels never requires
#    re-orienting.
#
#    Layout (matching the project's reference sketch):
#        Present harm        Future harm
#               Belief / concern
#        Trust science       Policy support
#        Social norms        Political orientation
#               Weather risk / prep
#        (Climate behavior, added below weather_risk_prep, Fig. 3 + supp only)
# -----------------------------------------------------------------------------
fixed_layout_coords <- tibble::tribble(
  ~name,                ~x,  ~y,
  "harm_present",       -1,   2,
  "harm_future",         1,   2,
  "belief_concern",      0,   1,
  "trust_science",      -1,   0,
  "policy_support",      1,   0,
  "social_norms",       -1,  -1,
  "politics",            1,  -1,
  "weather_risk_prep",   0,  -2,
  "climate_behavior",    0,  -3
)

# qgraph wants a plain numeric matrix (rows in the same order as its input
# correlation/weight matrix), not a data frame -- this helper builds that
# from fixed_layout_coords for whatever node subset/order you pass in.
fixed_layout_matrix <- function(node_order) {
  m <- as.matrix(fixed_layout_coords[match(node_order, fixed_layout_coords$name), c("x", "y")])
  rownames(m) <- node_order
  m
}

# -----------------------------------------------------------------------------
# 4. Final output dimensions, in mm (per the project's figure-size table) --
#    ggsave(..., units = "mm", device = cairo_pdf) in every figure script
#    now uses these.
# -----------------------------------------------------------------------------
FIG_DIMS_MM <- list(
  fig1 = c(width = 180, height = 120),
  # Taller/near-square, matching the original exploratory .qmd's own
  # fig.height=fig.width=7 (square) chunks for these two networks. The old
  # 124/155 values were sized for the shared fixed-grid layout (see
  # 05_figures2_3_ggm_redesign.R) and no longer apply -- fig2 and fig3 are
  # independent layouts now, so there's no ratio to preserve between them;
  # both are simply squared up.
  fig2 = c(width = 147, height = 147),
  fig3 = c(width = 147, height = 147),
  fig4 = c(width = 180, height = 105),
  # square and noticeably bigger than fig4 -- these are 8x8 stability
  # matrices with two-line cell labels, fig4's 180x105 landscape ratio was
  # too squat for that and made every cell cramped.
  fig_stability = c(width = 160, height = 160),
  fig5 = c(width = 180, height = 115),
  fig6 = c(width = 150, height = 100),
  fig7 = c(width = 160, height = 105),
  supp = c(width = 180, height = 180)  # 2x2 panel grid
)

save_ms_figure <- function(plot, filename, dims) {
  ggplot2::ggsave(
    filename = filename, plot = plot,
    width = dims["width"], height = dims["height"], units = "mm",
    device = grDevices::cairo_pdf
  )
}

# For base-graphics plotters (e.g. qgraph(), which draws directly to the
# open device rather than returning a ggplot object) -- takes a zero-arg
# function that performs the plotting call(s), opens a cairo_pdf device at
# the given mm dims, runs it, and closes the device. Used by
# 07_figure5_scm_hierarchical_v3.R (Figure 3's qgraph SCM diagram).
save_ms_basegraphics <- function(plot_fn, filename, dims) {
  grDevices::cairo_pdf(
    filename = filename,
    width = unname(dims["width"]) / 25.4,
    height = unname(dims["height"]) / 25.4
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_fn()
}

# Retained for any older code that still calls the old in-inch signature
# (width/height in inches) -- new figure scripts should call save_ms_figure
# with a FIG_DIMS_MM[["figN"]] entry instead.
FIG_WIDTH_IN <- 7.09       # 180mm, kept for reference only
NODE_LABEL_PT <- 9
AXIS_PT <- 9
TITLE_PT <- 10.5
