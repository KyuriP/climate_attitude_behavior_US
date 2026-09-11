# figure5_scm_hierarchical.R
# working SCM diagram (qgraph). reads scm_edges_finalized.csv from script 18.
#
# layout: qgraph's default layout kept crossing edges over labels and running
# edges through nodes. tried reusing the Rgraphviz layout from panel B but the
# coordinates weren't collision-free once fed into qgraph's vsize (harm_future
# and trust_science overlapped). ended up just building coords by hand: since
# the SCM is a DAG, layer each node by its longest path from a root so every
# edge points left->right, then space nodes >=1.3 apart within a layer and
# 1.6 apart between layers. checked by hand that no edge line comes within
# ~0.5-0.9 of a node it's not touching. this is SCM-only, doesn't touch
# fixed_layout_coords in 03_figure_style.R (that one's shared with fig 2 and
# the supp causal-discovery panels, don't want to break those).
#
# other stuff from the last styling pass: node colors pulled lighter/pastel
# to match fig 2, dropped "(worry)" from the weather-risk label everywhere,
# dropped the node-color legend (redundant with the node fill + label), bumped
# legend cex.
#
# two edges (belief_concern->weather_risk_prep, harm_future->trust_science)
# got added to the skeleton after adding a second predictor to those nodes in
# the lavaan model - coefficients below are from that refit (N=870, MLR).
# global fit for the 16-edge model, ORIGINAL edges: CFI=.979, TLI=.962,
# RMSEA=.086 [.074,.100], SRMR=.039. After the 2026-09-11 reversal of
# politics/belief_concern and policy_support/social_norms (see below and
# analysis_decisions_log.md Section 34): CFI=.978, TLI=.959, RMSEA=.089,
# SRMR=.046 -- essentially unchanged, no fit degradation.
#
# UPDATED 2026-09-11, two separate changes:
#  (1) belief_concern->politics and social_norms->policy_support now follow
#      the bootstrap asymmetry (was -.336 / -.325 for the old politics->
#      belief_concern / policy_support->social_norms direction, consistently
#      negative at both alphas) instead of the old theory-asserted direction.
#      These render SOLID now (data_aligned) -- BETA_TR below still has
#      placeholder NA for their coefficients pending a refit, see comment there.
#  (2) evidence_plot is read straight from scm_edges_finalized.csv's
#      final_tier, no per-edge override needed anymore: 00_config.R's
#      ORIENTATION_ASYMMETRY_EPS is now a .10 magnitude band (was .01,
#      sign-only), so belief_concern->weather_risk_prep (+.065),
#      harm_present->weather_risk_prep (+.078), and harm_future->harm_present
#      (+.076 pooled, but sign-flips per-alpha -- see 21/27 in r_patches)
#      all correctly fall out as "substantive" (dashed) from 04_scm_finalize.R's
#      audit directly, same as social_norms->climate_behavior always has.
#
# two line types - solid = bootstrap-supported (in the asserted direction),
# dashed = direction picked on substantive grounds or genuinely too weak/
# inconsistent to call.

suppressPackageStartupMessages({
  library(qgraph)
  library(dplyr)
  library(tibble)
})
source("r_patches/03_figure_style.R")

PREVIEW_MODE <- FALSE  # TRUE only for eyeballing layout before real betas are in

SCM_NODES <- c(
  "politics", "belief_concern", "harm_future", "trust_science",
  "harm_present", "policy_support", "weather_risk_prep", "social_norms",
  "climate_behavior"
)
n_scm <- length(SCM_NODES)

scm_layout_coords <- tibble::tribble(
  ~name,                 ~x,   ~y,
  "politics",            0.0, -0.6,
  "belief_concern",      1.6,  0.0,
  "harm_future",         3.2,  1.3,
  "harm_present",        3.2, -1.3,
  "trust_science",       4.8,  1.3,
  "policy_support",      4.8, -0.4,
  "weather_risk_prep",   6.4, -1.6,
  "social_norms",        6.4,  0.9,
  "climate_behavior",    8.0, -0.4
)

scm_layout_matrix <- function(node_order) {
  m <- as.matrix(scm_layout_coords[match(node_order, scm_layout_coords$name), c("x", "y")])
  rownames(m) <- node_order
  m
}

# path coefficients, standardized. first 14 rows straight from the supp
# table; last two are the refit values for the two newly-added edges.
BETA_TR <- tibble::tribble(
  ~from,                ~to,                   ~beta,
  "belief_concern",     "politics",             .530,  # refit 2026-09-11 (N=870, MLR); politics ~ belief_concern, single predictor -- see analysis_decisions_log.md Section 34
  "belief_concern",     "harm_future",          .856,
  "belief_concern",     "harm_present",         .291,
  "harm_future",        "harm_present",         .636,
  "belief_concern",     "trust_science",        .589,  # was .793, dropped once harm_future added as 2nd predictor
  "belief_concern",     "policy_support",       .335,
  "politics",           "policy_support",       .113,
  "trust_science",      "policy_support",       .491,
  "trust_science",      "social_norms",         .460,
  "social_norms",       "policy_support",       .135,  # refit 2026-09-11 (N=870, MLR); policy_support ~ belief_concern + trust_science + politics + social_norms
  "harm_present",       "weather_risk_prep",    .426,  # was .520, dropped once belief_concern added as 2nd predictor
  "harm_present",       "climate_behavior",     .311,
  "weather_risk_prep",  "climate_behavior",     .186,
  "social_norms",       "climate_behavior",     .114,
  "belief_concern",     "weather_risk_prep",    .113,  # SE=.050, Z=2.272, p=.023
  "harm_future",        "trust_science",        .238   # SE=.042, Z=5.720, p<.001
)
stopifnot(nrow(BETA_TR) == 16)

SCM_AUDIT_CSV <- "pipeline_outputs/scm_edges_finalized.csv"
if (!file.exists(SCM_AUDIT_CSV)) {
  stop("can't find ", SCM_AUDIT_CSV, " - run 18_finalize_scm_specification_v4.R first")
}
audited <- read.csv(SCM_AUDIT_CSV, stringsAsFactors = FALSE)

flagged <- audited |> filter(final_tier %in% c("FLAG_WEAK_EXISTENCE", "FLAG_UNCLASSIFIED"))
if (nrow(flagged) > 0) {
  print(as.data.frame(flagged[, c("from", "to", "p_adjacent", "final_tier")]), row.names = FALSE)
  stop(nrow(flagged), " edge(s) flagged by the audit, sort those before plotting")
}

scm_edges <- audited |>
  select(from, to, current_evidence, p_adjacent, final_tier, mismatch, note) |>
  inner_join(BETA_TR, by = c("from", "to")) |>
  mutate(evidence_plot = final_tier)
stopifnot(nrow(scm_edges) == 16)

na_pairs <- scm_edges |> filter(is.na(beta)) |> select(from, to)
na_beta  <- scm_edges |> filter(is.na(beta))
if (nrow(na_beta) > 0 && !PREVIEW_MODE) {
  print(as.data.frame(na_beta[, c("from", "to", "final_tier", "note")]), row.names = FALSE)
  stop(nrow(na_beta), " edge(s) still missing a real beta - refit and paste the value into BETA_TR, or set PREVIEW_MODE <- TRUE to just check the layout")
}
if (nrow(na_beta) > 0 && PREVIEW_MODE) {
  scm_edges$beta[is.na(scm_edges$beta)] <- 0.15  # placeholder width, not a real coefficient
}

if (any(scm_edges$mismatch %in% TRUE)) {
  cat("edges where the audit disagrees with what I'd tagged before (plotting with the audited tier):\n")
  print(as.data.frame(scm_edges[scm_edges$mismatch %in% TRUE,
                                 c("from", "to", "current_evidence", "final_tier")]),
        row.names = FALSE)
}

notes_present <- scm_edges |> filter(nzchar(note))
if (nrow(notes_present) > 0) {
  cat("edges with a caveat worth mentioning in the caption/text (not shown by line style alone):\n")
  for (i in seq_len(nrow(notes_present))) {
    cat(sprintf("  %s -> %s: %s\n", notes_present$from[i], notes_present$to[i], notes_present$note[i]))
  }
}

# edge width floored at .20 so the weak real edges (.113-.186) don't
# disappear now that we're not printing numeric labels anymore
W_scm <- matrix(0, n_scm, n_scm, dimnames = list(SCM_NODES, SCM_NODES))
for (k in seq_len(nrow(scm_edges))) {
  f <- scm_edges$from[k]; t <- scm_edges$to[k]; b <- scm_edges$beta[k]
  W_scm[f, t] <- max(abs(b), 0.20)
}

# two tiers: data_aligned = solid, substantive = dashed
evidence_lty_code <- c(data_aligned = 1L, substantive = 2L)
lty_scm <- matrix(1L, n_scm, n_scm, dimnames = list(SCM_NODES, SCM_NODES))
for (k in seq_len(nrow(scm_edges))) {
  lty_scm[scm_edges$from[k], scm_edges$to[k]] <-
    evidence_lty_code[[scm_edges$evidence_plot[k]]]
}

evidence_breaks <- c("data_aligned", "substantive")
evidence_legend_labels <- c(
  "Data-supported orientation",
  "Data-uncertain / theory-completed"
)
evidence_legend_lty <- unname(evidence_lty_code[evidence_breaks])

fig_title <- if (PREVIEW_MODE && nrow(na_beta) > 0) {
  "Working structural causal model [PREVIEW]"
} else {
  "Working structural causal model"
}

plot_scm_qgraph <- function() {
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.88, 0.12))

  par(mar = c(1, 2, 6, 2), family = "sans")
  qgraph::qgraph(
    W_scm, directed = TRUE, layout = scm_layout_matrix(SCM_NODES),
    labels = node_labels[SCM_NODES], color = node_fill_for(SCM_NODES),
    theme = "colorblind", edge.color = "#244E68", posCol = "#4A4A4A", negCol = "#4A4A4A",
    lty = lty_scm,
    curveAll = FALSE, curve = 1.5,
    fade = TRUE, mar = c(2, 2, 4, 2),
    esize = 6, asize = 4.4, cut = 0, vsize = 10, label.cex = .8,
    title = fig_title, title.cex = 1.2
  )

  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, legend = evidence_legend_labels, lty = evidence_legend_lty,
         lwd = 2, col = "#244E68", bty = "n", cex = 1, seg.len = 2.6)
}

out_path <- if (PREVIEW_MODE && nrow(na_beta) > 0) {
  "figures/fig5_scm_hierarchical_PREVIEW.pdf"
} else {
  "figures/fig5_scm_hierarchical.pdf"
}
save_ms_basegraphics(plot_scm_qgraph, out_path, c(width = 170, height = 129))
message("saved: ", out_path)
