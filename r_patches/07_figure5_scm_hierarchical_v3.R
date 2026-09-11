# figure5_scm_hierarchical.R
# Working SCM diagram (qgraph).
#
# Reads the finalized 16-edge SCM specification from:
#   pipeline_outputs/scm_edges_finalized.csv
#
# Layout:
# qgraph's automatic layout produced several edge/node collisions, so the
# coordinates below were set manually. The baseline working SCM is recursive,
# and the layout places nodes roughly from upstream to downstream while keeping
# crossing edges away from node labels.
#
# Current baseline SCM:
# N = 870, lavaan MLR
# CFI = .978
# TLI = .958
# RMSEA = .090, 90% CI [.077, .104]
# SRMR = .046
#
# Directional-evidence display:
#   solid  = direction met the data-support criterion
#   dashed = one of the four relationships with weak or inconsistent
#            directional evidence, varied in the directional-sensitivity
#            analysis
#
# `substantive` is retained below only because it is the legacy internal tier
# name in scm_edges_finalized.csv. It should not be interpreted as meaning
# that the direction was chosen from substantive theory.
#
# Edge widths use the current standardized path coefficients. A minimum
# plotted width of .20 is applied so that weaker retained paths remain visible.
#
# Current coefficients are from the clean baseline SCM refit on 2026-09-11.

suppressPackageStartupMessages({
  library(qgraph)
  library(dplyr)
  library(tibble)
})

source("r_patches/03_figure_style.R")


# =============================================================================
# Settings
# =============================================================================

PREVIEW_MODE <- FALSE  # TRUE only for checking layout with missing betas


# =============================================================================
# Nodes and layout
# =============================================================================

SCM_NODES <- c(
  "politics",
  "belief_concern",
  "harm_future",
  "trust_science",
  "harm_present",
  "policy_support",
  "weather_risk_prep",
  "social_norms",
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
  
  m <- as.matrix(
    scm_layout_coords[
      match(node_order, scm_layout_coords$name),
      c("x", "y")
    ]
  )
  
  rownames(m) <- node_order
  m
}


# =============================================================================
# Current baseline standardized path coefficients
# =============================================================================

BETA_TR <- tibble::tribble(
  ~from,                ~to,                   ~beta,
  
  "belief_concern",     "politics",             .530,
  
  "belief_concern",     "harm_future",          .856,
  
  "belief_concern",     "harm_present",         .291,
  "harm_future",        "harm_present",         .636,
  
  "belief_concern",     "trust_science",        .589,
  "harm_future",        "trust_science",        .238,
  
  "belief_concern",     "policy_support",       .306,
  "politics",           "policy_support",       .113,
  "trust_science",      "policy_support",       .424,
  "social_norms",       "policy_support",       .135,
  
  "trust_science",      "social_norms",         .690,
  
  "harm_present",       "weather_risk_prep",    .426,
  "belief_concern",     "weather_risk_prep",    .113,
  
  "harm_present",       "climate_behavior",     .295,
  "weather_risk_prep",  "climate_behavior",     .186,
  "social_norms",       "climate_behavior",     .127
)

stopifnot(nrow(BETA_TR) == 16)


# =============================================================================
# Read finalized edge audit
# =============================================================================

SCM_AUDIT_CSV <- "pipeline_outputs/scm_edges_finalized.csv"

if (!file.exists(SCM_AUDIT_CSV)) {
  stop(
    "Can't find ", SCM_AUDIT_CSV,
    " - run 18_finalize_scm_specification_v4.R first."
  )
}

audited <- read.csv(
  SCM_AUDIT_CSV,
  stringsAsFactors = FALSE
)


# =============================================================================
# Audit edge set before plotting
# =============================================================================

flagged <- audited |>
  filter(final_tier %in% c(
    "FLAG_WEAK_EXISTENCE",
    "FLAG_UNCLASSIFIED"
  ))

if (nrow(flagged) > 0) {
  
  print(
    as.data.frame(
      flagged[, c(
        "from",
        "to",
        "p_adjacent",
        "final_tier"
      )]
    ),
    row.names = FALSE
  )
  
  stop(
    nrow(flagged),
    " edge(s) flagged by the audit; resolve them before plotting."
  )
}


# Make sure audit and coefficient table contain exactly the same directed edges

missing_beta_edges <- audited |>
  anti_join(
    BETA_TR,
    by = c("from", "to")
  )

extra_beta_edges <- BETA_TR |>
  anti_join(
    audited,
    by = c("from", "to")
  )

if (
  nrow(missing_beta_edges) > 0 ||
  nrow(extra_beta_edges) > 0
) {
  
  cat("\nEdges in audit but not coefficient table:\n")
  print(
    as.data.frame(
      missing_beta_edges[, c("from", "to")],
      row.names = FALSE
    )
  )
  
  cat("\nEdges in coefficient table but not audit:\n")
  print(
    as.data.frame(
      extra_beta_edges[, c("from", "to")],
      row.names = FALSE
    )
  )
  
  stop(
    "SCM audit and coefficient table do not contain the same directed edge set."
  )
}


# =============================================================================
# Merge audit information with coefficients
# =============================================================================

scm_edges <- audited |>
  select(
    from,
    to,
    current_evidence,
    p_adjacent,
    final_tier,
    mismatch,
    note
  ) |>
  inner_join(
    BETA_TR,
    by = c("from", "to")
  ) |>
  mutate(
    evidence_plot = final_tier
  )

stopifnot(nrow(scm_edges) == 16)


# =============================================================================
# Missing-beta check
# =============================================================================

na_beta <- scm_edges |>
  filter(is.na(beta))

if (
  nrow(na_beta) > 0 &&
  !PREVIEW_MODE
) {
  
  print(
    as.data.frame(
      na_beta[, c(
        "from",
        "to",
        "final_tier",
        "note"
      )]
    ),
    row.names = FALSE
  )
  
  stop(
    nrow(na_beta),
    " edge(s) still missing a real beta. ",
    "Refit and add the coefficient, or set PREVIEW_MODE <- TRUE ",
    "only to inspect the layout."
  )
}

if (
  nrow(na_beta) > 0 &&
  PREVIEW_MODE
) {
  
  scm_edges$beta[
    is.na(scm_edges$beta)
  ] <- 0.15
}


# =============================================================================
# Print any audit discrepancies / notes
# =============================================================================

if (any(
  scm_edges$mismatch %in% TRUE,
  na.rm = TRUE
)) {
  
  cat(
    "Edges where the audit disagrees with an earlier tag ",
    "(plotting with the finalized audit tier):\n"
  )
  
  print(
    as.data.frame(
      scm_edges[
        scm_edges$mismatch %in% TRUE,
        c(
          "from",
          "to",
          "current_evidence",
          "final_tier"
        )
      ]
    ),
    row.names = FALSE
  )
}


notes_present <- scm_edges |>
  filter(
    !is.na(note),
    nzchar(note)
  )

if (nrow(notes_present) > 0) {
  
  cat(
    "Edges with audit notes worth checking ",
    "(not represented by line style alone):\n"
  )
  
  for (i in seq_len(nrow(notes_present))) {
    
    cat(
      sprintf(
        "  %s -> %s: %s\n",
        notes_present$from[i],
        notes_present$to[i],
        notes_present$note[i]
      )
    )
  }
}


# =============================================================================
# Edge weights
# =============================================================================

# A minimum width of .20 is used so weaker real paths
# (.113-.186) remain visually legible.

W_scm <- matrix(
  0,
  n_scm,
  n_scm,
  dimnames = list(
    SCM_NODES,
    SCM_NODES
  )
)

for (k in seq_len(nrow(scm_edges))) {
  
  f <- scm_edges$from[k]
  t <- scm_edges$to[k]
  b <- scm_edges$beta[k]
  
  W_scm[f, t] <- max(
    abs(b),
    0.20
  )
}


# =============================================================================
# Line types
# =============================================================================

# Legacy audit tier names:
#
# data_aligned = direction met the data-support criterion
# substantive  = weak/inconsistent direction varied in sensitivity analysis
#
# "substantive" is only an internal legacy label.

evidence_lty_code <- c(
  data_aligned = 1L,
  substantive  = 2L
)


lty_scm <- matrix(
  1L,
  n_scm,
  n_scm,
  dimnames = list(
    SCM_NODES,
    SCM_NODES
  )
)


for (k in seq_len(nrow(scm_edges))) {
  
  tier <- scm_edges$evidence_plot[k]
  
  if (!tier %in% names(evidence_lty_code)) {
    stop(
      "Unexpected final_tier for ",
      scm_edges$from[k],
      " -> ",
      scm_edges$to[k],
      ": ",
      tier
    )
  }
  
  lty_scm[
    scm_edges$from[k],
    scm_edges$to[k]
  ] <- evidence_lty_code[[tier]]
}


evidence_breaks <- c(
  "data_aligned",
  "substantive"
)

evidence_legend_labels <- c(
  "Data-supported direction",
  "Varied in directional sensitivity"
)

evidence_legend_lty <- unname(
  evidence_lty_code[
    evidence_breaks
  ]
)


# =============================================================================
# Figure title
# =============================================================================

fig_title <- if (
  PREVIEW_MODE &&
  nrow(na_beta) > 0
) {
  
  "Working structural causal model [PREVIEW]"
  
} else {
  
  "Working structural causal model"
}


# =============================================================================
# Plot
# =============================================================================

plot_scm_qgraph <- function() {
  
  layout(
    matrix(
      c(1, 2),
      nrow = 2
    ),
    heights = c(
      0.88,
      0.12
    )
  )
  
  
  # Main graph
  
  par(
    mar = c(1, 2, 6, 2),
    family = "sans"
  )
  
  qgraph::qgraph(
    W_scm,
    directed = TRUE,
    layout = scm_layout_matrix(SCM_NODES),
    
    labels = node_labels[SCM_NODES],
    color = node_fill_for(SCM_NODES),
    
    theme = "colorblind",
    
    edge.color = "#244E68",
    posCol = "#4A4A4A",
    negCol = "#4A4A4A",
    
    lty = lty_scm,
    
    curveAll = FALSE,
    curve = 1.5,
    
    fade = TRUE,
    
    mar = c(2, 2, 4, 2),
    
    esize = 6,
    asize = 4.4,
    cut = 0,
    
    vsize = 10,
    label.cex = .8,
    
    title = fig_title,
    title.cex = 1.2
  )
  
  
  # Legend
  
  par(
    mar = c(0, 0, 0, 0)
  )
  
  plot.new()
  
  legend(
    "center",
    horiz = TRUE,
    
    legend = evidence_legend_labels,
    lty = evidence_legend_lty,
    
    lwd = 2,
    col = "#244E68",
    
    bty = "n",
    cex = 1,
    seg.len = 2.6
  )
}


# =============================================================================
# Save
# =============================================================================

out_path <- if (
  PREVIEW_MODE &&
  nrow(na_beta) > 0
) {
  
  "figures/fig5_scm_hierarchical_PREVIEW.pdf"
  
} else {
  
  "figures/fig5_scm_hierarchical.pdf"
}


save_ms_basegraphics(
  plot_scm_qgraph,
  out_path,
  c(
    width = 170,
    height = 129
  )
)

message(
  "saved: ",
  out_path
)
