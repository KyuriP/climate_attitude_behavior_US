# =============================================================================
# 29_joint_pag_edgetype_audit.R  (v2, 2026-09-06)
#
# WHY THIS SCRIPT EXISTS
# -----------------------
# The scalar asymmetry statistic A_alpha = P(arrow at destination) - P(arrow
# at source), used everywhere so far (script 21's table, script 27's
# diagnostic, Table 2), collapses several genuinely different PAG relation
# types into one number. A_alpha near 0 can arise from two substantively
# different situations the MARGINAL per-endpoint mark proportions cannot
# tell apart:
#   - the edge is mostly bidirected (X <-> Y): arrowhead mass high at BOTH
#     ends -- consistent with latent confounding / neither endpoint an
#     ancestor of the other. A specific, informative PAG relation.
#   - the edge is mostly circle-circle (X o-o Y): mass low/undetermined at
#     both ends -- genuinely uninformative.
# Telling these apart needs the JOINT mark at both ends WITHIN THE SAME
# bootstrap resample. Confirmed (by reading the .qmd's Section 7.4 bootstrap
# chunk directly) that only the pooled/marginal sums are ever saved to disk
# (mark_props_ext_by_alpha.rds, fci_props_ext.rds) -- the per-resample list
# that would contain the joint pairing exists only transiently in memory
# during that one render. So this script reruns ONLY the FCI half of the
# bootstrap (PC is skipped -- Table 2's classification is FCI-based),
# classifying the joint edge type per resample immediately (never storing a
# per-resample array). NOTHING else reruns: not PC, not the SCM, not the
# intervention bootstrap, not behavior/attrition/measurement analyses.
#
# v2 CHANGES (three corrections made before running this for real)
# ------------------------------------------------------------------------
# 1. FAIL LOUDLY on any endpoint mark combination outside the 7 canonical
#    PAG relations (i_to_j, j_to_i, bidirected, i_ocirc_arrow_j,
#    j_ocirc_arrow_i, circle_circle, no_edge). v1 silently folded anything
#    unmatched into an "other" bucket -- changed so classify_joint() now
#    calls stop() with the resample index, alpha, node pair, and both raw
#    mark values if this ever happens, halting the run rather than quietly
#    treating an unanticipated encoding as "unresolved". The classification
#    step itself sits OUTSIDE the tryCatch that wraps the pcalg::fci() call
#    (which legitimately can fail to fit and should be skipped/counted, not
#    treated as a mark-encoding problem), so a stop() here is never
#    accidentally swallowed by that tryCatch and propagates all the way up
#    through furrr::future_map, aborting the whole alpha-loop.
# 2. n_boot: corrected after checking the actual rendered pipeline output,
#    which prints "1000 resamples per alpha threshold" -- contradicting an
#    earlier (WRONG) claim in this file that the primary analysis used
#    n_boot=2500 throughout. Independently corroborated via file timestamps:
#    pipeline_outputs/fci_props_ext.rds, which Table 2 (via script 18) and
#    Figure 5 trace back to, is dated 2026-03-28 and was never regenerated
#    by any rerun since -- it reflects whatever n_boot was set to back then
#    (1,000), not the live .qmd's current n_boot=2500. Script 18 only
#    recomputes fci_props_ext from a fresh in-memory object if one already
#    exists in the session; otherwise it falls back to this frozen on-disk
#    file (see 18_finalize_scm_specification_v4.R, ~line 61) -- which is
#    what actually happened, since scm_edges_finalized.csv (2026-09-05)
#    predates the fresh n_boot=2500 rerun (mark_props_ext_by_alpha.rds,
#    2026-09-06). So: Table 2/Figure 5's published numbers are genuinely
#    n_boot=1,000-based, matching main2.tex's stated methodology exactly --
#    NO manuscript-text change needed there. This joint audit is
#    deliberately a SEPARATE, higher-precision targeted analysis at
#    n_boot_joint=2,500 (see below -- decoupled from the shared `n_boot`
#    variable entirely, precisely so this script's resample count can never
#    again be confused with, or silently drift together with, the primary
#    analysis's). Manuscript wording (adopted as-is): "Primary
#    causal-discovery stability was assessed with 1,000 bootstrap resamples
#    per conditional-independence threshold. Because the joint
#    endpoint-type audit was used to make the final distinction between
#    directional, bidirected, and unresolved PAG relations, we repeated
#    that targeted FCI audit with 2,500 resamples per threshold to obtain a
#    more stable characterization of the joint endpoint frequencies."
# 3. NO derived classification. v1 added `likely_confounded`/
#    `likely_genuinely_unresolved` boolean columns using an ad hoc 20%
#    threshold on top of the raw proportions -- exactly the kind of
#    mechanical collapse this audit exists to avoid. Removed. The primary
#    output (`joint_summary` / joint_pag_edgetype_audit.csv) is now nothing
#    but from, to, alpha, and the 7 raw joint-type percentages -- full PAG
#    semantics preserved, no verdict imposed. A separate `cross_check` table
#    is still produced (sorted by a plain, non-thresholded sum of bidirected
#    + circle-circle mass, purely so 16 edges are easy to scan) but carries
#    no TRUE/FALSE classification either -- it's a sort order, not a rule.
#
# HOW TO USE: run within the same .qmd session, after Section 7.4 (so
# agg_ext, node_order_ext, context_idx, alphas are already defined -- same
# prerequisites as scripts 21-28, minus n_boot, which this script no longer
# depends on -- see note 2 above).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(furrr)
  library(future)
})

stopifnot(exists("agg_ext"), exists("node_order_ext"), exists("context_idx"),
          exists("alphas"))

# DECOUPLED from the shared `n_boot` session variable (2026-09-06, after a
# real mix-up: the primary FCI/PC stability analysis behind Table 2/Figure 3
# actually used n_boot=1,000 per alpha -- confirmed directly from the
# actual rendered pipeline output, and independently corroborated here via
# file timestamps: pipeline_outputs/fci_props_ext.rds, which Table 2 traces
# back to, is dated 2026-03-28 and was never regenerated by any subsequent
# rerun, so it reflects whatever n_boot was set to back then -- 1,000, not
# the current n_boot=2500). This joint audit is a deliberately separate,
# higher-precision targeted analysis and should never silently inherit
# whatever the primary pipeline's n_boot happens to be set to at the time
# it's run -- hence its own explicit constant here, independent in both
# directions.
n_boot_joint <- 2500L
cat("Joint PAG edge-type audit uses its own n_boot_joint =", n_boot_joint,
    "per alpha, independent of the primary pipeline's n_boot (whatever that",
    "is currently set to) -- see header note 2 for why these are",
    "deliberately different, separately-described analyses.\n")

set.seed(42)  # same seed as the .qmd's ext-bootstrap chunk, for continuity

p_ext <- length(node_order_ext)

# The 16 retained edges, in the SAME (from, to) direction as base_edges in
# clean_pipeline/05_scm_intervention_helpers.R / the working SCM -- this
# is just a reference list for the final summary table (which edges to
# report on, and in which direction to label "from_to" vs "to_from"). The
# audit below actually computes all p*(p-1)/2 pairs, not just these 16, so
# nothing here restricts what the bootstrap itself estimates.
# belief_concern->politics and social_norms->policy_support updated
# 2026-09-11 to track base_edges after the orientation-rule reversal (see
# analysis_decisions_log.md Section 34) -- r_patches/02 (referenced by the
# old comment here) still holds the pre-reversal direction and is stale/
# deliberately unmaintained, so this file now points at 05 instead.
retained_edges <- tibble::tribble(
  ~from,               ~to,
  "belief_concern",    "politics",
  "belief_concern",    "harm_future",
  "belief_concern",    "harm_present",
  "harm_future",       "harm_present",
  "belief_concern",    "trust_science",
  "harm_future",       "trust_science",
  "belief_concern",    "policy_support",
  "trust_science",     "policy_support",
  "politics",          "policy_support",
  "social_norms",      "policy_support",
  "trust_science",     "social_norms",
  "harm_present",      "weather_risk_prep",
  "belief_concern",    "weather_risk_prep",
  "harm_present",      "climate_behavior",
  "social_norms",      "climate_behavior",
  "weather_risk_prep", "climate_behavior"
)
stopifnot(nrow(retained_edges) == 16)

# --- Joint edge-type categories: exactly the 7 canonical PAG relations,
# nothing else. am[a,b] convention (established elsewhere in this pipeline,
# scripts 21/27, run_one_ext() in the .qmd): am[a,b] holds the mark AT NODE
# b. classify_joint(mark_i, mark_j, ...) below is called with mark_i = the
# TRUE mark at i and mark_j = the TRUE mark at j (i.e. the call site passes
# am[j,i] as mark_i and am[i,j] as mark_j -- see run_one_joint below; this
# ordering was the source of a real labeling bug in v1, caught via a
# temporal-impossibility check on harm_present->climate_behavior, see
# analysis_decisions_log.md Section 21).
joint_cats <- c("i_to_j", "j_to_i", "bidirected",
                "i_ocirc_arrow_j", "j_ocirc_arrow_i", "circle_circle",
                "no_edge")

classify_joint <- function(mark_i, mark_j, name_i, name_j, b, alpha) {
  if (mark_i == 0L && mark_j == 0L) return("no_edge")
  if (mark_i == 3L && mark_j == 2L) return("i_to_j")           # tail@i, arrow@j: i -> j
  if (mark_i == 2L && mark_j == 3L) return("j_to_i")           # arrow@i, tail@j: j -> i
  if (mark_i == 2L && mark_j == 2L) return("bidirected")       # i <-> j
  if (mark_i == 1L && mark_j == 2L) return("i_ocirc_arrow_j")  # circle@i, arrow@j: i o-> j
  if (mark_i == 2L && mark_j == 1L) return("j_ocirc_arrow_i")  # arrow@i, circle@j: j o-> i
  if (mark_i == 1L && mark_j == 1L) return("circle_circle")    # i o-o j
  stop(sprintf(
    paste0("UNANTICIPATED PAG mark combination -- resample b=%d, alpha=%s, ",
           "pair (%s, %s): mark_at_%s=%d, mark_at_%s=%d ",
           "(0=no edge, 1=circle, 2=arrow, 3=tail). This matches none of the ",
           "7 canonical categories. Per instruction: stopping here rather ",
           "than silently folding an unanticipated encoding into an ",
           "'other'/unresolved bucket. Investigate this exact combination ",
           "before rerunning."),
    b, alpha, name_i, name_j, name_i, mark_i, name_j, mark_j
  ))
}

run_one_joint <- function(b, data, nms, alpha, ctx_idx, p, pairs_i, pairs_j) {
  idx    <- sample(nrow(data), replace = TRUE)
  suff_b <- list(C = cor(as.matrix(data[idx, ])), n = length(idx))

  # Only the FCI FIT ITSELF is allowed to fail silently (a legitimate,
  # already-expected possibility, counted separately below as n_failed) --
  # the classification step happens OUTSIDE this tryCatch so a stop() from
  # classify_joint's safety net is never accidentally swallowed here.
  fci_b <- tryCatch(
    pcalg::fci(
      suffStat = suff_b, indepTest = pcalg::gaussCItest,
      alpha = alpha, labels = nms,
      contextVars = ctx_idx, jci = "1",
      selectionBias = FALSE, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(fci_b)) return(rep(NA_character_, length(pairs_i)))

  am <- fci_b@amat
  cats_this_resample <- character(length(pairs_i))
  for (k in seq_along(pairs_i)) {
    i <- pairs_i[k]; j <- pairs_j[k]
    cats_this_resample[k] <- classify_joint(am[j, i], am[i, j], nms[i], nms[j], b, alpha)
  }
  cats_this_resample
}

# All p*(p-1)/2 unordered pairs, i<j in node_order_ext's index order.
pair_grid <- t(utils::combn(seq_len(p_ext), 2))
pairs_i <- pair_grid[, 1]
pairs_j <- pair_grid[, 2]
n_pairs <- length(pairs_i)

n_cores <- max(1L, parallelly::availableCores() - 1L)
plan(multisession, workers = n_cores)

joint_counts <- list()
for (alph in names(alphas)) {
  cat("Joint edge-type bootstrap (FCI only), alpha =", alph, "-- ", n_boot_joint, "resamples\n")
  res <- furrr::future_map(
    seq_len(n_boot_joint), run_one_joint,
    data = agg_ext, nms = node_order_ext, alpha = alphas[[alph]],
    ctx_idx = context_idx, p = p_ext, pairs_i = pairs_i, pairs_j = pairs_j,
    .options = furrr::furrr_options(seed = TRUE), .progress = TRUE
  )
  # An uncaught stop() from classify_joint (an unanticipated mark
  # combination) surfaces here as an error from future_map itself -- this
  # for-loop, and the whole script, halts at this point. That is the
  # intended behavior, not a bug to catch.
  counts <- matrix(0L, nrow = n_pairs, ncol = length(joint_cats),
                    dimnames = list(NULL, joint_cats))
  n_failed <- 0L
  for (cats_vec in res) {
    if (all(is.na(cats_vec))) { n_failed <- n_failed + 1L; next }
    for (k in seq_along(cats_vec)) {
      cat_k <- cats_vec[k]
      if (!is.na(cat_k)) counts[k, cat_k] <- counts[k, cat_k] + 1L
    }
  }
  if (n_failed > 0) cat("  (", n_failed, "of", n_boot_joint, "resamples failed to fit and were skipped)\n")
  joint_counts[[alph]] <- counts / (n_boot_joint - n_failed)
}
plan(sequential)

# --- Assemble the PRIMARY summary table for the 16 retained edges -- raw
# proportions ONLY, no derived/thresholded columns. ---------------------------
build_summary_row <- function(from, to, alph) {
  i <- match(from, node_order_ext); j <- match(to, node_order_ext)
  stopifnot(!is.na(i), !is.na(j))
  # pair_grid stores i<j by node_order_ext INDEX, not by (from,to) order --
  # find which stored index corresponds to this pair and orient the reported
  # percentages to match the paper's claimed (from -> to) direction.
  lo <- min(i, j); hi <- max(i, j)
  row_idx <- which(pairs_i == lo & pairs_j == hi)
  stopifnot(length(row_idx) == 1)
  props <- joint_counts[[alph]][row_idx, ]
  # If the paper's "from" is the higher-index node (hi), the stored
  # i_to_j/j_to_i (etc.) labels are relative to (lo,hi), i.e. REVERSED
  # relative to (from,to) -- swap the directional pairs' labels accordingly.
  reversed <- (from == node_order_ext[hi])
  pct_from_to       <- if (reversed) props[["j_to_i"]]          else props[["i_to_j"]]
  pct_to_from       <- if (reversed) props[["i_to_j"]]          else props[["j_to_i"]]
  pct_from_ocirc_to <- if (reversed) props[["j_ocirc_arrow_i"]] else props[["i_ocirc_arrow_j"]]
  pct_to_ocirc_from <- if (reversed) props[["i_ocirc_arrow_j"]] else props[["j_ocirc_arrow_i"]]
  tibble::tibble(
    from = from, to = to, alpha = alph,
    pct_from_to        = round(100 * pct_from_to, 1),        # from -> to
    pct_to_from        = round(100 * pct_to_from, 1),        # to -> from
    pct_bidirected     = round(100 * props[["bidirected"]], 1),
    pct_from_ocirc_to  = round(100 * pct_from_ocirc_to, 1),  # from o-> to
    pct_to_ocirc_from  = round(100 * pct_to_ocirc_from, 1),  # to o-> from (from <-o to)
    pct_circle_circle  = round(100 * props[["circle_circle"]], 1),
    pct_absent         = round(100 * props[["no_edge"]], 1)
  )
}

joint_summary <- purrr::pmap_dfr(
  list(from = rep(retained_edges$from, 2), to = rep(retained_edges$to, 2),
       alph = rep(names(alphas), each = nrow(retained_edges))),
  build_summary_row
)

cat("\n=== Joint PAG edge-type proportions, 16 retained edges, both alpha",
    "(FCI only) -- raw proportions, no derived classification ===\n")
print(as.data.frame(joint_summary), row.names = FALSE)

# --- A plain sort order for scanning 16 edges, NOT a classification. Uses
# only a sum of two already-reported raw quantities (bidirected +
# circle-circle), no threshold, no TRUE/FALSE verdict. Computed on a local
# copy so the written primary CSV (joint_summary) carries no derived
# columns at all. --------------------------------------------------------
current_flip_set <- c("belief_concern->weather_risk_prep", "harm_present->weather_risk_prep",
                       "social_norms->climate_behavior", "harm_future->harm_present")
# updated 2026-09-11: out with politics->belief_concern / policy_support->
# social_norms (now reversed, data-aligned, no longer uncertain), in with
# the two weather_risk_prep edges -- see 05_scm_intervention_helpers.R and
# analysis_decisions_log.md Section 34.

cross_check <- joint_summary |>
  dplyr::mutate(
    edge_key = paste0(from, "->", to),
    in_current_flip_set = edge_key %in% current_flip_set,
    confound_or_unresolved_pct = pct_bidirected + pct_circle_circle
  ) |>
  dplyr::group_by(from, to, in_current_flip_set) |>
  dplyr::summarise(
    max_confound_or_unresolved_pct = max(confound_or_unresolved_pct),
    max_bidirected_pct = max(pct_bidirected),
    max_circle_circle_pct = max(pct_circle_circle),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(max_confound_or_unresolved_pct))

cat("\n=== Sort order only (not a classification) -- 16 edges by how much",
    "bidirected+circle-circle mass each has at its worse alpha; inspect the",
    "full joint_summary table above for the actual per-category breakdown",
    "before drawing any conclusion ===\n")
print(as.data.frame(cross_check), row.names = FALSE)

# --- Write outputs -----------------------------------------------------------
dir.create("pipeline_outputs", showWarnings = FALSE)
write.csv(joint_summary, "pipeline_outputs/joint_pag_edgetype_audit.csv", row.names = FALSE)
write.csv(cross_check,   "pipeline_outputs/joint_pag_edgetype_crosscheck.csv", row.names = FALSE)
cat("\nWrote pipeline_outputs/joint_pag_edgetype_audit.csv (32 rows = 16 edges x 2 alpha,",
    "raw proportions only) and joint_pag_edgetype_crosscheck.csv (16 rows, sort order only).\n")
cat("\nNo downstream script (SCM, intervention enumeration, Figure 5/7) is rerun by",
    "this script. Inspect the full joint-type table (not just the sort order) before",
    "deciding whether the 4-edge flip set changes.\n")
