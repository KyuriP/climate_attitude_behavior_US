# =============================================================================
# 18_finalize_scm_specification.R  (v4)
#
# Purpose: replace "the SCM's evidence tiers are asserted by hand in
# figscripts/07_figure5_scm_hierarchical.R's hardcoded scm_edges tribble"
# with "the SCM's evidence tiers are AUDITED against the real primary
# bootstrap array by explicit, documented rules, and any mismatch or weak
# spot is flagged instead of silently asserted." Mechanically applies one
# single, fixed existence rule (60%) with no edge-by-edge exceptions.
#
# Output: writes scm_edges_finalized.csv (from, to, current_evidence,
# p_adjacent, avg_n, asymmetry, existence_band, orientation_label,
# final_tier, current_tier_expected, mismatch, avg_n_published,
# asymmetry_published, avg_n_drift, asymmetry_drift, table2_stale, note)
# for figscripts/07_figure5_scm_hierarchical.R to read.
#
# v4 adds a 16th edge, harm_future -> trust_science, because causal discovery
# already surfaced it and there was no principled reason to leave it out.
# Surfaced by v3's own candidate-scan (the one that runs at the same 60%
# threshold as the skeleton itself) at p_adjacent = .657, asymmetry = +.2865.
# Two things distinguish this edge from the 15th (belief_concern ->
# weather_risk_prep, added in v3):
#   - Existence is close to the 60% floor (65.7%), not comfortably above it
#     like weather_risk_prep's 94.6% -- worth flagging as a genuinely
#     thinner margin, in case a bootstrap-seed sensitivity check is ever
#     wanted for this specific edge.
#   - Orientation, unlike weather_risk_prep, IS cleanly resolved: +.2865 is
#     well within the typical range for already-resolved edges (+.15 to
#     +.55), not a weak/both-ends-high-arrowhead pattern. So this edge is
#     simpler than weather_risk_prep in one respect (no
#     latent-confounder-shaped orientation signal) and more marginal in
#     another (existence closer to the floor).
# Same accounting as the 15th edge: current_evidence = "new_2026-09-03c",
# current_tier_expected = NA (nothing to compare a brand-new edge against),
# mismatch = FALSE.
#
# Everything else unchanged from v3 -- see that file's own header for the
# single-existence-threshold, weather_risk_prep-addition, and
# candidate-scan-threshold reasoning.
# ==============================================================================
#
# ============================ TESTING STATUS ================================
# Verified against the same real-numbers reconstruction of fci_props_ext used
# to validate v3 (built from every pairwise row of
# sensitivity_RCoT_vs_primary_ext.csv). Confirmed: 16 rows written;
# harm_future->trust_science lands in "data_aligned" (cleanly resolved,
# unlike weather_risk_prep's weak-orientation note); a direct networkx
# cycle-check (all 8 combinations of the still-3-edge flip set, now against
# this 16-edge fixed baseline) confirms all 8 remain acyclic -- i.e. adding
# this edge does not reopen any of the cyclicity questions resolved when TP
# was removed from the flip set.
# NOT yet re-run against the live fci_props_ext -- needs a real run
# before trusting the exact printed numbers below for the manuscript.
# ==============================================================================

suppressPackageStartupMessages(library(dplyr))

if (!exists("fci_props_ext")) {
  if (file.exists("pipeline_outputs/fci_props_ext.rds")) {
    fci_props_ext <- readRDS("pipeline_outputs/fci_props_ext.rds")
  } else {
    stop("fci_props_ext not found in session and fci_props_ext.rds doesn't exist -- ",
         "run the primary FCI extended-network bootstrap (Section 7.4) in a live ",
         "analysis session first.")
  }
}

EXT_NODES <- dimnames(fci_props_ext)[[1]]

# ---- The SCM's CURRENT 16 edges (14 original + belief_concern->weather_risk_prep
# [v3] + harm_future->trust_science [v4, see above]). This script does not
# invent or remove edges beyond what's listed here -- adding/removing one is
# still a modeling decision made outside this script, just recorded here as
# the audit target. ------------------------------------------------------------
CURRENT_SCM_EDGES <- tibble::tribble(
  ~from,                ~to,                   ~current_evidence,
  "politics",           "belief_concern",      "substantive",
  "belief_concern",     "harm_future",         "bootstrap_theory",
  "belief_concern",     "harm_present",        "bootstrap_theory",
  "harm_future",        "harm_present",        "bootstrap",
  "belief_concern",     "trust_science",       "substantive",
  "belief_concern",     "policy_support",      "bootstrap_theory",
  "politics",           "policy_support",      "bootstrap",
  "trust_science",      "policy_support",      "bootstrap_theory",
  "trust_science",      "social_norms",        "bootstrap",
  "policy_support",     "social_norms",        "substantive",
  "harm_present",       "weather_risk_prep",   "substantive",
  "harm_present",       "climate_behavior",    "bootstrap",
  "weather_risk_prep",  "climate_behavior",    "bootstrap",
  "social_norms",       "climate_behavior",    "weak",
  "belief_concern",     "weather_risk_prep",   "new_2026-09-03",
  "harm_future",        "trust_science",       "new_2026-09-03c"
)
stopifnot(nrow(CURRENT_SCM_EDGES) == 16)

# ---- Table 2's OWN published Avg.N / Asymmetry, transcribed verbatim from
# main.tex (\label{tab:orientation}) -- used ONLY for the staleness
# comparison below, never for classification. Still 14 rows -- neither of the
# 2 newly-added edges has a published entry to compare against, so their
# staleness columns come back NA and table2_stale is forced FALSE for those
# rows specifically (nothing to be stale relative to). ------------------------
PUBLISHED_TABLE2 <- tibble::tribble(
  ~from,                ~to,                   ~avg_n_published, ~asymmetry_published,
  "weather_risk_prep",  "climate_behavior",     .00,  .53,
  "harm_present",       "climate_behavior",     .08,  .52,
  "harm_future",        "harm_present",         .00,  .30,
  "trust_science",      "social_norms",         .00,  .18,
  "politics",           "policy_support",       .03,  .20,
  "belief_concern",     "harm_future",          .00,  .36,
  "belief_concern",     "policy_support",       .07,  .22,
  "trust_science",      "policy_support",       .00,  .22,
  "belief_concern",     "harm_present",         .00,  .12,
  "policy_support",     "social_norms",         .21, -.38,
  "politics",           "belief_concern",       .07, -.26,
  "belief_concern",     "trust_science",        .00, -.20,
  "harm_present",       "weather_risk_prep",    .36, -.17,
  "social_norms",       "climate_behavior",     .38,  .00
)
stopifnot(nrow(PUBLISHED_TABLE2) == 14)

# ---- Thresholds -- fixed here, before classifying anything (unchanged from v3) --
EXISTENCE_MIN <- 0.60
ORIENTATION_ASYMMETRY_EPS <- 0.01
TABLE2_DRIFT_FLAG_MIN <- 0.15

get_pair_stats <- function(from, to) {
  if (!(from %in% EXT_NODES) || !(to %in% EXT_NODES)) {
    stop(sprintf("Pair %s-%s: one or both nodes not found in fci_props_ext.", from, to))
  }
  mark_at_to   <- fci_props_ext[from, to, ]
  mark_at_from <- fci_props_ext[to, from, ]
  p_adjacent <- 1 - mark_at_from[["N"]]
  asymmetry <- mark_at_to[[">"]] - mark_at_from[[">"]]
  list(
    p_adjacent = p_adjacent,
    avg_n      = 1 - p_adjacent,
    asymmetry  = asymmetry
  )
}

classify_existence <- function(p_adjacent) {
  if (p_adjacent >= EXISTENCE_MIN) "included" else "weak"
}

classify_orientation <- function(asymmetry) {
  if (is.na(asymmetry)) return(NA)
  asymmetry > ORIENTATION_ASYMMETRY_EPS
}

# ---- Supplementary notes -- edges where a technically-clean tier label
# would understate real nuance. social_norms->climate_behavior and
# belief_concern->weather_risk_prep notes unchanged from v3. New:
# harm_future->trust_science (v4, see above). ---------------------------------
SUPPLEMENTARY_NOTES <- list(
  "social_norms->climate_behavior" = paste(
    "Unconditional RCoT p ~ 1e-16 (stable across 20 seeds); no single",
    "main-network variable, conditioned on individually, explains it away;",
    "conditional on the full network jointly, 60-65% support under two",
    "independent methods (npn primary: 0.6615; RCoT diagnostic: 12-13/20",
    "seeds significant at safe num_f, i.e. not resolved by adding more",
    "random features). See r_patches/17_diagnostic_rcot_social_norms_behavior.R."
  ),
  "belief_concern->weather_risk_prep" = paste(
    "Added 2026-09-03: existence is robust (p_adjacent=.946, higher than 10",
    "of the other 14 edges), clearing the same 60% rule as every other edge.",
    "Orientation is only weakly resolved (raw asymmetry +.065, versus +.15 to",
    "+.55 for other data_aligned edges), with high arrowhead mass at BOTH",
    "endpoints (.77 at belief_concern, .83 at weather_risk_prep) and little",
    "circle/undetermined mass at either -- a pattern at least as consistent",
    "with a shared latent upstream cause as with a clean directed effect.",
    "Rendered solid (data_aligned) per the same sign-only rule used",
    "everywhere else, but this note exists so the figure's solid line does",
    "not overstate how cleanly this one is resolved."
  ),
  "harm_future->trust_science" = paste(
    "Added 2026-09-03 (second candidate-scan pass, threshold matched to 60%):",
    "existence is close to the skeleton floor (p_adjacent=.657, versus .946",
    "for belief_concern->weather_risk_prep and 78%+ for every one of the",
    "original 14 edges) -- a genuinely thinner margin, worth a bootstrap-seed",
    "sensitivity check if one is wanted specifically for this edge. Unlike",
    "weather_risk_prep, orientation itself IS cleanly resolved here (raw",
    "asymmetry +.2865, squarely within the +.15 to +.55 range typical of",
    "already-resolved edges, no both-ends-high-arrowhead pattern). No",
    "established theoretical story for this direction has been discussed yet",
    "(unlike weather_risk_prep, where belief_concern as a second driver of",
    "weather-risk perception is an easy substantive story) -- worth a sentence in",
    "Methods/Discussion regardless."
  )
)

# ---- Audit the 16 current edges --------------------------------------------
audit_results <- CURRENT_SCM_EDGES |>
  left_join(PUBLISHED_TABLE2, by = c("from", "to")) |>
  rowwise() |>
  mutate(
    stats             = list(get_pair_stats(from, to)),
    p_adjacent        = stats$p_adjacent,
    avg_n             = stats$avg_n,
    asymmetry         = stats$asymmetry,
    existence_band    = classify_existence(p_adjacent),
    orientation_ok    = classify_orientation(asymmetry),
    orientation_label = if (isTRUE(orientation_ok)) "bootstrap-resolved" else "theory-needed",
    final_tier = case_when(
      existence_band == "weak" ~ "FLAG_WEAK_EXISTENCE",
      existence_band == "included" & isTRUE(orientation_ok) ~ "data_aligned",
      existence_band == "included" & !isTRUE(orientation_ok) ~ "substantive",
      TRUE ~ "FLAG_UNCLASSIFIED"
    ),
    current_tier_expected = case_when(
      current_evidence %in% c("bootstrap", "bootstrap_theory") ~ "data_aligned",
      current_evidence == "substantive" ~ "substantive",
      current_evidence == "weak" ~ "substantive",
      current_evidence %in% c("new_2026-09-03", "new_2026-09-03c") ~ NA_character_,
      TRUE ~ NA_character_
    ),
    mismatch = !is.na(current_tier_expected) & (final_tier != current_tier_expected),
    asymmetry_drift = asymmetry - asymmetry_published,
    avg_n_drift     = avg_n - avg_n_published,
    table2_stale = if (is.na(avg_n_published)) FALSE else
      ((abs(asymmetry_drift) >= TABLE2_DRIFT_FLAG_MIN) | (abs(avg_n_drift) >= TABLE2_DRIFT_FLAG_MIN)),
    note = coalesce(SUPPLEMENTARY_NOTES[[paste0(from, "->", to)]], ""),
  ) |>
  ungroup() |>
  select(from, to, current_evidence, p_adjacent, avg_n, asymmetry, existence_band,
         orientation_label, final_tier, current_tier_expected, mismatch,
         avg_n_published, asymmetry_published, avg_n_drift, asymmetry_drift,
         table2_stale, note)

cat("=============================================================\n")
cat("SCM edge audit -- current 16 edges vs. real bootstrap array (fci_props_ext)\n")
cat("=============================================================\n")
print(as.data.frame(audit_results[, c("from","to","p_adjacent","asymmetry","existence_band",
                                        "orientation_label","final_tier","mismatch")]),
      row.names = FALSE)

n_flagged <- sum(audit_results$final_tier %in% c("FLAG_WEAK_EXISTENCE", "FLAG_UNCLASSIFIED"))
n_mismatch <- sum(audit_results$mismatch, na.rm = TRUE)
if (n_flagged > 0) {
  message(sprintf(
    paste0("\n%d edge(s) flagged for existence below the skeleton threshold ",
           "(p_adjacent < %.2f) -- these are currently IN the SCM but wouldn't clear ",
           "the agreed rule. Needs a decision: drop, or keep with an explicit caveat."),
    n_flagged, EXISTENCE_MIN))
  print(as.data.frame(audit_results[audit_results$final_tier %in%
                                       c("FLAG_WEAK_EXISTENCE", "FLAG_UNCLASSIFIED"),
                                     c("from","to","p_adjacent","final_tier")]),
        row.names = FALSE)
}
if (n_mismatch > 0) {
  message(sprintf(
    paste0("\n%d edge(s) where the real data's classification doesn't match the ",
           "currently-asserted evidence tag -- worth a look before finalizing."),
    n_mismatch))
  print(as.data.frame(audit_results[audit_results$mismatch,
                                     c("from","to","current_evidence","final_tier")]),
        row.names = FALSE)
}

write.csv(audit_results, "pipeline_outputs/scm_edges_finalized.csv", row.names = FALSE)
message("\nWrote scm_edges_finalized.csv (16 edges) -- figscripts/07_figure5_scm_hierarchical.R ",
        "should read `final_tier` from this (data_aligned/substantive/FLAG_* only).")

# =============================================================================
# TABLE 2 STALENESS CHECK -- unchanged mechanics from v3. Neither new edge has
# a published row to compare against, so both are excluded from the "stale"
# printout below (table2_stale forced FALSE) but still shown in the full
# comparison for completeness.
# =============================================================================
stale_edges <- audit_results[audit_results$table2_stale, ]
cat("\n=============================================================\n")
cat("Table 2 staleness check (published vs. current array)\n")
cat("=============================================================\n")
print(as.data.frame(audit_results[, c("from","to","avg_n_published","avg_n",
                                        "asymmetry_published","asymmetry","table2_stale")]),
      row.names = FALSE)
if (nrow(stale_edges) > 0) {
  message(sprintf(
    paste0("\n%d edge(s) drifted by >= %.2f from Table 2's published Avg.N or ",
           "Asymmetry -- worth checking whether Table 2 needs regenerating from ",
           "the current bootstrap array (e.g. if this is downstream of the ",
           "weather_risk_prep redefinition)."),
    nrow(stale_edges), TABLE2_DRIFT_FLAG_MIN))
  print(as.data.frame(stale_edges[, c("from","to","avg_n_published","avg_n",
                                        "asymmetry_published","asymmetry")]),
        row.names = FALSE)
} else {
  message("\nNo edge drifted by >= ", TABLE2_DRIFT_FLAG_MIN,
          " from Table 2's published values -- current array is consistent with what's published.")
}

# ---- Candidate-scan, same 0.60 threshold as the skeleton itself. With 16
# edges now in current_pairs_set, both belief_concern-weather_risk_prep AND
# harm_future-trust_science correctly stop appearing here. -------------------
all_pairs <- combn(EXT_NODES, 2, simplify = FALSE)
current_pairs_set <- apply(CURRENT_SCM_EDGES[, c("from","to")], 1,
                            function(r) paste(sort(r), collapse = "|"))

candidates <- lapply(all_pairs, function(p) {
  key <- paste(sort(p), collapse = "|")
  if (key %in% current_pairs_set) return(NULL)
  st <- tryCatch(get_pair_stats(p[1], p[2]), error = function(e) NULL)
  if (is.null(st) || is.na(st$p_adjacent) || st$p_adjacent < EXISTENCE_MIN) return(NULL)
  data.frame(A = p[1], B = p[2], p_adjacent = st$p_adjacent, asymmetry = st$asymmetry)
}) |> Filter(Negate(is.null), x = _) |> bind_rows()

if (nrow(candidates) > 0) {
  candidates <- candidates |> arrange(desc(p_adjacent))
  cat("\n=============================================================\n")
  cat("Candidate pairs NOT currently in the SCM, p_adjacent >= ", EXISTENCE_MIN, "\n")
  cat("(informational -- NOT auto-added; each needs the same substantive\n")
  cat("judgment used to add the last two edges -- is this a real direct edge,\n")
  cat("or fully mediated through an existing SCM path? Read this list in full\n")
  cat("before assuming the 16-edge SCM is the complete picture.\n")
  cat("=============================================================\n")
  print(as.data.frame(candidates), row.names = FALSE)
  write.csv(candidates, "pipeline_outputs/scm_candidate_edges_for_review.csv", row.names = FALSE)
} else {
  message("\nNo candidate pairs found outside the current 16 clearing the skeleton threshold.")
}

message("\nDone.")
