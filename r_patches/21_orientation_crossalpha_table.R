# =============================================================================
# 21_orientation_crossalpha_table.R
#
# Purpose: descriptive input for a team decision on the orientation-classification
# rule (methodological review, item #2). Does NOT change which edges are in the
# working SCM (the 60% pooled existence rule is untouched) and does NOT rerun the
# orientation-uncertainty enumeration. It only reports, for each of the 16 SCM
# edges, the FCI bootstrap endpoint evidence separately at alpha=.05 and
# alpha=.01, alongside a *proposed* cross-alpha classification for review.
#
# Proposed rule under review (not yet adopted): an edge's direction counts as
# data-supported only if, at BOTH alpha=.05 and alpha=.01 separately, the
# dominant non-absence endpoint mark at the destination is an arrowhead (">")
# and the dominant non-absence endpoint mark at the source is not an arrowhead
# (i.e. a clean single-headed pattern in the stated direction at each alpha,
# not a circle/circle or arrow/arrow pattern). If either alpha fails this, or
# the two alphas disagree in sign, the edge is flagged unresolved -> the
# existing theory-completion process applies, and the edge would move into the
# orientation-sensitivity enumeration.
#
# Requires pipeline_outputs/mark_props_ext_by_alpha.rds (per-alpha FCI mark
# proportions, saved by the .qmd immediately after the extended bootstrap is
# pooled) and pipeline_outputs/fci_props_ext.rds (existing pooled array, used
# here only for the side-by-side "current pooled classification" column).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

mp_path <- "pipeline_outputs/mark_props_ext_by_alpha.rds"
pooled_path <- "pipeline_outputs/fci_props_ext.rds"

if (!file.exists(mp_path)) {
  stop("pipeline_outputs/mark_props_ext_by_alpha.rds not found. Re-render the .qmd's ",
       "extended-bootstrap section (7.4) first -- it now saves this object automatically.")
}
if (!file.exists(pooled_path)) {
  stop("pipeline_outputs/fci_props_ext.rds not found. Re-render Section 7.4 first.")
}

mark_props_ext <- readRDS(mp_path)
fci_props_ext  <- readRDS(pooled_path)

stopifnot(all(c("0.05", "0.01") %in% names(mark_props_ext)))
fci_non_n_e <- c("o", ">", "-")

# The current 16-edge working SCM (unchanged; this script does not alter
# skeleton membership). Kept as an explicit from/to list so results are
# reported in exactly the same edge set and order as Table 2 in the manuscript.
scm_edges <- tribble(
  ~from,             ~to,
  "politics",        "belief_concern",
  "harm_future",     "harm_present",
  "belief_concern",  "harm_present",
  "belief_concern",  "harm_future",
  "trust_science",   "policy_support",
  "belief_concern",  "policy_support",
  "politics",        "policy_support",
  "belief_concern",  "trust_science",
  "harm_future",     "trust_science",
  "trust_science",   "social_norms",
  "policy_support",  "social_norms",
  "harm_present",    "weather_risk_prep",
  "belief_concern",  "weather_risk_prep",
  "harm_present",    "climate_behavior",
  "weather_risk_prep", "climate_behavior",
  "social_norms",    "climate_behavior"
)

dominant_nonN_mark <- function(props, i, j) {
  fci_non_n_e[which.max(props[i, j, fci_non_n_e])]
}

alpha_evidence <- function(alph, from, to) {
  props <- mark_props_ext[[alph]]$fci
  arrow_at_to   <- props[from, to, ">"]
  arrow_at_from <- props[to, from, ">"]
  dom_to        <- dominant_nonN_mark(props, from, to)
  dom_from      <- dominant_nonN_mark(props, to, from)
  clean_directed <- (dom_to == ">") && (dom_from != ">")
  list(
    arrow_at_dest   = round(arrow_at_to, 3),
    arrow_at_source = round(arrow_at_from, 3),
    asymmetry       = round(arrow_at_to - arrow_at_from, 3),
    dom_mark_dest   = dom_to,
    dom_mark_source = dom_from,
    clean_directed  = clean_directed
  )
}

results <- scm_edges |>
  rowwise() |>
  mutate(
    ev05 = list(alpha_evidence("0.05", from, to)),
    ev01 = list(alpha_evidence("0.01", from, to)),
    asym_a05        = ev05$asymmetry,
    asym_a01        = ev01$asymmetry,
    dom_pattern_a05 = paste0(ev05$dom_mark_source, " @source / ", ev05$dom_mark_dest, " @dest"),
    dom_pattern_a01 = paste0(ev01$dom_mark_source, " @source / ", ev01$dom_mark_dest, " @dest"),
    clean_a05       = ev05$clean_directed,
    clean_a01       = ev01$clean_directed,
    same_sign       = sign(asym_a05) == sign(asym_a01) && asym_a05 != 0 && asym_a01 != 0,
    # Current pooled classification, for side-by-side comparison only.
    asymmetry_pooled = round(fci_props_ext[from, to, ">"] - fci_props_ext[to, from, ">"], 3),
    pooled_positive  = asymmetry_pooled > 0,
    proposed_status  = if (clean_a05 && clean_a01 && same_sign && asym_a05 > 0 && asym_a01 > 0) {
      "data-supported (proposed)"
    } else {
      "unresolved -> theory-completed (proposed)"
    },
    status_changes = proposed_status == "unresolved -> theory-completed (proposed)" && pooled_positive
  ) |>
  ungroup() |>
  select(from, to, asymmetry_pooled, pooled_positive,
         asym_a05, dom_pattern_a05, clean_a05,
         asym_a01, dom_pattern_a01, clean_a01,
         same_sign, proposed_status, status_changes)

dir.create("pipeline_outputs", showWarnings = FALSE)
write.csv(results, "pipeline_outputs/orientation_crossalpha_table.csv", row.names = FALSE)

cat("\n=== Cross-alpha orientation evidence, all 16 SCM edges ===\n")
print(as.data.frame(results), row.names = FALSE)

n_flip <- sum(results$status_changes)
cat("\n", n_flip, " of 16 edges currently classified as bootstrap-resolved (positive pooled ",
    "asymmetry) would move to unresolved/theory-completed under the proposed cross-alpha rule.\n",
    "Full table written to pipeline_outputs/orientation_crossalpha_table.csv -- review before ",
    "deciding whether to adopt this rule and, if so, before rerunning the orientation enumeration ",
    "with whichever edges end up in the unresolved set.\n", sep = "")
