# formalizes the 3-way edge classification specified 2026-09-08, instead
# of leaving it as a one-off hand count: reuses the paper's existing 60%
# adjacency rule (EXISTENCE_MIN) rather than inventing a new percentage, and
# requires the SAME plurality type at both alpha levels before calling
# anything resolved -- same cross-threshold-consistency logic the rest of
# this pipeline already uses for direction claims, just applied to
# bidirectedness too.
#
# three labels, no others:
#   - bidirected_dominant: passes EXISTENCE_MIN, and bidirected is the
#     plurality joint-endpoint type at BOTH .05 and .01.
#   - directionally_resolved: passes EXISTENCE_MIN, and the SAME directed
#     type (from_to or to_from) is the plurality at BOTH .05 and .01.
#   - mixed_threshold_sensitive: everything else -- plurality type differs
#     between alphas, or fails EXISTENCE_MIN. this bucket isn't one thing --
#     it covers genuine direction reversal (harm_future-harm_present),
#     directed/bidirected flips (politics-belief_concern,
#     harm_present-climate_behavior), and existence-thin edges
#     (harm_future-trust_science) -- read the plurality columns themselves
#     before assuming they're all the same kind of "mixed."
#
# reads r_patches/29's existing output, doesn't touch main9.tex, doesn't
# rerun any bootstrap.

source("clean_pipeline/00_config.R")

audit_path <- file.path(OUTPUT_DIR, "joint_pag_edgetype_audit.csv")
scm_path   <- file.path(OUTPUT_DIR, "scm_edges_finalized.csv")
stopifnot(file.exists(audit_path), file.exists(scm_path))

audit <- utils::read.csv(audit_path, stringsAsFactors = FALSE)
scm   <- utils::read.csv(scm_path, stringsAsFactors = FALSE)

pct_cols <- c(from_to = "pct_from_to", to_from = "pct_to_from",
              bidirected = "pct_bidirected", from_ocirc_to = "pct_from_ocirc_to",
              to_ocirc_from = "pct_to_ocirc_from", circle_circle = "pct_circle_circle",
              absent = "pct_absent")
stopifnot(all(pct_cols %in% names(audit)))

pct_mat <- as.matrix(audit[, unname(pct_cols)])
audit$plurality <- names(pct_cols)[apply(pct_mat, 1, which.max)]

wide <- reshape(
  audit[, c("from", "to", "alpha", "plurality")],
  timevar = "alpha", idvar = c("from", "to"), direction = "wide"
)
names(wide) <- gsub("^plurality\\.", "plurality_", names(wide))

classified <- merge(wide, scm[, c("from", "to", "p_adjacent")], by = c("from", "to"))
classified$passes_existence <- classified$p_adjacent >= EXISTENCE_MIN

classified$label <- with(classified, ifelse(
  passes_existence & plurality_0.05 == "bidirected" & plurality_0.01 == "bidirected",
  "bidirected_dominant",
  ifelse(
    passes_existence & plurality_0.05 == plurality_0.01 &
      plurality_0.05 %in% c("from_to", "to_from"),
    "directionally_resolved",
    "mixed_threshold_sensitive"
  )
))

cat("\n--- edge classification: existing 60% adjacency rule + same plurality",
    "type at both alpha levels ---\n")
print(classified[order(classified$label), c("from", "to", "p_adjacent",
                                             "plurality_0.05", "plurality_0.01", "label")],
      row.names = FALSE)
cat("\ncounts:\n")
print(table(classified$label))

write.csv(classified, file.path(OUTPUT_DIR, "tables", "scm_edge_confounding_classification.csv"),
          row.names = FALSE)
cat("\nWrote", file.path(OUTPUT_DIR, "tables", "scm_edge_confounding_classification.csv"), "\n")

# cross-check against 15's hardcoded confound_candidates -- same 8 edges,
# written out again here independently rather than parsed from 15's source,
# so this catches drift in EITHER direction if the audit data or the
# hardcoded list ever gets out of sync.
confound_candidates_in_15 <- data.frame(
  from = c("belief_concern", "politics", "policy_support", "trust_science",
           "harm_present", "belief_concern", "social_norms", "weather_risk_prep"),
  to   = c("harm_future", "policy_support", "social_norms", "social_norms",
           "weather_risk_prep", "weather_risk_prep", "climate_behavior", "climate_behavior")
)
key <- function(df) sort(paste(df$from, df$to, sep = "->"))
bidirected_here <- classified[classified$label == "bidirected_dominant", ]
if (!identical(key(confound_candidates_in_15), key(bidirected_here))) {
  warning("15_confound_sensitivity_diagnostic.R's hardcoded confound_candidates does NOT ",
          "match this script's bidirected_dominant set -- one of them is stale, don't ",
          "trust 15's output until this is resolved.")
} else {
  cat("\ncross-check OK: 15's hardcoded 8-edge confound_candidates matches this",
      "classification's bidirected_dominant set exactly.\n")
}
