# =============================================================================
# 12_allwave_sensitivity_check.R
#
# WHY THIS SCRIPT EXISTS
# -----------------------
# The Supplement's all-wave-vs-Waves-1-4 sensitivity paragraph (Section
# "Waves-1-4 vs. all-wave composite construction") currently reports
# coefficients (.309/.097/.276 for all-wave; .306/.113/.249 for Waves-1-4-only)
# that predate the ew5-alone weather_risk_prep redefinition and were never
# found anywhere in the current .qmd. This script actually runs that specific
# comparison regression under the CURRENT (ew5-alone) construction, using the
# objects the main .qmd already creates:
#   - df_main_allwave / df_extended_allwave : all-wave composites (built before
#     Section 4.4's `df_main <- df_main_pre5` reassignment)
#   - df_main / df_extended                 : Waves-1-4-only composites (the
#     primary analysis dataset, after that reassignment)
# Both already use weather_risk_prep = ew5 alone (confirmed by reading the
# construction code directly), so this rerun is a clean, apples-to-apples
# check under the current pipeline -- no other change needed.
#
# HOW TO USE
# ----------
# Run this in a live session AFTER Section 4.4 has executed (so all four
# objects above exist -- this is already true by the time you reach the later
# analysis sections, so if rerunning the whole .qmd top to bottom this
# just needs to be sourced anywhere after that point). It only reads existing
# objects and writes one small CSV; no bootstrap, no refitting of the SCM.
#
# OUTPUT
# ------
#   tables/allwave_sensitivity_check.csv -- one row per (construction, term),
#   with the standardized coefficient, so the manuscript sentence can be
#   updated with confirmed numbers.
#   Also prints the per-node all-wave vs. Waves-1-4 correlations (the r>=.992
#   claim earlier in that same paragraph) so that can be reconfirmed too.
# =============================================================================

stopifnot(exists("df_main_allwave"), exists("df_extended_allwave"),
          exists("df_main"), exists("df_extended"))

dir.create("tables", showWarnings = FALSE)

NODES <- c("belief_concern", "harm_present", "harm_future", "trust_science",
           "policy_support", "social_norms", "politics", "weather_risk_prep")

cat("=== Per-node correlation: all-wave vs. Waves-1-4-only composites ===\n")
common_ids <- intersect(df_main_allwave$participant_id, df_main$participant_id)
m_all <- df_main_allwave[match(common_ids, df_main_allwave$participant_id), ]
m_w14 <- df_main[match(common_ids, df_main$participant_id), ]
node_cors <- sapply(NODES, function(nd) {
  round(cor(m_all[[nd]], m_w14[[nd]], use = "pairwise.complete.obs"), 4)
})
print(node_cors)

cat("\n=== Sensitivity regression: climate_behavior ~ harm_present + social_norms + weather_risk_prep ===\n")

fit_one <- function(dat, label) {
  d <- dat[, c("climate_behavior", "harm_present", "social_norms", "weather_risk_prep")]
  d <- as.data.frame(scale(d))
  fit <- lm(climate_behavior ~ harm_present + social_norms + weather_risk_prep, data = d)
  co <- coef(fit)[-1]
  data.frame(
    construction = label,
    term = names(co),
    beta = round(unname(co), 3),
    row.names = NULL
  )
}

res_allwave <- fit_one(df_extended_allwave, "all_wave")
res_w14     <- fit_one(df_extended,         "waves_1_4_only")

result <- rbind(res_allwave, res_w14)
print(result)

write.csv(result, "tables/allwave_sensitivity_check.csv", row.names = FALSE)
message("Saved tables/allwave_sensitivity_check.csv")

cat("\n--- WHAT TO SEND BACK ---\n")
cat("1. The two 'node_cors' and 'result' console printouts above (or just the\n")
cat("   CSV) -- these give the confirmed beta_harm_present / beta_social_norms /\n")
cat("   beta_weather_risk for both constructions under the current ew5-alone\n")
cat("   pipeline, replacing the stale .309/.097/.276 and .306/.113/.249 numbers.\n")
