# 40_review_checks.R (2026-09-26) -- read-only checks requested in review.
# Nothing here changes a fitted model used in the paper.
setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
source("clean_pipeline/05_scm_intervention_helpers.R")
stopifnot(exists("df_extended"), nrow(df_extended) == 870)
cat("\n==== (2) free parameters in the baseline working SCM ====\n")
fit_b <- lavaan::sem(build_lavaan_syntax(base_edges), data = df_extended, estimator = "MLR", fixed.x = FALSE)
pt <- lavaan::parameterTable(fit_b)
print(pt[pt$op == "~~" & pt$lhs != pt$rhs, c("lhs", "op", "rhs", "free", "est")], row.names = FALSE)
print(lavaan::fitMeasures(fit_b, c("df", "npar", "aic", "bic")))
ss <- lavaan::standardizedSolution(fit_b)
print(ss[ss$op == "~~" & ss$lhs != ss$rhs, c("lhs", "rhs", "est.std", "se", "pvalue")], row.names = FALSE)
cat("n used:", lavaan::lavInspect(fit_b, "nobs"), "\n")
cat("lavaan version:", as.character(packageVersion("lavaan")), "\n")
# same check for the extra present-harm specification
e2 <- base_edges[!(base_edges$from == "harm_present" & base_edges$to == "climate_behavior"), ]
fit_2 <- lavaan::sem(paste(build_lavaan_syntax(e2), "harm_present ~~ climate_behavior", sep = "\n"),
                     data = df_extended, estimator = "MLR", fixed.x = FALSE)
pt2 <- lavaan::parameterTable(fit_2)
cat("\n-- hp_cb_cov spec off-diagonal covariances --\n")
print(pt2[pt2$op == "~~" & pt2$lhs != pt2$rhs, c("lhs", "op", "rhs", "free", "est")], row.names = FALSE)
print(lavaan::fitMeasures(fit_2, c("df", "npar")))

cat("\n==== (4c) nine-node PC-stable: retained SCM pairs, dominant marks ====\n")
latest <- read.csv(LATEST_POINTER_PATH, stringsAsFactors = FALSE)
pcx <- readRDS(latest$pc_props_ext_path)
cat("pc_props_ext dims:", paste(dim(pcx), collapse = "x"), " marks:", paste(dimnames(pcx)[[3]], collapse = ","), "\n")
for (k in seq_len(nrow(base_edges))) {
  a <- base_edges$from[k]; b <- base_edges$to[k]
  m_b <- pcx[a, b, ]; m_a <- pcx[b, a, ]
  cat(sprintf("%-18s -- %-18s  mark@%s: %s | mark@%s: %s\n", a, b,
              b, paste(sprintf("%s=%.2f", names(m_b), m_b), collapse = " "),
              a, paste(sprintf("%s=%.2f", names(m_a), m_a), collapse = " ")))
}
cat("40 DONE\n")
