# 43_review2_checks.R (2026-09-26): read-only documentation checks for review round 2.
# (A) which residual covariances lavaan's defaults added in every reported specification
# (B) per-intervention reduced-system diagnostics for the four cyclic specifications
setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
source("clean_pipeline/05_scm_intervention_helpers.R")
stopifnot(exists("df_extended"), nrow(df_extended) == 870)
lab <- c(belief_concern = "BC", harm_present = "HP", harm_future = "HF", policy_support = "PS",
         trust_science = "TS", social_norms = "SN", politics = "POL", weather_risk_prep = "WW",
         climate_behavior = "CB")
fit_spec <- function(edges, cf = character(0), ct = character(0)) {
  syn <- build_lavaan_syntax(edges)
  if (length(cf)) syn <- paste(c(syn, paste(cf, "~~", ct)), collapse = "\n")
  lavaan::sem(syn, data = df_extended, estimator = "MLR", fixed.x = FALSE)
}
covs_of <- function(fit, explicit = character(0)) {
  pt <- lavaan::parameterTable(fit)
  r <- pt[pt$op == "~~" & pt$lhs != pt$rhs & pt$free > 0, ]
  if (!nrow(r)) return("none")
  k <- paste(lab[r$lhs], lab[r$rhs], sep = "~~")
  paste(ifelse(paste(r$lhs, r$rhs) %in% explicit | paste(r$rhs, r$lhs) %in% explicit, paste0(k, "*"), k), collapse = ", ")
}
rows <- list()
# orientation enumeration (16 combos)
for (nm in names(scenario_list)) {
  e <- flip_edges(base_edges, scenario_list[[nm]])
  f <- fit_spec(e)
  fm <- lavaan::fitMeasures(f, c("df", "npar"))
  rows[[length(rows) + 1]] <- data.frame(family = "orientation", spec = nm, flipped = scenario_flip_labels[[nm]],
    acyclic = is_acyclic(e), df = fm[["df"]], npar = fm[["npar"]], auto_covariances = covs_of(f))
}
# confounding specs
conf8 <- data.frame(from = c("belief_concern","politics","social_norms","trust_science","harm_present","belief_concern","social_norms","weather_risk_prep"),
                    to   = c("harm_future","policy_support","policy_support","social_norms","weather_risk_prep","weather_risk_prep","climate_behavior","climate_behavior"))
drop_e <- function(ed, f, t) { k <- rep(TRUE, nrow(ed)); for (i in seq_along(f)) k <- k & !(ed$from == f[i] & ed$to == t[i]); ed[k, ] }
cspecs <- c(lapply(seq_len(8), function(i) list(nm = paste0("cov: ", lab[conf8$from[i]], "-", lab[conf8$to[i]]), f = conf8$from[i], t = conf8$to[i])),
            list(list(nm = "cov: all 8", f = conf8$from, t = conf8$to)),
            list(list(nm = "cov: HP-CB", f = "harm_present", t = "climate_behavior")))
for (s in cspecs) {
  f <- fit_spec(drop_e(base_edges, s$f, s$t), s$f, s$t)
  fm <- lavaan::fitMeasures(f, c("df", "npar"))
  rows[[length(rows) + 1]] <- data.frame(family = "confounding", spec = s$nm, flipped = "", acyclic = TRUE,
    df = fm[["df"]], npar = fm[["npar"]], auto_covariances = covs_of(f, paste(s$f, s$t)))
}
tab <- do.call(rbind, rows)
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
write.csv(tab, file.path(TABLES_DIR, "residual_covariances_by_spec.csv"), row.names = FALSE)
cat("\n==== (A) residual covariances per specification (* = explicitly specified replacement) ====\n")
print(tab[, c("family", "spec", "flipped", "acyclic", "df", "npar", "auto_covariances")], row.names = FALSE, right = FALSE)

# (B) cyclic combos: reduced-system diagnostics per intervention target
targets <- c("belief_concern","harm_present","harm_future","weather_risk_prep","social_norms","trust_science","policy_support","politics")
nodes <- all_nodes; n <- length(nodes)
bres <- list()
for (nm in names(scenario_list)) {
  e <- flip_edges(base_edges, scenario_list[[nm]]); if (is_acyclic(e)) next
  f <- fit_spec(e); ss <- lavaan::standardizedSolution(f); ss <- ss[ss$op == "~", ]
  B <- matrix(0, n, n, dimnames = list(nodes, nodes)); for (i in seq_len(nrow(ss))) B[ss$lhs[i], ss$rhs[i]] <- ss$est.std[i]
  rhoB <- max(Mod(eigen(B, only.values = TRUE)$values))
  for (k in targets) {
    o <- setdiff(nodes, k); Bo <- B[o, o]
    bres[[length(bres) + 1]] <- data.frame(combo = nm, rho_full = rhoB, target = lab[k],
      rho_reduced = max(Mod(eigen(Bo, only.values = TRUE)$values)), det_I_minus_Bred = det(diag(n - 1) - Bo),
      min_coef = min(B[B != 0]))
  }
}
btab <- do.call(rbind, bres)
write.csv(btab, file.path(TABLES_DIR, "cyclic_reduced_system_diagnostics.csv"), row.names = FALSE)
cat("\n==== (B) cyclic specifications: reduced-system diagnostics ====\n")
print(aggregate(cbind(rho_reduced, det_I_minus_Bred) ~ combo + rho_full, data = btab, FUN = function(x) c(min = min(x), max = max(x))))
cat("any negative coefficient in cyclic B:", any(btab$min_coef < 0), "\n")
print(btab[, c("combo", "target", "rho_reduced", "det_I_minus_Bred")], row.names = FALSE, digits = 3)
cat("43 DONE\n")
