# =============================================================================
# 37_presentharm_behavior_sensitivity_and_S8.R  (2026-09-26, review round)
#
# (A) Targeted confounding check for harm_present -> climate_behavior.
#     In the joint-endpoint FCI bootstrap this pair was directed-plurality at
#     alpha=.05 (45.8% vs 42.7% bidirected) but bidirected-plurality at
#     alpha=.01 (35.4% vs 44.7%), so it missed the "bidirected at BOTH alphas"
#     rule that defined the 8-edge set in clean_pipeline/15. The 8-edge rule
#     is NOT changed here; this is an additional, separately reported check.
#     Specs (same sample N=870, lavaan MLR, fixed.x=FALSE, do(X=+0.5) exact
#     mean propagation -- identical machinery to clean_pipeline/15):
#       baseline            16-edge working SCM
#       hp_cb_cov           harm_present -> climate_behavior replaced by
#                           harm_present ~~ climate_behavior
#       all8_plus_hp_cb     the 8 bidirected-at-both-alphas edges AND
#                           harm_present -> climate_behavior all replaced by
#                           residual covariances
# (B) Supplement S8 rebuilt: conditioning vs intervening on harm_present,
#     both computed from the covariance IMPLIED BY THE SAME fitted baseline
#     SCM (standardized), so scale/model are identical and the only
#     difference is conditioning vs intervention. belief_concern is shown as
#     the contrast case: it is an unconfounded root in the working model, so
#     the two quantities coincide there by construction.
# Writes pipeline_outputs/tables/presentharm_cb_sensitivity_*.csv and
# pipeline_outputs/tables/S8_conditioning_vs_intervention.csv
# =============================================================================
setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
if (!exists("df_extended") || nrow(df_extended) != 870) {
  source("clean_pipeline/00_config.R"); source("clean_pipeline/01_data_prep.R")
}
source("clean_pipeline/05_scm_intervention_helpers.R")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

drop_edges <- function(edges, f, t) {
  keep <- rep(TRUE, nrow(edges))
  for (i in seq_along(f)) keep <- keep & !(edges$from == f[i] & edges$to == t[i])
  edges[keep, ]
}
syntax_cov <- function(edges, cf = character(0), ct = character(0)) {
  s <- build_lavaan_syntax(edges)
  if (length(cf)) s <- paste(c(s, paste(cf, "~~", ct)), collapse = "\n")
  s
}
conf8 <- tibble::tribble(
  ~from,               ~to,
  "belief_concern",    "harm_future",
  "politics",          "policy_support",
  "social_norms",      "policy_support",
  "trust_science",     "social_norms",
  "harm_present",      "weather_risk_prep",
  "belief_concern",    "weather_risk_prep",
  "social_norms",      "climate_behavior",
  "weather_risk_prep", "climate_behavior"
)
nodes7 <- c("belief_concern", "harm_present", "harm_future", "weather_risk_prep",
            "social_norms", "trust_science", "policy_support")
nodes8 <- c(nodes7, "politics")
combos <- build_intervene_targets(nodes7, max_size = 3)
TOP_TRIPLE <- "harm_present+weather_risk_prep+social_norms"

specs <- list(
  baseline = list(edges = base_edges, cf = character(0), ct = character(0)),
  hp_cb_cov = list(edges = drop_edges(base_edges, "harm_present", "climate_behavior"),
                   cf = "harm_present", ct = "climate_behavior"),
  all8_plus_hp_cb = list(
    edges = drop_edges(base_edges, c(conf8$from, "harm_present"), c(conf8$to, "climate_behavior")),
    cf = c(conf8$from, "harm_present"), ct = c(conf8$to, "climate_behavior"))
)

fits <- list(); sn <- list(); cmb <- list(); fitstats <- list()
for (nm in names(specs)) {
  s <- specs[[nm]]
  fit <- lavaan::sem(syntax_cov(s$edges, s$cf, s$ct), data = df_extended,
                     estimator = "MLR", fixed.x = FALSE)
  stopifnot(lavaan::lavInspect(fit, "converged"))
  fits[[nm]] <- fit
  a <- single_node_ates(fit, s$edges, all_nodes, nodes8)
  sn[[nm]] <- tibble::tibble(spec = nm, node = nodes8, ate = unname(a[nodes8]))
  cmb[[nm]] <- purrr::imap_dfr(combos, function(cb, lab)
    tibble::tibble(spec = nm, target = lab, size = length(cb),
                   ate = target_ate(fit, s$edges, all_nodes, cb)))
  fm <- lavaan::fitMeasures(fit, c("chisq.scaled", "df", "cfi.robust", "rmsea.robust", "srmr", "aic", "bic"))
  fitstats[[nm]] <- tibble::tibble(spec = nm, !!!as.list(round(fm, 4)))
}
sn_all  <- dplyr::bind_rows(sn)
cmb_all <- dplyr::bind_rows(cmb) |> dplyr::group_by(spec) |>
  dplyr::mutate(rank = rank(-ate, ties.method = "min")) |> dplyr::ungroup()
fit_all <- dplyr::bind_rows(fitstats)

# residual correlation that replaces the path
rescor <- purrr::map_dfr(c("hp_cb_cov", "all8_plus_hp_cb"), function(nm) {
  ss <- lavaan::standardizedSolution(fits[[nm]])
  r <- ss[ss$op == "~~" & ((ss$lhs == "harm_present" & ss$rhs == "climate_behavior") |
                           (ss$rhs == "harm_present" & ss$lhs == "climate_behavior")), ]
  tibble::tibble(spec = nm, resid_cor = r$est.std, ci_lo = r$ci.lower, ci_hi = r$ci.upper, p = r$pvalue)
})

sn_wide <- tidyr::pivot_wider(sn_all, names_from = spec, values_from = ate)
top_by_spec <- cmb_all |> dplyr::filter(size == 3) |> dplyr::group_by(spec) |>
  dplyr::slice_max(ate, n = 1, with_ties = FALSE) |> dplyr::ungroup()
watch <- cmb_all |> dplyr::filter(target == TOP_TRIPLE)

write.csv(sn_wide, file.path(TABLES_DIR, "presentharm_cb_sensitivity_singlenode.csv"), row.names = FALSE)
write.csv(cmb_all, file.path(TABLES_DIR, "presentharm_cb_sensitivity_combos.csv"), row.names = FALSE)
write.csv(fit_all, file.path(TABLES_DIR, "presentharm_cb_sensitivity_fit.csv"), row.names = FALSE)
write.csv(rescor, file.path(TABLES_DIR, "presentharm_cb_sensitivity_rescor.csv"), row.names = FALSE)

cat("\n==== (A) single-node ATEs, do(X = +0.5) ====\n"); print(as.data.frame(sn_wide), digits = 3)
cat("\n==== baseline-best triple across specs ====\n"); print(as.data.frame(watch), digits = 3)
cat("\n==== best triple within each spec ====\n"); print(as.data.frame(top_by_spec), digits = 3)
cat("\n==== residual correlation HP~~CB ====\n"); print(as.data.frame(rescor), digits = 3)
cat("\n==== fit ====\n"); print(as.data.frame(fit_all))

# ---- (B) S8 -----------------------------------------------------------------
fb <- fits[["baseline"]]
R  <- stats::cov2cor(lavaan::fitted(fb)$cov)          # model-implied, standardized
cond <- function(x, y = "climate_behavior", v = 0.5) R[y, x] / R[x, x] * v
ates_b <- single_node_ates(fb, base_edges, all_nodes, c("harm_present", "belief_concern"))
std <- lavaan::standardizedSolution(fb)
b <- function(l, r) std$est.std[std$op == "~" & std$lhs == l & std$rhs == r]
s8 <- tibble::tibble(
  variable = c("harm_present", "belief_concern"),
  E_cb_given_X_eq_0.5 = c(cond("harm_present"), cond("belief_concern")),
  E_cb_do_X_eq_0.5    = unname(ates_b[c("harm_present", "belief_concern")])
) |> dplyr::mutate(difference = E_cb_given_X_eq_0.5 - E_cb_do_X_eq_0.5)
# pieces for the harm_present decomposition
pieces <- tibble::tibble(
  quantity = c("b_hp_cb", "b_hp_wr", "b_wr_cb", "b_bc_hp", "b_hf_hp", "b_bc_hf",
               "impl_cor_hp_cb", "impl_cor_bc_hp", "obs_cor_hp_cb"),
  value = c(b("climate_behavior", "harm_present"), b("weather_risk_prep", "harm_present"),
            b("climate_behavior", "weather_risk_prep"), b("harm_present", "belief_concern"),
            b("harm_present", "harm_future"), b("harm_future", "belief_concern"),
            R["harm_present", "climate_behavior"], R["belief_concern", "harm_present"],
            stats::cor(df_extended$harm_present, df_extended$climate_behavior)))
write.csv(s8, file.path(TABLES_DIR, "S8_conditioning_vs_intervention.csv"), row.names = FALSE)
write.csv(pieces, file.path(TABLES_DIR, "S8_components.csv"), row.names = FALSE)
cat("\n==== (B) S8: conditioning vs intervention (model-implied, standardized) ====\n")
print(as.data.frame(s8), digits = 4); print(as.data.frame(pieces), digits = 4)
cat("37 DONE\n")

# ---- (B2) S8 table: all variables, conditioning vs do(harm_present = 0.5) ----
do_hp <- scm_mean_propagate(fb, base_edges, all_nodes, list(harm_present = 0.5))
base0 <- scm_mean_propagate(fb, base_edges, all_nodes, list())
vars_s8 <- c("belief_concern", "harm_future", "weather_risk_prep", "climate_behavior",
             "trust_science", "social_norms", "politics", "policy_support")
s8_full <- tibble::tibble(
  variable = vars_s8,
  conditioning = vapply(vars_s8, function(v) R[v, "harm_present"] / R["harm_present", "harm_present"] * 0.5, numeric(1)),
  intervention = unname(do_hp[vars_s8] - base0[vars_s8]))
write.csv(s8_full, file.path(TABLES_DIR, "S8_presentharm_full.csv"), row.names = FALSE)
cat("\n==== (B2) S8 full table for harm_present ====\n"); print(as.data.frame(s8_full), digits = 4)
