# =============================================================================
# 27_diagnose_hf_hp_edge.R
#
# Purpose: diagnose whether harm_future -> harm_present's sign disagreement
# across alpha (found by script 21: asymmetry -.108 at alpha=.05, +.261 at
# alpha=.01) is a stable finding or a threshold/seed artifact, before deciding
# whether to promote this edge into the structural-orientation uncertainty
# set. This edge is far more consequential than the existing three (its
# coefficient, beta=.636, is the single largest in the 16-edge model), so the
# team wants it checked properly rather than mechanically promoted by a rule
# change. Five checks, run in order, matching the team's own specification:
#
#   1. Endpoint mark decomposition (adjacency, arrowhead/tail/circle
#      proportions at each endpoint, resulting asymmetry) separately at
#      alpha=.05 and alpha=.01 -- from the already-saved per-alpha bootstrap
#      array, no new computation.
#   2. Seed sensitivity: repeat the bootstrap (default 8 independent seeds x
#      1,000 resamples, at BOTH alpha levels each) and check whether the sign
#      pattern (negative at .05, positive at .01) recurs. This is the
#      expensive part of this script -- 8 x 2 x 1,000 = 16,000 FCI+PC fits,
#      parallelized the same way as scripts 02/23/the .qmd's own bootstrap.
#   3. Single-run FCI endpoint marks for this pair at both alpha levels, from
#      the already-fitted fci_ext_05/fci_ext_01 objects -- no new computation.
#   4. Acyclicity check: does reversing ONLY this edge in the 16-edge working
#      SCM still yield a DAG?
#   5. Substantive consequence: fit baseline (harm_future -> harm_present) vs.
#      reversed (harm_present -> harm_future) as the ONLY change to the
#      16-edge SCM; compare fit indices and the six single-node ATEs.
#
# Decision rule (team's own, stated in advance): if the sign reversal is
# stable across seeds and both alphas retain strong adjacency, treat the edge
# as unresolved under the cross-alpha criterion (4-edge/16-specification
# enumeration). If the alpha=.05 negative sign disappears across seeds or
# hovers near zero while alpha=.01 stays consistently positive, describe it as
# threshold-sensitive orientation evidence and do NOT promote it -- keep the
# rule from being driven by one bootstrap realization.
#
# HOW TO USE: run within the same .qmd session as scripts 21-23 (needs
# run_one_ext, agg_ext, context_idx, node_order_ext, alphas, fci_ext_05,
# fci_ext_01, df_extended all already in the session).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lavaan)
  library(igraph)
  library(purrr)
  library(tibble)
  library(furrr)
  library(future)
})

stopifnot(exists("run_one_ext"), exists("agg_ext"), exists("context_idx"),
          exists("node_order_ext"), exists("alphas"), exists("df_extended"),
          exists("fci_ext_05"), exists("fci_ext_01"))

HF <- "harm_future"; HP <- "harm_present"
mp_path <- "pipeline_outputs/mark_props_ext_by_alpha.rds"
if (!file.exists(mp_path)) stop("Run script 21's precondition first (re-render the ext-bootstrap chunk).")
mark_props_ext <- readRDS(mp_path)

# ============================================================================
# 1. Endpoint mark decomposition at each alpha (no new computation)
# ============================================================================
cat("\n============================================================\n")
cat("1. Endpoint mark decomposition, harm_future <-> harm_present\n")
cat("============================================================\n")

decompose_alpha <- function(alph) {
  props <- mark_props_ext[[alph]]$fci
  tibble::tibble(
    alpha              = alph,
    adjacency_pct_HFview = round(100 * (1 - props[HP, HF, "N"]), 1),
    adjacency_pct_HPview = round(100 * (1 - props[HF, HP, "N"]), 1),
    arrow_at_HF        = round(props[HP, HF, ">"], 3),
    arrow_at_HP        = round(props[HF, HP, ">"], 3),
    tail_at_HF         = round(props[HP, HF, "-"], 3),
    tail_at_HP         = round(props[HF, HP, "-"], 3),
    circle_at_HF       = round(props[HP, HF, "o"], 3),
    circle_at_HP       = round(props[HF, HP, "o"], 3),
    asymmetry_HP_minus_HF = round(props[HF, HP, ">"] - props[HP, HF, ">"], 3)
  )
}
decomp_table <- purrr::map_dfr(names(mark_props_ext), decompose_alpha)
print(as.data.frame(decomp_table), row.names = FALSE)
cat("(adjacency_pct reported from both readings as a sanity check -- they should match;\n",
    "a mismatch would indicate a data/indexing problem worth flagging.)\n", sep = "")

# ============================================================================
# 2. Seed sensitivity: N_REPLICATES independent 1,000-resample runs per alpha
# ============================================================================
N_REPLICATES <- 8L
N_BOOT_REPLICATE <- 1000L
cat("\n============================================================\n")
cat("2. Seed sensitivity:", N_REPLICATES, "replicate runs x", N_BOOT_REPLICATE,
    "resamples x 2 alpha levels\n")
cat("============================================================\n")

seed_results <- purrr::map_dfr(seq_len(N_REPLICATES), function(rep_i) {
  purrr::map_dfr(names(alphas), function(alph) {
    base_seed <- 9000L + rep_i * 100L
    plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))
    res <- furrr::future_map(
      seq_len(N_BOOT_REPLICATE), run_one_ext,
      data = agg_ext, nms = node_order_ext, alpha = alphas[[alph]],
      ctx_idx = context_idx,
      .options = furrr::furrr_options(seed = base_seed), .progress = FALSE
    )
    plan(sequential)
    p_ext <- length(node_order_ext)
    fci_c <- array(0L, dim = c(p_ext, p_ext, 4L),
                    dimnames = list(node_order_ext, node_order_ext, c("N", "o", ">", "-")))
    for (r in res) if (!is.null(r)) fci_c <- fci_c + r$fci
    props <- fci_c / N_BOOT_REPLICATE
    tibble::tibble(
      replicate = rep_i, alpha = alph,
      asymmetry = round(props[HF, HP, ">"] - props[HP, HF, ">"], 3),
      adjacency_pct = round(100 * (1 - props[HF, HP, "N"]), 1)
    )
  })
})

cat("\n-- Per-replicate asymmetry (harm_future -> harm_present), by alpha --\n")
print(as.data.frame(seed_results), row.names = FALSE)

seed_summary <- seed_results |>
  dplyr::group_by(alpha) |>
  dplyr::summarise(
    n_negative = sum(asymmetry < 0), n_positive = sum(asymmetry > 0),
    mean_asymmetry = round(mean(asymmetry), 3),
    min_asymmetry = round(min(asymmetry), 3), max_asymmetry = round(max(asymmetry), 3),
    .groups = "drop"
  )
cat("\n-- Summary across", N_REPLICATES, "replicates --\n")
print(as.data.frame(seed_summary), row.names = FALSE)
cat("\nDecision cue: if alpha=.05 is consistently negative (n_negative close to",
    N_REPLICATES, ") and alpha=.01 is consistently positive (n_positive close to",
    N_REPLICATES, "), the sign pattern is stable -- treat as a real cross-alpha\n",
    "disagreement. If alpha=.05's sign is inconsistent across replicates or hovers\n",
    "near zero, this looks threshold-sensitive rather than a stable finding.\n", sep = "")

# ============================================================================
# 3. Single-run FCI endpoint marks at both alpha (no new computation)
# ============================================================================
cat("\n============================================================\n")
cat("3. Single-run FCI endpoint marks, harm_future <-> harm_present\n")
cat("============================================================\n")

decode_fci_mark <- function(x) {
  dplyr::case_when(x == 1 ~ "circle", x == 2 ~ "arrow", x == 3 ~ "tail", TRUE ~ "none (no edge)")
}
hf_i <- which(node_order_ext == HF); hp_i <- which(node_order_ext == HP)

for (obj_name in c("fci_ext_05", "fci_ext_01")) {
  amat <- get(obj_name)@amat
  mark_at_HF <- decode_fci_mark(amat[hp_i, hf_i])  # mark at HF end, read from HP
  mark_at_HP <- decode_fci_mark(amat[hf_i, hp_i])  # mark at HP end, read from HF
  cat(obj_name, ": mark at harm_future end =", mark_at_HF,
      "| mark at harm_present end =", mark_at_HP, "\n")
}
cat("(A clean HF->HP relation single-run would show tail at HF, arrow at HP;\n",
    "circle at either end indicates the single-run PAG itself left this pair\n",
    "partially or fully unresolved -- context for whether the bootstrap sign\n",
    "disagreement reflects genuine instability or pooled mass around one PAG state.)\n", sep = "")

# ============================================================================
# 4. Acyclicity check on reversal
# ============================================================================
cat("\n============================================================\n")
cat("4. Acyclicity check: reverse ONLY harm_future -> harm_present\n")
cat("============================================================\n")

base_edges <- tibble::tribble(
  ~from,               ~to,
  "politics",          "belief_concern",
  "belief_concern",    "harm_future",
  "belief_concern",    "harm_present",
  "harm_future",       "harm_present",
  "belief_concern",    "trust_science",
  "harm_future",       "trust_science",
  "belief_concern",    "policy_support",
  "trust_science",     "policy_support",
  "politics",          "policy_support",
  "policy_support",    "social_norms",
  "trust_science",     "social_norms",
  "harm_present",      "weather_risk_prep",
  "belief_concern",    "weather_risk_prep",
  "harm_present",      "climate_behavior",
  "social_norms",      "climate_behavior",
  "weather_risk_prep", "climate_behavior"
)
reversed_edges <- base_edges |>
  dplyr::mutate(
    tmp_from = ifelse(from == HF & to == HP, HP, from),
    tmp_to   = ifelse(from == HF & to == HP, HF, to)
  ) |> dplyr::select(from = tmp_from, to = tmp_to)

is_acyclic <- function(edges) igraph::is_dag(igraph::graph_from_data_frame(edges, directed = TRUE))
cat("Baseline (harm_future -> harm_present) acyclic:", is_acyclic(base_edges), "\n")
cat("Reversed (harm_present -> harm_future) acyclic:", is_acyclic(reversed_edges), "\n")

# ============================================================================
# 5. Substantive consequence: baseline vs. reversed fit + six ATEs
# ============================================================================
cat("\n============================================================\n")
cat("5. Substantive consequence: baseline vs. reversed fit + six ATEs\n")
cat("============================================================\n")

all_nodes <- unique(c(base_edges$from, base_edges$to))
intervene_nodes <- c("belief_concern", "harm_present", "weather_risk_prep",
                      "social_norms", "trust_science", "policy_support")

build_lavaan_syntax <- function(edges) {
  edges |> dplyr::group_by(to) |>
    dplyr::summarise(rhs = paste(from, collapse = " + "), .groups = "drop") |>
    dplyr::mutate(line = paste(to, "~", rhs)) |> dplyr::pull(line) |> paste(collapse = "\n")
}
topo_order <- function(edges) {
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  igraph::topo_sort(g, mode = "out") |> names()
}
scm_mean_propagate <- function(fit, edges, all_nodes, intervene = list()) {
  ord <- union(topo_order(edges), all_nodes)
  std <- lavaan::standardizedSolution(fit) |> dplyr::filter(op == "~")
  means <- setNames(numeric(length(ord)), ord)
  for (node in ord) {
    if (!is.null(intervene[[node]])) { means[[node]] <- intervene[[node]]; next }
    parents <- edges$from[edges$to == node]
    if (length(parents) == 0) { means[[node]] <- 0; next }
    coefs <- std |> dplyr::filter(lhs == node, rhs %in% parents)
    beta <- setNames(coefs$est.std, coefs$rhs)
    means[[node]] <- sum(vapply(parents, function(p) beta[[p]] * means[[p]], numeric(1)))
  }
  means
}
single_node_ates <- function(fit, edges, all_nodes, nodes) {
  base_cb <- scm_mean_propagate(fit, edges, all_nodes, list())[["climate_behavior"]]
  purrr::map_dbl(nodes, function(nd) {
    do_cb <- scm_mean_propagate(fit, edges, all_nodes, setNames(list(0.5), nd))[["climate_behavior"]]
    do_cb - base_cb
  }) |> setNames(nodes)
}
fit_and_summarize <- function(edges, label) {
  fit <- tryCatch(
    lavaan::sem(build_lavaan_syntax(edges), data = df_extended, estimator = "MLR", fixed.x = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
    cat(label, ": FAILED/non-converged\n")
    return(NULL)
  }
  fm <- lavaan::fitmeasures(fit)
  get_fm <- function(nms) { hit <- nms[nms %in% names(fm)]; if (length(hit) == 0) NA_real_ else unname(fm[hit[1]]) }
  fit_row <- tibble::tibble(
    spec = label,
    cfi = round(get_fm(c("cfi.robust", "cfi")), 3), tli = round(get_fm(c("tli.robust", "tli")), 3),
    rmsea = round(get_fm(c("rmsea.robust", "rmsea")), 3), srmr = round(get_fm("srmr"), 3),
    aic = round(get_fm("aic"), 1), bic = round(get_fm("bic"), 1)
  )
  ates <- single_node_ates(fit, edges, all_nodes, intervene_nodes)
  list(fit_row = fit_row, ates = ates)
}

baseline_result <- fit_and_summarize(base_edges, "baseline (HF->HP)")
reversed_result  <- fit_and_summarize(reversed_edges, "reversed (HP->HF)")

cat("\n-- Fit comparison --\n")
print(as.data.frame(dplyr::bind_rows(baseline_result$fit_row, reversed_result$fit_row)), row.names = FALSE)

cat("\n-- Six single-node ATEs, baseline vs. reversed --\n")
ate_compare <- tibble::tibble(
  node = intervene_nodes,
  ate_baseline = round(unname(baseline_result$ates[intervene_nodes]), 4),
  ate_reversed = round(unname(reversed_result$ates[intervene_nodes]), 4)
) |> dplyr::mutate(diff = round(ate_reversed - ate_baseline, 4))
print(as.data.frame(ate_compare), row.names = FALSE)

dir.create("pipeline_outputs", showWarnings = FALSE)
write.csv(decomp_table, "pipeline_outputs/hf_hp_endpoint_decomposition.csv", row.names = FALSE)
write.csv(seed_results,  "pipeline_outputs/hf_hp_seed_sensitivity.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(baseline_result$fit_row, reversed_result$fit_row),
          "pipeline_outputs/hf_hp_reversal_fit_comparison.csv", row.names = FALSE)
write.csv(ate_compare, "pipeline_outputs/hf_hp_reversal_ate_comparison.csv", row.names = FALSE)
cat("\nWrote pipeline_outputs/hf_hp_{endpoint_decomposition,seed_sensitivity,",
    "reversal_fit_comparison,reversal_ate_comparison}.csv\n", sep = "")
