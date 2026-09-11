# =============================================================================
# 32_deterministic_8node_ates.R  (NEW, 2026-09-07)
#
# Figure 7 has always shown single-node intervention effects for exactly 6
# nodes (the ones in intervene_nodes across scripts 02/30/31: belief_concern,
# harm_present, weather_risk_prep, social_norms, trust_science,
# policy_support). But the working SCM has 8 non-outcome nodes total -- those
# 6 plus harm_future and politics -- and neither of the other two has ever had
# a single-node ATE computed. Extending the figure to all 8 needs those two
# additional numbers, computed the same exact-mean-propagation way as
# everything else in the intervention analysis (see script 30's header for why
# mean propagation is exact rather than simulated for this linear-recursive
# SCM).
#
# harm_future is included with no caveat -- it's a belief/attitude construct
# like present harm, and "do(harm_future=0.5)" has the same interpretation as
# "do(harm_present=0.5)": a hypothetical messaging intervention that shifts
# perceived future climate harm by 0.5 SD.
#
# politics is included too, but flagged: it's a demographic/ideological
# covariate, not something any real intervention manipulates. "do(politics=
# 0.5)" is not a policy lever the way the other 7 are -- it's included here
# only as a descriptive benchmark for comparison, and the figure/caption
# downstream need to say so explicitly rather than let it sit next to the
# other 7 rows implying equivalence.
#
# what stays the same: does not touch the 41-target combo analysis (scripts 02/31) or its outputs
# at all -- politics and harm_future stay out of intervene_nodes there, since
# that 41-target set is specifically the "policy-relevant nodes" combinatorial
# space referenced throughout the text and Table 3. This script only adds two
# more single-node rows for Figure 7.
#
# output:
# tables/orientation_enumeration_ate_deterministic_8node.csv -- same schema as
# tables/orientation_enumeration_ate_deterministic.csv (scenario, node,
# ate_climate_behavior), but with harm_future and politics rows appended (16
# scenarios x 8 nodes = 128 rows total, up from 96).
#
# run this after script 30 (or standalone -- it re-fits its own 16 lavaan
# models, doesn't depend on script 30 having been run in the same session).
# Needs df_extended already loaded in the session.
# Then re-run r_patches/10_figure7_uncertainty_pub.R, which has been updated
# to read the 8-node file and plot all 8 rows.
# =============================================================================

dir.create("tables", showWarnings = FALSE)

# --- 0: sanity check ---------------------------------------------------------
stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

# --- 1: baseline 16-edge working-SCM skeleton -- VERBATIM from script 30 ----
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
stopifnot(nrow(base_edges) == 16)

flip_candidates <- tibble::tribble(
  ~edge_label,                          ~from,            ~to,
  "politics -> belief_concern",          "politics",       "belief_concern",
  "policy_support -> social_norms",      "policy_support",  "social_norms",
  "social_norms -> climate_behavior",    "social_norms",    "climate_behavior",
  "harm_future -> harm_present",         "harm_future",     "harm_present"
)

n_flip <- nrow(flip_candidates)
stopifnot(n_flip == 4)

# --- 2: helpers -- VERBATIM from script 30/31 -------------------------------

flip_edges <- function(edges, flip_set_idx) {
  e <- edges
  if (length(flip_set_idx) > 0) {
    for (i in flip_set_idx) {
      row <- flip_candidates[i, ]
      match_idx <- which(e$from == row$from & e$to == row$to)
      if (length(match_idx) == 1) {
        tmp <- e$from[match_idx]
        e$from[match_idx] <- e$to[match_idx]
        e$to[match_idx]   <- tmp
      }
    }
  }
  e
}

is_acyclic <- function(edges) {
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  igraph::is_dag(g)
}

build_lavaan_syntax <- function(edges) {
  edges |>
    dplyr::group_by(to) |>
    dplyr::summarise(rhs = paste(from, collapse = " + "), .groups = "drop") |>
    dplyr::mutate(line = paste(to, "~", rhs)) |>
    dplyr::pull(line) |>
    paste(collapse = "\n")
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

# --- 3: enumerate all 16 orientation combinations -- VERBATIM from script 30

scenario_list <- setNames(
  lapply(0:(2^n_flip - 1), function(k) which(as.logical(intToBits(k)[1:n_flip]))),
  paste0("combo_", 0:(2^n_flip - 1))
)

all_nodes <- unique(c(base_edges$from, base_edges$to))

# --- 4: the 8 single-node targets -- the 6 already computed by script 30, ---
# plus harm_future and politics (see header for why those two were missing).
node_targets <- c(
  "belief_concern", "harm_present", "weather_risk_prep",
  "social_norms", "trust_science", "policy_support",
  "harm_future", "politics"
)
stopifnot(length(node_targets) == 8)

# --- 5: refit all 16 specifications, compute all 8 single-node ATEs ---------
det_8node_out <- scenario_list |>
  purrr::imap(function(flip_idx, scen_name) {
    e <- flip_edges(base_edges, flip_idx)
    if (!is_acyclic(e)) {
      message(scen_name, ": cyclic after flip, skipped.")
      return(NULL)
    }
    model_syntax <- build_lavaan_syntax(e)
    fit <- tryCatch(
      lavaan::sem(model_syntax, data = df_extended, estimator = "MLR", fixed.x = FALSE),
      error = function(err) NULL
    )
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      message(scen_name, ": fit failed or did not converge, skipped.")
      return(NULL)
    }
    base_means <- scm_mean_propagate(fit, e, all_nodes, list())
    base_cb    <- base_means[["climate_behavior"]]
    purrr::map_dfr(node_targets, function(node) {
      do_means <- scm_mean_propagate(fit, e, all_nodes, setNames(list(0.5), node))
      tibble::tibble(
        scenario = scen_name,
        node = node,
        ate_climate_behavior = do_means[["climate_behavior"]] - base_cb
      )
    })
  }) |>
  dplyr::bind_rows()

stopifnot(nrow(det_8node_out) == 16 * 8)  # 128 rows, all 16 specs converged

# --- 6: write output ---------------------------------------------------------
if (requireNamespace("readr", quietly = TRUE)) {
  readr::write_csv(det_8node_out, "tables/orientation_enumeration_ate_deterministic_8node.csv")
} else {
  write.csv(det_8node_out, "tables/orientation_enumeration_ate_deterministic_8node.csv", row.names = FALSE)
}

message("Saved tables/orientation_enumeration_ate_deterministic_8node.csv (",
        nrow(det_8node_out), " rows: 16 specs x 8 nodes). ",
        "Re-run r_patches/10_figure7_uncertainty_pub.R next to regenerate Figure 7 with all 8 rows.")

# Quick sanity print of the two new nodes' baseline (combo_0) values, so it's
# obvious at a glance whether this ran correctly before moving on to the figure.
det_8node_out |>
  dplyr::filter(scenario == "combo_0", node %in% c("harm_future", "politics")) |>
  print()
