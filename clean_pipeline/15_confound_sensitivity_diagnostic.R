# confounding/existence sensitivity diagnostic -- prompted by the joint pag
# endpoint audit (r_patches/29 + pipeline_outputs/joint_pag_edgetype_audit.csv).
# 8 of the 16 scm edges were bidirected-plurality in the fci bootstrap at
# BOTH alpha=.05 and alpha=.01 -- more consistent with a shared latent cause
# than with the directed arrow the scm currently assumes. this script checks
# what happens to the intervention ates/rankings if those 8 edges are
# represented as residual covariances (x ~~ y) instead of directed paths,
# one at a time and all together, plus two more targeted checks: reversing
# harm_future->harm_present (the one edge that's genuinely direction-
# uncertain, not bidirected-dominated) and dropping harm_future->trust_science
# entirely (its existence itself is thin -- p_adjacent=.657, near the .60
# floor, and 45.6% "absent" in the joint audit at alpha=.01).
#
# this does NOT touch base_edges in 05 or the published scm -- it builds its
# own edge sets locally, on top of the same baseline. nothing here changes
# scm_edges_finalized.csv or any manuscript number by itself.

source("clean_pipeline/05_scm_intervention_helpers.R")

stopifnot(exists("df_extended"))
stopifnot(nrow(df_extended) == 870)

TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- the 8 edges that were bidirected-plurality at BOTH alpha levels ------
# UPDATED 2026-09-11: policy_support/social_norms direction flipped to match
# base_edges in 05 (social_norms->policy_support now, not policy_support->
# social_norms) -- this list is matched against base_edges by exact from/to
# in drop_edges() below, so it has to track base_edges's current direction or
# that edge silently fails to drop (the directed path AND the covariance
# term would both end up in the model instead of one replacing the other).
confound_candidates <- tibble::tribble(
  ~edge_label,                            ~from,               ~to,
  "belief_concern -> harm_future",        "belief_concern",    "harm_future",
  "politics -> policy_support",           "politics",          "policy_support",
  "social_norms -> policy_support",       "social_norms",      "policy_support",
  "trust_science -> social_norms",        "trust_science",     "social_norms",
  "harm_present -> weather_risk_prep",    "harm_present",      "weather_risk_prep",
  "belief_concern -> weather_risk_prep",  "belief_concern",    "weather_risk_prep",
  "social_norms -> climate_behavior",     "social_norms",      "climate_behavior",
  "weather_risk_prep -> climate_behavior","weather_risk_prep", "climate_behavior"
)
stopifnot(nrow(confound_candidates) == 8)

# note: two nodes lose ALL their edges in the all-8-together spec, not just
# one:
#  - social_norms: its one incoming edge (trust_science->) AND both outgoing
#    edges (->policy_support, ->climate_behavior) are all in this list --
#    same conclusion as before the direction flip, just 1-in/2-out now
#    instead of 2-in/1-out.
#  - weather_risk_prep: both incoming (harm_present->, belief_concern->) and
#    its one outgoing edge (->climate_behavior) are all in this list,
#    unaffected by the direction flip above.
# both end up fully disconnected from the mean-propagation graph in that one
# spec -- exogenous, mean 0, and incapable of reaching climate_behavior
# through any ~ path. their single-node ates there are exactly 0 by
# construction, not an empirical result, and the top triple's ate in that
# spec is driven entirely by harm_present's own direct edge -- flagged again
# below where it prints.

# ---- helpers ----------------------------------------------------------

# edges, minus the given from/to pairs
drop_edges <- function(edges, pairs_from, pairs_to) {
  drop_flag <- rep(FALSE, nrow(edges))
  for (i in seq_along(pairs_from)) {
    drop_flag <- drop_flag | (edges$from == pairs_from[i] & edges$to == pairs_to[i])
  }
  edges[!drop_flag, ]
}

# regression syntax from the (already-reduced) edges, plus explicit residual
# covariance lines for whichever pairs were dropped as confound-replacements
build_syntax_with_cov <- function(edges_kept, cov_from = character(0), cov_to = character(0)) {
  reg_syntax <- build_lavaan_syntax(edges_kept)
  if (length(cov_from) == 0) return(reg_syntax)
  cov_lines <- paste(cov_from, "~~", cov_to)
  paste(c(reg_syntax, cov_lines), collapse = "\n")
}

intervene_nodes_7 <- c("belief_concern", "harm_present", "harm_future",
                        "weather_risk_prep", "social_norms", "trust_science",
                        "policy_support")
combo_targets <- build_intervene_targets(intervene_nodes_7, max_size = 3)
stopifnot(length(combo_targets) == 63)

# single-node reporting uses all 8 (7 non-political + politics), same set fig
# 7 / 07_intervention_ates_8node.R uses -- politics is a benchmark covariate
# here, not a real intervention target, same caveat 07 already carries.
# combos above stay at 7 (politics excluded), matching 08.
intervene_nodes_8 <- c(intervene_nodes_7, "politics")

# fit one spec and pull everything we want to compare, in one place so no
# two specs can accidentally use different node sets/estimator settings
run_spec <- function(spec_name, description, edges_kept, cov_from = character(0), cov_to = character(0)) {
  if (!is_acyclic(edges_kept)) {
    message(spec_name, ": cyclic, skipped -- ", description)
    return(list(spec = spec_name, description = description, status = "cyclic_skipped",
                singlenode = NULL, combo = NULL))
  }
  model_syntax <- build_syntax_with_cov(edges_kept, cov_from, cov_to)
  fit <- tryCatch(
    lavaan::sem(model_syntax, data = df_extended, estimator = "MLR", fixed.x = FALSE),
    error = function(err) NULL
  )
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
    message(spec_name, ": fit failed/non-converged -- ", description)
    return(list(spec = spec_name, description = description, status = "fit_failed",
                singlenode = NULL, combo = NULL))
  }

  sn_ates <- single_node_ates(fit, edges_kept, all_nodes, intervene_nodes_8)
  singlenode <- tibble::tibble(
    spec = spec_name, node = intervene_nodes_8,
    ate_climate_behavior = unname(sn_ates[intervene_nodes_8])
  )

  combo_rows <- purrr::imap(combo_targets, function(combo, label) {
    tibble::tibble(
      spec = spec_name, target_label = label, target_size = length(combo),
      ate_climate_behavior = target_ate(fit, edges_kept, all_nodes, combo)
    )
  }) |> dplyr::bind_rows()

  list(spec = spec_name, description = description, status = "ok",
       singlenode = singlenode, combo = combo_rows)
}

# ---- the specs ----------------------------------------------------------
# 1 baseline + 8 one-at-a-time confound-replacements + 1 all-8-together +
# 1 direction-reversal + 1 existence-drop = 12

specs <- list()

specs[["baseline"]] <- list(
  description = "current 16-edge working SCM, no changes",
  edges = base_edges, cov_from = character(0), cov_to = character(0)
)

for (i in seq_len(nrow(confound_candidates))) {
  nm <- paste0("confound_1of8_", confound_candidates$from[i], "_", confound_candidates$to[i])
  specs[[nm]] <- list(
    description = paste0(confound_candidates$edge_label[i], " replaced with residual covariance"),
    edges = drop_edges(base_edges, confound_candidates$from[i], confound_candidates$to[i]),
    cov_from = confound_candidates$from[i], cov_to = confound_candidates$to[i]
  )
}

specs[["confound_all8"]] <- list(
  description = "all 8 bidirected-plurality edges replaced with residual covariances at once",
  edges = drop_edges(base_edges, confound_candidates$from, confound_candidates$to),
  cov_from = confound_candidates$from, cov_to = confound_candidates$to
)

specs[["reverse_harmfuture_harmpresent"]] <- list(
  description = "harm_future <-> harm_present reversed (direction-uncertain edge, not bidirected-dominated)",
  edges = flip_edges(base_edges, which(flip_candidates$from == "harm_future" & flip_candidates$to == "harm_present")),
  cov_from = character(0), cov_to = character(0)
)

specs[["drop_harmfuture_trustscience"]] <- list(
  description = "harm_future -> trust_science dropped entirely, no replacement (existence sensitivity: p_adjacent=.657, 45.6% absent at alpha=.01)",
  edges = drop_edges(base_edges, "harm_future", "trust_science"),
  cov_from = character(0), cov_to = character(0)
)

stopifnot(length(specs) == 12)

# ---- run them all ---------------------------------------------------------

results <- purrr::imap(specs, function(s, nm) {
  run_spec(nm, s$description, s$edges, s$cov_from, s$cov_to)
})

status_tbl <- tibble::tibble(
  spec = names(results),
  description = purrr::map_chr(results, "description"),
  status = purrr::map_chr(results, "status")
)
cat("\n--- confound/existence sensitivity: 12 specs ---\n")
print(as.data.frame(status_tbl), row.names = FALSE)

singlenode_all <- results |> purrr::map("singlenode") |> purrr::compact() |> dplyr::bind_rows()
combo_all <- results |> purrr::map("combo") |> purrr::compact() |> dplyr::bind_rows()

write.csv(singlenode_all, file.path(TABLES_DIR, "confound_sensitivity_singlenode_ate.csv"), row.names = FALSE)
write.csv(combo_all, file.path(TABLES_DIR, "confound_sensitivity_combo_ate.csv"), row.names = FALSE)
write.csv(status_tbl, file.path(TABLES_DIR, "confound_sensitivity_status.csv"), row.names = FALSE)

# ---- quick readout: single-node ates, baseline vs each spec ---------------
sn_wide <- singlenode_all |>
  tidyr::pivot_wider(names_from = spec, values_from = ate_climate_behavior)
cat("\n--- single-node ates on climate_behavior, baseline vs each spec ---\n")
print(as.data.frame(sn_wide), row.names = FALSE)

# ---- quick readout: does the top triple/pair change rank? -----------------
baseline_top_triple <- "harm_present+weather_risk_prep+social_norms"
baseline_top_pair    <- "belief_concern+weather_risk_prep"

top_watch <- combo_all |>
  dplyr::filter(target_label %in% c(baseline_top_triple, baseline_top_pair,
                                     "harm_present+social_norms")) |>
  tidyr::pivot_wider(names_from = spec, values_from = ate_climate_behavior, id_cols = target_label)
cat("\n--- baseline's top triple/pair ate, tracked across every spec ---\n")
print(as.data.frame(top_watch), row.names = FALSE)

cat("\n--- for each spec, where does the baseline top triple (",
    baseline_top_triple, ") rank among all 63 combos? ---\n", sep = "")
rank_watch <- combo_all |>
  dplyr::group_by(spec) |>
  dplyr::mutate(rank = rank(-ate_climate_behavior, ties.method = "min")) |>
  dplyr::filter(target_label == baseline_top_triple) |>
  dplyr::ungroup() |>
  dplyr::select(spec, ate_climate_behavior, rank)
print(as.data.frame(rank_watch), row.names = FALSE)

cat("\nWrote", file.path(TABLES_DIR, "confound_sensitivity_singlenode_ate.csv"), ",",
    file.path(TABLES_DIR, "confound_sensitivity_combo_ate.csv"), "and",
    file.path(TABLES_DIR, "confound_sensitivity_status.csv"), "\n")
cat("(social_norms's AND weather_risk_prep's ates in confound_all8 are\n",
    "mechanically 0 -- both lose every edge touching them in that one spec,\n",
    "see the note above confound_candidates. not a finding, just the shape\n",
    "of that particular spec -- the one-at-a-time specs are the ones that\n",
    "isolate each edge's own contribution.)\n", sep = "")
