# the missing piece for the pag-level intervention question: there's no
# saved single-run pag for the EXTENDED 9-node network anywhere in the repo.
# fci_05/fci_01 in 03 are main-network-only (that's figure S4, 8 nodes, no
# climate_behavior) -- the extended network only has the bootstrap-aggregated
# fci_props_ext (proportions across resamples), never a discrete single pag.
# pag2mag/lv-ida both need a real amat.pag (0/1/2/3 codes), not a proportion
# array, so step one is computing that fit -- same fci() call as run_one_ext
# inside 03's bootstrap loop, just on the real (non-resampled) data once, at
# BOTH alpha levels (not just .05 -- the rest of this paper has already shown
# structural conclusions can depend on alpha, no reason this would be exempt).
#
# on the contextVars/jci="1" args below: these are NOT saying the 8 wave1-4
# attitude variables are "context" in JCI's usual environmental-indicator
# sense. all jci="1" does here is encode an exogeneity restriction -- the
# system variable (climate_behavior, wave 5) cannot cause a context variable.
# passing the 8 earlier-wave variables as contextVars is just the mechanism
# for the temporal constraint we actually want (later behavior can't cause
# earlier attitudes), nothing substantive about "context" beyond that.
#
# needs 01+02+03 already run live in this session (agg_ext, node_order_ext,
# context_idx come straight from 03 -- same requirement 14 already has).
#
# the lv-ida part below (listMags/lv.ida) needs
# https://github.com/dmalinsk/lv-ida downloaded separately -- it's Malinsky's
# own standalone code, not part of pcalg, so it's not something this repo can
# vendor automatically. grab lvida.R from that repo and point LVIDA_PATH at
# it. if it's not there, this script still computes and saves both pags and
# reports the edge-type table, then stops before the mag/lv-ida part with a
# message instead of erroring.

source("clean_pipeline/00_config.R")
suppressPackageStartupMessages(library(pcalg))

stopifnot(exists("agg_ext"), exists("node_order_ext"), exists("context_idx"))

LVIDA_PATH <- "lv-ida/lvida.R"  # change this if you put it somewhere else

# ---- step 1: the single-run extended pag, both alphas ---------------------
# same suffStat + same contextVars/jci/selectionBias args as run_one_ext's
# fci_b call in 03 -- if these don't match, this isn't the same pag the rest
# of the pipeline (joint audit, scm finalize) is built on.
suffStat_ext_gauss <- list(C = cor(as.matrix(agg_ext)), n = nrow(agg_ext))

fci_ext_05 <- pcalg::fci(suffStat = suffStat_ext_gauss, indepTest = pcalg::gaussCItest,
                          alpha = .05, labels = node_order_ext,
                          contextVars = context_idx, jci = "1",
                          selectionBias = FALSE, verbose = FALSE)
fci_ext_01 <- pcalg::fci(suffStat = suffStat_ext_gauss, indepTest = pcalg::gaussCItest,
                          alpha = .01, labels = node_order_ext,
                          contextVars = context_idx, jci = "1",
                          selectionBias = FALSE, verbose = FALSE)

saveRDS(fci_ext_05, file.path(OUTPUT_DIR, "fci_ext_05_singlerun.rds"))
saveRDS(fci_ext_01, file.path(OUTPUT_DIR, "fci_ext_01_singlerun.rds"))
cat("Saved the extended-network single-run pags -- fci_ext_05_singlerun.rds,",
    "fci_ext_01_singlerun.rds, both in", OUTPUT_DIR, "\n")

# ---- step 2: the actual edge type for all 16 scm edges, both pags ---------
# same 7 categories as r_patches/29's joint bootstrap audit (from_to, to_from,
# bidirected, from_ocirc_to, to_ocirc_from, circle_circle, absent) so this
# table lines up column-for-column against joint_pag_edgetype_audit.csv --
# this is the discrete single-fit answer, that csv is the bootstrap-resample
# answer, and the point is to compare them directly, not just count circles.
mark_name <- function(code) c("no_edge", "circle", "arrowhead", "tail")[code + 1]

classify_relation <- function(mark_at_from, mark_at_to) {
  dplyr::case_when(
    mark_at_from == "no_edge" | mark_at_to == "no_edge" ~ "absent",
    mark_at_from == "tail"      & mark_at_to == "arrowhead" ~ "from_to",
    mark_at_from == "arrowhead" & mark_at_to == "tail"      ~ "to_from",
    mark_at_from == "arrowhead" & mark_at_to == "arrowhead" ~ "bidirected",
    mark_at_from == "circle"    & mark_at_to == "arrowhead" ~ "from_ocirc_to",
    mark_at_from == "arrowhead" & mark_at_to == "circle"    ~ "to_ocirc_from",
    mark_at_from == "circle"    & mark_at_to == "circle"    ~ "circle_circle",
    TRUE ~ paste0("other(", mark_at_from, "/", mark_at_to, ")")  # shouldn't
    # happen in a valid pag (e.g. tail/tail or tail/circle) -- if this shows
    # up, don't ignore it, it means one of the two mark reads above is wrong.
  )
}

relation_arrow <- c(absent = "(absent)", from_to = "->", to_from = "<-",
                     bidirected = "<->", from_ocirc_to = "o->",
                     to_ocirc_from = "<-o", circle_circle = "o-o")

edge_mark_report <- function(fci_fit, alpha_label) {
  am <- fci_fit@amat
  rows <- lapply(seq_len(nrow(base_edges)), function(i) {
    f <- base_edges$from[i]; t <- base_edges$to[i]
    # amat.pag convention: the mark code is keyed by the COLUMN index, so
    # am[i,j] is the mark AT NODE j, not at i (checked against pcalg's own
    # amatType doc -- easy to get backwards, so it's spelled out here)
    m_from <- mark_name(am[t, f])  # column = from -> mark at "from"
    m_to   <- mark_name(am[f, t])  # column = to   -> mark at "to"
    rel <- classify_relation(m_from, m_to)
    tibble::tibble(
      alpha = alpha_label, from = f, to = t,
      mark_at_from = m_from, mark_at_to = m_to,
      relation = rel,
      arrow = paste(f, if (rel %in% names(relation_arrow)) relation_arrow[[rel]] else rel, t)
    )
  })
  dplyr::bind_rows(rows)
}

# base_edges comes from 05_scm_intervention_helpers.R -- source it if not
# already live (only need the tibble, not the whole lavaan/igraph setup)
if (!exists("base_edges")) source("clean_pipeline/05_scm_intervention_helpers.R")

marks_05 <- edge_mark_report(fci_ext_05, "0.05")
marks_01 <- edge_mark_report(fci_ext_01, "0.01")
marks_all <- dplyr::bind_rows(marks_05, marks_01)

cat("\n--- edge type for all 16 scm edges, single-run extended pag, both alphas ---\n")
print(as.data.frame(marks_all[, c("alpha", "arrow", "relation")]), row.names = FALSE)

n_bidirected_05 <- sum(marks_all$relation[marks_all$alpha == "0.05"] == "bidirected")
n_bidirected_01 <- sum(marks_all$relation[marks_all$alpha == "0.01"] == "bidirected")
n_circlecircle_05 <- sum(marks_all$relation[marks_all$alpha == "0.05"] == "circle_circle")
n_circlecircle_01 <- sum(marks_all$relation[marks_all$alpha == "0.01"] == "circle_circle")
cat("\nbidirected in this one fit: alpha=.05:", n_bidirected_05, "/ 16   alpha=.01:", n_bidirected_01, "/ 16\n")
cat("circle-circle (fully unresolved) in this one fit: alpha=.05:", n_circlecircle_05,
    "/ 16   alpha=.01:", n_circlecircle_01, "/ 16\n")
cat("(this is the real single-fit answer, not the bootstrap's resample-variability\n",
    "percentages in joint_pag_edgetype_audit.csv -- compare the two directly, they're\n",
    "answering different questions and don't have to agree.)\n", sep = "")

write.csv(marks_all, file.path(OUTPUT_DIR, "tables", "extended_pag_scm_edge_marks.csv"), row.names = FALSE)

# ---- step 3: lv-ida, if it's available -------------------------------------
if (!file.exists(LVIDA_PATH)) {
  cat("\n", LVIDA_PATH, " not found -- stopping here. download lvida.R from\n",
      "https://github.com/dmalinsk/lv-ida, put it at that path (or edit\n",
      "LVIDA_PATH above), and rerun for the mag-count + effect-bounds part.\n", sep = "")
} else {
  source(LVIDA_PATH)

  # correlation matrix, not raw covariance -- agg_ext's variables aren't
  # already standardized, and passing the correlation matrix as lv.ida's
  # mcov is the standard way to get standardized-scale regression
  # coefficients out of it directly (same idea as lavaan's
  # standardizedSolution() used everywhere else in 05's helpers), instead of
  # a raw-covariance-scale number that would need separate rescaling.
  mcov_ext <- cor(agg_ext)

  # x_pos/y_pos below are matched against node_order_ext, but mcov_ext's row/
  # column order comes from agg_ext -- almost certainly the same order
  # already, but worth proving rather than assuming, since a silent mismatch
  # here would misassign every effect to the wrong node without erroring.
  stopifnot(identical(colnames(mcov_ext), node_order_ext))

  NMAGS_CAP <- 500
  lvida_targets <- setdiff(node_order_ext, "climate_behavior")  # all 8 non-outcome
  # nodes -- 7 non-political + politics, exactly the set fig 7 / 07's 8-node
  # ate table uses, so this maps onto the same comparison.
  y_pos <- match("climate_behavior", node_order_ext)

  mag_count_report <- function(fci_fit, alpha_label) {
    am <- fci_fit@amat
    mags <- tryCatch(listMags(am, nMags = NMAGS_CAP, method = "global"),
                      error = function(e) NULL)
    n <- if (is.null(mags)) NA_integer_ else length(mags)
    tibble::tibble(alpha = alpha_label, nMags_cap = NMAGS_CAP,
                   n_mags_found = n, hit_cap = !is.na(n) && n == NMAGS_CAP)
  }
  mag_counts <- dplyr::bind_rows(mag_count_report(fci_ext_05, "0.05"),
                                  mag_count_report(fci_ext_01, "0.01"))
  cat("\n--- listMags() mag-completion count per alpha ---\n")
  print(as.data.frame(mag_counts), row.names = FALSE)
  if (any(mag_counts$hit_cap, na.rm = TRUE)) {
    cat("at least one alpha hit the", NMAGS_CAP, "cap -- that count is a LOWER bound,",
        "not the true number of completions. raise nMags and rerun before treating",
        "either count as exhaustive.\n")
  }

  # doing `effects_std[!is.na(effects_std)]` before saving anything would
  # silently drop the unidentified mags instead of recording them. per
  # lv.ida()'s own code (lv-ida/lvida.R, global method, read directly rather
  # than assumed): it returns exactly one value per mag, and that value is
  # one of three genuinely different things --
  #   0   = x is NOT a possible ancestor of y in that mag -> a real,
  #         identified answer: "no effect under this structure"
  #   NA  = x COULD be an ancestor but the adjustment lv.ida uses can't
  #         identify the effect from that mag -> genuinely unidentified,
  #         not "no effect" and not "no data"
  #   number = identified nonzero effect under that mag
  # collapsing NA into "dropped" would make a 30-mag equivalence class look
  # like a 22-mag one with no explanation for the other 8. this version keeps
  # every mag, tags each with its status, and only computes min/median/max
  # over the mags lv.ida() actually identified (zero + nonzero) -- never over
  # an already-filtered vector that lost the unidentified count along the way.
  run_lvida_report <- function(fci_fit, alpha_label) {
    am <- fci_fit@amat
    summaries <- list()
    long_rows <- list()
    # per-alpha mag total from listMags() itself (mag_counts, computed
    # above) -- lv.ida() should return exactly one entry per compatible mag,
    # so check that explicitly rather than let a length mismatch surface
    # indirectly, later, as a confusing count discrepancy somewhere
    # downstream.
    expected_n_mags <- mag_counts$n_mags_found[mag_counts$alpha == alpha_label]
    for (nd in lvida_targets) {
      x_pos <- match(nd, node_order_ext)
      # fail loudly on a real lv.ida() error instead of silently recording
      # it as "one unidentified mag" -- a length-1 NA fallback would hide a
      # genuine failure behind what looks like ordinary structural
      # uncertainty. didn't happen in the runs so far (every node/alpha
      # completed and its count matched expected_n_mags below), but this is
      # exactly the kind of silent-failure shape this project keeps finding,
      # so don't leave the door open for it here.
      effects_std <- tryCatch(
        lv.ida(x_pos, y_pos, mcov_ext, am, method = "global", nMags = NMAGS_CAP),
        error = function(e) stop("lv.ida failed for ", nd, " at alpha=", alpha_label,
                                  ": ", e$message)
      )
      stopifnot(length(effects_std) == expected_n_mags)
      # rescale to the same "+0.5 sd shift" scale the rest of clean_pipeline
      # reports ates on: dY = 0.5 * tau, for standardized tau, baseline mean
      # 0. conceptually straightforward under the same linearity assumption
      # already used throughout 05's helpers -- still worth checking this
      # equivalence empirically once this actually runs, not just assuming it.
      effects_ate_scale <- 0.5 * effects_std

      status <- dplyr::case_when(
        is.na(effects_std) ~ "unidentified",
        effects_std == 0 ~ "zero",
        TRUE ~ "identified_nonzero"
      )
      n_mags_total <- length(effects_std)
      n_zero <- sum(status == "zero")
      n_identified_nonzero <- sum(status == "identified_nonzero")
      n_unidentified <- sum(status == "unidentified")

      # everything lv.ida() actually identified -- zero AND nonzero, never
      # just nonzero -- is what min/median/max should summarize over.
      identified_vals <- effects_ate_scale[status != "unidentified"]
      n_identified_total <- length(identified_vals)

      summaries[[nd]] <- tibble::tibble(
        alpha = alpha_label, node = nd,
        n_mags_total = n_mags_total,
        n_zero = n_zero,
        n_identified_nonzero = n_identified_nonzero,
        n_unidentified = n_unidentified,
        ate_min = if (n_identified_total) min(identified_vals) else NA_real_,
        ate_q25 = if (n_identified_total) stats::quantile(identified_vals, .25, names = FALSE) else NA_real_,
        ate_median = if (n_identified_total) stats::median(identified_vals) else NA_real_,
        ate_q75 = if (n_identified_total) stats::quantile(identified_vals, .75, names = FALSE) else NA_real_,
        ate_max = if (n_identified_total) max(identified_vals) else NA_real_
      )

      long_rows[[nd]] <- tibble::tibble(
        alpha = alpha_label, node = nd,
        mag_index = seq_along(effects_std),
        effect_std_scale = effects_std,
        effect_ate_scale = effects_ate_scale,
        status = status
      )
    }
    list(summary = dplyr::bind_rows(summaries), long = dplyr::bind_rows(long_rows))
  }

  lvida_05_full <- run_lvida_report(fci_ext_05, "0.05")
  lvida_01_full <- run_lvida_report(fci_ext_01, "0.01")
  lvida_05 <- lvida_05_full$summary
  lvida_01 <- lvida_01_full$summary
  lvida_all <- dplyr::bind_rows(lvida_05, lvida_01)
  lvida_long_all <- dplyr::bind_rows(lvida_05_full$long, lvida_01_full$long)

  cat("\n--- lv-ida possible-effect summary on climate_behavior, rescaled to the +0.5sd ate scale ---\n")
  cat("(dY = 0.5 * tau, tau from lv.ida() run on the correlation matrix -- comparable in\n",
      "principle to scm_edges_finalized.csv / the deterministic ate tables' .206/.195/etc,\n",
      "but that equivalence hasn't been checked empirically yet. n_mags_total splits into\n",
      "n_zero (x not a possible ancestor of y in that mag -- a real identified zero),\n",
      "n_identified_nonzero, and n_unidentified (x could be an ancestor but lv.ida can't\n",
      "identify the effect from that mag). ate_min/median/max are computed over n_zero +\n",
      "n_identified_nonzero only -- i.e. every mag lv.ida() actually identified, never over\n",
      "an already-filtered vector. THIS TABLE IS A SUMMARY -- see the full per-mag status\n",
      "breakdown printed next before drawing any conclusion from ranges alone -- don't\n",
      "collapse this down without seeing the shape, and don't silently drop unidentified\n",
      "mags before summarizing.)\n", sep = "")
  print(as.data.frame(lvida_all), row.names = FALSE)

  cat("\n--- full per-mag lv-ida results (+0.5sd ate scale), ALL mags incl. unidentified, NOT summarized ---\n")
  for (al in c("0.05", "0.01")) {
    cat("\nalpha =", al, ":\n")
    for (nd in lvida_targets) {
      sub <- lvida_long_all[lvida_long_all$alpha == al & lvida_long_all$node == nd, ]
      n_total <- nrow(sub)
      n_zero <- sum(sub$status == "zero")
      n_na <- sum(sub$status == "unidentified")
      nonzero_vals <- sub$effect_ate_scale[sub$status == "identified_nonzero"]
      cat("  ", nd, ": ", n_total, " mags total -- ", n_zero, " zero, ",
          length(nonzero_vals), " identified nonzero, ", n_na, " unidentified(NA)",
          sep = "")
      if (length(nonzero_vals) > 0) {
        cat(" -- nonzero values: ", paste(sprintf("%.4f", sort(nonzero_vals)), collapse = ", "), "\n", sep = "")
      } else {
        cat("\n")
      }
    }
  }

  write.csv(lvida_all, file.path(OUTPUT_DIR, "tables", "extended_pag_lvida_singlenode.csv"), row.names = FALSE)
  write.csv(lvida_long_all, file.path(OUTPUT_DIR, "tables", "extended_pag_lvida_full_vectors.csv"), row.names = FALSE)
  write.csv(mag_counts, file.path(OUTPUT_DIR, "tables", "extended_pag_mag_counts.csv"), row.names = FALSE)
  cat("\nWrote", file.path(OUTPUT_DIR, "tables", "extended_pag_lvida_singlenode.csv"), "(summary, incl.\n",
      "n_zero/n_identified_nonzero/n_unidentified counts), ",
      file.path(OUTPUT_DIR, "tables", "extended_pag_lvida_full_vectors.csv"),
      "(every mag, tagged with status -- zero/identified_nonzero/unidentified, nothing dropped), and\n",
      file.path(OUTPUT_DIR, "tables", "extended_pag_mag_counts.csv"), "\n")
  cat("\n(not attempting pairs/triples through lv-ida yet -- setting two nodes at once can\n",
      "block pathways between them, so single-node effects don't just add. singles first.)\n", sep = "")
}
