# Bootstrap PAG -> MAG -> LV-IDA sensitivity
# Purpose: propagate sampling variability in the FCI PAG into PAG-compatible
# single-node total-effect conclusions. Each bootstrap resample is one unit of
# the outer bootstrap summary. MAGs are NOT pooled across resamples as if they
# were probability draws.
source("clean_pipeline/00_config.R")
suppressPackageStartupMessages({
  library(pcalg)
  library(dplyr)
  library(tibble)
  library(furrr)
  library(future)
})
# Requires the current extended analysis objects to already exist in-session.
# These are created in clean_pipeline/03_bootstrap_causal_discovery.R after
# df_net_ext exists.
stopifnot(exists("agg_ext"), exists("node_order_ext"), exists("context_idx"))
stopifnot(identical(node_order_ext, NODE_ORDER_EXT))
stopifnot(identical(colnames(agg_ext), node_order_ext))
stopifnot(nrow(agg_ext) == 870L)
stopifnot(!anyNA(agg_ext))
stopifnot(all(vapply(agg_ext, function(x) is.finite(stats::sd(x)) && stats::sd(x) > 0,
                     logical(1))))
LVIDA_PATH <- "lv-ida/lvida.R"
if (!file.exists(LVIDA_PATH)) {
  stop(LVIDA_PATH,
       " not found. Download the current lvida.R from https://github.com/dmalinsk/lv-ida ",
       "and put it at that path (or edit LVIDA_PATH).")
}
source(LVIDA_PATH)
# ------------------------- RUN SETTINGS ------------------------------------
# Pilot (20L) confirmed clean: n_cap_hit = 0, n_fci_failed = 0,
# n_mag_failed = 0, n_invalid_correlation = 0 at both alphas. The only
# failures were n_effect_computation_failed (1/20 at alpha=.01, 2/20 at
# alpha=.05) -- a macOS main-thread C-stack limit hit inside the recursive
# graph-traversal step (dsepset.reach()/is.ancestor()), confirmed via
# Cstack_info() to be an OS-level limitation that ulimit -s cannot raise on
# this machine (unlike Linux, macOS does not apply RLIMIT_STACK to the main
# thread). These failures are caught and labeled cleanly (no wrong values --
# the computation just aborts), and drop out for all 8 target nodes at once,
# so they cannot selectively bias which effects look positive vs. null.
# Moving to the manuscript run at 1000L; report the excluded fraction
# (expect roughly the same ~7.5% rate as the pilot) the same way
# n_cap_hit is already reported.
N_BOOT_LVIDA <- 1000L
# listMags() truncates when nMags is reached. A resample returning exactly this
# many MAGs is therefore treated as cap-ambiguous and EXCLUDED, not analyzed
# as though the enumeration were exhaustive. Raise this if the pilot hits it.
# Raised 2026-09-10: the 1000L manuscript run hit the cap for exactly 2
# resamples (boot_id 666, 934, both alpha=0.01; batches 4 and 5) out of the
# other 998 resamples' real max of 372 compatible MAGs -- so 500 was cutting
# it close but not for most of the distribution. Raising to 2000L for real
# headroom above that observed max. If these two (or others) still hit even
# 2000, that's a sign the true count may be unboundedly large for that
# specific resampled PAG structure -- at that point treat a small
# residual n_cap_hit as a legitimate, reported exclusion (like
# effect_computation_failed) rather than keep chasing an ever-higher cap.
NMAGS_CAP <- 2000L
# A single pathological resample (a PAG whose MAG enumeration or
# graph-traversal step blows up combinatorially) can otherwise block an
# entire batch indefinitely -- observed 2026-09-10: one worker sat at ~99%
# CPU far longer than the other 199 resamples in its batch combined. This is
# a judgment call, not measured from real per-resample timing data: normal
# resamples in the pilot finished in well under a minute each, so 300s is a
# generous margin for genuine cases while still bounding the worst case.
# Tighten or loosen based on what normal resamples actually take on your
# machine once you have a feel for it.
PER_BOOT_TIMEOUT_SECS <- 300L
lvida_targets <- setdiff(node_order_ext, "climate_behavior")
y_pos <- match("climate_behavior", node_order_ext)
# Suppress the very verbose unconditional cat() output inside listMags().
list_mags_quiet <- function(am, n_mags_cap) {
  out <- NULL
  invisible(utils::capture.output(
    out <- listMags(am, nMags = n_mags_cap, method = "global")
  ))
  out
}
# Exact global LV-IDA calculation for a PRECOMPUTED MAG list.
# This mirrors the global branch of Daniel Malinsky's lv.ida() implementation,
# but reuses the same MAG list for all eight exposure nodes so MAG enumeration
# is done ONCE per bootstrap PAG rather than once per target.
lvida_from_mags <- function(x_pos, y_pos, mcov, mags) {
  beta_hat <- rep(NA_real_, length(mags))
  for (i in seq_along(mags)) {
    gMag <- mags[[i]]
    # If X is not an ancestor of Y in this MAG, total effect is identified zero.
    if (!is.ancestor(x_pos, y_pos, gMag)) {
      beta_hat[i] <- 0
      next
    }
    gMag_x <- remove.visible.edges(x_pos, gMag)
    dsepset <- setdiff(dsepset.reach(x_pos, y_pos, -1, gMag_x), x_pos)
    descendants_x <- which(vapply(seq_len(nrow(gMag)), function(k) {
      is.ancestor(x_pos, k, gMag)
    }, logical(1)))
    if (gMag_x[x_pos, y_pos] != 0) {
      beta_hat[i] <- NA_real_
    } else if (length(intersect(dsepset, descendants_x)) != 0) {
      beta_hat[i] <- NA_real_
    } else {
      beta_hat[i] <- lm.cov(mcov, y_pos, c(x_pos, dsepset))
    }
  }
  beta_hat
}
# ------------------- ONE-TIME IMPLEMENTATION CHECK -------------------------
# Before bootstrapping, prove that the precomputed-MAG helper reproduces the
# official lv.ida() output on the full sample at both alpha thresholds.
check_precomputed_helper <- function(alpha_value) {
  suff <- list(C = cor(as.matrix(agg_ext)), n = nrow(agg_ext))
  fit <- pcalg::fci(
    suffStat = suff, indepTest = pcalg::gaussCItest,
    alpha = alpha_value, labels = node_order_ext,
    contextVars = context_idx, jci = "1",
    selectionBias = FALSE, verbose = FALSE
  )
  mags <- list_mags_quiet(fit@amat, NMAGS_CAP)
  if (length(mags) == NMAGS_CAP) {
    stop("Full-sample PAG hit NMAGS_CAP=", NMAGS_CAP,
         " during helper validation. Raise NMAGS_CAP before continuing.")
  }
  mcov <- suff$C
  for (nd in lvida_targets) {
    x_pos <- match(nd, node_order_ext)
    custom <- lvida_from_mags(x_pos, y_pos, mcov, mags)
    official <- NULL
    invisible(utils::capture.output(
      official <- lv.ida(x_pos, y_pos, mcov, fit@amat,
                         method = "global", nMags = NMAGS_CAP)
    ))
    # lv.ida()'s own beta.hat starts life as rep(NA, n.mags) -- plain logical
    # NA, not NA_real_ -- and if every MAG for this node falls into an NA
    # branch (never hits beta.hat[i] <- 0 or a numeric lm.cov() result), the
    # whole vector never gets coerced to double. lvida_from_mags() always
    # returns double (rep(NA_real_, ...)). all.equal() treats a logical-NA
    # vector and a double-NA_real_ vector as a mode mismatch even when every
    # value genuinely agrees -- confirmed directly (2026-09-09, belief_concern
    # at alpha=.05, 1 MAG): both sides print "NA", typeof(custom)="double",
    # typeof(official)="logical", all.equal() on the raw values is FALSE, but
    # coercing both to numeric first makes it TRUE. as.numeric() below fixes
    # the comparison only -- lvida_from_mags() itself was never wrong, this
    # never affects the actual bootstrap effects computed later in this
    # script, which only ever reads from lvida_from_mags()'s (always-double)
    # output, never from lv.ida() directly.
    if (!isTRUE(all.equal(as.numeric(custom), as.numeric(official), tolerance = 1e-10,
                          check.attributes = FALSE))) {
      stop("Precomputed-MAG helper does not reproduce lv.ida() for ", nd,
           " at alpha=", alpha_value)
    }
  }
  cat("Helper check passed at alpha =", alpha_value,
      "with", length(mags), "compatible MAG(s).\n")
  invisible(TRUE)
}
check_precomputed_helper(.05)
check_precomputed_helper(.01)
# ------------------------ ONE BOOTSTRAP RESAMPLE ----------------------------
run_one_boot_lvida <- function(b, data, nms, alpha, ctx_idx,
                               y_pos, targets, nmags_cap) {
  # Bound this resample's total wall-clock time. If exceeded, R raises an
  # error wherever execution currently is -- inside the FCI-fit tryCatch,
  # the MAG-enumeration tryCatch, or the effect-computation tryCatch below
  # -- so it gets caught and labeled exactly like any other failure at that
  # stage, with "reached elapsed time limit" in the message column
  # (filterable afterward to tell a genuine timeout apart from a real
  # failure of that stage). transient = TRUE plus the on.exit reset keeps
  # this scoped to this one call and this one worker process.
  setTimeLimit(cpu = PER_BOOT_TIMEOUT_SECS, elapsed = PER_BOOT_TIMEOUT_SECS,
               transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  idx <- sample(nrow(data), replace = TRUE)
  data_b <- as.matrix(data[idx, , drop = FALSE])
  C_b <- cor(data_b)
  if (any(!is.finite(C_b))) {
    return(list(
      meta = tibble::tibble(boot_id = b, run_status = "invalid_correlation",
                            n_mags = NA_integer_),
      effects = NULL
    ))
  }
  suff_b <- list(C = C_b, n = nrow(data_b))
  fci_result <- tryCatch(
    list(
      fit = pcalg::fci(
        suffStat = suff_b, indepTest = pcalg::gaussCItest,
        alpha = alpha, labels = nms,
        contextVars = ctx_idx, jci = "1",
        selectionBias = FALSE, verbose = FALSE
      ),
      error = NULL
    ),
    error = function(e) list(fit = NULL, error = conditionMessage(e))
  )
  if (is.null(fci_result$fit)) {
    return(list(
      meta = tibble::tibble(boot_id = b, run_status = "fci_failed",
                            n_mags = NA_integer_, message = fci_result$error),
      effects = NULL
    ))
  }
  am <- fci_result$fit@amat
  mags_result <- tryCatch(
    list(mags = list_mags_quiet(am, nmags_cap), error = NULL),
    error = function(e) list(mags = NULL, error = conditionMessage(e))
  )
  if (is.null(mags_result$mags)) {
    return(list(
      meta = tibble::tibble(boot_id = b, run_status = "mag_enumeration_failed",
                            n_mags = NA_integer_, message = mags_result$error),
      effects = NULL
    ))
  }
  mags <- mags_result$mags
  n_mags_b <- length(mags)
  # Equality to the cap is ambiguous: the true class may contain exactly the
  # cap or may have been truncated. Do not use it as exhaustive.
  if (n_mags_b == nmags_cap) {
    return(list(
      meta = tibble::tibble(boot_id = b, run_status = "mag_cap_hit",
                            n_mags = n_mags_b,
                            message = paste0("Reached NMAGS_CAP=", nmags_cap)),
      effects = NULL
    ))
  }
  mcov_b <- C_b
  dimnames(mcov_b) <- list(nms, nms)
  # Everything above this point (resampling, FCI, MAG enumeration) already has
  # its own tryCatch with a labeled failure status. This is the one remaining
  # uncaught step -- if lvida_from_mags() throws for some rare resampled PAG
  # structure, an uncaught error here is why a boot_id could vanish from
  # run_status entirely instead of showing up under any known failure count.
  # Give it the same treatment: catch it, label it, keep the message.
  effect_result <- tryCatch(
    {
      effect_rows <- lapply(targets, function(nd) {
        x_pos <- match(nd, nms)
        effects_std <- lvida_from_mags(x_pos, y_pos, mcov_b, mags)
        stopifnot(length(effects_std) == n_mags_b)
        effects_ate_scale <- 0.5 * effects_std
        status <- dplyr::case_when(
          is.na(effects_std) ~ "unidentified",
          effects_std == 0   ~ "zero",
          effects_std > 0    ~ "identified_positive",
          effects_std < 0    ~ "identified_negative"
        )
        tibble::tibble(
          boot_id = b,
          node = nd,
          mag_index = seq_along(effects_std),
          effect_std_scale = effects_std,
          effect_ate_scale = effects_ate_scale,
          status = status
        )
      })
      list(rows = dplyr::bind_rows(effect_rows), error = NULL)
    },
    error = function(e) list(rows = NULL, error = conditionMessage(e))
  )
  if (is.null(effect_result$rows)) {
    return(list(
      meta = tibble::tibble(boot_id = b, run_status = "effect_computation_failed",
                            n_mags = n_mags_b, message = effect_result$error),
      effects = NULL
    ))
  }
  list(
    meta = tibble::tibble(boot_id = b, run_status = "ok",
                          n_mags = n_mags_b, message = NA_character_),
    effects = effect_result$rows
  )
}
# ------------------------------ RUN -----------------------------------------
# Batched to bound peak memory. Running all N_BOOT_LVIDA resamples in one
# future_map() call keeps that many worker processes' accumulated memory
# alive for the whole call -- this is what exhausted memory on the first
# attempt at 1000L. Splitting into batches of BATCH_SIZE, tearing the
# multisession plan down and rebuilding it between batches (releasing each
# worker's memory), and checkpointing every batch to disk immediately bounds
# peak memory to roughly one batch's worth of work. It also means a crash
# partway through loses at most the in-progress batch, not the whole run --
# rerunning this script skips any batch whose file already exists on disk
# rather than recomputing it.
BATCH_SIZE <- 200L
stopifnot(N_BOOT_LVIDA %% BATCH_SIZE == 0)
n_batches <- N_BOOT_LVIDA %/% BATCH_SIZE
# Cap worker count directly rather than relying solely on
# availableCores() - 1, which can ask for more concurrent FCI/MAG-enumeration
# workers than this machine's free RAM supports. Lower this further (e.g. 2L)
# if a batch still runs out of memory; raise it if batches finish comfortably
# and you want more speed.
N_WORKERS_CAP <- 4L
n_cores <- min(N_WORKERS_CAP, max(1L, parallelly::availableCores() - 1L))
cat("Using", n_cores, "parallel worker(s),", n_batches, "batch(es) of",
    BATCH_SIZE, "resamples each, per alpha.\n")
batch_dir <- file.path(OUTPUT_DIR, "tables", "lvida_batches")
dir.create(batch_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)
for (alph in names(ALPHAS)) {
  for (batch in seq_len(n_batches)) {
    batch_file <- file.path(batch_dir,
      paste0("batch_alpha_", alph, "_", sprintf("%03d", batch), ".rds"))
    if (file.exists(batch_file)) {
      cat("Alpha =", alph, "batch", batch, "of", n_batches,
          "already on disk -- skipping.\n")
      next
    }
    boot_ids <- ((batch - 1L) * BATCH_SIZE + 1L):(batch * BATCH_SIZE)
    cat("Bootstrap PAG/LV-IDA, alpha =", alph, "-- batch", batch, "of",
        n_batches, "(boot_id", min(boot_ids), "-", max(boot_ids), ")\n")
    future::plan(multisession, workers = n_cores)
    res <- furrr::future_map(
      boot_ids,
      run_one_boot_lvida,
      data = agg_ext,
      nms = node_order_ext,
      alpha = ALPHAS[[alph]],
      ctx_idx = context_idx,
      y_pos = y_pos,
      targets = lvida_targets,
      nmags_cap = NMAGS_CAP,
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
    future::plan(sequential)
    meta <- dplyr::bind_rows(lapply(res, `[[`, "meta")) |>
      dplyr::mutate(alpha = alph, .before = 1)
    # Same loud check as before, scoped to this batch's boot_ids.
    if (nrow(meta) != length(boot_ids)) {
      present_ids <- sort(unique(meta$boot_id))
      missing_ids <- setdiff(boot_ids, present_ids)
      stop("Alpha = ", alph, ", batch ", batch, ": expected ", length(boot_ids),
           " meta rows, got ", nrow(meta), ". Missing boot_id(s): ",
           paste(missing_ids, collapse = ", "),
           ". run_one_boot_lvida() did not return a result for some b -- ",
           "investigate before trusting any bootstrap summary.")
    }
    effects <- dplyr::bind_rows(lapply(res, `[[`, "effects"))
    if (nrow(effects) > 0) effects$alpha <- alph
    saveRDS(list(meta = meta, effects = effects), batch_file)
    cat("Run status, alpha =", alph, "batch", batch, ":\n")
    print(table(meta$run_status, useNA = "ifany"))
    if (any(meta$run_status == "ok")) {
      ok_mags <- meta$n_mags[meta$run_status == "ok"]
      cat("MAG counts among successful resamples this batch: median =",
          stats::median(ok_mags), "max =", max(ok_mags), "\n")
    }
    rm(res, meta, effects)
    gc()
  }
}
# ------------------------- COMBINE BATCHES -----------------------------------
batch_files <- list.files(batch_dir, pattern = "^batch_alpha_.*\\.rds$",
                          full.names = TRUE)
stopifnot(length(batch_files) == length(names(ALPHAS)) * n_batches)
all_batches <- lapply(batch_files, readRDS)
run_status <- dplyr::bind_rows(lapply(all_batches, `[[`, "meta"))
boot_lvida_full <- dplyr::bind_rows(lapply(all_batches, `[[`, "effects")) |>
  dplyr::select(alpha, boot_id, node, mag_index,
                effect_std_scale, effect_ate_scale, status)
stopifnot(nrow(boot_lvida_full) > 0)
for (alph in names(ALPHAS)) {
  n_here <- sum(run_status$alpha == alph)
  if (n_here != N_BOOT_LVIDA) {
    stop("Alpha = ", alph, ": combined run_status has ", n_here,
         " rows across all batches, expected ", N_BOOT_LVIDA,
         ". Check the .rds files in ", batch_dir, " before trusting results.")
  }
}
# ----------------------- RESAMPLE-LEVEL SUMMARY -----------------------------
# One bootstrap resample is the outer unit. MAGs are summarized WITHIN each
# PAG first. We do not pool all MAGs across resamples into one pseudo-sample.
resample_summary <- boot_lvida_full |>
  dplyr::group_by(alpha, boot_id, node) |>
  dplyr::summarise(
    n_mags_b = dplyr::n(),
    n_zero = sum(status == "zero"),
    n_positive = sum(status == "identified_positive"),
    n_negative = sum(status == "identified_negative"),
    n_unidentified = sum(status == "unidentified"),
    any_positive = n_positive > 0,
    any_negative = n_negative > 0,
    all_zero = n_zero == n_mags_b,
    all_positive = n_positive == n_mags_b,
    any_unidentified = n_unidentified > 0,
    min_positive = if (n_positive > 0)
      min(effect_ate_scale[status == "identified_positive"]) else NA_real_,
    max_positive = if (n_positive > 0)
      max(effect_ate_scale[status == "identified_positive"]) else NA_real_,
    .groups = "drop"
  )
run_counts <- run_status |>
  dplyr::group_by(alpha) |>
  dplyr::summarise(
    n_attempted = dplyr::n(),
    n_success = sum(run_status == "ok"),
    n_fci_failed = sum(run_status == "fci_failed"),
    n_mag_failed = sum(run_status == "mag_enumeration_failed"),
    n_cap_hit = sum(run_status == "mag_cap_hit"),
    n_invalid_correlation = sum(run_status == "invalid_correlation"),
    n_effect_computation_failed = sum(run_status == "effect_computation_failed"),
    .groups = "drop"
  )
node_summary <- resample_summary |>
  dplyr::group_by(alpha, node) |>
  dplyr::summarise(
    n_resamples_with_data = dplyr::n(),
    pct_resamples_any_positive = 100 * mean(any_positive),
    pct_resamples_all_positive = 100 * mean(all_positive),
    pct_resamples_all_zero = 100 * mean(all_zero),
    pct_resamples_any_unidentified = 100 * mean(any_unidentified),
    pct_resamples_any_negative = 100 * mean(any_negative),
    median_n_mags = stats::median(n_mags_b),
    max_n_mags = max(n_mags_b),
    median_resample_min_positive = if (any(any_positive))
      stats::median(min_positive, na.rm = TRUE) else NA_real_,
    median_resample_max_positive = if (any(any_positive))
      stats::median(max_positive, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) |>
  dplyr::left_join(run_counts, by = "alpha")
cat("\n--- BOOTSTRAP PAG/LV-IDA NODE SUMMARY ---\n")
print(as.data.frame(node_summary), row.names = FALSE)
cat("\nInterpretation:\n",
    "- pct_resamples_any_positive: successful bootstrap PAGs in which at least one compatible MAG gave an identified positive effect.\n",
    "- pct_resamples_all_positive: successful bootstrap PAGs in which EVERY compatible MAG gave an identified positive effect (strong criterion).\n",
    "- pct_resamples_all_zero: successful bootstrap PAGs in which every compatible MAG implied zero total effect.\n",
    "- pct_resamples_any_unidentified: successful bootstrap PAGs containing at least one MAG for which LV-IDA could not identify the effect.\n",
    "These are bootstrap stability summaries, not posterior probabilities and not formal confidence intervals.\n",
    sep = "")
# ------------------------------ SAVE ----------------------------------------
dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
write.csv(run_status,
          file.path(OUTPUT_DIR, "tables", "bootstrap_lvida_run_status.csv"),
          row.names = FALSE)
write.csv(boot_lvida_full,
          file.path(OUTPUT_DIR, "tables", "bootstrap_lvida_full_vectors.csv"),
          row.names = FALSE)
write.csv(resample_summary,
          file.path(OUTPUT_DIR, "tables", "bootstrap_lvida_resample_summary.csv"),
          row.names = FALSE)
write.csv(node_summary,
          file.path(OUTPUT_DIR, "tables", "bootstrap_lvida_node_summary.csv"),
          row.names = FALSE)
cat("\nWrote bootstrap LV-IDA outputs to ", file.path(OUTPUT_DIR, "tables"), ".\n", sep = "")
cat("IMPORTANT: for the manuscript run, require n_cap_hit = 0. If the pilot hits the MAG cap, raise NMAGS_CAP and rerun.\n")
