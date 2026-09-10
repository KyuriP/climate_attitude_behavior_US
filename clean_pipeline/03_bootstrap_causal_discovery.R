# fci/pc-stable bootstrap, straight from the qmd's causal-setup/bootstrap
# chunks. only real change: everything writes to a dated filename + manifest
# row instead of overwriting one bare .rds, so 04 can't read a stale array
# without knowing it.


source("clean_pipeline/00_config.R")

stopifnot(exists("df_net_main"), exists("df_net_ext"))

suppressPackageStartupMessages({
  library(pcalg); library(graph); library(dplyr)
  library(furrr); library(future)
})

RUN_DATE <- format(Sys.Date(), "%Y%m%d")
tag <- function(name) file.path(OUTPUT_DIR, sprintf("%s_nboot%d_%s.rds", name, N_BOOT, RUN_DATE))

n_cores <- max(1L, parallelly::availableCores() - 1L)

# ---------------------------------------------------------------------------
# Main (8-node) network -- verbatim from `causal-setup` + `bootstrap`
# ---------------------------------------------------------------------------
node_order_cd  <- NODE_ORDER_MAIN
p              <- length(node_order_cd)
agg_data       <- df_net_main[, node_order_cd]
suffStat_gauss <- list(C = cor(agg_data), n = nrow(agg_data))

fci_marks <- c("N", "o", ">", "-")
pc_marks  <- c("N", "-", ">")

run_one <- function(b, data, nms, alpha) {
  p      <- length(nms)
  idx    <- sample(nrow(data), replace = TRUE)
  suff_b <- list(C = cor(as.matrix(data[idx, ])), n = length(idx))
  fci_arr <- array(0L, dim = c(p, p, 4L))
  pc_arr  <- array(0L, dim = c(p, p, 3L))
  tryCatch({
    fci_b <- pcalg::fci(suffStat = suff_b, indepTest = pcalg::gaussCItest,
                         alpha = alpha, labels = nms, verbose = FALSE)
    am <- fci_b@amat
    for (i in seq_len(p)) for (j in seq_len(p)) {
      if (i == j) next
      slot <- am[i, j] + 1L
      fci_arr[i, j, slot] <- fci_arr[i, j, slot] + 1L
    }
  }, error = function(e) NULL)
  tryCatch({
    pc_b <- pcalg::pc(suffStat = suff_b, indepTest = pcalg::gaussCItest,
                       alpha = alpha, labels = nms, skel.method = "stable", verbose = FALSE)
    am <- as(pc_b@graph, "matrix")
    for (i in seq_len(p)) for (j in seq_len(p)) {
      if (i == j) next
      slot <- if      (am[i,j]==0L && am[j,i]==0L) 1L
              else if (am[i,j]==1L && am[j,i]==1L) 2L
              else if (am[i,j]==1L && am[j,i]==0L) 3L
              else if (am[i,j]==0L && am[j,i]==1L) 2L
              else                                  1L
      pc_arr[i, j, slot] <- pc_arr[i, j, slot] + 1L
    }
  }, error = function(e) NULL)
  list(fci = fci_arr, pc = pc_arr)
}

plan(multisession, workers = n_cores)
set.seed(SEED)
results <- list()
for (alph in names(ALPHAS)) {
  cat("Main-network bootstrap alpha =", alph, "\n")
  results[[alph]] <- furrr::future_map(
    seq_len(N_BOOT), run_one, data = agg_data, nms = node_order_cd, alpha = ALPHAS[[alph]],
    .options = furrr::furrr_options(seed = TRUE), .progress = TRUE)
}
plan(sequential)

mark_props <- list()
for (alph in names(ALPHAS)) {
  fci_c <- array(0L, dim = c(p, p, 4L), dimnames = list(node_order_cd, node_order_cd, fci_marks))
  pc_c  <- array(0L, dim = c(p, p, 3L), dimnames = list(node_order_cd, node_order_cd, pc_marks))
  for (res in results[[alph]]) {
    if (is.null(res)) next
    fci_c <- fci_c + res$fci; pc_c <- pc_c + res$pc
  }
  fci_sums <- apply(fci_c, c(1,2), sum); diag(fci_sums) <- N_BOOT
  pc_sums  <- apply(pc_c,  c(1,2), sum); diag(pc_sums)  <- N_BOOT
  if (any(abs(fci_sums - N_BOOT) > 1)) warning(sprintf("Alpha %s: FCI cell sums deviate from %d", alph, N_BOOT))
  if (any(abs(pc_sums  - N_BOOT) > 1)) warning(sprintf("Alpha %s: PC cell sums deviate from %d",  alph, N_BOOT))
  mark_props[[alph]] <- list(fci = fci_c / N_BOOT, pc = pc_c / N_BOOT, fci_counts = fci_c, pc_counts = pc_c)
}

fci_props_main <- (mark_props[["0.05"]]$fci_counts + mark_props[["0.01"]]$fci_counts) / (2 * N_BOOT)
pc_props_main  <- (mark_props[["0.05"]]$pc_counts  + mark_props[["0.01"]]$pc_counts)  / (2 * N_BOOT)
dimnames(fci_props_main) <- list(node_order_cd, node_order_cd, fci_marks)
dimnames(pc_props_main)  <- list(node_order_cd, node_order_cd, pc_marks)

saveRDS(fci_props_main, tag("fci_props_main"))
saveRDS(pc_props_main,  tag("pc_props_main"))
cat("Main-network bootstrap complete:", N_BOOT, "resamples per alpha.\n")

# ---------------------------------------------------------------------------
# Extended (9-node, +climate_behavior) network -- verbatim from `ext-setup`
# + `ext-bootstrap`
# ---------------------------------------------------------------------------
node_order_ext <- NODE_ORDER_EXT
p_ext          <- length(node_order_ext)
agg_ext        <- df_net_ext[, node_order_ext]
context_idx    <- which(node_order_ext %in% node_order_cd)
abbr_ext       <- ABBR_EXT

fci_marks_ext <- c("N","o",">","-")
pc_marks_ext  <- c("N","-",">")

run_one_ext <- function(b, data, nms, alpha, ctx_idx) {
  p      <- length(nms)
  idx    <- sample(nrow(data), replace = TRUE)
  suff_b <- list(C = cor(as.matrix(data[idx, ])), n = length(idx))
  fci_arr <- array(0L, dim = c(p, p, 4L))
  pc_arr  <- array(0L, dim = c(p, p, 3L))
  tryCatch({
    fci_b <- pcalg::fci(suffStat = suff_b, indepTest = pcalg::gaussCItest,
                         alpha = alpha, labels = nms, contextVars = ctx_idx, jci = "1",
                         selectionBias = FALSE, verbose = FALSE)
    am <- fci_b@amat
    for (i in seq_len(p)) for (j in seq_len(p)) {
      if (i == j) next
      slot <- am[i, j] + 1L
      fci_arr[i, j, slot] <- fci_arr[i, j, slot] + 1L
    }
  }, error = function(e) NULL)
  tryCatch({
    pc_b <- pcalg::pc(suffStat = suff_b, indepTest = pcalg::gaussCItest,
                       alpha = alpha, labels = nms, skel.method = "stable", verbose = FALSE)
    am <- as(pc_b@graph, "matrix")
    att_nms <- nms[nms != "climate_behavior"]
    bg_pairs <- Filter(function(a) am[a,"climate_behavior"]==1 && am["climate_behavior",a]==1, att_nms)
    if (length(bg_pairs) > 0) {
      bg_result <- pcalg::addBgKnowledge(gInput = pc_b@graph, x = bg_pairs,
                                          y = rep("climate_behavior", length(bg_pairs)))
      am <- as(bg_result, "matrix")
    }
    for (i in seq_len(p)) for (j in seq_len(p)) {
      if (i == j) next
      slot <- if      (am[i,j]==0L && am[j,i]==0L) 1L
              else if (am[i,j]==1L && am[j,i]==1L) 2L
              else if (am[i,j]==1L && am[j,i]==0L) 3L
              else if (am[i,j]==0L && am[j,i]==1L) 2L
              else                                  1L
      pc_arr[i, j, slot] <- pc_arr[i, j, slot] + 1L
    }
  }, error = function(e) NULL)
  list(fci = fci_arr, pc = pc_arr)
}

plan(multisession, workers = n_cores)
set.seed(SEED)
results_ext <- list()
for (alph in names(ALPHAS)) {
  cat("Extended-network bootstrap alpha =", alph, "\n")
  results_ext[[alph]] <- furrr::future_map(
    seq_len(N_BOOT), run_one_ext, data = agg_ext, nms = node_order_ext, alpha = ALPHAS[[alph]],
    ctx_idx = context_idx, .options = furrr::furrr_options(seed = TRUE), .progress = TRUE)
}
plan(sequential)

mark_props_ext <- list()
for (alph in names(ALPHAS)) {
  fci_c <- array(0L, dim = c(p_ext, p_ext, 4L), dimnames = list(node_order_ext, node_order_ext, fci_marks_ext))
  pc_c  <- array(0L, dim = c(p_ext, p_ext, 3L), dimnames = list(node_order_ext, node_order_ext, pc_marks_ext))
  for (res in results_ext[[alph]]) {
    if (is.null(res)) next
    fci_c <- fci_c + res$fci; pc_c <- pc_c + res$pc
  }
  mark_props_ext[[alph]] <- list(fci = fci_c / N_BOOT, pc = pc_c / N_BOOT, fci_counts = fci_c, pc_counts = pc_c)
}

fci_props_ext <- (mark_props_ext[["0.05"]]$fci_counts + mark_props_ext[["0.01"]]$fci_counts) / (2 * N_BOOT)
pc_props_ext  <- (mark_props_ext[["0.05"]]$pc_counts  + mark_props_ext[["0.01"]]$pc_counts)  / (2 * N_BOOT)
dimnames(fci_props_ext) <- list(node_order_ext, node_order_ext, fci_marks_ext)
dimnames(pc_props_ext)  <- list(node_order_ext, node_order_ext, pc_marks_ext)

fci_ext_path  <- tag("fci_props_ext")
pc_ext_path   <- tag("pc_props_ext")
mark_ext_path <- tag("mark_props_ext_by_alpha")
saveRDS(fci_props_ext,  fci_ext_path)
saveRDS(pc_props_ext,   pc_ext_path)
saveRDS(mark_props_ext, mark_ext_path)
cat("Extended-network bootstrap complete:", N_BOOT, "resamples per alpha.\n")

# ---------------------------------------------------------------------------
# Provenance: append to the manifest, then overwrite the "latest" pointer.
# This is the actual fix for the bug that started this refactor -- every
# array this script writes is dated and versioned, and 04_scm_finalize.R
# reads ONLY through LATEST_bootstrap_run.csv, never a bare hardcoded
# filename, so a superseded run can never again be read silently.
# ---------------------------------------------------------------------------
manifest_row <- data.frame(
  timestamp = as.character(Sys.time()), n_boot = N_BOOT,
  alphas = paste(names(ALPHAS), collapse = ";"),
  fci_props_ext_path = fci_ext_path, pc_props_ext_path = pc_ext_path,
  mark_props_ext_path = mark_ext_path,
  node_order_ext = paste(node_order_ext, collapse = ";")
)
write.table(manifest_row, MANIFEST_PATH, sep = ",", row.names = FALSE,
            col.names = !file.exists(MANIFEST_PATH), append = file.exists(MANIFEST_PATH))
write.csv(manifest_row, LATEST_POINTER_PATH, row.names = FALSE)
cat("\nManifest updated:", MANIFEST_PATH, "\nLatest pointer:", LATEST_POINTER_PATH, "\n")
cat("*** fci_props_ext for THIS run: n_boot =", N_BOOT, "-- generated", RUN_DATE, "***\n")

# single-run (non-bootstrap) FCI/PC-stable fits -- fig 2 panel B and supp fig 8
# need one fit to the real sample at each alpha, not a bootstrap resample.
# nothing else in 00-04 computes these. from the qmd's fci-05/fci-01/pc-05/
# pc-01 chunks, main network only (these are illustrative, not a quantitative
# claim).

fci_05 <- pcalg::fci(suffStat = suffStat_gauss, indepTest = pcalg::gaussCItest,
                      alpha = .05, labels = node_order_cd, verbose = FALSE)
fci_01 <- pcalg::fci(suffStat = suffStat_gauss, indepTest = pcalg::gaussCItest,
                      alpha = .01, labels = node_order_cd, verbose = FALSE)
pc_05  <- pcalg::pc(suffStat = suffStat_gauss, indepTest = pcalg::gaussCItest,
                     alpha = .05, labels = node_order_cd, skel.method = "stable", verbose = FALSE)
pc_01  <- pcalg::pc(suffStat = suffStat_gauss, indepTest = pcalg::gaussCItest,
                     alpha = .01, labels = node_order_cd, skel.method = "stable", verbose = FALSE)

cat("Single-run FCI/PC-stable fits (fci_05, fci_01, pc_05, pc_01) computed",
    "-- needed for Figure 2B and Supplementary Figure 8, not for any bootstrap table.\n")
