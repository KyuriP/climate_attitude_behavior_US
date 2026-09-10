# =============================================================================
# 11c: two follow-ups from 11b's output
#
# (1) pc_amat_05 is NOT in the same 0/1/2/3 mark convention as fci_amat_05 /
#     cci_amat_05 -- it's pcalg's plain adjacency matrix (1/0 asymmetric =
#     directed, 1/1 symmetric = undirected). Running it through
#     summarize_marks()/edge_type_counts() from 11b mislabels every PC edge
#     as "circle_involved", which is an artifact of the wrong encoding, not
#     a real finding about PC's uncertainty. This block gives PC's real
#     edge-type counts in its own native convention.
#
# (2) CCI came back bit-for-bit identical at alpha=.05 and alpha=.01 -- same
#     17 edges, 100% arrowhead-arrowhead, no circles or tails at either
#     threshold. That's worth a sanity check before trusting it: does
#     fit_cci() ever produce a non-arrowhead mark at all, even on a trivial,
#     confounder-free, acyclic synthetic chain it should recover easily?
#     If yes -> the all-bidirected result on the real data is a genuine,
#     citable finding about this dataset. If the synthetic case ALSO comes
#     back all-bidirected -> something in the fit_cci() call itself needs
#     fixing before CCI results go anywhere near the paper.
#
# Run this after 11b (needs fit_cci(), summarize_marks(), pc_amat_05,
# cci_amat_05, cci_amat_01 already in memory).
# =============================================================================

# ---- (1) PC's real edge types, in its own native convention ---------------
pc_edge_types <- function(amat) {
  p <- nrow(amat)
  directed <- 0; undirected <- 0
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      a <- amat[i, j]; b <- amat[j, i]
      if (a == 0 && b == 0) next
      if (a == 1 && b == 1) undirected <- undirected + 1
      else directed <- directed + 1
    }
  }
  c(directed = directed, undirected = undirected, total = directed + undirected)
}
cat("=== PC-stable, corrected (native 0/1 convention, NOT PAG marks) ===\n")
print(pc_edge_types(pc_amat_05))

# ---- (2) Does CCI ever produce a tail/circle on an easy synthetic case? ---
set.seed(1)
n_sim <- 3000
x1 <- rnorm(n_sim)
x2 <- 0.6 * x1 + rnorm(n_sim)
x3 <- 0.6 * x2 + rnorm(n_sim)
x4 <- 0.6 * x3 + rnorm(n_sim)
sim_dat    <- data.frame(x1, x2, x3, x4)
sim_labels <- names(sim_dat)
sim_suffStat <- list(C = cor(sim_dat), n = n_sim)

sim_amat <- fit_cci(sim_suffStat, alpha = 0.05, labels = sim_labels, p = 4)
cat("\n=== CCI on a clean synthetic chain (x1 -> x2 -> x3 -> x4, no latents, no cycles) ===\n")
print(sim_amat)
cat("\n--- per-edge marks ---\n")
print(summarize_marks(sim_amat))
