# =============================================================================
# 11b: CCI mark-composition diagnostic
#
# Checks whether cci_amat_05's "every edge is arrowhead-arrowhead" pattern
# is alpha-specific, and puts FCI/PC/CCI edges on the SAME 8-node main
# network so the comparison is apples-to-apples (the earlier snippet
# mixed the 9-node extended fci_amat_ext against the 8-node main-network
# pc_amat_05/cci_amat_05 -- that mismatch is why the "FCI: 25" count isn't
# comparable to "PC: 9" / "CCI: 26").
#
# Run this AFTER r_patches/09_supp_causal_graphs_redesign.R and
# r_patches/10_cci_addition_and_figure4_pc_appendix.R, so that
# fci_amat_05, pc_amat_05, cci_amat_05, cci_amat_01 all exist.
# =============================================================================

mark_name <- function(v) c("none", "circle", "arrowhead", "tail")[v + 1]

summarize_marks <- function(amat) {
  p <- nrow(amat)
  rows <- list()
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      mij <- amat[i, j]; mji <- amat[j, i]
      if (mij == 0 && mji == 0) next
      rows[[length(rows) + 1]] <- data.frame(
        from = rownames(amat)[i], to = rownames(amat)[j],
        mark_at_to = mark_name(mij), mark_at_from = mark_name(mji)
      )
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

edge_type_counts <- function(amat) {
  s <- summarize_marks(amat)
  if (nrow(s) == 0) return(c(total = 0))
  type <- ifelse(s$mark_at_to == "arrowhead" & s$mark_at_from == "arrowhead", "bidirected",
           ifelse(s$mark_at_to == "circle" | s$mark_at_from == "circle", "circle_involved",
           ifelse((s$mark_at_to == "arrowhead" & s$mark_at_from == "tail") |
                  (s$mark_at_to == "tail" & s$mark_at_from == "arrowhead"), "directed", "other")))
  c(table(type), total = nrow(s))
}

cat("=== CCI alpha=.05: per-edge marks ===\n"); print(summarize_marks(cci_amat_05))
cat("\n=== CCI alpha=.01: per-edge marks ===\n"); print(summarize_marks(cci_amat_01))

cat("\n=== Edge-type counts, SAME 8-node main network, alpha=.05 ===\n")
cat("FCI:\n"); print(edge_type_counts(fci_amat_05))
cat("PC:\n");  print(edge_type_counts(pc_amat_05))
cat("CCI:\n"); print(edge_type_counts(cci_amat_05))

cat("\n=== CCI alpha=.01 edge-type counts (compare to CCI alpha=.05 above) ===\n")
print(edge_type_counts(cci_amat_01))
