# =============================================================================
# 11d: is CCI's failure on x1->x2 in 11c a root-cause limitation, or an
#      artifact of x1 being FIRST in the variable/label order?
#
# 11c's synthetic chain (x1 -> x2 -> x3 -> x4, no latents, no cycles) came
# back with x2->x3 and x3->x4 correctly resolved as directed, but x1-x2
# came back bidirected -- which is wrong (x1 is a clean root cause with no
# confounding). Before concluding anything about CCI's real behavior on
# the actual data, we need to know: does the SAME true edge (x1's edge)
# fail again once it's no longer first in the column order? If yes, this
# is a genuine root-cause/orientation limitation. If instead a DIFFERENT
# edge fails now (whichever variable is first this time), that points to
# an ordering-dependent bug in fit_cci()/CCI.KP rather than anything about
# root causes specifically.
#
# Run this after 11c (reuses sim_dat, fit_cci(), summarize_marks()).
# =============================================================================

sim_dat2     <- sim_dat[, c("x2", "x3", "x4", "x1")]
sim_labels2  <- names(sim_dat2)
sim_suffStat2 <- list(C = cor(sim_dat2), n = n_sim)

sim_amat2 <- fit_cci(sim_suffStat2, alpha = 0.05, labels = sim_labels2, p = 4)
cat("=== Same synthetic chain, column order permuted (x1 now LAST) ===\n")
print(sim_amat2)
cat("\n--- per-edge marks ---\n")
print(summarize_marks(sim_amat2))

cat("\nCompare: in 11c (x1 first), the x1-x2 edge came back bidirected while\n")
cat("x2-x3 and x3-x4 were correctly directed. Here, check which edge (if any)\n")
cat("comes back bidirected -- same true edge (x1's) or whichever variable is\n")
cat("first in THIS ordering (x2's edges)?\n")
