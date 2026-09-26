# re-fit the single-run nine-node FCI-JCI PAG at alpha=.05 on the current data
# and report the trust_science--policy_support endpoints (circle vs tail)
setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
suppressPackageStartupMessages(library(pcalg))
source("clean_pipeline/00_config.R")
stopifnot(exists("df_net_ext"))
ne <- NODE_ORDER_EXT; agg <- df_net_ext[, ne]; ctx <- which(ne %in% NODE_ORDER_MAIN)
ss <- list(C = cor(as.matrix(agg)), n = nrow(agg))
f05 <- pcalg::fci(ss, pcalg::gaussCItest, alpha = .05, labels = ne, contextVars = ctx,
                  jci = "1", selectionBias = FALSE, verbose = FALSE)
saved <- readRDS("pipeline_outputs/fci_ext_05_singlerun.rds")@amat
code <- c("0" = "none", "1" = "circle", "2" = "arrowhead", "3" = "tail")
cat("refit: mark at trust_science =", code[as.character(f05@amat["policy_support", "trust_science"])],
    "| mark at policy_support =", code[as.character(f05@amat["trust_science", "policy_support"])], "\n")
cat("refit identical to saved rds:", identical(unname(f05@amat), unname(saved[ne, ne])), "\n")
cat("any circle marks anywhere in refit PAG:", any(f05@amat == 1), "\n")
# the plotAG() exploratory version in the qmd used pairwise-complete cor
ss2 <- list(C = cor(agg, use = "pairwise.complete.obs"), n = nrow(agg))
f05b <- pcalg::fci(ss2, pcalg::gaussCItest, alpha = .05, labels = ne, contextVars = ctx,
                   jci = "1", selectionBias = FALSE, verbose = FALSE)
cat("qmd-style (pairwise cor): mark at trust_science =", code[as.character(f05b@amat["policy_support", "trust_science"])], "\n")
cat("39 DONE\n")
