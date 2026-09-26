# 46_ggm_connectedness_vs_effect.R (2026-09-26): read-only check for the framing
# "the most connected nodes are not necessarily the most effective intervention
# targets". Re-estimates the GGMs exactly as 02_ggm.R, then tabulates for each
# attitude/context node: (i) strength centrality in the 8-node GGM, (ii) strength
# in the 9-node GGM, (iii) its partial correlation with climate behavior, and
# (iv) the baseline working-SCM single-node effect on behavior.
setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
source("clean_pipeline/00_config.R")
source("clean_pipeline/01_data_prep.R")
source("clean_pipeline/02_ggm.R")
nodes <- rownames(W_main)
str_main <- rowSums(abs(W_main))
str_ext  <- rowSums(abs(W_ext))[nodes]
pcor_cb  <- W_ext[nodes, "climate_behavior"]
ate <- read.csv("pipeline_outputs/tables/confound_sensitivity_singlenode_ate.csv")
ate <- ate[ate$spec == "baseline", ]
eff <- setNames(ate$ate_climate_behavior, ate$node)[nodes]
out <- data.frame(node = nodes, strength_8node = round(str_main, 3), strength_9node = round(str_ext, 3),
                  pcor_with_behavior = round(pcor_cb, 3), scm_effect = round(eff, 3))
out$rank_strength_8 <- rank(-out$strength_8node)
out$rank_pcor_cb    <- rank(-abs(out$pcor_with_behavior))
out$rank_effect     <- rank(-out$scm_effect)
out <- out[order(-out$scm_effect), ]
print(out, row.names = FALSE)
write.csv(out, "pipeline_outputs/tables/ggm_connectedness_vs_scm_effect.csv", row.names = FALSE)
cat("46 DONE\n")

# cross-check against numbers quoted in the manuscript text
cat("\nW_main HP-HF:", round(W_main["harm_present","harm_future"],3),
    " PS-TS:", round(W_main["policy_support","trust_science"],3), "\n")
cat("W_ext CB row:\n"); print(round(W_ext["climate_behavior", ], 3))
cat("n rows df_net_ext:", nrow(df_net_ext), " nodes:", paste(colnames(df_net_ext), collapse=","), "\n")
cat("46b DONE\n")
