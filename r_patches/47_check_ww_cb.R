# 47_check_ww_cb.R: where does the manuscript's "partial r = .229" (weather worry -- behavior) come from?
setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
if (!exists("df_net_ext")) { source("clean_pipeline/00_config.R"); source("clean_pipeline/01_data_prep.R"); source("clean_pipeline/02_ggm.R") }
S <- cor(df_net_ext); P <- -cov2cor(solve(S)); diag(P) <- 1
cat("zero-order r(WW,CB) npn:", round(S["weather_risk_prep","climate_behavior"],3), "\n")
cat("unregularized partial r(WW,CB):", round(P["weather_risk_prep","climate_behavior"],3), "\n")
cat("unregularized partial r row CB:\n"); print(round(P["climate_behavior",],3))
cat("EBICglasso W_ext(WW,CB):", round(W_ext["weather_risk_prep","climate_behavior"],3), "\n")
Sr <- cor(df_net_ext_raw, use="pairwise.complete.obs"); cat("zero-order r(WW,CB) raw:", round(Sr["weather_risk_prep","climate_behavior"],3), "\n")
cat("47 DONE\n")
