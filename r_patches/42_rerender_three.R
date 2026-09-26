setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
for (f in c("r_patches/36_revision_figures_1_2A.R", "r_patches/38_figure_causal_combined.R", "r_patches/07_figure5_scm_hierarchical_v3.R"))
  cat(f, ":", tryCatch({ source(f, local = new.env(parent = globalenv())); "OK" }, error = function(e) conditionMessage(e)), "\n")
