# =============================================================================
# 36_revision_figures_1_2A.R  (2026-09-26, reviewer figure comments)
#
# Regenerates only the manuscript figures touched by the reviewer's plotting
# comments, in a fresh session, from the repo root:
#   Rscript r_patches/36_revision_figures_1_2A.R
#
#   Figure 1  (fig2_ggm_main_fixed.pdf / fig3_ggm_extended_fixed.pdf):
#             one shared layout, identical positions for the 8 shared nodes.
#   Figure 2A (fig4_stability_dotmatrix.pdf): row labels "Name (ABBR)",
#             WR -> WW.
#
# Plotting only: networks are re-estimated exactly as in 02_ggm.R, and the
# FCI bootstrap proportions are READ from the frozen 1,000-resample run named
# in LATEST_bootstrap_run.csv -- nothing is re-bootstrapped.
# =============================================================================
setwd("/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US")
source("clean_pipeline/00_config.R")
source("clean_pipeline/01_data_prep.R")
source("clean_pipeline/02_ggm.R")

cat("Retained edges, main network:",
    sum(network_main$graph[upper.tri(network_main$graph)] != 0), "\n")

source("r_patches/05_figures2_3_ggm_redesign.R")

latest         <- read.csv(LATEST_POINTER_PATH, stringsAsFactors = FALSE)
fci_props_ext  <- readRDS(latest$fci_props_ext_path)
node_order_ext <- NODE_ORDER_EXT
abbr_ext       <- ABBR_EXT
cat("FCI proportions read from:", latest$fci_props_ext_path, "\n")

source("r_patches/06_figure4_stability_redesign.R")
cat("DONE\n")
