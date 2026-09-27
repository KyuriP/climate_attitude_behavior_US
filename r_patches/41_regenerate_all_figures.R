# 41_regenerate_all_figures.R (2026-09-26): re-export the manuscript figures
# with the final reader-facing labels (weather worry / political orientation /
# trust in science / belief and concern). Plotting only -- no model is refit
# except the GGM re-estimation that 36 already does (identical to 02_ggm.R).
# Not included: the item dendrogram and PCA scree figures
# (r_patches/20_supplementary_measurement_figures.R), which need objects from
# the exploratory .qmd in the same session, so run that one separately.
REPO <- "/Users/Kyuri1/Documents/Kyuri_P/PhD_UvA/Climate_Attitude_Behavior/climate_attitude_behavior_US"
if (dir.exists(REPO)) setwd(REPO)  # otherwise run from the repository root
run <- function(label, file, pre = NULL) {
  cat("\n######## ", label, " (", file, ")\n", sep = "")
  res <- tryCatch({ if (!is.null(pre)) pre(); source(file, local = new.env(parent = globalenv())); "OK" },
                  error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("######## ", label, ": ", res, "\n", sep = ""); invisible(res)
}
status <- c(
  fig1_and_2A_inputs = run("Figure 1 (GGM panels)", "r_patches/36_revision_figures_1_2A.R"),
  fig2 = run("Figure 2 (combined)", "r_patches/38_figure_causal_combined.R"),
  fig3 = run("Figure 3 (SCM)", "r_patches/07_figure5_scm_hierarchical_v3.R"),
  fig4_S7 = run("Figure 4 + S7", "r_patches/34_figure_bootstrap_lvida_ridge.R"),
  fig5_ncc = run("Figure 5, NCC version (shared-cause sensitivity)", "r_patches/51_figure_confound_sensitivity.R"),
  fig5_long = run("Figure 5, longer version (pathway dependence)", "r_patches/44_figure_pathway_dependence.R"),
  fig45_ncc_option = run("Optional NCC Figure 4 (Figures 4 and 5 combined)", "r_patches/52_figure_sensitivity_combined.R"),
  S1 = run("Figure S1 (heatmap)", "r_patches/04_figure1_heatmap.R"),
  S4 = run("Figure S4 (stability matrices)", "clean_pipeline/20_stability_matrix_figures.R",
           pre = function() { assign("fci_marks", c("N", "o", ">", "-"), envir = globalenv())
                              assign("pc_marks", c("N", "-", ">"), envir = globalenv())
                              assign("node_order_cd", NODE_ORDER_MAIN, envir = globalenv()) }),
  S5 = run("Figure S5 (all combos)", "r_patches/08b_supp_figure_all_combo_interventions.R"),
  S6 = run("Figure S6 (cyclic)", "r_patches/35_figure_cyclic_feedback_equilibrium.R"))
cat("\n==== SUMMARY ====\n"); print(status)
