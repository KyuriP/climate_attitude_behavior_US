# ==============================================================================
# 20_supplementary_measurement_figures.R
#
# Purpose: export clean, standalone Supplement figures for the exploratory
# measurement-reduction stage -- an item dendrogram, a PCA scree plot, and
# the Wave-5 behavior-item correlation structure. These reuse the exact
# objects already computed in climate_analysis_avg_v2_altweather.qmd (no new
# analysis is run here); this script only re-renders existing in-session
# results as standalone, publication-sized PDFs for main2.tex's Supplementary
# Materials S2, replacing the "% INSERT NEW FIGURE" comment placeholders
# there with real files.
#
# Run climate_analysis_avg_v2_altweather.qmd first (interactively, or knit
# through the relevant sections):
#   - Section 3.5 (Hierarchical clustering) and 3.6 (PCA) for `cor_kept`,
#     `hc_ward`, and `eig`.
#   - Section 3.8 (Wave-5 behavioral EDA) for `cor_beh`.
# ==============================================================================

stopifnot(exists("hc_ward"), exists("eig"), exists("cor_beh"))

dir.create("figures", showWarnings = FALSE)
dir.create("pipeline_outputs", showWarnings = FALSE)

# ---- S2 figure: attitude-item dendrogram (Ward.D2) --------------------------
# Same clustering already reported in Section 3.5 of the qmd and summarized
# in Supplementary S2; a clean single-panel version (Ward.D2 only, k = 8),
# sized for print, replacing the 4-panel exploratory linkage comparison.
pdf("figures/figS_item_dendrogram.pdf", width = 9, height = 5)
par(mar = c(2, 4, 3, 1))
plot(hc_ward, main = "Candidate climate-attitude items (Ward.D2 clustering)",
     cex = 0.8, xlab = "", sub = "", ylab = "Height (1 - r)")
rect.hclust(hc_ward, k = 8, border = 2:6)
dev.off()

# ---- S2 figure: PCA scree plot ----------------------------------------------
# Same eigen-decomposition already reported in Section 3.6 of the qmd
# (PC1 = 14.95, PC2 = 1.73); re-rendered as a standalone PDF with a cleaner,
# print-appropriate theme. Full eigenvalue vector also written out in case a
# reviewer wants the tail beyond PC1/PC2.
scree_df <- data.frame(component = seq_along(eig$values),
                        eigenvalue = eig$values)
write.csv(scree_df, "pipeline_outputs/pca_eigenvalues.csv", row.names = FALSE)

library(ggplot2)
p_scree <- ggplot(scree_df, aes(x = component, y = eigenvalue)) +
  geom_col(fill = "#4C72B0", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = scree_df$component) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "Component", y = "Eigenvalue", title = NULL)
ggsave("figures/figS_pca_scree.pdf", p_scree, width = 7, height = 4.5)

# ---- S2 figure: Wave-5 behavior-item correlation structure ------------------
# Same correlation matrix already reported in Section 3.8 of the qmd; saved
# as its own standalone PDF (via pheatmap's own file-output argument) rather
# than only appearing inline in the knitted notebook.
pheatmap::pheatmap(cor_beh, clustering_method = "ward.D2",
                    color = colorRampPalette(c("white", "#F1948A"))(100),
                    display_numbers = TRUE, number_format = "%.2f",
                    number_color = "black", fontsize_number = 9,
                    fontsize_row = 9, fontsize_col = 9,
                    main = "Wave-5 behavioral item correlations",
                    border_color = NA,
                    filename = "figures/figS_behavior_corr.pdf",
                    width = 6, height = 5)

cat("Saved:\n",
    " figures/figS_item_dendrogram.pdf\n",
    " figures/figS_pca_scree.pdf\n",
    " figures/figS_behavior_corr.pdf\n",
    " pipeline_outputs/pca_eigenvalues.csv\n")
