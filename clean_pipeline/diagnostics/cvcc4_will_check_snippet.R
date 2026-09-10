# not part of the numbered pipeline, run by hand. checks whether
# social_norms built from cvcc4_will (descriptive norm) instead of
# cvcc4_should (injunctive, what's actually used) looks meaningfully
# different for climate_behavior.
#
# run after 01 + 02 in the same session (needs df_extended,
# df_att_person_pre5, network_ext already loaded).


stopifnot(exists("df_extended"), exists("df_att_person_pre5"), exists("network_ext"))

cvcc4_will_pre5 <- df_att_person_pre5 |>
  dplyr::select(participant_id, cvcc4_will)

check_df <- df_extended |>
  dplyr::inner_join(cvcc4_will_pre5, by = "participant_id")
stopifnot(nrow(check_df) == nrow(df_extended), !anyNA(check_df$cvcc4_will))

cat("N =", nrow(check_df), "\n\n")

cat("--- (1) Simple correlations with climate_behavior ---\n")
cat("  social_norms (cvcc4_should, RETAINED)  :",
    round(cor(check_df$social_norms, check_df$climate_behavior, use = "complete.obs"), 3), "\n")
cat("  cvcc4_will (descriptive norm, EXCLUDED):",
    round(cor(check_df$cvcc4_will, check_df$climate_behavior, use = "complete.obs"), 3), "\n")
cat("  cor(cvcc4_should, cvcc4_will) themselves:",
    round(cor(check_df$social_norms, check_df$cvcc4_will, use = "complete.obs"), 3), "\n\n")

cat("--- (2) Swap-and-refit: EBICglasso edge weight into climate_behavior ---\n")
cat("(same spec as 02_ggm.R's extended network -- tuning=.5, nonparanormal, EBICglasso)\n")
alt_net_dat <- check_df |>
  dplyr::mutate(social_norms = cvcc4_will) |>
  dplyr::select(dplyr::all_of(NODE_ORDER_EXT)) |>
  huge::huge.npn(npn.func = "truncation") |>
  as.data.frame()
colnames(alt_net_dat) <- NODE_ORDER_EXT
alt_net  <- bootnet::estimateNetwork(alt_net_dat, default = "EBICglasso", tuning = .5, corMethod = "cor")
alt_wmat <- qgraph::getWmat(alt_net)
orig_wmat <- qgraph::getWmat(network_ext)

cat("  using cvcc4_should (current, published):",
    round(orig_wmat["social_norms", "climate_behavior"], 3), "\n")
cat("  using cvcc4_will instead                :",
    round(alt_wmat["social_norms", "climate_behavior"], 3), "\n\n")

cat("If the cvcc4_will number is close to zero (or the edge drops out of the\n")
cat("graph entirely) while cvcc4_should's isn't, that's a quantitative version\n")
cat("of the same conclusion the qmd's construct-validity argument already\n")
cat("reached -- worth one sentence in the Supplement's construct-selection\n")
cat("discussion if you want the number to back up the qualitative argument.\n")
