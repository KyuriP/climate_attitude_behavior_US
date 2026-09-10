# ggm setup -- network-setup/network-main chunks from the qmd. same
# construction 05_figures2_3_ggm_redesign.R already uses, so running this
# first just means that script skips re-estimating the network.
#
# needs df_main/df_extended from 01. produces network_main/network_ext.


source("clean_pipeline_alt_will/00_config.R")

stopifnot(exists("df_main"), exists("df_extended"))

suppressPackageStartupMessages({
  library(huge)
  library(qgraph)
  library(bootnet)
  library(dplyr)
})

main_nodes     <- NODE_ORDER_MAIN
extended_nodes <- NODE_ORDER_EXT

# ---- Nonparanormal transform (huge::huge.npn, truncation) -- verbatim -----
df_net_main_raw <- df_main |> dplyr::select(dplyr::all_of(main_nodes))
df_net_main <- df_net_main_raw |> huge::huge.npn(npn.func = "truncation") |> as.data.frame()
colnames(df_net_main) <- main_nodes

df_net_ext_raw <- df_extended |> dplyr::select(dplyr::all_of(extended_nodes))
df_net_ext <- df_net_ext_raw |> huge::huge.npn(npn.func = "truncation") |> as.data.frame()
colnames(df_net_ext) <- extended_nodes

cat("Main nodes:", paste(main_nodes, collapse = ", "), "\n")
cat("Extended nodes:", paste(extended_nodes, collapse = ", "), "\n")

# ---- EBICglasso network estimation (tuning=.5) -- verbatim -----------------
network_main <- bootnet::estimateNetwork(
  df_net_main, default = "EBICglasso", tuning = .5, corMethod = "cor"
)
network_ext <- bootnet::estimateNetwork(
  df_net_ext, default = "EBICglasso", tuning = .5, corMethod = "cor"
)

W_main <- network_main$graph
W_ext  <- network_ext$graph

n_possible_main <- choose(length(main_nodes), 2)
n_retained_main <- sum(W_main[upper.tri(W_main)] != 0)
cat(sprintf("\nMain GGM: %d of %d possible edges retained.\n", n_retained_main, n_possible_main))
cat("(figures/fig2_ggm_main_fixed.pdf has exactly 23 -- checked mechanically\n")
cat(" against its vector line objects, not eyeballed. the qmd's own prose\n")
cat(" elsewhere says 24 and is the stale one -- this printed count is the one\n")
cat(" to trust, it comes straight from W_main.)\n")
