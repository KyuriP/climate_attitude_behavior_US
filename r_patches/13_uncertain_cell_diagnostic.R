# =============================================================================
# 13: why did each "uncertain" (o) cell in Figure 4 / figS_stability_pc land
#     there? plot_edge_matrix()'s classifier (r_patches/10, Part 4) treats
#     "uncertain" as a catch-all for three different situations:
#       (1) either end's dominant mark is a circle
#       (2) the two ends disagree about whether the edge exists at all
#           (one end's dominant mark is N, the other's isn't)
#       (3) both ends are dominantly tail ("mutual" pattern)
#     This prints the actual dominant mark + proportion at each end for
#     every cell currently classified "uncertain", for both FCI and PC, so
#     you can see which of the three situations applies to which cell.
#
# Run after r_patches/10 (needs fci_props_combined, pc_props_combined,
# fci_marks, pc_marks, node_order_cd already in memory).
# =============================================================================
inspect_uncertain <- function(props, marks, nodes = node_order_cd) {
  p_local <- length(nodes)
  cells <- expand.grid(i = seq_len(p_local), j = seq_len(p_local),
                        stringsAsFactors = FALSE) |>
    dplyr::filter(i > j) |>
    dplyr::mutate(row_node = nodes[i], col_node = nodes[j]) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      mark_at_col = marks[which.max(props[i, j, ])],
      mark_at_row = marks[which.max(props[j, i, ])],
      prop_at_col = max(props[i, j, ]),
      prop_at_row = max(props[j, i, ]),
      edge_cat = dplyr::case_when(
        mark_at_row == "N" & mark_at_col == "N" ~ "absent",
        mark_at_row == "-" & mark_at_col == ">" ~ "directed",
        mark_at_row == ">" & mark_at_col == "-" ~ "directed",
        mark_at_row == ">" & mark_at_col == ">" ~ "bidirected",
        TRUE ~ "uncertain"
      )
    ) |> dplyr::ungroup() |>
    dplyr::filter(edge_cat == "uncertain") |>
    dplyr::select(row_node, col_node, mark_at_row, prop_at_row, mark_at_col, prop_at_col)
  cells
}

cat("=== FCI: cells classified 'uncertain' -- why? ===\n")
print(inspect_uncertain(fci_props_combined, fci_marks))
cat("\n=== PC-stable: cells classified 'uncertain' -- why? ===\n")
print(inspect_uncertain(pc_props_combined, pc_marks))
