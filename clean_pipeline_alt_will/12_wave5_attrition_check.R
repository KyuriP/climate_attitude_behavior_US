# wave5 attrition check, from 24_wave5_attrition_check.R. feeds supp S11's
# opening attrition numbers (age SMD .29 etc). compares the N=870 analytic
# sample against the 1,117 who have waves 1-4 but no wave-5 behavior, on the
# 8 nodes + demographics, using SMD instead of p-values.
#
# demographic var names checked against data_henry's own codebook/raw data --
# the Core Wave 1 codebook pdf is a different (european) survey, doesn't
# apply here, don't use it for these variable names.


source("clean_pipeline_alt_will/00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

stopifnot(exists("df_main"), exists("df_extended"), exists("df_long"))
stopifnot(all(c("dem_age", "dem_educ", "dem_income", "dem_male",
                "dem_race_cat", "dem_urban_cat") %in% names(df_long)))

retained_ids <- df_extended$participant_id
attrited_ids <- setdiff(df_main$participant_id, retained_ids)

cat("Retained (has Wave-5 behavior), N =", length(retained_ids), "\n")
cat("Attrited (Waves 1-4 only, no Wave-5 behavior), N =", length(attrited_ids), "\n")
stopifnot(length(retained_ids) == 870)
stopifnot(length(retained_ids) + length(attrited_ids) == nrow(df_main))

group_of <- function(id) ifelse(id %in% retained_ids, "retained_N870", "attrited")

# --- (a) Waves-1-4 node comparison: standardized mean differences ----------
node_cols <- NODE_ORDER_MAIN  # the 8 retained attitude/context nodes, already in df_main
smd_numeric <- function(x, g) {
  m1 <- mean(x[g == "retained_N870"], na.rm = TRUE); s1 <- sd(x[g == "retained_N870"], na.rm = TRUE)
  m2 <- mean(x[g == "attrited"],      na.rm = TRUE); s2 <- sd(x[g == "attrited"],      na.rm = TRUE)
  pooled_sd <- sqrt((s1^2 + s2^2) / 2)
  list(mean_retained = m1, mean_attrited = m2, smd = (m1 - m2) / pooled_sd)
}

df_main_grp <- df_main |> dplyr::mutate(group = group_of(participant_id))

node_comparison <- purrr::map_dfr(node_cols, function(nd) {
  r <- smd_numeric(df_main_grp[[nd]], df_main_grp$group)
  tibble::tibble(variable = nd, mean_retained = round(r$mean_retained, 3),
                 mean_attrited = round(r$mean_attrited, 3), smd = round(r$smd, 3))
}) |> dplyr::mutate(flag_ge_.10 = abs(smd) >= .10)

cat("\n=== (a) Waves-1-4 node comparison: retained (N=870) vs. attrited (N=",
    length(attrited_ids), ") ===\n", sep = "")
print(as.data.frame(node_comparison), row.names = FALSE)

# --- (b) Demographics: one row per participant (values repeat identically
# across waves in this dataset; dplyr::first(na.omit(.)) is defensive in
# case any participant has a wave-specific gap). ----------------------------
dem_participant <- df_long |>
  dplyr::group_by(participant_id) |>
  dplyr::summarise(
    age       = dplyr::first(na.omit(dem_age)),
    educ      = dplyr::first(na.omit(dem_educ)),
    income    = dplyr::first(na.omit(dem_income)),
    male      = dplyr::first(na.omit(dem_male)),
    race_cat  = dplyr::first(na.omit(as.character(dem_race_cat))),
    urban_cat = dplyr::first(na.omit(as.character(dem_urban_cat))),
    .groups = "drop"
  ) |>
  dplyr::filter(participant_id %in% c(retained_ids, attrited_ids)) |>
  dplyr::mutate(group = group_of(participant_id))

cat("\nDemographic coverage: N =", nrow(dem_participant),
    "of", length(retained_ids) + length(attrited_ids), "analytic-sample participants.\n")

dem_numeric_comparison <- purrr::map_dfr(c("age", "educ", "income"), function(v) {
  r <- smd_numeric(dem_participant[[v]], dem_participant$group)
  tibble::tibble(variable = v, mean_retained = round(r$mean_retained, 3),
                 mean_attrited = round(r$mean_attrited, 3), smd = round(r$smd, 3))
}) |> dplyr::mutate(flag_ge_.10 = abs(smd) >= .10)

prop_diff_table <- function(data, var, group_col = "group") {
  tab <- data |>
    dplyr::count(.data[[group_col]], .data[[var]]) |>
    tidyr::pivot_wider(names_from = all_of(group_col), values_from = n, values_fill = 0) |>
    dplyr::mutate(
      p_retained = retained_N870 / sum(retained_N870),
      p_attrited = attrited / sum(attrited),
      std_prop_diff = round(
        (p_retained - p_attrited) / sqrt((p_retained * (1 - p_retained) + p_attrited * (1 - p_attrited)) / 2),
        3
      ),
      p_retained = round(p_retained, 3), p_attrited = round(p_attrited, 3)
    )
  names(tab)[names(tab) == var] <- "category"
  tab |> dplyr::mutate(variable = var, flag_ge_.10 = abs(std_prop_diff) >= .10) |>
    dplyr::select(variable, category, p_retained, p_attrited, std_prop_diff, flag_ge_.10)
}

dem_gender_comparison <- dem_participant |>
  dplyr::mutate(male_cat = dplyr::case_when(male == 1 ~ "male", male == 0 ~ "female", TRUE ~ "self_described")) |>
  prop_diff_table(var = "male_cat")

dem_race_comparison  <- prop_diff_table(dem_participant, "race_cat")
dem_urban_comparison <- prop_diff_table(dem_participant, "urban_cat")

cat("\n=== (b) Demographics: continuous/ordinal ===\n")
print(as.data.frame(dem_numeric_comparison), row.names = FALSE)
cat("\n=== (b) Demographics: gender ===\n")
print(as.data.frame(dem_gender_comparison), row.names = FALSE)
cat("\n=== (b) Demographics: race ===\n")
print(as.data.frame(dem_race_comparison), row.names = FALSE)
cat("\n=== (b) Demographics: urbanicity ===\n")
print(as.data.frame(dem_urban_comparison), row.names = FALSE)

all_flags <- c(node_comparison$flag_ge_.10, dem_numeric_comparison$flag_ge_.10,
                dem_gender_comparison$flag_ge_.10, dem_race_comparison$flag_ge_.10,
                dem_urban_comparison$flag_ge_.10)
n_flagged <- sum(all_flags)
summary_msg <- if (n_flagged == 0) {
  paste0("Nearly all differences are small -- report this table in the Supplement ",
         "and skip the IPW sensitivity.")
} else {
  paste0("Some differences exceed .10 -- review which variables before deciding whether an ",
         "IPW sensitivity for the behavior-facing analyses is warranted.")
}
cat("\n", n_flagged, " of ", length(all_flags),
    " comparisons have |SMD| (or standardized proportion difference) >= .10.\n",
    summary_msg, "\n", sep = "")
cat("\nManuscript (Supplementary S11) reports age SMD=.29 (55.4 vs 50.9 years),\n")
cat("income SMD=.19, and urbanicity (suburban +.15, rural -.10) as the flagged\n")
cat("variables -- compare directly against the tables above.\n")

write.csv(node_comparison,        file.path(OUTPUT_DIR, "attrition_waves1to4_comparison.csv"), row.names = FALSE)
write.csv(dem_numeric_comparison, file.path(OUTPUT_DIR, "attrition_demographics_numeric.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(dem_gender_comparison, dem_race_comparison, dem_urban_comparison),
          file.path(OUTPUT_DIR, "attrition_demographics_categorical.csv"), row.names = FALSE)
cat("\nWrote attrition_waves1to4_comparison.csv, attrition_demographics_numeric.csv,",
    "attrition_demographics_categorical.csv to", OUTPUT_DIR, "\n")
