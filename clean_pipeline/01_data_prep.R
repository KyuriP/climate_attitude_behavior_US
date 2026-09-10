# data prep -- pulled straight from the qmd (libraries/load-data/cleaning/
# composites/wave restriction chunks). left out the pure EDA stuff (dist
# plots, correlation heatmaps, PCA/EFA) since that was just item screening,
# not part of the actual pipeline -- still in the qmd if needed.
#
# watch out: df_main/df_extended get reassigned partway through to the
# waves-1-4-only version. everything downstream should use those, not
# df_main_allwave/df_extended_allwave (those are just for the sensitivity
# check against the all-wave version).
#
# needs data_henry/*.parquet in the working dir.


suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(psych)
})

select   <- dplyr::select
filter   <- dplyr::filter
mutate   <- dplyr::mutate
arrange  <- dplyr::arrange
describe <- psych::describe
alpha    <- psych::alpha

# ---------------------------------------------------------------------------
# load-data
# ---------------------------------------------------------------------------
codebook    <- read_parquet("data_henry/codebook.parquet")
participant <- read_parquet("data_henry/participant.parquet")
rawdat      <- read_parquet("data_henry/w1w2w3w4w5_indices_weights_jul12_2022.parquet")

df_long <- rawdat |>
  dplyr::filter(PID %in% participant$participant_id,
                WAVE %in% c("1","2","3","4","5")) |>
  dplyr::rename(participant_id = PID) |>
  dplyr::mutate(wave = as.integer(as.character(WAVE)))

df_w5 <- df_long |> dplyr::filter(wave == 5)

screen_items_raw <- c(
  "cc1", "cc2", "cc3", "cc6",
  "cc4_world", "cc4_wealthUS", "cc4_poorUS", "cc4_comm",
  "cc5_world", "cc5_wealthUS", "cc5_poorUS", "cc5_comm",
  "cc_pol_tax", "cc_pol_car",
  "cvcc9_cc", "cvcc6",
  "cvcc4_should", "cvcc4_will", "cvcc4_personal",
  "cvcc_worryothers",
  "pol_ideology", "pol7", "pol7_pi",
  "ew5", "ew6",
  "cc10", "cc11", "cc12"
)

beh_items <- c(
  "cc_behavior_meat", "cc_behavior_travel", "cc_behavior_activ",
  "cc_behavior_discuss", "cc_behavior_evacuate", "cc_behavior_move"
)

cat("Rows:", nrow(df_long),
    "| Participants:", dplyr::n_distinct(df_long$participant_id),
    "| Waves:", paste(sort(unique(df_long$wave)), collapse = ","), "\n")
cat("Wave 5 N:", nrow(df_w5), "\n")

# ---------------------------------------------------------------------------
# cleaning-helpers
# ---------------------------------------------------------------------------
rescale01 <- function(x, from) (x - from[1]) / (from[2] - from[1])

na_special <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  dplyr::case_when(x %in% c(98, 99, -98, -99) ~ NA_real_, TRUE ~ x)
}

recode_binary01 <- function(x) {
  x <- na_special(x)
  dplyr::case_when(x == 1 ~ 0, x == 2 ~ 1, TRUE ~ NA_real_)
}

recode_cc1_midpoint <- function(x) {
  x  <- suppressWarnings(as.numeric(x))
  ux <- sort(unique(x[!is.na(x)]))
  if (all(ux %in% c(0, 1, 2))) {
    dplyr::case_when(x == 0 ~ 0, x == 1 ~ 0.5, x == 2 ~ 1, TRUE ~ NA_real_)
  } else if (all(ux %in% c(0, 1, 99))) {
    dplyr::case_when(x == 0 ~ 0, x == 1 ~ 1, x == 99 ~ 0.5, TRUE ~ NA_real_)
  } else {
    dplyr::case_when(x == 0 ~ 0, x == 1 ~ 0.5, x == 2 ~ 1, x == 99 ~ 0.5, TRUE ~ NA_real_)
  }
}

recode_cc1_dropdk <- function(x) {
  x  <- suppressWarnings(as.numeric(x))
  ux <- sort(unique(x[!is.na(x)]))
  if (all(ux %in% c(0, 1, 2))) {
    dplyr::case_when(x == 0 ~ 0, x == 2 ~ 1, TRUE ~ NA_real_)
  } else {
    dplyr::case_when(x == 0 ~ 0, x == 1 ~ 1, TRUE ~ NA_real_)
  }
}

recode_cc2_ord <- function(x) {
  x <- na_special(x)
  out <- dplyr::case_when(x == 4 ~ 0, x == 2 ~ 1, x == 3 ~ 2, x == 1 ~ 3, TRUE ~ NA_real_)
  rescale01(out, from = c(0, 3))
}

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# ---------------------------------------------------------------------------
# clean-long-data
# ---------------------------------------------------------------------------
df_long_h <- df_long |>
  dplyr::mutate(
    cc1_raw = suppressWarnings(as.numeric(cc1)),
    cc1     = recode_cc1_midpoint(cc1),

    cc3     = rescale01(na_special(cc3), from = c(1, 4)),
    cc6     = rescale01(na_special(cc6), from = c(1, 4)),

    cc4_world    = rescale01(na_special(cc4_world),    from = c(1, 4)),
    cc4_wealthUS = rescale01(na_special(cc4_wealthUS), from = c(1, 4)),
    cc4_poorUS   = rescale01(na_special(cc4_poorUS),   from = c(1, 4)),
    cc4_comm     = rescale01(na_special(cc4_comm),     from = c(1, 4)),

    cc5_world    = rescale01(na_special(cc5_world),    from = c(1, 4)),
    cc5_wealthUS = rescale01(na_special(cc5_wealthUS), from = c(1, 4)),
    cc5_poorUS   = rescale01(na_special(cc5_poorUS),   from = c(1, 4)),
    cc5_comm     = rescale01(na_special(cc5_comm),     from = c(1, 4)),

    cc_pol_tax   = rescale01(na_special(cc_pol_tax), from = c(1, 5)),
    cc_pol_car   = rescale01(na_special(cc_pol_car), from = c(1, 5)),

    cvcc9_cc         = rescale01(na_special(cvcc9_cc), from = c(1, 5)),
    cvcc4_should     = rescale01(na_special(cvcc4_should), from = c(1, 5)),
    cvcc4_will       = rescale01(na_special(cvcc4_will), from = c(1, 5)),
    cvcc4_personal   = rescale01(na_special(cvcc4_personal), from = c(1, 5)),
    cvcc_worryothers = rescale01(na_special(cvcc_worryothers), from = c(1, 4)),

    pol_ideology = rescale01(na_special(pol_ideology), from = c(1, 5)),
    ew5          = rescale01(na_special(ew5), from = c(1, 4)),
    ew6          = rescale01(na_special(ew6), from = c(1, 4)),

    cc2     = recode_cc2_ord(cc2),
    cvcc6   = rescale01(na_special(cvcc6), from = c(1, 5)),
    pol7    = recode_binary01(pol7),
    pol7_pi = recode_binary01(pol7_pi),
    cc10    = rescale01(na_special(cc10), from = c(1, 7)),
    cc11    = rescale01(na_special(cc11), from = c(1, 7)),
    cc12    = rescale01(na_special(cc12), from = c(1, 7))
  )

# ---------------------------------------------------------------------------
# clean-wave5-behavior
# ---------------------------------------------------------------------------
df_w5_h <- df_w5 |>
  dplyr::mutate(across(any_of(beh_items), ~rescale01(na_special(.x), from = c(1, 5))))

# ---------------------------------------------------------------------------
# participant-averages (all-wave version)
# ---------------------------------------------------------------------------
att_items_h <- screen_items_raw

df_att_person <- df_long_h |>
  dplyr::group_by(participant_id) |>
  dplyr::summarise(across(any_of(att_items_h), mean_or_na), .groups = "drop")

cat("Person-level averaged item dataset:", nrow(df_att_person), "x", ncol(df_att_person), "\n")

# ---------------------------------------------------------------------------
# composite-construction (all-wave version -- gets superseded below)
# ---------------------------------------------------------------------------
df_main <- df_att_person |>
  dplyr::transmute(
    participant_id,
    belief_concern = rowMeans(dplyr::pick(cc1, cc6), na.rm = TRUE),
    harm_present   = rowMeans(dplyr::pick(cc4_world, cc4_wealthUS, cc4_poorUS, cc4_comm), na.rm = TRUE),
    harm_future    = rowMeans(dplyr::pick(cc5_world, cc5_wealthUS, cc5_poorUS, cc5_comm), na.rm = TRUE),
    policy_support = rowMeans(dplyr::pick(cc_pol_tax, cc_pol_car), na.rm = TRUE),
    trust_science  = cvcc9_cc,
    social_norms   = cvcc4_should,
    politics       = pol_ideology,
    # ALT-SPEC: weather_risk_prep = ew5 alone. ew6 (preparedness behavior) is
    # excluded -- different construct (self-reported protective behavior, not
    # worry/risk perception) and empirically far more weakly connected to the
    # rest of the system. Variable name kept as weather_risk_prep so nothing
    # downstream needs to change.
    weather_risk_prep = ew5
  )

cat("Main composite dataset:", nrow(df_main), "x", ncol(df_main), "\n")
cat("Any missing values in df_main:", any(is.na(df_main)), "\n")
print(colSums(is.na(df_main)))

# ---------------------------------------------------------------------------
# cc1-sensitivity (diagnostic only -- does not feed anything downstream)
# ---------------------------------------------------------------------------
df_long_cc1 <- df_long |>
  dplyr::mutate(
    cc1_mid = recode_cc1_midpoint(cc1),
    cc1_na  = recode_cc1_dropdk(cc1),
    cc6_h   = rescale01(na_special(cc6), from = c(1, 4))
  )

df_cc1_sens <- df_long_cc1 |>
  dplyr::group_by(participant_id) |>
  dplyr::summarise(
    bel_mid = mean_or_na(rowMeans(cbind(cc1_mid, cc6_h), na.rm = TRUE)),
    bel_na  = mean_or_na(rowMeans(cbind(cc1_na,  cc6_h), na.rm = FALSE)),
    .groups = "drop"
  )

cat("[cc1 sensitivity] corr(midpoint, drop-DK) =",
    round(cor(df_cc1_sens$bel_mid, df_cc1_sens$bel_na, use = "complete.obs"), 3),
    "| missing in drop-DK version:", sum(is.na(df_cc1_sens$bel_na)), "\n")

# ---------------------------------------------------------------------------
# climate-behavior-composite
# ---------------------------------------------------------------------------
df_behavior_w5 <- df_w5_h |>
  dplyr::transmute(
    participant_id,
    beh_meat     = cc_behavior_meat,
    beh_travel   = cc_behavior_travel,
    beh_activ    = cc_behavior_activ,
    beh_discuss  = cc_behavior_discuss,
    beh_evacuate = cc_behavior_evacuate,
    beh_move     = cc_behavior_move,
    climate_behavior = rowMeans(
      dplyr::pick(beh_meat, beh_travel, beh_activ, beh_discuss, beh_evacuate, beh_move),
      na.rm = TRUE)
  )

cat("Wave-5 behavior dataset:", nrow(df_behavior_w5), "\n")
cat("Missing climate_behavior:", sum(is.na(df_behavior_w5$climate_behavior)), "\n")
cat("climate_behavior alpha =", round(psych::alpha(
  df_behavior_w5 |> dplyr::select(beh_meat, beh_travel, beh_activ, beh_discuss, beh_evacuate, beh_move)
)$total$raw_alpha, 3), "\n")

# ---------------------------------------------------------------------------
# analysis-datasets (all-wave version)
# ---------------------------------------------------------------------------
df_extended <- df_main |>
  dplyr::inner_join(df_behavior_w5 |> dplyr::select(participant_id, climate_behavior),
                     by = "participant_id")

cat("Main belief-network dataset (all-wave): N =", nrow(df_main), "\n")
cat("Extended (+ behavior) dataset (all-wave): N =", nrow(df_extended), "\n")

# ---------------------------------------------------------------------------
# Waves-1-4-only restriction -- THIS is the version actually used downstream.
# Preserves the qmd's own reassignment order exactly: all-wave df_main/
# df_extended are stashed as *_allwave, then df_main/df_extended are
# OVERWRITTEN with the wave<=4 version.
# ---------------------------------------------------------------------------
df_long_h_pre5 <- df_long_h |> dplyr::filter(wave <= 4)

df_att_person_pre5 <- df_long_h_pre5 |>
  dplyr::group_by(participant_id) |>
  dplyr::summarise(across(any_of(att_items_h), mean_or_na), .groups = "drop")

df_main_pre5 <- df_att_person_pre5 |>
  dplyr::transmute(
    participant_id,
    belief_concern    = rowMeans(dplyr::pick(cc1, cc6), na.rm = TRUE),
    harm_present      = rowMeans(dplyr::pick(cc4_world, cc4_wealthUS, cc4_poorUS, cc4_comm), na.rm = TRUE),
    harm_future       = rowMeans(dplyr::pick(cc5_world, cc5_wealthUS, cc5_poorUS, cc5_comm), na.rm = TRUE),
    policy_support    = rowMeans(dplyr::pick(cc_pol_tax, cc_pol_car), na.rm = TRUE),
    trust_science     = cvcc9_cc,
    social_norms      = cvcc4_should,
    politics          = pol_ideology,
    weather_risk_prep = ew5  # ALT-SPEC: ew5 alone
  )

df_main_allwave     <- df_main
df_extended_allwave <- df_extended

df_main     <- df_main_pre5
df_extended <- df_main_pre5 |>
  dplyr::inner_join(df_behavior_w5 |> dplyr::select(participant_id, climate_behavior),
                     by = "participant_id")

cat("\n*** FINAL df_main (Waves 1-4 only): N =", nrow(df_main), "***\n")
cat("*** FINAL df_extended (Waves 1-4 + Wave-5 behavior): N =", nrow(df_extended), "***\n")
cat("(Compare against the manuscript's stated N=1,987 / N=870 before trusting anything downstream.)\n")

# ---------------------------------------------------------------------------
# composite-reliability (diagnostic only, all-wave df_att_person -- matches
# original qmd chunk exactly, does not affect df_main/df_extended)
# ---------------------------------------------------------------------------
r2 <- function(x, y) round(cor(x, y, use = "complete.obs"), 3)

cat("\n=== Main attitude composites (reliability, all-wave df_att_person) ===\n")
cat("belief_concern r(cc1, cc6) =", r2(df_att_person$cc1, df_att_person$cc6), "\n")
cat("harm_present alpha =", round(psych::alpha(
  df_att_person |> dplyr::select(cc4_world, cc4_wealthUS, cc4_poorUS, cc4_comm))$total$raw_alpha, 3), "\n")
cat("harm_future alpha =", round(psych::alpha(
  df_att_person |> dplyr::select(cc5_world, cc5_wealthUS, cc5_poorUS, cc5_comm))$total$raw_alpha, 3), "\n")
cat("policy_support r(cc_pol_tax, cc_pol_car) =", r2(df_att_person$cc_pol_tax, df_att_person$cc_pol_car), "\n")
cat("weather_risk_prep (ALT-SPEC) = ew5 alone -- [reference only] r(ew5, ew6) =",
    r2(df_att_person$ew5, df_att_person$ew6), "-- NOT used to construct weather_risk_prep\n")
cat("trust_science = single item (cvcc9_cc); social_norms = single item (cvcc4_should); ",
    "politics = single item (pol_ideology)\n")
