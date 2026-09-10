# =============================================================================
# 17_diagnostic_rcot_social_norms_behavior.R
#
# Targeted diagnostic, NOT a full sensitivity-check script like 14/15/16 --
# this doesn't run fci()/pc() at all. Purpose: script 16's bootstrap (N_BOOT=10)
# came back p_adjacent_RCoT = 0.0 for social_norms-climate_behavior, but that
# number is ambiguous on its own -- FCI's skeleton search drops an edge the
# moment ANY ONE tested conditioning set fails to reject independence (you
# need every test to reject to keep an edge, only one non-rejection to lose
# it), so a lower-power test can look "sparser" even when a real association
# exists, especially compounded across the many conditioning sets FCI tries
# internally and quantized to 0.1 increments at N_BOOT=10. This script
# decomposes the question directly instead of guessing from the aggregate:
#   1. Is there ANY unconditional association between social_norms and
#      climate_behavior at all, under a fully nonparametric test?
#   2. If yes, does conditioning on any ONE other main-network variable kill
#      it -- and if so, which one? (That would be a real, reportable
#      mediation/confounding story, not a null result.)
#   3. Does it survive conditioning on ALL other main-network variables at
#      once (closest to what a saturated regression/FCI's most demanding
#      test would ask)?
# Also checks p-value STABILITY across several random seeds -- RCoT's default
# num_f2=5 (features for the non-conditioning variables) is a small number,
# and its null approximation depends on the specific random Fourier feature
# draw, so a single seed's p-value could be noisier than it looks. If the
# verdict flips across seeds, that itself is informative (says num_f2=5 is
# too low for a stable answer here, not that the true relationship is
# ambiguous).
#
# Requires the same `library(RCIT)` fix as script 16 -- see that script's
# header for the full diagnosis of the hbe/momentchi2 packaging bug this
# works around.
#
# Not tested against the live data -- tested against synthetic data with a
# known true structure
# (see TESTING NOTE below) to confirm the script itself runs correctly and
# that its printed interpretation matches a known ground truth.
# =============================================================================

if (!requireNamespace("RCIT", quietly = TRUE)) {
  stop("RCIT is required and isn't installed -- see script 16's header for ",
       "install instructions.")
}
suppressPackageStartupMessages(library(RCIT))  # attaches RCIT + momentchi2 + MASS

stopifnot(exists("df_extended"))

if (exists("node_order_cd")) {
  MAIN_NODES <- node_order_cd
} else {
  MAIN_NODES <- c(
    "belief_concern", "harm_present", "harm_future", "policy_support",
    "trust_science", "social_norms", "politics", "weather_risk_prep"
  )
  message("node_order_cd not found in the session -- falling back to the ",
          "hardcoded 8-node list (see scripts 14/15/16 headers). Confirm ",
          "this still matches Section 7.1's node_order_cd.")
}
OTHER_NODES <- setdiff(MAIN_NODES, "social_norms")

sn <- as.numeric(df_extended[["social_norms"]])
cb <- as.numeric(df_extended[["climate_behavior"]])

SEEDS <- 1:20  # cheap (RCoT is ~ms per call) -- checks p-value stability,
               # not just a single draw's answer

run_across_seeds <- function(z = NULL, label) {
  ps <- vapply(SEEDS, function(s) {
    out <- if (is.null(z)) {
      RCIT::RCoT(sn, cb, seed = s)
    } else {
      RCIT::RCoT(sn, cb, z, seed = s)
    }
    out$p
  }, numeric(1))
  cat(sprintf("%-45s  median p=%.4f  [min=%.4f, max=%.4f]  seeds<.05: %d/%d\n",
              label, median(ps), min(ps), max(ps), sum(ps < 0.05), length(SEEDS)))
  invisible(ps)
}

cat("=============================================================\n")
cat("social_norms <-> climate_behavior: RCoT across", length(SEEDS), "seeds\n")
cat("=============================================================\n\n")

cat("---- Step 1: UNCONDITIONAL (no other variables at all) ----\n")
p_uncond <- run_across_seeds(NULL, "unconditional")
cat("\nInterpretation: if this is consistently significant (small p, stable\n")
cat("across seeds), there IS a real raw association -- any later 'independence'\n")
cat("found below is about what EXPLAINS it, not whether it exists at all. If\n")
cat("this is consistently non-significant, that's a clean, assumption-free\n")
cat("null result on its own, independent of anything FCI's search does.\n\n")

cat("---- Step 2: conditioning on ONE other main-network variable at a time ----\n")
p_single <- list()
for (v in OTHER_NODES) {
  p_single[[v]] <- run_across_seeds(as.numeric(df_extended[[v]]), sprintf("conditioning on %s", v))
}
cat("\nInterpretation: a variable that flips this from significant to\n")
cat("non-significant is a candidate mediator/confounder of the social_norms-\n")
cat("climate_behavior link -- worth naming specifically if so, rather than\n")
cat("reporting a blanket 'edge not supported.'\n\n")

cat("---- Step 3: conditioning on ALL other main-network variables at once ----\n")
p_all <- run_across_seeds(as.matrix(df_extended[, OTHER_NODES, drop = FALSE]),
                           "conditioning on all others")
cat("\nInterpretation: this is the closest match to what FCI's skeleton search\n")
cat("actually demands (surviving every tested conditioning set) -- if step 1\n")
cat("was significant but this isn't, the bootstrap's p_adjacent_RCoT=0.0 is\n")
cat("more likely a real 'explained away by the rest of the network' result\n")
cat("than a power artifact.\n")

cat("\n=============================================================\n")
cat("Summary table\n")
cat("=============================================================\n")
summary_df <- data.frame(
  test = c("unconditional", paste0("| ", OTHER_NODES), "| all others"),
  median_p = c(median(p_uncond), vapply(p_single, median, numeric(1)), median(p_all)),
  frac_seeds_sig = c(mean(p_uncond < 0.05),
                      vapply(p_single, function(p) mean(p < 0.05), numeric(1)),
                      mean(p_all < 0.05))
)
print(summary_df, row.names = FALSE)

# =============================================================================
# Step 4 (added after seeing step 3's real instability on the dataset): is the
# seed-to-seed noise in "conditioning on all others" just num_f=100 (the
# default number of random Fourier features used for the CONDITIONING set)
# being too low for a 7-dimensional z, or is it a real property of the
# relationship? Tested this on synthetic data first,
# before adding it here -- and found a real trap worth knowing about before
# just cranking num_f up: at num_f approaching or exceeding the sample size n
# (tested n=500, num_f=500-1000), a GENUINELY NULL relationship started
# returning spurious near-zero p-values in most seeds -- a real false-positive
# problem from the random-feature dimension outgrowing the data, not evidence
# of a wrongly "unstable" true test. At a larger n (2000, closer to the live
# main-network size) the same null case stayed clean (p mostly .6-.9, 0/20
# significant) even at num_f=1000 -- so the danger zone is specifically
# num_f approaching n, not num_f in absolute terms. This step therefore caps
# the num_f values it tries at n/4 (a conservative margin below where the
# false-positive problem started appearing in testing) rather than trying
# arbitrarily large values -- if you want to push higher than that cap,
# increase MAX_NUM_F_FRACTION below deliberately, but re-run the TRUE-NULL
# sanity check in this script's own testing notes first, since the safe
# threshold depends on n and hasn't been mapped precisely.
# =============================================================================
cat("\n=============================================================\n")
cat("Step 4: is 'conditioning on all others'' instability a num_f artifact?\n")
cat("=============================================================\n")
n_data <- nrow(df_extended)
MAX_NUM_F_FRACTION <- 0.25  # conservative cap -- see comment above
num_f_candidates <- unique(pmin(c(100, 200, 400, 800), floor(n_data * MAX_NUM_F_FRACTION)))
num_f_candidates <- num_f_candidates[num_f_candidates >= 100]
cat(sprintf("n = %d -- trying num_f in {%s} (capped at %.0f%% of n = %d)\n\n",
            n_data, paste(num_f_candidates, collapse = ", "),
            MAX_NUM_F_FRACTION * 100, floor(n_data * MAX_NUM_F_FRACTION)))

z_all <- as.matrix(df_extended[, OTHER_NODES, drop = FALSE])
for (nf in num_f_candidates) {
  ps_nf <- vapply(SEEDS, function(s) RCIT::RCoT(sn, cb, z_all, num_f = nf, seed = s)$p, numeric(1))
  cat(sprintf("num_f=%-5d  median p=%.4f  [min=%.4f, max=%.4f]  seeds<.05: %d/%d\n",
              nf, median(ps_nf), min(ps_nf), max(ps_nf), sum(ps_nf < 0.05), length(SEEDS)))
}
cat("\nInterpretation: if the seeds<.05 count climbs toward 20/20 and the p-value\n")
cat("range tightens as num_f increases (without approaching the false-positive\n")
cat("zone above), that supports 'step 3's instability was a feature-count power\n")
cat("issue, and the true joint-conditional relationship is more likely real than\n")
cat("the raw num_f=100 result suggested. If it stays unstable even at the ")
cat("largest\nsafe num_f tried, that's better evidence the instability is real,\n")
cat("not an artifact.\n")

message("\nDone. Paste this whole console block (including the per-line seed ranges ",
        "above and Step 4's num_f table) back for interpretation -- the summary ",
        "table alone loses the seed-stability information, which matters here.")
