# clean_pipeline/ -- confirmatory analysis, out of the qmd

Pulled the confirmatory chain out of `climate_analysis_avg_v2_altweather.qmd`
into scripts that run top to bottom on their own, instead of depending on
whatever's still sitting in an R session. The qmd stays as the exploratory
record (EDA, clustering, EFA, the construct-validity calls) -- that part
doesn't need to be a clean script, it's a one-time decision trail.

## raw data not included

`01_data_prep.R` reads the raw participant-level survey files from
`data_henry/` (`w1w2w3w4w5_indices_weights_jul12_2022.parquet`,
`codebook.parquet`, `participant.parquet`). That folder is gitignored and
not part of this repo, since it's individual-level survey data with
demographic fields, not something to publish alongside the code. Running
the pipeline from scratch requires getting those files separately; get in
touch for access to the raw data.

## run order

```
run_all.R                                 (optional -- runs 01-14 fresh and freezes the result, see "reference runs" below)
01_data_prep.R
02_ggm.R
03_bootstrap_causal_discovery.R           (also gives single-run fci_05/fci_01/pc_05/pc_01 for fig 2B/supp fig 8)
04_scm_finalize.R
05_scm_intervention_helpers.R             (sourced automatically by 06-16, don't run directly)
06_intervention_ates_singlenode.R
07_intervention_ates_8node.R
08_intervention_ates_combo.R
09_intervention_bootstrap_ci.R
10_orientation_enumeration_fit.R
11_interaction_moderation.R
12_wave5_attrition_check.R
13_ipw_attrition_sensitivity.R            (needs 12 conceptually, not a code dependency)
14_behavior_outcome_sensitivity.R         (needs 03's run_one_ext/node_order_ext/context_idx still live in session -- see its header)
15_confound_sensitivity_diagnostic.R      (optional -- needs 01+02 live, see its header)
16_extended_pag_lvida.R                   (optional -- needs 03's agg_ext/node_order_ext/context_idx live, and lv-ida downloaded separately, see its header)
17_edge_confounding_classification.R      (optional -- just reads pipeline_outputs csvs, no live session objects needed)
```

Everything's pulled with `sed` from the qmd's / r_patches' actual text, not
retyped, so the numerical logic matches what already produced the current
numbers, except where a script says otherwise (11's exact cross-check;
05/13 sharing one `base_edges` instead of separate copies). 15, 16, and 17
are new, not migrated from anywhere -- see "still open" below for why they
exist.

## what each one does

- `00_config.R` -- N_BOOT (1,000, the causal-discovery bootstrap),
  N_BOOT_INTERVENTION (1,000, the separate participant-level one), alpha
  thresholds, node lists, EXISTENCE_MIN, output paths.
- `run_all.R` -- sources 01-14 in order in one fresh session, then freezes
  whatever that run actually wrote (by mtime, not a hand-typed filename
  list) into `pipeline_outputs/reference_runs/<timestamp>/`, alongside a
  manifest (sessionInfo(), the exact scm_edges_finalized.csv used + its
  md5, elapsed time) and a check of the 8 single-node baseline ATEs against
  the values locked in on 2026-09-08. refuses to run if 15/16/17 objects
  are already in the session -- restart R first, this is meant to be a
  clean run, not a continuation of whatever was open before. see
  "reference runs" below.
- `01_data_prep.R` -- raw parquet through df_main/df_extended (waves 1-4
  only; all-wave versions kept as df_main_allwave/df_extended_allwave for
  the sensitivity check). also produces df_long, df_behavior_w5 for 12/13/14.
- `02_ggm.R` -- npn transform + EBICglasso (tuning=.5), both networks.
  prints the retained-edge count off the actual fitted graph (23 of 28 is
  right -- checked against fig2's actual vector line objects).
- `03_bootstrap_causal_discovery.R` -- FCI/PC-stable bootstrap, both
  networks, dated + manifested output so nothing downstream can silently
  read a stale array again. also computes single-run fci_05/fci_01/pc_05/
  pc_01 for fig 2 panel B / supp fig 8. 14 needs this script's run_one_ext/
  node_order_ext/context_idx live in the session, not just saved output.
- `04_scm_finalize.R` -- adapted from 18_finalize_scm_specification_v4.R.
  reads the bootstrap array through LATEST_bootstrap_run.csv, stops on an
  n_boot mismatch instead of silently using whatever's on disk.
- `05_scm_intervention_helpers.R` -- one shared base_edges (16 edges) +
  the deterministic mean-propagation helpers, used by 06 onward instead of
  each script keeping its own copy. cross-checks against 04's
  scm_edges_finalized.csv on load.
- `06/07/08` (from scripts 30/32/31) -- exact ATEs: single-node (6 nodes),
  8-node (+harm_future, politics as benchmark only), and all 63 combo
  targets, all x 16 orientation scenarios.
- `09` (from script 22) -- participant bootstrap CIs, the "1,000 participant
  bootstrap resamples" in the manuscript.
- `10` (from 02_full_orientation_enumeration_v4.R's fit half) -- CFI/TLI/
  RMSEA/SRMR/AIC/BIC across the 16 scenarios. doesn't touch ATEs, that's
  06/08's job.
- `11` (from script 11 parts B/C/D) -- politics x belief_concern interaction
  + politics-conditional shift intervention, plus an exact closed-form
  cross-check of the Monte Carlo shift means.
- `12` (from script 24) -- retained (N=870) vs attrited (N=1,117) SMD
  comparison, feeds supp S11's opening numbers.
- `13` (from script 28) -- IPW weights, overlap/ESS/truncation, weighted vs
  unweighted behavior equation + single-node ATEs. feeds the rest of S11.
- `14` (from script 23) -- mitig4/evacuate/move alternative outcomes: GGM
  weights, FCI adjacency, behavior refit, single-node ATEs. feeds S12
  (mitig4 only -- evacuate/move aren't cited with specific numbers yet).
- `15` -- confound/existence sensitivity diagnostic, prompted by the joint
  pag endpoint audit (`r_patches/29_joint_pag_edgetype_audit.R`): the 8 scm
  edges that were bidirected-plurality in the fci bootstrap at BOTH alpha
  levels, replaced one at a time (then all together) with a residual
  covariance instead of a directed path; plus reversing
  harm_future->harm_present on its own and dropping
  harm_future->trust_science entirely. reuses 05's base_edges/helpers as-is.
- `16` -- the single-run extended (9-node) pag doesn't exist anywhere else
  in the repo -- fci_05/fci_01 in 03 are main-network-only, that's figure
  S4. computes it at BOTH alpha levels with the same jci/contextVars args
  as 03's extended bootstrap (that's an exogeneity restriction -- wave-5
  behavior can't cause the 8 earlier-wave variables -- not a substantive
  "context variable" claim), prints the actual edge type (->, <-, <->,
  o->, <-o, o-o, or absent) for all 16 scm edges in that one discrete fit,
  same 7 categories r_patches/29's bootstrap audit uses so the two compare
  directly. if lv-ida is downloaded separately (see its header): reports
  listMags()'s mag-completion count per alpha AND whether it hit the 500
  cap (a capped count is a floor, not the true number), then runs lv.ida()
  for all 8 non-outcome nodes (matching fig 7's node set), rescaled to the
  same +0.5sd ate scale the rest of the pipeline reports on (not yet
  checked empirically that this rescaling is right). doesn't attempt
  pairs/triples -- single-node effects don't just add once you allow
  latent confounding, singles first.
- `17` -- formalizes the bidirected/directed/mixed classification rule
  (2026-09-08): an edge is bidirected_dominant only if it clears the
  existing 60% adjacency threshold AND bidirected is the plurality joint
  type at BOTH alpha levels (no new percentage invented); directionally_
  resolved needs the same directed type at both alphas; everything else is
  mixed_threshold_sensitive (that bucket isn't one thing -- covers genuine
  direction reversal, directed/bidirected flips, and thin-existence edges,
  see its header). cross-checks its own output against 15's hardcoded
  8-edge list and warns if they've drifted apart.

## reviewed, not migrated

checked each of these against main9.tex directly (grepped for method
names/numbers, didn't just trust the script's stated purpose):

- not cited anywhere in the manuscript currently: `14_sensitivity_rank_copula_fci.R`,
  `15_sensitivity_mixed_ci_fci.R`, `16_sensitivity_rcot_fci.R` (no "copula"/
  "RCoT"/"kernel" anywhere in main9.tex), `17_diagnostic_rcot_social_norms_behavior.R`,
  `19_systematic_interaction_screen.R` (screens 10 candidates, ms only
  reports the one 11 already covers), `21_orientation_crossalpha_table.R`
  (a proposed alt rule, not adopted -- its own header says so),
  `26_intervention_ci_followup.R` (extra pair CIs main9.tex doesn't cite).
- already covered elsewhere: `27_diagnose_hf_hp_edge.R` -- its finding
  (harm_future->harm_present threshold sensitivity) is already the 4th flip
  candidate in 05/06-10, nothing separate to migrate.
- CCI chain, decision already made: `11b/11c/11d_cci_*`, `10_cci_addition_and_figure4_pc_appendix.R` --
  established CCI's order-dependent orientation issue, supports keeping CCI
  at appendix level only. not number-producing for main text.
- pure diagnostic, not number-producing: `13_uncertain_cell_diagnostic.R`.
- dead: `02_full_orientation_enumeration_v3.R` (superseded by v4),
  `cache_fig4_fig9_objects.R` (interactive-session caching helper, not
  needed here).
- figure scripts, not touched beyond path fixes: `03_figure_style.R`,
  `04_figure1_heatmap.R`, `05_figures2_3_ggm_redesign.R`,
  `06_figure4_stability_redesign.R`, `07_figure5_scm_hierarchical_v3.R`,
  `08b_supp_figure_all_combo_interventions.R`, `09_supp_causal_graphs_redesign.R`,
  `20_supplementary_measurement_figures.R`.
- fig 7 -- keeping all 5 variants (`10_figure7_uncertainty_pub.R` plus
  `_hump`/`_jitter`/`_v2`/`_violin`). the unsuffixed one is what main9.tex's
  `\includegraphics` actually points to; the rest are alternate visual
  encodings worth keeping around, not dead code.

## fixed along the way

- `r_patches/08b_supp_figure_all_combo_interventions.R` and
  `r_patches/10_figure7_uncertainty_pub.R` were reading their input CSVs
  from a bare `tables/...` path; clean_pipeline writes those under
  `pipeline_outputs/tables/...`. fixed both read paths in place, nothing
  else touched. everything reading a bare `pipeline_outputs/...` path was
  already fine since OUTPUT_DIR literally is that string.
- `_hump`/`_jitter`/`_violin` were reading `tables/orientation_enumeration_ate.csv` --
  the superseded Monte Carlo file (old script 02 v4, n=20,000, the one with
  the spurious-negative-floor noise bug the deterministic chain replaced).
  each had an `if (!exists("full_results"))` guard, so sourcing one right
  after the correct script in the same session would silently paper over
  it -- run standalone they'd plot stale numbers. fixed all three to read
  `pipeline_outputs/tables/orientation_enumeration_ate_deterministic_8node.csv`
  (same schema, their 6-node filter just ignores the extra rows). `_v2`'s
  bare path was fixed too but it's still deliberately reading the 6-node
  deterministic file, not the 8-node one the main figure uses -- a real
  scope difference, worth a look, not silently resolved here.
- Table S9's climate_behavior shift numbers were stale (published .202/.205,
  current 16-edge model gives .205/.207, confirmed two independent ways --
  the Monte Carlo run and the exact closed-form cross-check in 11 agree to
  3 decimals). fixed directly in Overleaf, not here.
- checked whether any of the above path fixes moved a number that's
  actually in the manuscript: rounded/keyed diff on every affected file,
  zero real differences -- pure path corrections, no statistical content
  moved.

## reference runs (2026-09-08)

After finding three separate stale-artifact bugs from mixing outputs across
different pipeline states (the S9 interaction numbers, the "24 of 28 edges"
prose, 06's own hardcoded baseline comment), I wanted one clean
`01->14` run, frozen with its own provenance, before adding 15/16/17 on top
-- rather than trusting whichever pipeline_outputs files happen to be
sitting around from whatever session last touched them.

`run_all.R` does that: fresh session, 01 through 14 in order, then copies
every file that run actually produced (by mtime) into
`pipeline_outputs/reference_runs/<run_tag>/`, writes
`reference_run_manifest.txt` (timestamps, sessionInfo(), the finalized
16-edge file's path/mtime/md5, the baseline check table), and checks the 8
single-node baseline ATEs against the values already confirmed twice this
session (06 and 15 agreeing to 9 decimals): belief_concern=.2022,
harm_present=.2157, harm_future=.1399, weather_risk_prep=.0921,
politics=.1080, social_norms=.0481, trust_science=.0115, policy_support=.0077.

If all 8 match: that run is the frozen reference bundle -- don't rerun
01-14 again without a reason, and treat that `reference_runs/<run_tag>/`
folder (plus the live pipeline_outputs/ state, which is what it's a copy
of) as authoritative until Results/Table 3 get rewritten from it.

**Important:** run_all.R does NOT clear the R session when it finishes --
`df_extended`, `agg_ext`, `node_order_ext`, `context_idx`, etc. are all
still live. 15 and 16 both need those objects live (`stopifnot(exists(...))`
checks in 05 and 16), so run 15 and 16 in that *same* session, right after
run_all.R, rather than restarting R in between -- otherwise there's no
guarantee they're seeing the same fitted model/data run_all.R just froze.
17 is the exception -- it only reads pipeline_outputs csvs, no live session
objects needed, so it can run anytime after 04 (or after run_all.R) in a
fresh session too.

## still open

- the big one: whether the working scm should keep being a plain DAG built
  by hand-picking one orientation per edge, or move to a pag-consistent
  semi-Markovian representation (directed mag edges = observed-level
  causal/ancestral paths, bidirected mag edges = correlated disturbances /
  latent common causes), with lv-ida for single-node effects instead of a
  point-estimate ATE. 15 is the quick diagnostic version of this (manual
  covariance replacement on the 8 edges the joint audit flagged as
  bidirected-plurality at both alpha levels). 16 is the more principled
  version (an actual pag object + listMags()/lv.ida()), but it's brand new
  and unrun -- don't treat either one's output as settled until it's been
  run and read. not touching main9.tex until this is resolved.
- `run_all.R` is written -- not yet run. once it's run clean (all 8
  baseline ATEs matching), that frozen `reference_runs/<timestamp>/`
  bundle is what 15/16/17 should be read against, and what Table 3/Results
  gets rewritten from -- see "reference runs" below.
- resyncing N_BOOT to 2500 -- would need rerunning 03+04 and diffing the new
  scm_edges_finalized.csv against the current one before touching any
  manuscript text. not done, not urgent.
- `_v2`'s 6-node vs 8-node scope difference (see fig 7 note above).
