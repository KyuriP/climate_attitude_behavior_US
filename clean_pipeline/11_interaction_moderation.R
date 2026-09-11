# 11_interaction_moderation.R -- REMOVED 2026-09-11.
#
# This script formerly tested a belief_concern x politics interaction on
# policy_support (Parts A-D; reported in the main text as an "additional
# check" and in Supplementary Section S9). It has been removed entirely,
# per the 2026-09-11 decision to drop the interaction/moderation analysis
# from the paper -- see analysis_decisions_log.md Section 38 for the full
# rationale. In short:
#
#   - With belief_concern now upstream of politics (2026-09-11 edge
#     reversal, see Section 34/35), the original motivation for testing
#     whether the belief_concern -> policy_support effect varies by
#     political orientation no longer holds: politics is not exogenous to
#     belief_concern anymore, so the interaction's clean asymmetric
#     interpretation is gone.
#   - policy_support is a sink node in the working SCM (zero outgoing
#     edges) -- an interaction that lives only in policy_support's own
#     equation has no downstream route to climate_behavior, so it cannot
#     materially affect the paper's intervention conclusions either way.
#   - The analysis was a side question relative to the paper's main
#     sequence (conditional-dependence network -> causal discovery ->
#     working SCM -> intervention analysis -> directional uncertainty ->
#     latent-confounding sensitivity -> PAG/LV-IDA robustness -> feedback/
#     cyclic extension) and no longer earns its place in it.
#
# Nothing downstream depends on this script's outputs
# (interaction_coefficient.csv, interaction_shift_results_montecarlo.csv,
# interaction_shift_results_exact.csv) -- confirmed by grepping
# clean_pipeline/ and r_patches/ for references to those filenames and to
# 11_interaction_moderation.R itself; the only hits were this file and
# run_all.R's own (now permanent) exclusion comment.
#
# Permanently excluded from clean_pipeline/run_all.R's scripts_in_order
# (see that file). Left in place as a stub rather than deleted from the
# repo so the removal rationale stays attached to a real file at its old
# path; git-rm it instead if you'd rather it disappear entirely.
