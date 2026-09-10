# clean_pipeline_alt_will/ -- side-by-side test: social_norms = cvcc4_will

Sara raised a concern about the injunctive framing of social_norms
(currently cvcc4_should, "people should permanently shift..."). This is a
full isolated copy of clean_pipeline/ testing cvcc4_will (descriptive norm,
"people will permanently shift...") instead, so it could actually be checked
without touching anything validated.

Same as clean_pipeline/00-14 except:

1. `01_data_prep.R` builds social_norms from cvcc4_will instead of
   cvcc4_should -- the only content change.
2. `00_config.R`'s OUTPUT_DIR is `pipeline_outputs_alt_will` instead of
   `pipeline_outputs`.
3. every internal source() call points here instead of clean_pipeline/.

writes only into pipeline_outputs_alt_will/, can't touch pipeline_outputs/
or anything validated against main9.tex.

## outcome (ran end to end)

cvcc4_will collapses the evidence for two of social_norms' three edges --
trust_science->social_norms drops from p_adjacent 1.00 to 0.019,
policy_support->social_norms from 0.79 to 0.072, both well under the 0.60
inclusion threshold. social_norms->climate_behavior stays solid. 04's audit
flags both (FLAG_WEAK_EXISTENCE) but doesn't drop them automatically -- its
edge list is fixed, the bootstrap only audits it -- so 06-14 ran through
fine on the unchanged topology. taking the audit's flag seriously though,
cvcc4_will would leave social_norms hanging off climate_behavior alone,
disconnected from the belief/trust/policy part of the model.

downstream: every ATE moved by something in the .001-.03 range, the top-
ranked intervention pair flipped (belief_concern+weather_risk_prep drops out
of the ranking entirely, harm_present+weather_risk_prep takes over), and
social_norms' attrition SMD crossed the .10 flag threshold (-.001 -> -.115).

raw correlation with climate_behavior also favors cvcc4_should (.377 vs
.248 for cvcc4_will; the two items correlate .448 with each other), and the
isolated network-edge check alone (.104 vs .103) looked neutral -- it took
the full refit to see the actual structural cost.

**decision: keeping cvcc4_should.** this folder stays around in case it's
worth revisiting.

## how to rerun it

fresh R session (don't mix objects from a clean_pipeline/ session -- names
like df_extended, network_ext, base_edges are reused and mixing sessions
would quietly combine both variants):

```
source("clean_pipeline_alt_will/01_data_prep.R")
source("clean_pipeline_alt_will/02_ggm.R")
source("clean_pipeline_alt_will/03_bootstrap_causal_discovery.R")   # slow, ~2000 pcalg fits
source("clean_pipeline_alt_will/04_scm_finalize.R")
source("clean_pipeline_alt_will/06_intervention_ates_singlenode.R")
source("clean_pipeline_alt_will/07_intervention_ates_8node.R")
source("clean_pipeline_alt_will/08_intervention_ates_combo.R")
source("clean_pipeline_alt_will/09_intervention_bootstrap_ci.R")   # slow, ~1000 lavaan refits
source("clean_pipeline_alt_will/10_orientation_enumeration_fit.R")
source("clean_pipeline_alt_will/11_interaction_moderation.R")
source("clean_pipeline_alt_will/12_wave5_attrition_check.R")
source("clean_pipeline_alt_will/13_ipw_attrition_sensitivity.R")
source("clean_pipeline_alt_will/14_behavior_outcome_sensitivity.R")  # needs 03's run_one_ext live, same as primary
source("clean_pipeline_alt_will/compare_to_primary.R")
```

05_scm_intervention_helpers.R (auto-sourced by 06 onward) hardcodes the
current 16-edge structure and only checks that 04's fresh output still
lists the same edges -- it won't stop just because an edge's evidence got
weak, only if an edge outright disappears or a new one shows up.
compare_to_primary.R prints the value-level diff against pipeline_outputs/
(rounded to 4 decimals, flags anything past rounding noise) -- edges first,
then every ATE/CI/interaction/sensitivity number.
