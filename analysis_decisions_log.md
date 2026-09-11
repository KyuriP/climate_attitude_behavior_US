# Analysis decisions log — Study 1 (climate attitudes → action)

Purpose: a shareable record of a set of specific investigations carried out
during the final analysis/figure pass, organized as what was tried, what was
found, and what was decided. This complements the methodological narrative in
`study1-methodology-confirmed-story.md` rather than repeating it; the items
below are the ones not already covered there.

## 1. Completing the intervention-effect simulation (single nodes → all 41 targets)

**Why this came up.** The intervention barplot originally reported only two
hand-picked joint-intervention scenarios (harm_present+social_norms;
belief_concern+harm_present) out of every possible combination of the six
intervenable nodes, and even those two were never re-simulated across the
eight structural-uncertainty specifications the way the six single-node
effects were.

**What was tried.** Every non-empty subset of size 1-3 of the six intervenable
nodes was enumerated (6 singles + 15 pairs + 20 triples = 41 targets) and
simulated under each of the same 8 orientation specifications already being
fit for the single-node results -- no new model refits, just additional
`simulate_scm_generic()` calls against models that were already fitted.
Combined-target results were written to a separate table
(`orientation_enumeration_combo_ate.csv` /
`orientation_uncertainty_band_combos_full.csv`) so the existing single-node
table and figure were left untouched.

**What was found.** Ranking all 41 targets by their baseline (combo_0) effect
on climate_behavior:
- Best single node: belief_concern (+0.202 SD), narrowly ahead of harm_present
  (+0.192 SD).
- Best pair: belief_concern + weather_risk_prep (+0.253 SD, range 0.247-0.266
  across the 8 specifications).
- Best triple, and best target overall: harm_present + weather_risk_prep +
  social_norms (+0.311 SD, range 0.286-0.312).

As expected, larger combinations dominate the top of the list simply by
having more levers; the value of the full enumeration is being able to name
the actual best pair and best triple with a number, rather than reporting
effect sizes only for the two combinations that happened to be chosen by
hand.

**Decision.** Report the best single, best pair, and best triple (with their
structural-uncertainty ranges) in the manuscript's intervention table/figure,
replacing the two hand-picked combined scenarios.

## 2. Structural uncertainty in social_norms -> climate_behavior: a fit improvement that isn't usable

**Why this came up.** Among the three edges whose direction is not resolved by
the FCI bootstrap asymmetry statistic and is instead assigned by theory
(politics->belief_concern, policy_support->social_norms,
social_norms->climate_behavior), the orientation-enumeration procedure fits
all 2^3 = 8 sign combinations and compares them by AIC/BIC. That comparison
showed the specification that reverses social_norms->climate_behavior fits
decisively better than the working model: AIC and BIC each improve by
approximately 11 points (baseline AIC = -4790.01, BIC = -4670.80; reversed-edge
AIC = -4801.05, BIC = -4681.83), a larger margin than any of the other seven
alternative specifications produce.

**What was checked.** Before treating this as evidence for a genuine feedback
effect (climate_behavior -> social_norms rather than the reverse), the two
variables' construction in the data was checked directly. social_norms is
built as an average across Waves 1-4; climate_behavior is measured only at
Wave 5. That ordering makes the reversed direction temporally impossible as a
causal claim -- climate_behavior cannot cause a variable that was already
measured, on average, before it existed. A statistically better-fitting
reversed edge under a linear-Gaussian path model is therefore almost
certainly a model-comparison artifact (a reversed edge can improve global fit
for reasons that have nothing to do with the true causal direction), not
evidence of genuine bidirectional or reversed influence.

**Decision.** Retain social_norms -> climate_behavior at its original,
theory-assigned direction in the working SCM and in `flip_candidates`, rather
than removing or reversing it on the strength of the fit comparison. The
fit-based finding is noted here for the record, in case a reviewer raises the
same AIC/BIC comparison independently, but it does not change the model. This
is a separate question from the edge's *existence* uncertainty (~60-65%
conditional support; see Section 11 of `study1-methodology-confirmed-story.md`),
which remains an open, genuine finding -- this section concerns only its
*direction*, which is not in genuine doubt.

## 3. Encoding structural-orientation uncertainty in the merged intervention/uncertainty figure

**Why this came up.** The intervention-effect and orientation-uncertainty
figures were merged into one main-text figure (each candidate node's
baseline effect, plus how that effect moves across the seven non-baseline
structural specifications). The open design question was how to draw those
seven alternative-orientation values per node -- a small, complete enumeration
rather than a sample, which rules out some otherwise-natural chart choices.

**What was tried, in order:**
1. A beeswarm of the seven points per node -- rejected as visually noisy for
   only seven values per row.
2. A half-density "hump" silhouette with an unconstrained kernel-density
   bandwidth and padded tails -- the tails extended past each node's actual
   data range, and for nodes with a tight cluster plus one outlier (e.g.
   trust_science), the default bandwidth smeared the two into a single
   lopsided wedge rather than a legible two-part shape.
3. Fixed-width-bin histograms, at two different bin counts -- both read as
   visually busy relative to only seven data points.
4. A full/symmetric violin plot (`geom_violin`, which clips automatically to
   the real data range) and a plain jittered point cloud, built as two
   alternate comparison versions once it was clear the hump needed further
   work.
5. A revised hump with the density curve clamped to start and end exactly at
   each node's real minimum and maximum (removing the padding problem from
   attempt 2) -- this fixed the overshoot but the lopsided-wedge shape
   persisted for clustered-plus-outlier nodes, since the bandwidth itself was
   still too wide.
6. The bandwidth was then explicitly capped to one-sixth of each node's own
   value range (floored at one-twentieth, to avoid a degenerate near-zero
   bandwidth), on top of the range-clamping from attempt 5. This let a tight
   cluster and a distant outlier show as a distinct, visually separated bump
   rather than one long ramp.

**What was found.** Attempt 6 -- range-clamped density with a range-scaled
bandwidth cap -- was the only version that avoided both the overshoot problem
and the lopsided-wedge problem simultaneously, while still reading as a
smooth, legible shape rather than a bar chart.

**Decision.** The final figure uses the capped-bandwidth, range-clamped hump
as the main encoding, with a short tick mark under each node's row for each
of the seven real alternative-specification values (so the exact numbers
remain visible underneath the smoothed shape), and the baseline specification
marked separately with a diamond and its own "+X.XXX SD" label. The full
violin and plain-jitter versions were kept as alternate scripts
(`10_figure7_uncertainty_pub_violin.R`, `10_figure7_uncertainty_pub_jitter.R`)
in case a co-author prefers a different encoding when reviewing the figure.

## 4. Repository cleanup before sharing

Ahead of pushing `belief_network_US` to a repository, the following cleanup
was done:
- Removed script versions and duplicate files already confirmed superseded
  and not referenced anywhere in the pipeline (archived under
  `r_patches/_archive_2026-09-05/`, with a README documenting what replaced
  what), plus one additional confirmed-superseded script
  (`08_figure6_intervention_relabel.R`, replaced by the merged
  `10_figure7_uncertainty_pub.R`).
- Removed five outdated analysis notebooks (`.qmd`) and their HTML knits from
  the repository root, none of which were referenced by any current script --
  only `climate_analysis_avg_v2_altweather.qmd` is the live, authoritative
  analysis file.
- Removed a scratch folder of old script copies and superseded draft
  manuscript PDFs from earlier in the collaboration.
- Rewrote script comments and headers throughout `r_patches/` into
  professional third-person engineering documentation, preserving all
  technical content and rationale while removing conversational and
  collaborator-identifying phrasing.
- Moved the 13 CSV/RDS files that had accumulated at the repository root
  (bootstrap outputs, sensitivity-check comparisons, the finalized SCM edge
  list) into `pipeline_outputs/`, and updated every script and the `.qmd`
  that reads or writes them. One additional file, a Table 2 regeneration
  comparison that had already served its purpose (see Section 6 below), was
  deleted outright rather than moved, since nothing in the pipeline reads it.

## 4b. Why the social-norms construct is an injunctive-norm item, specifically

**Why this came up.** The `social_norms` node in the working SCM is a single
item (`cvcc4_should`), and the manuscript repeatedly flags that any
downstream social-norms finding should be read narrowly, as a possible
injunctive-norm pathway rather than evidence about social norms in general.
This choice is worth stating plainly for a collaborator asking why that
narrowing is necessary.

**What was checked.** Three candidate norm-related items were available in
the item pool: `cvcc4_should` ("people should permanently shift some
lower-carbon behaviors"), `cvcc4_will` ("people will..."), and
`cvcc4_personal` (personal intention). Descriptively, these three do not
behave as interchangeable measures of one construct: `cvcc4_should` has the
highest mean, `cvcc4_personal` is next, and `cvcc4_will` is lowest --
respondents are more willing to endorse what people *should* do than to say
they personally will do it, and least willing to say others actually will.

**What was found.** The three items map onto three different constructs, not
three noisy measures of the same one: `cvcc4_should` is an **injunctive
norm** (what people ought to do); `cvcc4_will` is a **descriptive
expectation** about others' likely behavior; `cvcc4_personal` is a
**personal intention**, not a social norm at all.

**Decision.** Retain `cvcc4_should` as the sole `social_norms` node, because
it is the clearest and most defensible standalone indicator of perceived
normative obligation among the three candidates, and exclude the other two
as measuring different constructs rather than as redundant alternatives.
This is also why the manuscript's discussion of `social_norms` is
consistently qualified as an injunctive-norm pathway specifically: the item
retained cannot speak to descriptive norms (what others actually do) or to
personal intention, only to what respondents believe people ought to do.
The item's COVID-specific framing ("permanently shift... following COVID")
is a separate, acknowledged limitation of this measure, independent of the
injunctive-vs-descriptive distinction above.

## 5. Auditing which edges are genuinely orientation-uncertain, and why "trust in science is fragile" didn't survive

**Why this came up.** The manuscript's headline uncertainty finding was that
trust in science's intervention effect was fragile to structural orientation:
across a five-edge, 32-combination enumeration, its effect ranged from .020
to .171 SD -- more than an eightfold swing -- driven mainly by whether
belief/concern -> trust science or the reverse was assumed. Regenerating
Table 2 against the current pooled FCI bootstrap array (a routine audit, not
originally aimed at this finding) surfaced something that called this
directly into question.

**What was found.** Of the five edges originally treated as requiring a
substantive (theory-assigned) direction, two -- belief/concern -> trust
science and present harm -> weather risk (worry) -- now show *positive*
bootstrap asymmetry under the current data (previously they were negative,
which is why they were assigned by theory in the first place). Both
therefore now agree with the direction already assumed in the model, which
means the causal-discovery evidence no longer actually leaves their
direction unresolved. Only three edges remain genuinely uncertain: politics
-> belief/concern, policy support -> social norms, and social norms ->
climate behavior.

**What was checked.** Re-running the orientation-sensitivity enumeration
restricted to only these three genuinely uncertain edges (2^3 = 8
combinations instead of 2^5 = 32) collapsed trust in science's range from
.020-.171 SD down to .0035-.0383 SD -- it is now one of the *more* stable
nodes in the model, not the least. Every other node's range was similarly
narrow, and model fit varied only slightly across all 8 specifications (CFI
.978-.981), confirming the earlier wide fit range (CFI .954-.976) was
specifically driven by specifications that reversed one of the two edges
that should not have been in the flip set to begin with.

**Decision.** Treat this as a correction, not a new finding to layer on top
of the old one: the three-edge, eight-combination enumeration is now the
model's sensitivity analysis, and the previously reported eightfold
trust-in-science range is retired as an artifact of including two edges that
current evidence no longer supports treating as uncertain. This reverses the
paper's stated headline conclusion (structural orientation uncertainty
materially changes which nodes look effective) to its opposite (once the
uncertain edges are correctly identified, this system's intervention
conclusions are comparatively robust to them) -- the abstract, Discussion,
and Conclusion were rewritten accordingly (Section 6).

## 6. Syncing the manuscript text to the corrected 16-edge pipeline

**Why this came up.** The live manuscript (`main2.tex`) had fallen out of
sync with the analysis pipeline: it still described the original 14-edge,
five-flippable-edge, 32-combination model, and Figure 5 (interventions) and
Figure 6 (uncertainty) were still two separate figures with a broken
reference to a deleted file, even though the pipeline had already merged
them into one (`fig7_uncertainty_pub.pdf`) and added two edges (16 total).

**What was done.** The Results section, Table 2 (edge-orientation evidence),
Table 3 (interventions), the merged intervention/uncertainty figure, the
abstract, Discussion, Conclusion, and the full Supplementary
structural-orientation-sensitivity section (methods, both tables, and prose)
were all rewritten against the real current numbers: 16 edges with their
current bootstrap avg-N/asymmetry values, the real fit statistics
(CFI=.979, TLI=.962, RMSEA=.086, SRMR=.039), the actual best intervention
pair and triple from the full 41-target search (Section 1 above) in place of
the two hand-picked combinations, and the corrected three-edge orientation
enumeration (Section 5 above) in place of the five-edge version. Every
change was checked against the real pipeline output tables directly (not
estimated), and the two edges whose classification changed since original
publication (belief/concern -> trust science; present harm -> weather risk)
are flagged with footnotes in Table 2 rather than silently reclassified.

**Verification.** Rather than only reviewing the LaTeX by eye, the full
manuscript was compiled end to end with `pdflatex` (two passes, to resolve
cross-references) after each round of edits. This caught real problems a
read-through would have missed: a genuinely broken figure reference (see
Section 7), and a set of six Supplement cross-references that had silently
been broken for some time under a labeling-convention mismatch
(`\ref{sec:ccipcappendix}` in the main text vs. `\label{supp:algorithms}`
in the Supplement, and five more like it) -- these are now aliased so both
names resolve, and the document compiles cleanly with zero errors and zero
undefined references.

**Decision.** Adopt this as a standing check: recompile `main2.tex` after
any substantive edit, rather than relying on a visual read-through, since
LaTeX cross-reference and file-path errors are otherwise silent until
someone happens to notice a "??" in a rendered PDF.

**Still open, deliberately not touched:** the Supplement's full standardized
path-coefficient table (`tab:pathcoef` / `tab:supp_pathcoef`) still has the
original 14-edge coefficients. Only 5 of the 16 real values were available
without a fresh pull from the live session, so the table was left as-is
rather than partially updated with a guess for the rest. A single line was
added to the `.qmd` (right after the existing `std_tr` computation) to
write the full current table to `pipeline_outputs/std_path_coefficients_16edge.csv`
the next time that chunk is run.

## 7. Figure 3/4: resolving a broken main-text reference and reviving two orphaned Supplement figures

**Why this came up.** `main2.tex`'s Figure 3 (`fig:causal`) pointed at
`figures/fig4_stability_matrix.pdf`, a file that does not exist -- this was
the one fatal error blocking full compilation of the manuscript.

**What was found.** The actual current main-text design (a continuous
dot-matrix: dot size = adjacency stability, dot fill = orientation
asymmetry, in `06_figure4_stability_redesign.R`) was already built,
FCI-only as required, and already rendered -- just under a different,
draft-stage filename (`fig4_stability_dotmatrix.pdf`). Separately, the
script that had originally been a candidate for the main-text figure
(`10_cci_addition_and_figure4_pc_appendix.R`) turned out not to be dead code
needing archival: its own "Figure 4 rebuild" section had already been
redirected to produce two Supplement-only figures instead (an FCI detail
view and the PC-stable matrix), both already rendered, but never actually
wired into the manuscript -- the relevant Supplement figure blocks were
still commented-out placeholders.

**Decision.** Rather than build anything new, three small connections were
made: Figure 3's image reference and caption were corrected to match the
real dot-matrix output and its actual encoding (the old caption still
described a retired categorical absent/directed/bidirected scheme); the two
commented-out Supplement figure placeholders were activated with the
already-rendered PDFs; and script 10's header comment, which still claimed
to feed the main text, was corrected to describe its actual current
Supplement-only role. No script logic changed -- this was a wiring and
documentation fix, verified by a full `pdflatex` recompile (Section 6).

## 8. Weather-risk terminology cleanup

**Why this came up.** A pass through the manuscript for consistent terminology
found that the six-word phrase "weather risk (worry)" appeared throughout the
Working SCM and intervention sections, while every earlier section (Introduction,
Methods, GGM Results) and the Supplement's own tables used plain "weather risk."

**What was found.** All 17 occurrences of the "(worry)" qualifier were confined
to the sections rewritten during this session's 16-edge/3-flippable-edge sync
(roughly the Working SCM description, Table 2's footnotes, the intervention
results, and the merged uncertainty figure's caption) -- it had been added there
for extra precision but never propagated to, or removed from, anywhere else.

**Decision.** Standardized on plain "weather risk" everywhere, matching the rest
of the manuscript and the Supplement's own tables. The construct's precise
definition -- a single item about worry/perceived personal risk, explicitly not
a weather-risk/preparedness composite -- remains documented once, where it
belongs: the Measurement section's "Construct-specific decisions" (S1/S2).

## 9. Supplementary measurement-reduction figures (item dendrogram, PCA scree plot, behavior-item correlation structure)

**Why this came up.** Wanted a proper, "neat" supplementary treatment
of the exploratory measurement-reduction stage -- something more than the couple
of summary numbers currently in S2.

**What was found.** The candidate-item hierarchical clustering, the PCA
eigen-decomposition, and the Wave-5 behavior-item correlation matrix are all
already computed in `climate_analysis_avg_v2_altweather.qmd` (Sections 3.5, 3.6,
3.8) -- they just were never exported as standalone figures or wired into the
manuscript, and the qmd's own EDA narrative (which items cluster together, which
of the six EFA factors each item loads on) was never carried over into S2's much
terser prose.

**Decision.** No new analysis was run. `r_patches/20_supplementary_measurement_figures.R`
was written to re-render the already-computed objects (`hc_ward`, `eig`, `cor_beh`)
as three standalone, print-sized PDFs (`figS_item_dendrogram.pdf`,
`figS_pca_scree.pdf`, `figS_behavior_corr.pdf`) once run (it assumes the
qmd has already been run through Section 3.8, the same convention every other
`r_patches` script uses). Separately, S2's prose was substantially expanded using
findings that were already sitting in the qmd's own EDA write-up but had not
reached the manuscript: which items group together in the Ward.D2 dendrogram
(future-harm cluster; present-harm+belief/concern cluster; a broader
belief-policy-response cluster; the norm-item triplet; a peripheral
weather-risk/risk-perception/public-worry cluster), and the substantive content
of each of the six EFA factors. The three new figure blocks are wired in as
commented-out placeholders -- captions and labels already written -- following
the same convention already used elsewhere in this Supplement for not-yet-rendered
figures; they need only be uncommented once the script above has been run.

## 10. Compile-check gap found: missing bibliography entries

**Why this came up.** Recompiling `main2.tex` after the Section 9 changes
surfaced two issues that the earlier "verified clean" compile (Section 6) had
not caught.

**What was found.** Both are pre-existing and unrelated to today's edits. First,
the verification environment does not have `main2.tex` and its `figures/` folder
in the same place (they live in two different connected folders on my
machine), so pdflatex fell back to draft placeholders for 8 figures and reported
them as errors; copying `figures/` alongside `main2.tex` for a one-off check
confirmed the manuscript itself compiles cleanly (32 pages, exit 0) once the
images are reachable -- this was a gap in the check, not in the manuscript.
Second, and genuinely unresolved: eight citation keys used in-text have no
matching entry in `references.bib` -- `pearl2009causality` and
`peters2017elements` (Methods, Pearl's $do$-operator), and `drews2016public`,
`fielding2026intention`, `bergquist2022meta`, `helferich2023norms`,
`lind2024comparing`, `zia2024machine` (Introduction).

**Decision.** Left these as-is rather than fabricate bibliography entries --
they need the actual dataset source details (same treatment as the existing
dataset-citation TODO in Methods). Flagged directly rather than guessed at.

## 11. Figure polish from a live review pass (legend, Figure 3 panel, Figure 5 table, orphaned combo figure)

**Why this came up.** Reviewing the rendered Figure 7 (intervention/uncertainty) and
Figure 3 (FCI stability) turned up four separate small gaps.

**What was found and decided, one by one:**

- *Figure 7's legend was wrong.* The "Alternative orientation" key showed an open
  circle, but no circle is ever actually drawn on the chart -- those values are
  encoded only by the density hump and the short tick marks under each row. The
  circle came from a leftover ggplot trick (an invisible, zero-size point layer
  added only to force a legend entry) whose legend glyph was never updated when
  the encoding moved from dots to hump+ticks. Fixed in both
  `10_figure7_uncertainty_pub.R` (the deployed script) and
  `10_figure7_uncertainty_pub_hump.R` by changing the legend shape for
  "Alternative orientation" to a vertical-tick glyph (pch 124, "|"), matching what
  is actually drawn.

- *Figure 3 was missing its second panel.* `06_figure4_stability_redesign.R`'s own
  header comment says the dot-matrix panel was always meant to be paired with a
  single-run FCI example graph via LaTeX minipages (the same two-panel technique
  already used for Figure 2's GGM main/extended panels) -- but when the broken
  image-path error was fixed earlier this session (Section 7), only the dot-matrix
  panel was wired in, dropping the example-graph panel that this design always
  called for. Restored the two-panel minipage layout, using the existing
`figures/fci_pag_alpha05_example.png` as panel B. Already flagged
  that this particular rendering is rough (arrowhead placement, title overlap)
  and plan to redraw it myself; no change needed on the LaTeX side once I
  do, since it will replace the same filename.

- *Figure 5 (SCM) had no coefficient table.* Added a compact standardized-path-
  coefficient table (edge, $\beta$) in a minipage beside the diagram, using the
  real, already-audited `BETA_TR` values from `07_figure5_scm_hierarchical_v3.R`
  itself (the same table the script uses to size the diagram's own edge widths) --
  no new computation, no rerun required. The full standard errors and
  significance tests for these same 16 edges remain in
  Table~\ref{tab:supp_pathcoef}, which is still waiting on the `.qmd`'s
  `std_tr` export (Section 6's still-open item).

- *A 41-target combo figure existed but was never wired in.* `figures/figS_all_combo_interventions.pdf`
  (from `08b_supp_figure_all_combo_interventions.R`) was already fully rendered
  -- showing all 35 pair/triple combinations not already covered by the
  main-text best-pair/best-triple callouts -- but was never referenced anywhere
  in `main2.tex`. Added it to Supplementary S7 with a short paragraph explaining
  its relationship to the main-text summary.

**Decision.** All four are wiring/consistency fixes to designs and results that
already existed, not new analysis. Recompiled (33 pages, exit 0, no new errors)
to confirm.

## 12. Real 16-edge path coefficients landed; Supplement path-coefficient table completed

**Why this came up.** This was the one item still blocking full completion:
full standardized-path-coefficient table (with SEs and significance, not just point
estimates) needed the `.qmd`'s `std_tr` object, which required rerunning the
notebook with the `write.csv` line added earlier this session.

**What was found.** `pipeline_outputs/std_path_coefficients_16edge.csv` now exists
(16 rows: From, To, Beta, SE, Z, p, sig). Cross-checked against the `BETA_TR` table
already hard-coded in `07_figure5_scm_hierarchical_v3.R` -- all 16 point estimates
match exactly, confirming both sources agree. This also caught a real, small error
in the main text: the intervention-implications paragraph attributed $\beta=.113$ to
the social-norms $\rightarrow$ climate-behavior path, but that value actually belongs
to belief/concern $\rightarrow$ weather risk; the social-norms path is $\beta=.114$.

**Decision.** Replaced the stale 14-edge `tab:supp_pathcoef` table with the real
16-edge version (adding SE and significance-star columns, with a footnote flagging
the two edges added after the original discovery-stage backbone). Updated the two
structural equations in S6 that were still missing the corresponding predictor terms
($T$'s equation was missing $\beta_{FT}F$; $W$'s equation was missing $\beta_{BW}B$).
Updated the S6 model-fit paragraph from the old 14-edge fit statistics to the real
16-edge ones already used in the main text (CFI $=.979$, TLI $=.962$, RMSEA $=.086$
[.074, .100], SRMR $=.039$), and corrected "14 selected relationships" to "16."
Fixed the $\beta=.113$/$.114$ mix-up in the main text. Also confirmed
`pipeline_outputs/pca_eigenvalues.csv` (from rerunning `20_supplementary_measurement_figures.R`)
matches the eigenvalues already reported in text (14.95, 1.73) -- no discrepancy.
Activated the three Supplement figure placeholders from Section 9 now that the
corresponding PDFs exist (`figS_item_dendrogram.pdf`, `figS_pca_scree.pdf`,
`figS_behavior_corr.pdf`). Recompiled with figures reachable: 35 pages, exit 0, no
new errors -- only the same 8 pre-existing missing-bibliography-entry warnings
remain (Section 10).

## 13. Responding to the 18-point methodological review — Section 1 (SCM construction rule) and independent text fixes

**Why this came up.** A detailed co-author review of the full manuscript identified 5
[BLOCKER]-level issues to resolve before further prose/figure polish, 9 [HIGH] issues,
and 4 [POLISH] issues. Per the reviewer's own explicit priority, work started on the
blockers rather than the polish items.

**Blocker 1 (make the 16-edge SCM construction rule explicit and reproducible).**
The rule was never actually undocumented — it was already decided and locked in
`study1-methodology-confirmed-story.md` (existence: pooled FCI bootstrap adjacency
proportion ≥ 60%, fixed in advance, mechanical, no exceptions; orientation: sign of
the pooled bootstrap asymmetry statistic, the same one already reported in Table 2) —
but `main2.tex` never stated it. Both the Methods subsection
(`sec:workingscmmethods`) and the Results subsection (`sec:workingscm`) were rewritten
to state the two-part rule explicitly and to describe the final graph as one unified
16-edge working SCM, replacing the "14 original edges + 2 added later" chronology.
The Table 2 footnote marking the two most recently qualified edges was reworded from
"added after the 14-edge specification... not part of the originally published model"
to a statement of their actual existence-proportion margin (65.7% and 94.6%,
respectively, versus ≥78% for every other edge) — factual and reproducible rather than
narrating debugging history.

**Independent factual fixes applied at the same time** (safe regardless of how
Blocker 2 — the orientation-criterion dispute — is eventually resolved; see the
the open question recorded further down in this entry
- Fixed a stale figure cross-reference in S2: prose pointed to "Supplementary Figure
  S2" for the full item-level correlation heatmap, but that heatmap was never
  actually rendered/activated (still a commented placeholder, `fig:supp_fullcor`),
  and the newly-added dendrogram/scree figures shifted the real numbering. Reworded
  to point to the actual rendered dendrogram figure via `\ref{}` instead of a
  hand-typed number.
- Fixed two mechanism-description errors in the intervention-implications paragraph
  (Results, ~line 319) that no longer matched the 16-edge graph: belief/concern's
  mechanism chain was described as only belief→future_harm→present_harm→behavior,
  omitting the direct belief→present_harm and belief→weather_risk edges also present
  in the 16-edge model; and the belief/concern + weather-risk joint-intervention
  result was described as reaching behavior through "largely separate paths," which
  contradicts the explicit belief→weather_risk edge (a joint hard intervention on
  both nodes actually severs that edge, which is the correct explanation for why the
  joint effect is modestly sub-additive rather than evidence of independence).
- Fixed a self-contradiction around trust-science's leverage claim (~line 363): one
  paragraph said the belief→trust-science edge is now bootstrap-resolved and
  excluded from the orientation-uncertainty enumeration, while a later paragraph said
  the same negligible-leverage claim "depends primarily on a single weakly identified
  edge" — leftover phrasing from the earlier five-edge version of the analysis.
  Reworded to state clearly that this was true previously, not now.
- Added one explicit limitation paragraph (Discussion) stating that the intervention
  simulations condition on the completed observed-variable SCM, and that the
  orientation-sensitivity analysis varies specific unresolved edge directions but
  does not and cannot enumerate every possible latent-confounding structure
  consistent with the FCI skeleton (reviewer item #7).

**Verified.** Full scratch-directory compile (`pdflatex → bibtex → pdflatex →
pdflatex`) succeeded, exit 0 at every step, 0 fatal errors, 35 pages.

**Still open / not yet done, per the reviewer's own explicit priority (blockers before
polish):**
- Blocker 2 (orientation-classification criterion) — genuine tension between the
  already-locked sign-based rule and the reviewer's request for a CI-based or
  cross-α-consistency criterion; needs a team decision, not a unilateral rewrite.
- Blocker 3 (bootstrap-participants CIs for intervention effects) — not started.
- Blocker 4 (4-item mitigation-only behavior-outcome sensitivity) — not started.
- Blocker 5 (Wave-5 attrition/selection analysis, N=870 of 1,987) — not started.
- All HIGH and POLISH items deferred per reviewer's stated priority.

## 14. Supplementary single-run causal-graph figure (fig_supp_causal_redesign.pdf) — illegible arrowheads

**Why this came up.** Flagged the rendered supplementary 2x2 single-run
causal-graph figure (FCI/PC-stable at α=.05/.01) as illegible — arrowheads not
visible on the edges. Staged the current PDF into the cloud workspace and viewed it
directly to confirm: the directional marks render as small dark ticks/notches, not
recognizable triangular arrowheads.

**Diagnosis.** In `r_patches/09_supp_causal_graphs_redesign.R`, the directional
portion of each edge was drawn as only the innermost `arrow_fraction = .17` of the
already-trimmed line, with a `grid::arrow(length = unit(2.8, "mm"))` closed
arrowhead. At that fraction, the colored/arrow segment is too short relative to its
own arrowhead size for the triangular head to render as a distinct, legible shape —
it collapses into a short dark smudge. Also confirmed, separately, that the rendered
PDF still shows "Weather risk (worry)" node labels, while the shared label map in
`03_figure_style.R` already has the corrected "Weather risk" (no qualifier) — this
PDF simply predates that terminology fix, same stale-render pattern as earlier
figures this session.

**Fix applied** (not yet rerun): increased `arrow_fraction` from .17 to
.42 (arrow segment now covers a visibly larger, unambiguous portion of the edge);
increased arrowhead length from 2.8mm to 3.6mm and angle to 22° for a fuller,
clearer triangle; increased arrow-segment linewidth from .95 to 1.15 and added
`lineend="round", linejoin="mitre"` for cleaner rendering. A rerun of script 09
will also pick up the already-correct "Weather risk" label.

**Decision.** Rerun script 09 to regenerate `figures/fig_supp_causal_redesign.pdf`
with both fixes; re-inspect the rendered PDF before considering this closed.

## 15. Team decisions on review items #2-#6, and four new scripts written

**Context.** Following the 18-point review, made explicit calls
on the remaining open items rather than leaving them for unilateral resolution:

- **#2 (orientation criterion) -- changing before publication.** Rejected a
  naive CI on the pooled asymmetry statistic (would mostly quantify Monte Carlo
  precision of the bootstrap frequency estimate, not real sampling uncertainty,
  short of an expensive nested bootstrap not judged worth it). Final rule:
  keep the existing >=60% pooled adjacency (existence) criterion unchanged;
  classify a direction as data-supported only if the dominant non-absence FCI
  endpoint pattern agrees on the same direction SEPARATELY at alpha=.05 and
  alpha=.01 (not just in the pooled combination). Disagreement between the two
  alphas, or a persistently circle/circle (bidirected/partially-oriented)
  pattern at either alpha, moves the edge into the orientation-sensitivity set.
  Explicitly willing to let this change the 3-edge/8-specification story.
  **Immediate next action (not yet decided further): generate the edge-by-edge
  cross-alpha evidence table first, review it, THEN decide the final unresolved
  set and only then rerun the orientation enumeration.**
- **#3 (intervention CIs) -- participant bootstrap, baseline orientation only,
  not nested inside the 8 orientation specifications** (kept as two separate
  uncertainties: sampling/parameter vs. structural). N=1,000 resamples, refit
  the baseline 16-edge SCM each time, recompute all 41 intervention targets.
  Report 95% percentile CIs for the six singles + the current best pair/triple,
  plus how often each pair/triple ranks first across resamples (more honest
  than a CI on a winner selected once). Deterministic mean propagation through
  the fitted structural equations instead of Monte Carlo inside the bootstrap
  loop (linear-recursive model -> no simulation noise, ~1000x fewer draws
  needed); the existing 20,000-draw Monte Carlo simulation stays as an
  independent validation check, untouched.
- **#4 (behavior-outcome sensitivity) -- targeted, not a second full
  pipeline.** Four-item mitigation/engagement outcome from food/travel/
  activism/discussion (`beh_meat`, `beh_travel`, `beh_activ`, `beh_discuss`).
  Emergency preparedness (`beh_evacuate`) and relocation (`beh_move`) analyzed
  SEPARATELY, not pre-averaged into a two-item "adaptation" composite (that
  can be added later if the data support it). Rerun, for each of the three
  alternative outcomes: extended GGM edge weights into the outcome; FCI/PC
  bootstrap existence proportions for the 8 attitude-to-behavior adjacencies;
  the behavior equation of the working SCM refit on the alternative outcome;
  the six single-node ATEs on the alternative outcome (full 16-edge SCM
  refit). Key question: does weather_risk_prep survive once the two most
  weather-proximal behaviors are removed from the outcome.
- **#5 (Wave-5 attrition) -- codebook checked directly, not assumed.** The
  "Core Wave 1_Codebook.pdf" attached to the project turned out to describe a
  DIFFERENT (European multi-country) survey and does not apply to this
  dataset -- confirmed by cross-checking item wording against manuscript text
  (this dataset's items explicitly reference "the United States," "most
  Americans," etc.). The actual demographic variable names were confirmed
  directly against this dataset's own `data_henry/codebook.parquet` and raw
  parquet file: `dem_age`, `dem_educ` (1-6 ordinal), `dem_income` (1-6
  ordinal), `dem_male` (1/0/77-self-described), `dem_race_cat` (7 categories),
  `dem_urban_cat` (rural/suburban/urban) -- all non-missing at every wave for
  every participant, one stable value per participant. Compare N=870
  (retained) vs. N=1,117 (attrited) on the 8 Waves-1-4 nodes and these
  demographics using standardized mean/proportion differences (not p-values,
  per team decision -- N~2,000 makes trivial differences "significant"). No
  IPW built preemptively; only consider it if |SMD| exceeds .10 broadly.
- **#6 (non-Gaussian CI robustness)**: use the existing rank-copula/mixed-CI/
  RCoT sensitivity outputs already on record rather than opening a new
  methodological branch; validate against final data objects, one Supplement
  table + short paragraph. Not yet started this round.
- **After #2-#5**: freeze analysis, then do the publication cleanup (remove
  version-history language, finish figure redesigns, Table 2 simplification,
  RMSEA language, shift-intervention sensitivity, reporting/ethics/authorship,
  final prose pass) -- deferred per team's stated priority, unchanged from the
  original review response.

**Four new scripts written this round** (none yet run -- all still need to
run in the live R session):

- `r_patches/21_orientation_crossalpha_table.R` -- the immediate next action
  above. Requires one addition already made to the .qmd (a `saveRDS()` of the
  per-alpha `mark_props_ext` object, previously only pooled before saving);
  re-render Section 7.4, then run this script. Produces
  `pipeline_outputs/orientation_crossalpha_table.csv`, a proposed (not yet
  adopted) cross-alpha classification for review -- does not touch the
  existing skeleton or rerun the enumeration.
- `r_patches/22_intervention_bootstrap_ci.R` -- Blocker 3. Self-contained;
  requires `df_extended` in the session. ~1,000 lavaan refits, deterministic
  propagation (no Monte Carlo). Writes
  `pipeline_outputs/intervention_bootstrap_ci.csv` and pair/triple rank-1
  frequency tables.
- `r_patches/23_behavior_outcome_sensitivity.R` -- Blocker 4. Requires
  `df_behavior_w5` (Section 4.3) and the Section 7.4 bootstrap machinery
  (`run_one_ext`, `context_idx`, `node_order_ext`, `alphas`, `n_boot`) already
  in the session; reuses that machinery unmodified rather than reimplementing
  it. This is the most expensive of the four (three outcome variants x the
  full FCI/PC bootstrap). Writes four
  `pipeline_outputs/behavior_sensitivity_*.csv` files.
- `r_patches/24_wave5_attrition_check.R` -- Blocker 5. Requires `df_main`,
  `df_extended`, `df_long` in the session. Writes three
  `pipeline_outputs/attrition_*.csv` files. No IPW implemented.

**Not yet done**: none of these four scripts have been run; all output
described above is prospective (the exact CSV columns each script produces),
not yet real numbers. Per my own stated sequencing, #21 (the cross-alpha
table) should be reviewed before deciding the final unresolved-edge set, but
#22-24 do not depend on that decision and can run in parallel.

## 16. Script 24 syntax bug (caught before running) + script 22 parallelized

My editor flagged a real syntax error in `24_wave5_attrition_check.R`
before running it: the summary message used an `if/else` block as one
argument to `cat()`, with each branch containing two comma-separated string
literals -- invalid, since commas are not statement separators inside `{ }`
blocks (only inside a function call's argument list). Fixed by building the
message with `paste0()` first and passing the single resulting string to
`cat()`. Checked all four new scripts (21-24) for the same anti-pattern
(grepped every `if (...) {` / `} else {` block) and for brace/paren/bracket
balance (no R interpreter available on this bridge, so this and manual review
are the available checks) -- only script 24 had the bug; the rest were clean.

Separately, `22_intervention_bootstrap_ci.R`'s 1,000-resample loop was
originally sequential; parallelized it with the same `furrr`/`multisession`
pattern already used in script 23 and in the .qmd's own bootstrap, since
`furrr` is already a project dependency. Each resample sets its own seed
(`2026 + b`) inside the worker function itself, so results are identical
regardless of execution order or worker count -- parallelizing does not
change what the script computes, only how long it takes.

## 17. Real results from scripts 21-24 (first run)

**Script 21 (cross-alpha orientation table).** The literal "clean single-headed
dominant pattern at both alpha" operationalization flags 14 of 16 edges as
unresolved -- this reproduces the exact false-alarm problem already documented
in `study1-methodology-confirmed-story.md` Section 6 (most real FCI-resolved
edges show meaningful arrowhead mass at BOTH endpoints because FCI tolerates
latent confounding; a dominance-based rule mistakes this for non-resolution).
A weaker, more usable reading of the same table -- same SIGN of asymmetry at
both alpha=.05 and alpha=.01 -- gives a far more actionable result: only 2 of
16 edges disagree in sign across alpha: social_norms->climate_behavior (both
~0, already known as existence- not direction-uncertain) and, newly,
harm_future->harm_present (asymmetry -.108 at alpha=.05 vs. +.261 at
alpha=.01), despite this edge's strong pooled asymmetry (+.526) and its status
as one of the two structurally strongest edges in the model. Still an open
decision on which operationalization to actually adopt, and the
harm_future-harm_present sign disagreement flagged as worth investigating
regardless of that decision, given how central this edge is to the model's
main narrative ("future harm -> present harm" is one of the two edges the
Results section calls out as the strongest fitted paths).

**Script 22 (intervention bootstrap CIs).** Two substantive findings:
(1) belief_concern (.206, 95% CI [.182,.227]) and harm_present (.195, CI
[.159,.230]) have heavily overlapping CIs -- confirms reviewer item #10's
concern that the "belief/concern narrowly ahead of present harm" ranking claim
is not statistically well-supported at this sample size. (2) The "best pair"
claim is NOT robust: across 1,000 resamples, harm_present+social_norms wins
38.3%, belief_concern+weather_risk_prep wins 37.8%, harm_present+weather_
risk_prep wins 23.9% -- essentially a three-way statistical tie, not the clear
single winner currently reported in the manuscript. By contrast, the "best
triple" claim IS robust: harm_present+weather_risk_prep+social_norms wins
100% of resamples.

**Script 23 (behavior-outcome sensitivity).** Addresses the key question
directly: weather_risk_prep's connection to behavior is not purely an
artifact of the two weather-proximal outcome items. For the 4-item mitigation/
engagement-only outcome (mitig4), weather_risk_prep remains a significant
predictor (SCM beta=.135, p=1e-4; ATE=.068), though smaller than harm_present
(beta=.348, ATE=.204) and below the skeleton's own 60% FCI-existence
threshold (38.5%) on its own in that outcome's bootstrap. For the two
adaptation items analyzed separately: weather_risk_prep dominates `evacuate`
(as expected -- existence 100%, beta=.263), while harm_present actually edges
out weather_risk_prep for `move`/relocation (beta=.152 vs. .119) -- a mildly
counterintuitive but plausible result (general harm perception, not
weather-specific worry, drives relocation intent more).

**Script 24 (Wave-5 attrition).** The 8 Waves-1-4 attitude nodes are all well
below |SMD|=.10 except weather_risk_prep (-.125, borderline). Demographics
tell a different story: age has SMD=.29 (retained participants notably older,
mean 55.4 vs. 50.9), income SMD=.19, and urbanicity (suburban +.15, rural
education is borderline (.10). Per my ownstated decision rule ("if nearly everything is |SMD|<.10, report and stop; if
clear differences exist, consider IPW"), this does NOT clear the "nearly
everything is small" bar -- age in particular is a clear, non-trivial
difference. Worth noting that an IPW sensitivity for the behavior-facing
analyses looks warranted by my own stated criterion, even though the
attitude-node covariates most relevant to the SCM itself are largely fine.

All four CSVs are in `pipeline_outputs/`; full column detail there rather than
repeated here.

## 18. Manuscript text lock (precision-edited wording) + scripts 26-27 for the harm_future/harm_present edge

**Why this came up.** Given the real results from scripts 21-24 (Section 17),
**Why this came up.** Given the real results from scripts 21-24 (Section 17),
made two explicit decisions: lock in manuscript wording for the
bootstrap-CI findings now, and diagnose the `harm_future -> harm_present`
edge rigorously before deciding whether to expand the 3-edge structural-
uncertainty set to 4 edges. Refined the wording, rejecting "statistically
indistinguishable" (too strong
without a difference-CI) and "three-way statistical tie" (too strong without
pairwise-difference testing), in favor of the exact replacement language
below.

**Manuscript edit (`main2.tex`, Results/intervention section).** Locked in,
using this exact wording:
- Belief/concern vs. present harm: "Belief/concern and present harm produced
  similarly large single-node effects, with substantially overlapping
  bootstrap 95% confidence intervals (.182-.227 and .159-.230 SD,
  respectively, from 1,000 participant bootstrap resamples of the baseline
  working SCM)."
- Best pair: removed the "largest joint effect" / unique-best-pair framing.
  Replaced with: "no single pair was consistently best across bootstrap
  resamples: harm present + social norms ranked first in 38.3% of 1,000
  participant bootstrap resamples, belief/concern + weather risk in 37.8%,
  and harm present + weather risk in 23.9%. Belief/concern + weather risk's
  point estimate (Delta Y=.253) is modestly below the sum of the two nodes'
  individual effects (.202+.087=.289), consistent with partial pathway
  overlap through the working SCM's direct belief/concern -> weather-risk
  edge rather than independence."
- Best triple: kept (it IS robust -- 100% of resamples), added: "...the only
  top-tier target whose ranking was not sensitive to sampling uncertainty:
  this triple ranked first in all 1,000 bootstrap resamples."

Full scratch-directory recompile (figures pulled from
`belief_network_US/figures/`, which is where the real 13 referenced PDFs/PNGs
actually live -- not the Downloads `figures/` subfolders, which are stale
partial copies from earlier passes): `pdflatex -> bibtex -> pdflatex ->
pdflatex`, all exit 0, 35 pages, only the same 8 pre-existing
missing-bib-entry warnings (drews2016public, fielding2026intention,
bergquist2022meta, helferich2023norms, pearl2009causality,
peters2017elements, lind2024comparing, zia2024machine -- unchanged from prior
passes, still need to double check). Extracted the three
rendered paragraphs via `pdftotext` and confirmed the wording matches exactly
with no LaTeX artifacts.

**Script 26 (`r_patches/26_intervention_ci_followup.R`, NEW).** Optionally
computing the bootstrap CI on the difference
Delta Y_belief - Delta Y_harm if a stronger claim was wanted later, plus CIs
for the two additional near-tied pairs. Reuses the already-saved
`pipeline_outputs/intervention_bootstrap_ate_matrix.rds` from script 22 --
NO rerun of the 1,000-resample bootstrap needed. Computes: (1) point
estimates for `harm_present+social_norms` and `harm_present+weather_risk_prep`
via one deterministic refit on the real sample; (2) 95% CIs for all three
tied pairs from the saved matrix; (3) the paired bootstrap CI on
belief_concern - harm_present via row-wise subtraction within each of the
1,000 resamples (preserves the resample-level correlation between the two
targets, which a naive independent-CI subtraction would not). Writes
`pipeline_outputs/intervention_bootstrap_tied_pairs_ci.csv` and
`intervention_bootstrap_belief_vs_harm_diff_ci.csv`. Balance-checked and
namespace-audited. **Not yet run.**

**Script 27 (`r_patches/27_diagnose_hf_hp_edge.R`, NEW).** Implements
my own 5-step diagnostic protocol for `harm_future -> harm_present`
(beta=.636, the single largest coefficient in the model, flagged by script 21
as the one edge -- besides the known social_norms -> climate_behavior case --
whose FCI bootstrap asymmetry changes sign across alpha=.05 vs .01):
1. Endpoint mark decomposition at each alpha separately (adjacency,
   arrowhead/tail/circle proportions at both endpoints, resulting asymmetry)
   -- from the already-saved per-alpha object, no new computation.
2. Seed sensitivity: 8 independent replicates x 1,000 resamples each at both
   alpha levels (16,000 FCI+PC fits total), parallelized via
   `furrr`/`multisession`, extracting only the HF-HP cell's asymmetry each
   time, to check whether the sign pattern (negative at .05, positive at .01)
   recurs or is a single-realization artifact. This is the expensive step.
3. Single-run FCI endpoint marks at both thresholds, from the already-fit
   `fci_ext_05`/`fci_ext_01` objects.
4. Acyclicity check: `igraph::is_dag()` on the baseline 16-edge SCM vs. a
   version with only harm_future/harm_present reversed.
5. Substantive consequence: fits baseline vs. reversed SCM, compares
   CFI/TLI/RMSEA/SRMR/AIC/BIC and the six single-node ATEs.
Writes 4 CSVs: `hf_hp_endpoint_decomposition.csv`, `hf_hp_seed_sensitivity.csv`,
`hf_hp_reversal_fit_comparison.csv`, `hf_hp_reversal_ate_comparison.csv`.
Balance-checked and namespace-audited. **Not yet run.**

**Decision rule for after script 27 runs (my own criterion, recorded
here so it isn't re-litigated later):** if the sign reversal is stable across
seeds and both thresholds retain the adjacency strongly, treat
harm_future-harm_present as unresolved under the cross-alpha robustness
criterion and move to a 4-edge / 16-specification structural-uncertainty
enumeration. If the alpha=.05 negative sign disappears across seeds or hovers
near zero while alpha=.01 stays consistently positive, do NOT promote the
edge automatically -- describe it instead as threshold-sensitive orientation
evidence, and do not let the final rule be driven by one bootstrap
realization.

**Still open:** whether the age-driven attrition imbalance (SMD=.29, Section
17 / script 24) warrants an IPW sensitivity for the behavior-facing analyses
-- not yet decided.

## 19. Script 27 real results: harm_future -> harm_present diagnostic

**Step 1 (per-alpha endpoint decomposition).** Adjacency is 100% at both
alpha=.05 and alpha=.01 -- never in question. Endpoint marks: at alpha=.05,
arrowhead-at-HF=.573 vs. arrowhead-at-HP=.466 (asymmetry_HP_minus_HF=-.108,
mild lean toward HP->HF); at alpha=.01, arrowhead-at-HF=.428 vs.
arrowhead-at-HP=.689 (asymmetry=+.261, clearer lean toward the originally
claimed HF->HP). Confirms script 21's flag exactly.

**Step 2 (seed sensitivity, 8 replicates x 1,000 resamples x 2 alpha =
16,000 FCI+PC fits).** The sign pattern is NOT a single-realization
artifact: asymmetry at alpha=.05 is negative in all 8 replicates (range
-.039 to -.127, mean approx -.068) and at alpha=.01 is positive in all 8
(range .210 to .289, mean approx .248). Adjacency stayed at 100% in every
replicate at both alpha. The sign never crosses zero at either threshold
across any replicate.

**Step 5 (substantive consequence).** Fit statistics for baseline (HF->HP)
vs. reversed (HP->HF) 16-edge SCMs are IDENTICAL to reported precision (CFI
.979, TLI .963, RMSEA .086, SRMR .039, AIC -4790, BIC -4670.8 for both) --
consistent with the two specifications being Markov-equivalent (HF and HP
share belief_concern as a common parent; reversing the one edge between them
does not change the implied covariance structure). Of the six policy-relevant
single-node ATEs: belief_concern, weather_risk_prep, social_norms,
trust_science, and policy_support are UNCHANGED (diff=0 in every case);
harm_present's own ATE moves by only +.0046 SD (.1951 -> .1997, from the new
indirect harm_present->harm_future->trust_science->social_norms path created
under the reversal) -- negligible relative to its .159-.230 bootstrap CI.

Applying my own decision rule: if the
reversal IS stable across seeds and both thresholds retain the adjacency
strongly (100%) -- this is the branch of that rule that says treat
harm_future-harm_present as unresolved under the cross-alpha robustness
criterion and move to a 4-edge / 16-specification structural-uncertainty
enumeration. However, Step 5 shows this uncertainty is close to
inconsequential for the paper's actual conclusions -- the reversal looks
like genuine Markov-equivalence, not a substantively different causal story.
Leaning toward:statistical criterion is satisfied), but report explicitly in the
Supplement that all six policy-relevant intervention effects are essentially
identical under both orientations of this specific edge -- this makes the
4-edge enumeration a robustness demonstration rather than a source of
substantive ambiguity, which is a stronger, cleaner story for reviewers than
either quietly leaving it as 3 edges or promoting it without comment.

Steps 3 (single-run FCI endpoint marks) and 4 (acyclicity check) print to
console only (no CSV) -- not yet reviewed directly, but the identical fit
statistics in Step 5 already imply the reversed edge list fit without
issue (a genuine unresolved cycle would have produced a different, likely
My own console output for these two steps would complete
the record but is not decision-relevant given Step 5.

**Next:** still deciding whether to proceed with the 4-edge/
16-specification enumeration given this read, and separately, script 28
(IPW attrition sensitivity, per my most recent note) is in progress.

## 20. Real 4-edge/16-spec enumeration script delivered; script 28 (IPW) real results; joint PAG edge-type audit (script 29) written in response to a methodological refinement

**Orientation-enumeration v4.** Per my own decision after reviewing script
27, wrote `r_patches/02_full_orientation_enumeration_v4.R`: adds
harm_future->harm_present as the 4th flippable edge (flip_candidates grows
from 3 to 4 rows), giving 2^4=16 structural specifications instead of 8.
Mechanically almost everything in the v3 script generalizes automatically
(scenario generation, the combo-ATE loop, the uncertainty-band summaries were
already written in terms of `nrow(flip_candidates)` / `n_flip`) -- only the
flip_candidates table and the header commentary needed real changes. Kept
explicit per-scenario acyclicity + convergence checks (already present
mechanically) but removed v3's "UNEXPECTED, STOP" framing for a cyclic
result, since with 4 flippable edges we no longer have a strong prior that
all 16 will be acyclic -- per my own explicit note, this is now
checked, not assumed. Added a new `hp_check` diagnostic (section 6b)
specifically for harm_present's ATE range across all 16 specs, since
harm_present is both an intervenable node and an endpoint of the new
flippable edge. Output file paths are unchanged from v2/v3 so
`10_figure7_uncertainty_pub.R` and `08b_supp_figure_all_combo_interventions.R`
need no code changes -- only their "8 specifications" comments/captions were
Not yet run** as of this section.

Ran it.
Overlap: retained median P(retained)=.467 vs. attrited .419, reasonable
overlap, no extreme separation. Weights are very well-behaved: ESS=810.9 of
N=870 untruncated (93% efficiency -- attrition here is only mildly related
to the observed covariates), only 18 of 870 weights truncated at the
[1st,99th] percentile, ESS=816.3 after truncation, max weight only 2.90
untruncated. Balance: age's SMD (retained vs. the full N=1,987 sample) drops
from .167 unweighted to .030 IPW-weighted -- exactly the variable I was
most concerned about, and the correction works as intended; the 8
attitude/context nodes were already small (<=.07 unweighted) and stay small.
Behavior equation refit under weights: harm_present .309->.313,
weather_risk_prep .184->.185, social_norms .114->.115 -- trivial movement,
well within SE. Six single-node ATEs: max difference is belief_concern
(+.0033 SD), all others <=.0006 -- negligible. **Conclusion: the one-sentence
Supplement note I'd proposed is now justified by the actual numbers** --
"Because Wave-5 respondents differed from nonrespondents on age and some
demographic characteristics, we repeated the behavior-facing analyses using
inverse-probability-of-observation weights based on pre-Wave-5 variables;
the substantive conclusions were unchanged (Supplementary Materials)."

**Methodological refinement raised (important, addressed before any
further downstream work):** the scalar asymmetry statistic A_alpha
(P(arrow at destination) - P(arrow at source)) used throughout scripts
21/27 and the v4 flip-set decision collapses several different PAG relation
types into one number. A_alpha near 0 can mean the edge is mostly
BIDIRECTED (X<->Y, arrowhead mass high at BOTH ends -- consistent with
latent confounding, informative) OR mostly CIRCLE-CIRCLE (Xo-oY, low/
undetermined mass at both ends -- genuinely uninformative) -- both produce
similar marginal asymmetry but mean very different things. Recovering which
case holds needs the JOINT mark at both ends WITHIN THE SAME bootstrap
resample, which the pipeline's existing saved objects
(mark_props_ext_by_alpha.rds, fci_props_ext.rds) cannot provide -- these are
already-pooled/marginal sums; the per-resample `results_ext[[alph]]` list
(which DOES contain the joint pairing, transiently) is never written to
disk, confirmed by reading the .qmd's Section 7.4 bootstrap chunk directly
before writing anything. Per my own explicit fallback note, a
targeted rerun of ONLY the FCI half of the bootstrap (not PC, not SCM, not
interventions, not measurement/attrition analyses) is needed to capture the
joint edge type.

**Script 29 (`r_patches/29_joint_pag_edgetype_audit.R`, NEW).** Reruns
exactly the FCI half of run_one_ext() (same suffStat resampling, same
pcalg::fci() call/args, same n_boot=2500 and alphas=c(.05,.01) already
defined in the session) but classifies the joint edge type PER RESAMPLE,
immediately, into 7 canonical categories (i_to_j, j_to_i, bidirected,
i_ocirc_arrow_j [i o-> j], j_ocirc_arrow_i [i <-o j], circle_circle,
no_edge, plus an "other" safety-net bucket for any mark combination that
should not occur from a valid FCI fit) -- never storing a per-resample
array, so memory use matches the original bootstrap. Computes this for
ALL p*(p-1)/2 pairs but reports only the 16 retained edges, at both alpha,
oriented to match the paper's claimed (from, to) direction. Also produces a
`cross_check` table comparing combined bidirected+circle-circle mass
against the CURRENT 4-edge flip set, to surface any disagreement in either
direction (an edge outside the flip set with high confound/unresolved mass,
or an edge inside it that the joint marks actually resolve cleanly). Writes
`pipeline_outputs/joint_pag_edgetype_audit.csv` (32 rows = 16 edges x 2
alpha) and `joint_pag_edgetype_crosscheck.csv` (16 rows). Deliberately does
NOT apply any new classification rule automatically or rerun anything
downstream -- purely descriptive, per my own explicit note that this
audit should settle the question before Table 2 is frozen, not trigger
another automatic pipeline cascade. Balance-checked (Python brace/bracket
checker) and reviewed for the comma-in-block anti-pattern -- clean.
Not yet run.**

**Decision sequencing (my own plan, recorded so it isn't re-litigated):**
run script 29 first. If it confirms the same 4 edges are the genuinely
uncertain ones, keep the 4-edge/16-spec enumeration (v4, already written)
and proceed to run it, then freeze Table 2/Figure 5/§3.5/Abstract/
Discussion/Conclusion around that result. If script 29 changes which edges
are genuinely unresolved, the v4 orientation-enumeration script's
flip_candidates table needs updating to match BEFORE it is run -- do not
run v4 and script 29 as independent, order-agnostic steps.

## 21. Script 29 real results -- caught and fixed a direction-labeling bug; corrected findings are substantively important

**The bug.** Script 29's `run_one_joint()` called
`classify_joint(am[i, j], am[j, i])`. This pipeline's OWN established
convention (already documented and used correctly in scripts 21/27 and in
`run_one_ext()` in the .qmd) is that `am[a,b]` holds the mark AT NODE b, not
at a. The call above passed `am[i,j]` as if it were the mark at i -- backward.
This swapped the labels of exactly the two single-direction categories
(`i_to_j`/`j_to_i`) and the two circle-arrow categories
(`i_ocirc_arrow_j`/`j_ocirc_arrow_i`) for every one of the 16 edges, in both
alpha columns. The four symmetric categories (bidirected, circle_circle,
no_edge, other) were NOT affected (they don't depend on which side is
labeled i vs. j), and `joint_pag_edgetype_crosscheck.csv` was entirely
unaffected (it only uses the symmetric categories).

**How it was caught.** Before accepting the first run's output, checked it
against a temporal impossibility: `harm_present -> climate_behavior` and
`weather_risk_prep -> climate_behavior` both showed 0% mass for the
claimed forward direction and non-trivial mass (45.8% and 21.7% at
alpha=.05) for the reverse. climate_behavior is a Wave-5 outcome, strictly
after the Waves-1-4 harm_present/weather_risk_prep measurements -- reverse
causation is impossible, so the labels had to be swapped. Confirmed this is
a REPORTING bug, not a computation bug (the tallies themselves are
correct, just mislabeled), so **no rerun was needed** -- the corrected
values are just the two swappable column-pairs (`pct_from_to` <->
`pct_to_from`, `pct_from_ocirc_to` <-> `pct_to_ocirc_from`) exchanged.
Fixed the script's call site (`r_patches/29_joint_pag_edgetype_audit.R`,
swapped the two `am[]` arguments) for any future rerun; documented the fix
in-line at the call site.

**Corrected findings (real, substantively important).** Of the 16 retained
edges, 12 show non-trivial-to-dominant BIDIRECTED mass (not circle-circle --
specifically arrowhead-at-both-ends, the PAG relation consistent with
latent confounding). Five edges are more bidirected than either single
direction combined: `politics->policy_support` (64.8%/47.5% bidirected vs.
only 11.5%/16.9% total directed either way), `trust_science->social_norms`
(76.9%/66% bidirected vs. 18.4%/21.3% directed), `belief_concern->
weather_risk_prep` (73.5%/59.5% bidirected vs. 17.9%/22.3% directed),
`harm_present->weather_risk_prep` (66.7%/50.1% bidirected vs. 17.1%/20.5%
directed), `weather_risk_prep->climate_behavior` (62.4%/62.5% bidirected vs.
21.7%/14.5% directed). `social_norms->climate_behavior` (already the
weakest-existence edge, already in the flip set) is essentially never
cleanly directed either way (0%/0%) -- mostly bidirected or absent.

**Critically, in every one of these 12 "fixed" edges, when a clean single
direction IS resolved, it favors the SAME direction already claimed in the
working SCM -- no edge shows a flipped clean-directed majority.** So no
edge's stated DIRECTION looks wrong. What the audit changes is confidence
framing: several edges (especially the five listed above, plus
`harm_present->climate_behavior` at 42.7-44.7% bidirected) are more
honestly described as "adjacency robust, claimed direction favored when
resolvable, but frequently compatible with latent confounding" rather than
"bootstrap-resolved" as Table 2 currently labels them.

**harm_future->harm_present, re-examined under this more granular joint
metric, independently reproduces script 27's threshold-sensitive-direction
finding**: corrected dominant direction is reverse-favored at alpha=.05
(47.2% harm_present->harm_future vs. 29.4% harm_future->harm_present) and
forward-favored at alpha=.01 (27.1% vs. 18%) -- consistent with script 27,
now confirmed via joint marks rather than only the marginal asymmetry
statistic. Its bidirected mass (9.9%/24.7%) is much lower than the five
edges above -- a qualitatively DIFFERENT profile: "genuinely direction-
flippy but low-confounding" (HF-HP) vs. "direction is fine but often
confounded" (the five edges above). Both are "uncertain" in some sense, but
for different reasons that the marginal asymmetry statistic alone cannot
exactly my original point.

**What this means for the decision posed (the either/or: "if audit
confirms the same 4 edges, freeze; if it changes which edges are
unresolved, rerun the enumeration with the new set"):** neither branch
applies cleanly. The 4-edge flip set is not contradicted -- no currently-
fixed edge's direction is actually reversed by this audit, so there is no
NEW edge that needs to become flippable in the orientation-enumeration
sense. But the audit surfaces a different, real issue: several fixed edges
carry substantial latent-confounding-compatible mass that Table 2's current
"bootstrap-resolved" language does not convey. This is an argument for
adding a confounding-compatibility column/note to Table 2 for the relevant
edges (methods-level transparency), not for expanding the structural-
orientation flip set (which is specifically about DIRECTION ambiguity).
Decision on how to reflect this in Table 2/Methods is still open -- reported
back, not decided unilaterally here.

## 22. Three pre-run corrections to script 29 (v2); n_boot manuscript-text discrepancy found; age SMD (.29 vs .167) reconciled

**Three corrections applied before running v2:**
1. **Fail loudly on unanticipated marks.** `classify_joint()` now `stop()`s
   with the resample index, alpha, node-pair names, and both raw mark
   values if a mark combination doesn't match one of the 7 canonical PAG
   relations, instead of silently folding it into an "other" bucket (v1's
   behavior). The classification step was moved OUTSIDE the tryCatch that
   wraps the `pcalg::fci()` call itself (which legitimately can fail to fit
   and should be counted/skipped, not conflated with a mark-encoding
   problem), so a stop() here is never accidentally swallowed and
   propagates all the way up through `furrr::future_map()`, halting the
   entire run.
2. **n_boot reconciled, not changed.** Checked: script 29 inherits
   n_boot=2500 from the live session, the SAME value the primary Section
   7.4 bootstrap already uses. Cross-checked against main2.tex directly --
   the manuscript states "1,000 ... resamples ... at each ... threshold"
   (and "divided by 2,000" when pooling) in several places (Methods, Fig. 4
   caption, Table S caption, Supplement), but the ACTUAL, already-completed
   Section 7.4 analysis (whose own .qmd narrative text elsewhere already
   says "N = 2500") has used n_boot=2500 per alpha (5,000 pooled) all
   along -- confirmed this is what scripts 21/22/23/27/28 all inherited too,
   and therefore what every already-reported bootstrap number in the
   current manuscript (Table 2, Figure 4, existence percentages) is
   actually based on. This is a PRE-EXISTING manuscript-text/code
   discrepancy, not something script 29 introduces -- and not something
   script 29 should "fix" by switching to 1,000 (that would make this
   audit inconsistent with the primary procedure it's meant to check, not
   consistent with it). **Recommended fix: correct main2.tex's several
   "1,000"/"2,000" mentions to "2,500"/"5,000"** (the multiple locations:
   lines ~111, 118, 197, 214, 228, 325, 393, 825, 837) rather than rerunning
   the entire causal-discovery-dependent portion of the paper at 1,000 to
Still deciding   places in the manuscript.
3. **No derived classification.** Removed v1's `likely_confounded`/
   `likely_genuinely_unresolved` boolean columns (a 20% threshold rule).
   The primary output (`joint_pag_edgetype_audit.csv`) is now from, to,
   alpha, and the 7 raw joint-type percentages only -- nothing else. A
   separate `cross_check`/`joint_pag_edgetype_crosscheck.csv` table still
   exists purely as a scan-order convenience (sorted by a plain,
   non-thresholded sum of bidirected + circle-circle mass), explicitly
   labeled as "not a classification."

Also fixed the v1 direction-labeling bug's call site properly in v2 (already
fixed in Section 21, carried forward): `classify_joint(am[j,i], am[i,j], ...)`.
Balance-checked and reviewed for the comma-in-block anti-pattern -- clean.
**v2 not yet run** (v1's results, already corrected for the labeling
bug, remain the most recent real numbers on hand -- see Section 21 for the
corrected reading. v2 additionally re-verifies via the loud-failure check
that no 8th mark pattern exists; if it finishes cleanly, the corrected v1
numbers stand confirmed rather than superseded).

**Age SMD reconciliation (.29 in script 24 vs. .167 in script 28) --
verified arithmetically, not a bug.** Script 24 compares retained (N=870,
mean age 55.363) directly against attrited (N=1,117, mean age 50.914):
SMD = (55.363-50.914)/pooled_sd = .29, implying pooled_sd~15.34. Script 28's
balance check instead compares the (IPW-weighted) retained sample against
the FULL N=1,987 Waves-1-4 sample -- the population IPW is meant to
restore representativeness of, which is the standard target for an IPW
balance check, not the attrited-only group. The full-sample mean age is a
weighted blend of both groups: (870*55.363 + 1117*50.914)/1987 = 52.86.
Retained-vs-full mean difference = 55.363-52.86 = 2.50, roughly 56% of the
retained-vs-attrited difference (4.45) -- exactly the expected shrinkage
from comparing against a population that already contains the retained
group itself as part of the blend. 2.50/15.34 (approx. same pooled SD)
predicts SMD~.163, matching script 28's actual reported .167 (unweighted)
closely. **Conclusion: both numbers are correct; they answer different
questions** -- .29 (script 24) characterizes how different the two groups
are from each other (the standard "who is missing" attrition diagnostic);
.167->.030 (script 28) checks whether IPW successfully restores the
weighted-retained sample's resemblance to the full pre-attrition cohort
(the standard IPW-performance diagnostic). Recommendation: keep both, each
in its own place (Methods/attrition paragraph uses .29; the IPW Supplement
paragraph uses .167->.030), with the comparison population named explicitly
each time so a reader never has reason to expect the two figures to match.

## 23. Script 29 v2 real run: audit CONFIRMS the 4-edge uncertainty set -- cleared to run the v4 enumeration

Ran clean: no `stop()` fired across 2 alpha x 2,500 resamples x 120 pairs
(~600,000 classifications) -- the 7 canonical categories are exhaustive for
this data/algorithm, verified rather than assumed. Output values match the
hand-corrected v1 reading exactly (cross-validation of both the bug fix and
the earlier manual correction).

**Checked every one of the 12 currently-FIXED edges at both alpha: none has
its claimed direction contradicted.** Where any single-direction mass
exists at all, it favors the paper's stated direction at BOTH alpha for
every fixed edge, with no new cross-alpha sign disagreement beyond the
already-promoted harm_future/harm_present. Closest margin among fixed
edges: trust_science->policy_support (43.0/17.7 at .05, 35.4/29.0 at .01) --
tighter than most, but consistently forward-favoring at both alpha, so it
does not meet the sign-disagreement criterion and stays fixed.

**The 3 already-theory-completed edges are unchanged in status, with one
minor phrasing nuance surfaced:** politics->belief_concern and
policy_support->social_norms both show their (small) single-direction mass
leaning toward the REVERSE of the currently-assigned direction (e.g.
politics->belief_concern: 0.5-0.6% forward vs. 34.2-42.7% reverse, dominated
by ~38-58% bidirected either way) -- slightly sharper than Table 2's current
"negative or effectively zero asymmetry" phrasing, which reads as neutral
rather than mildly reverse-leaning. Doesn't change these edges' status
(already flagged as theory-assigned, not bootstrap-resolved), but worth a
still open.
social_norms->climate_behavior remains 0%/0% both directions at both
alpha -- no informative directed mass either way, consistent with its
already-known "weak"/lowest-existence status.

**Decision, per my own stated rule:** the joint audit CONFIRMS the same
four-edge uncertainty set. Cleared to run the v4 orientation enumeration
(`r_patches/02_full_orientation_enumeration_v4.R`, already on the machine)
as the next and, per my own framing, the only remaining substantive
rerun. The separate confounding-transparency observation from Section 21/22
(5+ fixed edges dominated by bidirected mass even though direction is
unambiguous) remains open as a Table 2 transparency question -- does not
block or change the enumeration.

## 24. n_boot provenance resolved: I was right, the earlier "fix main2.tex to 2,500" idea is RETRACTED

Independently checked my own rendered copy of the primary pipeline
and found its actual printed completion message reads "1000 resamples per
alpha threshold" -- contradicting my earlier claim (Section 20/22) that the
CURRENT primary analysis uses n_boot=2500 throughout. Cross-checked this
myself via file timestamps before accepting either claim:
- `pipeline_outputs/fci_props_ext.rds`: dated 2026-03-28 -- OLD, untouched
  by anything this session.
- `pipeline_outputs/mark_props_ext_by_alpha.rds`: dated 2026-09-06 12:05 --
  the FRESH save from this session's Section 7.4 rerun (n_boot=2500, the
  live .qmd's current setting), added this session per Section (earlier).
- `pipeline_outputs/scm_edges_finalized.csv` (which directly produces
  Table 2 and feeds Figure 5): dated 2026-09-05 19:52 -- generated BEFORE
  the Sept-6 fresh 2500-based rerun.

Checked script 18's own loading logic (`18_finalize_scm_specification_v4.R`,
line 61-63): `if (!exists("fci_props_ext")) { ... readRDS("pipeline_outputs/
fci_props_ext.rds") }` -- it uses an in-memory `fci_props_ext` if the
session already has one, and ONLY falls back to the on-disk .rds file
otherwise. Since scm_edges_finalized.csv (Sept 5, 19:52) predates the
fresh n_boot=2500 mark_props_ext rerun (Sept 6, 12:05), no freshly-computed
2500-based `fci_props_ext` could have been in memory when script 18 last
ran -- it must have fallen back to the on-disk file, which is the frozen
March-28 (n_boot=1000-vintage) version. **This confirms my finding
independently: Table 2 and Figure 5's currently-published numbers are
genuinely based on n_boot=1000 per alpha, matching main2.tex's stated
methodology exactly.** My earlier recommendation (Section 20/22) to change
main2.tex's "1,000"/"2,000" language to "2,500"/"5,000" was WRONG and is
retracted -- I had checked the live .qmd's `n_boot` variable value (2500)
without checking whether the object that ACTUALLY produced Table 2's
numbers (scm_edges_finalized.csv, and the frozen fci_props_ext.rds behind
it) reflected that same value or an older one. The correct account,
matching my own independent verification exactly:
- Primary FCI/PC stability analysis (Table 2, Figure 3/4): n_boot=1,000
  per alpha, matching current manuscript text -- NO CHANGE NEEDED.
- Joint PAG edge-type audit (script 29): n_boot=2,500 per alpha, a
  separate, deliberately higher-precision targeted audit -- to be
  described as such in the Supplement, not conflated with the primary
  bootstrap count.
- Intervention bootstrap (script 22/26): 1,000 participant resamples --
  already correctly described in the manuscript, unrelated to either of
  the above (a different resampling unit -- participants, not
  conditional-independence test decisions).

**Process fix (no rerun needed, results unchanged):** to prevent this exact
ambiguity from recurring, decoupled script 29 from the shared `n_boot`
session variable -- it now defines its own explicit `n_boot_joint <- 2500`
constant rather than inheriting whatever the primary pipeline's `n_boot` is
currently set to. This makes the joint audit's resample count independent
of the primary analysis's setting by construction, so a future change to
the primary `n_boot` (in either direction) can never silently change the
joint audit's precision or vice versa. Also flagged as a standing
reproducibility risk: `fci_props_ext.rds` is no longer re-saved by the
current .qmd's Section 7.4 chunk (only `mark_props_ext_by_alpha.rds` is,
since this session's earlier addition) -- it is a static, frozen artifact
from March 28 that script 18 will keep reading via disk fallback
indefinitely unless someone either re-adds a `saveRDS(fci_props_ext, ...)`
call (which would then silently pick up whatever n_boot is live at the time
Section 7.4 is next rendered) or manually regenerates it deliberately at
n_boot=1,000 for full future reproducibility. No action taken on this
still open whether it needs hardening.

**Manuscript-text decision (adopted as-is):** describe the two
bootstrap analyses separately rather than unifying the resample count:
"Primary causal-discovery stability was assessed with 1,000 bootstrap
resamples per conditional-independence threshold. Because the joint
endpoint-type audit was used to make the final distinction between
directional, bidirected, and unresolved PAG relations, we repeated that
targeted FCI audit with 2,500 resamples per threshold to obtain a more
stable characterization of the joint endpoint frequencies." No rerun of
anything required by this decision.

## Section 25: Manuscript integration begins — deterministic 16-spec ATE rerun (script 30)

**Date**: 2026-09-06.

**Context**: With the 16-specification structural-orientation enumeration
(script 02 v4) confirmed clean (16/16 acyclic+converged) and the joint
PAG edge-type audit (script 29 v2) confirmed clean, moved to drafting
the actual manuscript replacement text for Section 3.5 (orientation
sensitivity), Table 2's harm_future<->harm_present row, and the relevant
figure caption. Before finalizing exact numbers, one loose end remained:
script 02 v4 computes single-node ATEs via Monte Carlo simulation
(`simulate_scm_generic()`, n=20000 draws with `rnorm()` noise), which
produced a spurious negative lower bound for policy_support (-.00707) in
scenarios where its only outgoing edge is flipped, making its true
structural effect on climate_behavior exactly 0. This is simulation noise,
not a real negative effect.

**Decision**: rerun the 16-spec single-node ATEs deterministically
via exact mean propagation, rather than reporting "0 (up to Monte Carlo
noise) to .011" in the manuscript text.

**Action**: wrote `r_patches/30_deterministic_16spec_ates.R`. Reuses
`base_edges`, `flip_candidates`, `flip_edges()`, `is_acyclic()`,
`build_lavaan_syntax()`, `topo_order()`, and the scenario-enumeration logic
VERBATIM from script 02 v4 (so the 16 fitted lavaan models, and therefore
all fit indices, are identical to the already-reported v4 fit table — this
script does not touch model fitting, only how ATEs are computed from each
fit). Replaces `simulate_scm_generic()` with `scm_mean_propagate()` /
`single_node_ates()`, the exact deterministic helpers already used and
verified in scripts 22/23/26/27 (in particular, script 27's single
baseline-vs-reversed HF/HP comparison already used exactly these two
functions). Output: `tables/orientation_uncertainty_band_full_deterministic.csv`
(same 6-row shape as v4's Monte Carlo version) plus a side-by-side
`tables/orientation_ate_deterministic_vs_montecarlo.csv` diff table.
Verification: Python balance checker passed; manual review confirmed no
comma-in-block issues; package-namespace calls audited against
already-verified usages elsewhere in the pipeline (dplyr, lavaan, igraph,
purrr, tibble, readr — no new functions introduced). Pushed to
not yet run.

**Also resolved this round**: wondered directly why not just set
n_boot=2,500 everywhere given the confusion over 1,000 vs. 2,500 (see
Section 24). Answer given: doing so would require rerunning the entire
causal-discovery-dependent portion of the paper (Table 2, Figures 3/4,
script 18's edge audit, Figure 5) at real time/risk cost, for zero
scientific benefit — 1,000 resamples per threshold is already a
methodologically standard and adequate count, and the joint-audit's 2,500
is a deliberate, separately-labeled higher-precision choice for that one
targeted analysis, not evidence the primary analysis is undercounted.
Keeping the two explicitly distinct and documented (per Section 24) is
correct and requires no further pipeline changes.

**Also flagged, still pending confirmation**: my drafted "Figure 5
caption" replacement text (diamonds, ranges, "all 15 alternative acyclic
specifications") textually matches the CURRENT content/structure of
Figure 7's caption in main2.tex (line ~350, `fig7_uncertainty_pub.pdf`),
not Figure 5's actual current caption (line ~274, the SCM diagram with
solid/dash-dot edge styling). Not yet applied to main2.tex pending my own
confirmation of which figure I intend the new text for, to avoid
misapplying it to the wrong figure.

**Not yet done**: any main2.tex edits (Section 3.5 replacement, Table 2's
harm_future<->harm_present row, the figure caption(s), the "narrow relative
to baseline" wording fix, the two Supplementary data tables at lines
~1325/~1372 needing 8->16-row data refreshes). These will be applied in one
coordinated pass once script 30's deterministic numbers are back, followed
by a full compile-check per the project's standard verification discipline.

## Section 26: main2.tex updated for the 4-edge structural-orientation story (Results + Supplement)

**Date**: 2026-09-06.

Applied the full coordinated edit pass to `main2.tex` now that script 30's
deterministic numbers were confirmed and I answered the two open
questions (figure = `fig:interventions`, i.e. the current Figure 7/diamonds
plot; trust_science gets the same structural-zero treatment as
policy_support). Backup of the pre-edit file kept at
`Downloads/main2.tex.bak_before_orientation_edits`. Edits applied via five
python `do()`-pattern patch scripts (exact-match replace, assert-unique),
run in sequence on the device:

1. **Orientation-methods paragraph** (main text, edge-orientation-rule
   description): "three edges" -> "four edges," noting the added
   cross-threshold consistency check.
2. **Table 2 intro paragraph**: explains that future_harm->harm_present is
   the "thirteenth" edge with positive *pooled* asymmetry that nonetheless
   reverses sign per-alpha, so it joins the three theory-completed edges as
   unresolved.
3. **Table 2 body**: moved future_harm->harm_present out of the "Bootstrap
   orientation evidence aligned with selected direction" category into a
   new "Threshold-sensitive orientation" category, showing both per-alpha
   asymmetry values ($-.11$ / $+.26$) instead of a single pooled number, with
   a new footnote ($^\S$) explaining the reversal and referencing the pooled
   value for continuity with the table's existing convention.
4. **Section 3.5 (`sec:orientresults`) full body replacement**: my own
   drafted text, adapted to reference `\ref{fig:interventions}` by label
   (not figure number, per my own instruction), using the deterministic
   (script 30) numbers throughout, and including the new trust_science/
   policy_support structural-zero paragraph exactly as I phrased it. Ends
   with my final wording fix ("The overall intervention ranking and the
   absolute magnitude of the effects were broadly stable across
   specifications") replacing the old "narrow relative to baseline" claim.
5. **`fig:interventions` caption** (the actual Figure 7 diamonds/silhouette
   plot): updated to 16 specifications / four unresolved edges / 15
   alternatives, deterministic numbers, and the trust_science + policy_support
   structural-zero explanation. NOTE: kept the existing "smoothed silhouette
   + tick marks" descriptive language because that is what
   `10_figure7_uncertainty_pub.R` actually renders (verified by reading the
script itself), rather than replaced with my own draft caption's "combined
range bars + points," which is a reasonable lay gloss of the same
   silhouette+tick-mark encoding (the silhouette's bounds ARE the min/max
   range; the tick marks ARE the individual alternative-specification
   points), so no figure regeneration was judged necessary. Noted here
   in case I pictured a literally different chart type.
6. **Supplement S10 (`supp:orientation`)**: full rewrite -- the "three
   relationships" intro (now four, with a new paragraph explaining HP's
   different, threshold-sensitivity-based reason for being unresolved,
   distinct from PB/PS/SC's substantive-completion reason); the "$2^3=8$"
   enumeration language (now $2^4=16$); a new paragraph documenting that
   single-node ATEs are now computed by exact deterministic mean
   propagation rather than 20,000-draw Monte Carlo simulation (the 41-target
   combo-ATEs remain Monte Carlo, unchanged); `tab:supp_effectranges`
   (6-row single-node range table) updated to the deterministic numbers;
   `tab:supp_orientationfit` expanded from 8 to the full 16 rows (added HP
   and its six new pairwise/triple/quadruple combinations, using
   `tables/fit_supplement_table.csv` from script 02 v4 -- confirmed the 8
   rows shared with the old table are numerically identical, a useful
   consistency check that v4's refits didn't change anything already
   published); closing paragraphs updated from "narrow relative to
   baseline" framing to "broadly stable ranking/magnitude," with the
   trust_science/policy_support structural-zero point folded in.

**NOT changed, deliberately**: Figure 5's own caption (`fig:scm`, the SCM
diagram with solid/dash-dot line styles) still says "three edges" require
substantive completion -- left as is because future_harm->harm_present's
*pooled* asymmetry is still positive, so it legitimately renders solid under
the figure's existing sign-only line-style rule; recasting it as a 4th
dash-dot edge would require a real figure-regeneration decision that I
did not ask for. Also not touched: Discussion and Conclusion sections
(still reference "three edges," "eightfold range" comparisons) -- explicitly
deferred until Results/Supplement are frozen; the version-history/
five-edge/32-combination comparisons throughout, which remain accurate
descriptions of the project's own history and were not asked to change.

**Verification**: full `pdflatex -> bibtex -> pdflatex -> pdflatex` compile
in a scratch directory (figures copied fresh from
`belief_network_US/figures/`, not the stale Downloads copies). Zero fatal
(`^!`) errors beyond the two pre-existing, already-known issues (missing
`fci_pag_alpha05_example.png`; 8 missing BibTeX entries -- both flagged
since a much earlier segment, unrelated to this edit). Zero broken `\ref`/
`\label` (`??`) anywhere in the compiled output. Confirmed via `pdftotext`
that the new Table 2 category, the new §3.5 numbers, the 16-row fit table,
and the structural-zero language all render as intended. A cosmetic
hyperref "Hfootnote" warning was checked and confirmed pre-existing
(identical `\footnote{}` count, 2, in the pre-edit backup and the edited
file) -- not introduced by this round's edits.

**Outstanding after this round**: Figure 5's own caption/line-style question
above (flagged, not resolved); Discussion/Conclusion rewrite (explicitly
deferred); the pre-existing missing-figure and missing-citations issues;
IPW sensitivity paragraph still needs adding to the Supplement (script 28
wording already drafted in Section 24/25); all HIGH/
POLISH items from the original 18-point review remain deferred per my ownsequencing.

## Section 27: Figure 5 caption fix, IPW Supplement section, Discussion/Conclusion/Abstract rewrite

**Date**: 2026-09-06.

Confirmed Results are now essentially frozen and settled the remaining
publication sequence: caption consistency fix -> IPW supplement ->
Discussion/Conclusion rewrite -> abstract update -> final figure polish ->
provenance/ethics/BibTeX/missing-file cleanup -> submission-format pass.
This round completed the first four.

**Figure 5 (`fig:scm`) caption**: fixed the apparent "3 vs 4 edges"
contradiction flagged earlier. Per that framing -- solid/dash-dot in the SCM
diagram shows evidential basis for the *baseline* orientation, not
membership in the later sensitivity set -- the caption now says so
explicitly and adds a sentence clarifying that future_harm->harm_present
renders solid (its pooled asymmetry is positive) but is still varied in
the orientation-sensitivity analysis because its per-alpha signs disagree.
No figure regeneration needed.

**New Supplement section S11** (`supp:ipw`, "Attrition sensitivity:
inverse-probability weighting for Wave-5 retention"): written from
my drafted paragraph, with two numbers checked against the actual
script-28 output files before inserting (per my own explicit note not
to insert the age-SMD number without reconciling .29 vs .167 first):
- Confirmed via `pipeline_outputs/ipw_balance_check.csv`: age SMD .167
  (unweighted) -> .030 (weighted), both vs. the full N=1,987 sample. The
  separate .29 SMD (from script 24) is the retained-vs-attrited comparison,
  used only as motivating context, and is labeled as such.
- Caught and corrected a factual slip in my own draft: it read "the  maximum final weight was 2.9," but `pipeline_outputs/ipw_weight_
  distribution.csv` shows 2.896 (~2.9) is the *untruncated* max; the
  truncated weights actually used in the analysis top out at 1.965
  (~1.97). Reported both explicitly (2.90 untruncated, 1.97 truncated/
  final) rather than silently using the higher, incorrect number.
- Confirmed n_truncated=18 and the single-node ATE diffs (largest
  belief/concern +.0033 SD, all others <=.0006 SD) exactly match
  `pipeline_outputs/ipw_single_node_ate_compare.csv`.

**Discussion**: rewrote the "Third, ..." finding paragraph and the
latent-confounding limitations paragraph. Removed all "three edges" /
"eightfold range" / "an earlier version of this analysis" history per
my own instruction, and folded in results that were previously only in
the main Results section but not yet summarized in Discussion: the
best-pair/best-triple findings, the four-edge/16-specification robustness
result, and (in the limitations paragraph) explicit mention of parameter/
sampling uncertainty, measurement limitations, and the new IPW check
(cross-referenced by `\nameref{supp:ipw}`).

**Conclusion**: fully rewritten on the same principle -- states the final
4-edge/16-specification story directly, with zero reference to the
project's superseded 3-edge/5-edge/eightfold history.

**Abstract**: paragraph 4 updated from "three edge directions" / (implicit
8 specs) to "four edge directions" / 16 specifications, naming the
future_harm/present_harm threshold-sensitivity explicitly and updating the
closing claim to "broadly stable" (matching the Results/Discussion
wording) rather than "narrow range."

**Deliberately left unresolved**: whether the "five-edge,
32-combination" / "eightfold" historical narrative should also be removed
from Supplement S10 -- the instruction to self was explicit for "Discussion/
Conclusion" but also said more broadly "it should disappear from the
paper." Supplement S10 currently still documents this history (as a
methods audit trail), which has real transparency value; did not
unilaterally strip it pending my own confirmation of scope.

**Verification**: full compile check after each patch and again after all
four edits together (fresh `pdflatex -> bibtex -> pdflatex -> pdflatex`,
figures copied fresh from `belief_network_US/figures/`). Zero new fatal
errors, zero broken `\ref`/`\label`, only the two long-standing pre-existing
issues (missing PNG, missing BibTeX entries). `pdftotext` used to confirm
the new Discussion/Conclusion/Abstract/Figure-5-caption/S11 text all
render as intended and that no old-analysis phrases ("eightfold," "earlier
version of this analysis," "five-edge, 32-combination") remain anywhere
outside Supplement S10.

**Remaining in my stated sequence**: final figure polish;
provenance/ethics/BibTeX/missing-file cleanup (the long-deferred 8 missing
citations + missing PNG); submission-format pass.


## Section 28: Response to my own 16-point manuscript review

Wrote up a detailed, itemized review of the rendered PDF (2026-09-06,
afternoon) identifying 16 remaining issues before the paper could be
called submission-ready, with a requested sequence: numbers/figures
synchronization -> 16-spec Supplement update -> Methods orientation-rule
fix -> Figure 3/4 cleanup -> IPW/behavior sensitivity appendices ->
provenance/admin -> final prose trim. This section records what was done
in response, item by item, plus a full re-verification compile.

**Item 1 (numbers/figures sync) -- done for text; figure regeneration
pending.** Section 3.4 prose and Table 3 (`tab:interventions`)
updated from the stale Monte Carlo baselines (.202/.192/.087/.051/.029/
.011) to the deterministic script-30 values (.206/.195/.093/.057/.034/
.016), including the arithmetic-overlap sentence (.202+.087=.289 ->
.206+.093=.299). Table 3's "(best pair)"/"(best triple)" labels reworded
to "(highest point estimate; not the most stable pair)" and "(most robust
triple)" per my own exact suggested wording, with the caption clarified
that no single pair was consistently ranked best across bootstrap
resamples. **Deliberately NOT touched**: a separate sentence in the
belief/concern x politics interaction paragraph reading "$\Delta=.205$
vs.\ $.202$ SD" -- this may be a different computed quantity (from the
interaction/moderation simulation, not the single-node baseline) that
happens to numerically coincide with the old baseline. Flaggedrather than guessed at.

**Item 2 (regenerate Figure 5/`fig:interventions`, not just its caption)
-- script written, not yet run.** Confirmed the root cause: `r_patches/
10_figure7_uncertainty_pub.R` (v1) loads `tables/
orientation_enumeration_ate.csv` (the Monte Carlo table) by default via an
`if (!exists("full_results"))` fallback, so the rendered PDF still shows
the old baseline diamonds/ranges even though the surrounding text has
moved to the deterministic numbers. Wrote `r_patches/
10_figure7_uncertainty_pub_v2.R`: identical plotting logic (confirmed via
grep that nothing in the actual plotting code hardcodes the old
alternative-specification count -- only comments did), but the data-
loading block now reads `tables/orientation_enumeration_ate_deterministic.
csv` (script 30's output) unconditionally, with two new `stopifnot()`
checks (6 baseline rows, 6x15=90 alternative rows) that fail loudly rather
than silently plotting something wrong. Verified via `Rscript -e
'parse(...)'` that the script is syntactically valid R. Already written
and committed to `r_patches/10_figure7_uncertainty_pub_v2.R`;
not run yet (no R interpreter on this side of the device bridge --
still need to run it; it's fast, no bootstrap, just a data load + a
ggplot render).

**Item 3 (Figure 5/SCM caption "3 vs 4 edges") -- done, previously
completed** (Section 27): resolved via the solid/dash-dot line-style
distinction plus an explicit note that the sensitivity analysis
additionally varies future_harm->harm_present for its own,
threshold-sensitivity reason.

**Item 4 (Figure S6/`figS_all_combo_interventions.pdf` still 8-spec) --
diagnosed as a stale render, not a code problem; rerun needed.** Traced
the figure's data source: `08b_supp_figure_all_combo_interventions.R`
reads `tables/orientation_uncertainty_band_combos_full.csv`, written by
`02_full_orientation_enumeration_v4.R`. Confirmed that CSV's `n_specs`
column is uniformly `16` in the current file (mtime 13:32 today), i.e. the
underlying data is already correct. But `figS_all_combo_interventions.pdf`
itself has an earlier mtime (09:29 today) than that CSV -- the image
predates the last rerun of the CSV and is therefore stale, not wrong-code.
08b's actual plotting logic has no hardcoded "8" (only a header comment
does, cosmetic only). **No code change needed here** -- just needs
to rerun `08b_supp_figure_all_combo_interventions.R` after the current
`v4` CSV, which may already effectively be done, since v4
was already rerun today.

**Items 5-7 (Methods cross-threshold orientation rule; Table 2 bidirected-
mass caveat; fit-description reword) -- done and compile-verified.**
Applied via `patch_review_batch1.py`: (a) new opening paragraph in
`sec:orientmethod` stating explicitly that an edge counts as unresolved
either because pooled asymmetry is negative/zero OR because its preferred
direction differs across alpha=.05 vs .01 thresholds even when pooled
asymmetry is positive; (b) a general caveat sentence inserted immediately
before Table 2 that a positive asymmetry value does not exclude
substantial bidirected bootstrap mass and should not be read as evidence
against latent confounding; (c) the fit sentence reworded from "showed
good comparative and residual fit" to "Comparative fit and SRMR were
favorable, while RMSEA indicated some remaining misfit" (same numbers,
more honest framing given RMSEA=.086).

**Item 8 (Figure 3B roughness) -- caption fixed; image regeneration
flagged, not executed.** Caption text fixed (via patch_review_batch1.py):
removed "one bootstrap draw's endpoint marks before pooling" (factually
wrong -- it's a single full-sample fit, not a bootstrap resample),
replaced with an accurate description. **Investigated the actual image**:
Figure 3B is `figures/fci_pag_alpha05_example.png` -- confirmed this file
DOES exist on disk (145KB, PNG, 1152x1152, dated Sep 5 17:10) and DOES
copy correctly and compile correctly in a fresh scratch-directory check
(the "missing PNG" issue logged as a recurring pre-existing problem in
earlier compile checks this project did not reproduce this time -- worth
double-checking my own local compile if I still see that
error, but it is not currently reproducible from this file). **Also
found**: a second, newer file `figures/fci_pag_alpha05.png` (1500x1500,
generated today at 12:01, alongside `scm_layout_from_pag.csv` and
`pipeline_outputs/std_path_coefficients_16edge.csv` -- the latter's
16 coefficients match exactly what's already in Supplement Table S4, so
no manuscript numbers need to change there). This looks like my own
in-progress redesign attempt (color-coded nodes matching the SCM figure's
palette, apparently aiming for the "shared layout" idea), generated
outside the tracked `r_patches/` scripts (likely from a personal notebook,
`climate_analysis_avg_v2_altweather.html`, modified at the same timestamp).
Viewed both images directly: the new one is **not yet usable as a
replacement** -- its title text ("FCI — PAG (α = 0.05)") overlaps the
Belief Concern node, and the layout is otherwise the same messy
force-directed graph as the original, just recolored. **Left untouched
and left flagged** rather than swapped in, since it reads as my own
unfinished draft rather than a finished figure.

**Items 9-11 (S11 SMD reconciliation; S10 old-history removal; S7/Figure
S6 prose "8"->"16") -- done and compile-verified**, via
`patch_review_batch2.py` and `patch_review_batch3.py`: S11 now explicitly
states that the .29 (retained-vs-attrited) and .167/.030 (retained-vs-full-
sample, unweighted/weighted) SMDs are different comparisons, not
contradictory findings. S10 had three separate historical paragraphs
rewritten to remove all "five-edge," "32-combination," "eightfold,"
"earlier version" language, keeping only the substantive methodological
content (why BT/HW are excluded; why trust_science isn't fragile; the fit
range summary). S7's prose and Figure S6's own caption both updated from
"8" to "16" structural specifications (text only; the image itself is the
separate Item 4 issue above).

**Item 12 (behavior-outcome sensitivity write-up, mitigation/engagement
items only) -- new Supplement subsection written.** Discovered that
`pipeline_outputs/behavior_sensitivity_{ate,fci_adjacency,ggm,
scm_equation}.csv` already contained a complete, previously undocumented
analysis restricting the outcome to the four mitigation/engagement items
(excluding emergency-preparedness and relocation) -- contradicting a
stale decision-log note that this analysis "not yet started." Added new
`S12. Behavior-outcome sensitivity: mitigation/engagement items only`
(`supp:mitig4`) via `patch_review_S12.py`, reporting: FCI adjacency to the
narrower outcome (present harm 94.3%, social norms 71.5%, policy support
60.4%, weather risk dropping to 38.5%, belief/concern 52.9%); the
refitted behavior equation (weather risk beta=.135, SE=.035, p=.0001;
present harm beta=.348, SE=.035, p<.001; social norms beta=.127, SE=.034,
p=.0002); and the corresponding intervention effects (present harm .204,
weather risk .068, social norms .064 SD). Interpretation: weather risk's
behavioral role is real but markedly less bootstrap-stable once the
outcome is narrowed, unlike present harm's, which stays >90% stable
regardless of outcome definition.

**Items 13-16 -- not started, each for a specific reason:**
- *Item 13 (non-Gaussian CI sensitivity, RCoT/mixedCI vs. primary)*:
  inspected the raw data (`pipeline_outputs/sensitivity_{mixedCI,RCoT}_
  vs_primary_{main,ext}.csv`) and found substantial, non-trivial adjacency
  differences for several edges under RCoT specifically (e.g.
  belief_concern-politics adjacency existence .968 primary vs .3 RCoT;
  trust_science-politics .987 vs .5). This needs careful, honest framing
  rather than a quick automated table -- deliberately deferred rather than
  rushed into the Supplement.
- *Item 14 (dataset citation)*: cannot be resolved yet -- the
`\footnote{TODO: cite the source of the dataset.}` placeholder  remains untouched in main2.tex.
- *Item 15 (figure reorder -- GGM as Figure 1, correlation matrix to
Supplementary S1)*: not executed pending confirmation, given the  risk of breaking figure numbering and cross-references throughout the
  document.
- *Item 16 (final prose trim)*: not started.

**Full re-verification compile** (fresh scratch directory, `pdflatex ->
bibtex -> pdflatex -> pdflatex`, figures copied fresh from
`belief_network_US/figures/`): zero fatal (`^!`) errors; the only
undefined-citation warnings are the same 8 known, pre-existing missing
BibTeX keys (unchanged from every earlier check this project); `pdftotext`
confirms the new S12 section, the relabeled Table 3, the Table 2 caveat
sentence, the Methods cross-threshold paragraph, and the Figure 4 caption
pointing to "Table S4 in the Supplementary Materials" all render as
intended. The previously-logged "missing PNG" error
(`fci_pag_alpha05_example.png`) did **not** reproduce in this check -- the
file exists on disk and copies/compiles cleanly; if my own local
compile still shows it, that would point to a stale local copy rather than
a repo-side problem.

**Remaining in my stated sequence**: still need to (a) run
`10_figure7_uncertainty_pub_v2.R` to regenerate Figure 5 with the
deterministic data, (b) rerun `08b_supp_figure_all_combo_interventions.R`
to refresh Figure S6 with the already-correct 16-spec data, (c) decide
whether to finish my own Figure 3B redesign attempt
(`fci_pag_alpha05.png`) or have a clean version written from scratch, (d)
supply the dataset citation (Item 14), (e) confirm or decline the figure
reorder (Item 15), and (f) confirm or correct the flagged $.205$ vs.\
$.202$ SD sentence in the belief/concern x politics paragraph. Item 13
(RCoT/mixedCI sensitivity) and Item 16 (prose trim) remain open tasks for
a future session.


## Section 29: Figure 5 / Figure S6 regeneration confirmed; new S13 (RCoT/mixed-CI sensitivity) written

Ran `10_figure7_uncertainty_pub_v2.R` and `08b_supp_figure_all_combo_
interventions.R`. Verified both outputs directly (staged and viewed both
PDFs): Figure 5 now shows the deterministic baselines (+.206/+.195/+.093/
+.057/+.034/+.016 SD) with trust_science's and policy_support's silhouettes
correctly bottoming at zero; Figure S6 now spans the 16-spec ranges (top
row, present harm + weather risk + social norms, tops out at ~.311,
matching Table 3) and its own caption already read "16." Full compile
check with both fresh figures: 38 pages, zero fatal errors, same 8
pre-existing missing BibTeX keys, no new issues.

**Item 13 (RCoT/mixed-CI sensitivity) -- new Supplement S13 written**
(`supp:rcot_mixedci`, via `patch_review_S13.py`). Read both sensitivity
scripts' headers closely before writing anything, since the raw comparison
CSVs alone would have been misleading:
- Both `16_sensitivity_rcot_fci.R` and `15_sensitivity_mixed_ci_fci.R` use
  `N_BOOT <- 10` (a deliberate small pilot), not the 1,000 used for the
  primary bootstrap -- made this explicit in the S13 text so the reported
  proportions aren't read as comparably precise.
- Did a clean Python join of `sensitivity_RCoT_vs_primary_main.csv`
  against the working SCM's 13 non-outcome edges: 9/13 show perfect
  agreement (both 1.0); the two edges already marked as theory-added
  rather than purely discovery-supported (future harm->trust science,
  belief/concern->weather risk, both daggered in Table S4) show moderate
  drops (.986->.8, 1.0->.7) -- consistent with their already-weaker
  status; one core edge, political orientation->belief/concern, drops
  from .968 to .3 with no comparable excuse -- reported as an open
  caveat rather than resolved or hidden.
- For the three behavior-facing edges: present harm and weather risk to
  climate behavior both stay adjacent but reduced (.939->.4, .923->.6);
  social norms->climate behavior drops to exactly 0 under RCoT --
  reinforces, rather than newly raises, the existing Section-3
  acknowledgment that this edge's existence is the model's least certain
  feature.
- For mixed-CI: read `15_sensitivity_mixed_ci_fci.R`'s header in full --
  it explicitly documents a sparse-discrete-cell risk (cells as small as
  ~3 observations in preparatory checks) for the four ordinal single-item
  nodes, and explains the small N_BOOT=10 as a deliberate cautious pilot
  because the test's real-data reliability had never been established.
  Given that even the single strongest, best-supported path in the whole
  model (belief/concern->future harm, beta=.856) collapses to adjacency 0
  under this test, characterized the near-uniform collapse as evidence
  the test lacks power/stability on this data, not as evidence against
  the discovered structure -- explicitly declined to treat mixed-CI's
  numbers as comparable evidence to RCoT's.
- Caught and fixed my own drafting error before sending: first draft cited
  RCoT to a fabricated key (`zhang2012kernel`) that does not exist in
  references.bib. Checked references.bib directly (`grep -in "rcot|
  strobl|zhang"`), found no matching entry, and removed the citation
  entirely rather than leave an unverifiable reference -- described RCoT
  and the RCIT package by name instead, uncited.
- Rank-copula (the third sensitivity check, `14_sensitivity_rank_copula_
  fci.R`) has no vs-primary comparison CSV on disk (only its own
  standalone `bootstrap_fci_stability_rankcopula_main.csv`, main network
  only, no ext version) -- mentioned in S13's closing sentence as on
  record but not yet joined into a comparable table, flagged for a future
round if I want it added.

**Process note for future compile checks**: `references.bib` lives only in
`~/Downloads/`, not in `~/belief_network_US/` (confirmed via `ls` --
no such file there at all). A scratch-check copy command that falls back
to the `belief_network_US` copy on a Downloads-copy failure will silently
compile with ZERO real citations resolved (bibtex reports "I couldn't
open database file references.bib" and every one of the ~30 in-text
citations shows as undefined) -- this produced a false alarm mid-session
(looked like 16 missing keys instead of the real, stable 8) before being
caught by checking bibtex's own log line ("I couldn't open database
file...") rather than trusting the citation-count alone. Always copy
`references.bib` from `~/Downloads/` specifically; do not rely on the
belief_network_US fallback.

Full re-verification after S13 (fresh `pdflatex -> bibtex -> pdflatex ->
pdflatex`, both references.bib and figures copied correctly this time):
40 pages, zero fatal errors, same 8 known missing keys, only the same 2
pre-existing `??` marks (`pearl2009causality`/`peters2017elements`).


## Section 30: Final-cleanup pass per my own ordered list -- S13 rewrite, .205/.202 fix, prose trim, figure reorder

Wrote up an ordered finishing list: (1) fix/remove Figure 3B, (2) dataset
citation, (3) figure reorder, (4) .205/.202 sentence, (5) prose trim, (6)
consistency sweep -- explicitly stating I want no further new analysis
opened, S13 reconsidered (10-resample RCoT/mixed-CI should not be
presented as a formal edge-by-edge robustness table), and Figure 3
reworked so the single-run PAG moves to the Supplement and the main
causal-discovery figure is the bootstrap-stability matrix alone.

**S13 rewritten short and exploratory** (superseding the version added in
Section 29, via `patch_finalcleanup1.py`, span-replaced since the
subsection's on-disk line-wrapping didn't match an exact multi-line string
match): removed all edge-by-edge percentage claims (no more "9 of 13
perfectly stable," no more precise per-edge RCoT numbers); kept RCoT as a
qualitative exploratory note (broad pattern reassuring, small number of
edges incl. politics->belief/concern and social_norms->climate_behavior
showing visible reduction, flagged as needing a fully-powered rerun rather
than resolved here); dropped mixed-CI's numerical results entirely per
my own explicit preference, keeping only one sentence noting it collapsed
near-uniformly and that we did not directly diagnose why, so we don't
report its numbers. No new analysis was run -- purely a framing/reporting
change on data already on record.

**Interaction sentence (.205/.202) fixed in main text**, using my own
exact provided replacement sentence; the full moderation analysis
(Table~S-numbered, the $N_{\mathrm{sim}}=50{,}000$ shift-intervention
comparison) stays intact and untouched in Supplement S9.

**Stale phrase removed**: "this corrected enumeration" -> "this
enumeration" in S10 (missed in the earlier history-language cleanup
pass).

**Abstract trimmed**: paragraph 4's future/present-harm mechanics clause
removed (already explained in Results); closing paragraph's meta-comment
("They also show why it matters to identify, specifically, which edges
are genuinely unresolved...") merged into one tighter sentence.

**Discussion's "Third" finding paragraph condensed**: the three sentences
re-deriving the four-edge/16-specification enumeration mechanics (already
fully explained in Results Section 3.4) collapsed into one sentence that
states the robustness conclusion and points to
Section~\ref{sec:orientresults} for the mechanics, rather than repeating
them.

**Conclusion's parallel sentence tightened** the same way (dropped the
repeated "central relationship between future and present harm" clause,
kept the substantive claim).

**Full consistency sweep** run before and after this batch (grepped for:
1000, 2500, "8 specifications", "16 specifications", weather_risk_prep,
belief_concern, the old stale baselines .192/.087/.051/.029/.011, "best
pair", "earlier version", "corrected", TODO, "eightfold", "five-edge",
"32-combination", "three edges"). Only two real hits, both fixed above
("corrected enumeration"; the interaction sentence). The `TODO
: cite the source of the dataset.` placeholder remains --
confirmed this is the one genuine remaining item still needing real input,
not something fixable from the pipeline data.

**Figure 3B resolved without a redraw.** Investigated whether the messy
main-text single-run PAG panel (`fci_pag_alpha05_example.png`, Figure 2
panel B) could simply be dropped rather than fixed, per my own stated
preference to move it to the Supplement. Checked the existing Supplement
Figure `fig:causal_supp` (`figures/fig_supp_causal_redesign.pdf`, in S4)
first -- viewed it directly and confirmed it already shows exactly this
content, and more of it: single-run FCI AND PC-stable PAGs, at BOTH
$\alpha=.05$ and $\alpha=.01$, using a clean, shared, publication-quality
node layout (filled circles, consistent styling) -- i.e., already the
"clean shared-layout PAG" already wanted, already in the
Supplement, just not yet cross-referenced from the main-text figure.
**This means my own in-progress `fci_pag_alpha05.png` redraw is no
longer needed for this purpose** -- noting this so time isn't spent on it
again, unless it turns out to be needed for some other reason.

**Figure reorder executed** via `patch_figreorder.py`, three edits:
1. Removed the correlation-matrix figure (`fig:correlations`,
   `fig1_heatmap_construct.pdf`) from its main-text position (previously
   opening the Results before the GGM figure).
2. Re-inserted it into the Supplement as the very first supplementary
   figure (right after the `\setcounter{figure}{0}` reset, before S1),
   with a caption update noting it's now supplementary background for
   the main-text GGM figure. The single in-text `\ref{fig:correlations}`
   mention (Results, "were all positively correlated...") was left
   untouched -- it now correctly resolves to "Figure S1" automatically,
   no manual renumbering needed anywhere since every figure reference in
   this document uses `\ref`/`\label`, never a hardcoded numeral
   (confirmed via a `grep` sweep for bare "Figure N" text -- zero hits).
3. Rewrote the `fig:causal` figure environment: dropped the two-panel
   minipage structure and panel B's PAG image entirely, widened panel A
   (the bootstrap-stability dot matrix) to stand alone, and rewrote the
   caption to point to the existing `fig:causal_supp` (single-run PAGs)
   and `fig:pc_stability_supp` (PC-stable matrix) in the Supplement
   instead of describing a panel B that no longer exists in the main
   text.

Turned out no other main-text figure needed physically moving: the
existing order (GGM, causal, SCM, interventions) already matched the
requested Figure 1-4 sequence once the correlation matrix was removed
from the front -- confirmed by checking there were no other
`\begin{figure}` blocks between them.

**Full re-verification compile** (fresh `pdflatex -> bibtex -> pdflatex ->
pdflatex`, references.bib from `~/Downloads/` specifically -- see the
process note in Section 29): 39 pages, zero fatal errors, zero undefined
or multiply-defined references/labels, same 8 known missing BibTeX keys.
`pdftotext` confirms the new numbering exactly as intended: Figure 1 =
GGM main + behavior extension, Figure 2 = FCI bootstrap stability (single
panel), Figure 3 = working SCM, Figure 4 = intervention/uncertainty,
Figure S1 = the relocated correlation matrix (immediately before "S1.
Item coding..."), and Figure 2's caption correctly cross-references
Figure S5 (`fig:causal_supp`) and Figure S6 (`fig:pc_stability_supp`) for
the single-run PAGs and PC-stable matrix.

**Remaining before submission**: only the dataset citation/provenance
statement (still needs real input -- exact publication/data source,
recruitment description, ethics/consent wording, data-use statement) and
a final read-through; everything else on my ordered finishing list
is now done and compile-verified.


## Section 31: Drew my own clean FCI PAG; Figure 2 and Figure 3 layout/styling fixes

Decided to keep the single-run PAG in the main text after all and drew
my own clean, SCM-matching version (`figures/FCI_Graph.png`, same
`node_family_colors` palette as the SCM/GGM figures) rather than using the
existing Supplement figure. She asked for two layout fixes.

**Figure 2 (fig:causal): panel B restored, vertically misaligned titles
fixed.** Reverted the single-panel version from Section 30 back to a
two-panel minipage layout using my new `FCI_Graph.png` in place of the old
`fci_pag_alpha05_example.png`. Kept the corrected caption wording from
earlier this session ("fit once to the full sample, not one bootstrap
resample") rather than reverting to the factually-loose "one bootstrap
draw's endpoint marks" phrasing I'd pasted back into my own note -- the two
images' baked-in "A"/"B" titles rendered at different heights when both
panels were simply top-aligned (`[t]` minipages), since the two images have
different aspect ratios and different internal top-margins before their
title text starts. Rather than guess a fix, built the exact two-panel
layout in a standalone test document in the cloud sandbox (which has
`pdflatex`, `pdftoppm`, and Python/PIL available, unlike the device
bridge), rasterized it, and measured the actual pixel row where each
panel's title starts: 242px vs 288px at 200 DPI, a 46px (~16.5pt)
mismatch. A `\vspace*{-16.5pt}` at the top of panel B's minipage was tried
first and made things drastically worse (large negative vspace at the very
top of a minipage does not behave as naive intuition suggests); wrapping
the image in `\raisebox{16.5pt}{...}` instead worked cleanly, reducing the
misalignment to 1px. Applied to `main2.tex`, verified visually (rendered
the real compiled page) that "A FCI bootstrap stability" and "B Example
FCI PAG" now sit on exactly the same line.

**Figure 3 (fig:scm): the 16-edge coefficient table (removed in Section 30
in favor of pointing to Supplement Table S4) restored**,
alongside the diagram as before, but with three improvements I wanted
for: (1) the two minipages (graph, table) switched from `[t]` to `[c]`
vertical alignment so the shorter graph+legend block centers against the
taller 16-row table rather than leaving a large empty gap below the graph;
(2) the header row given a colored background (`#244E68`, "blue_dark" --
the same accent color already used throughout the R figures, confirmed by
reading `r_patches/03_figure_style.R`'s `COL` list) with white bold text;
(3) the whole table set in `\sffamily` (via the `helvet` package,
Helvetica as LaTeX's standard metric-compatible stand-in for Arial) to
match `FIG_FONT <- "Arial"`, the font actually used in all the R-generated
figures, confirmed from the same style file. Added `\usepackage[table]
{xcolor}` and `\usepackage{helvet}` to the preamble for this (not
previously loaded). First attempt at the original .60/.36 textwidth split
produced an 8.23pt overfull hbox in the table (a few of the longer edge
labels, e.g. "Present harm -> Climate behavior", didn't fit in a plain
`tabular`'s un-wrapped `l` column at that width); rebalanced to .57/.40
and confirmed zero overfull warnings in that figure afterward. The
caption's closing sentence was also reverted to describe the coefficients
as "listed alongside the diagram" (matching the restored inline table)
while still pointing to Table S4 for full SEs and significance tests.

**Both fixes iterated and pixel-verified before touching the real file**:
built and rendered standalone test documents in the cloud sandbox with the
actual production image files (staged via the device bridge) rather than
guessing spacing values and re-checking through slow device-bridge
round-trips. Full compile-check on the real `main2.tex` afterward (fresh
`pdflatex -> bibtex -> pdflatex -> pdflatex`): 39 pages, zero fatal
errors, zero overfull/underfull warnings introduced by either change (the
only two overfull hboxes remaining, at lines 927-931 and 1320-1323, are
pre-existing prose paragraphs unrelated to this edit), same 8 known
missing BibTeX keys. Rendered pages 8 and 10 directly and visually
confirmed both figures look correct at real document scale.


## Section 32: Figure 2 width recalibration, Figure 3 error fix / smaller font / bigger graph

Changed Figure 2's panel A from `.45\textwidth` to `.50\textwidth` and wanted to
recalibrate the alignment; separately, in a screenshot of my own editor, I saw
a live "Undefined control sequence" error on `\rowcolor{scmblue}` in my copy of
Figure 3, and asked for a smaller table font and a graph slightly bigger than the
table.

**Important divergence discovered**: the actual `main2.tex` on `~/Downloads/` (the
one this session edits) did NOT match either of my two most recent pasted
snippets -- it still had the last state pushed in Section 31 (Figure 2 at
`.45/.45\textwidth`, raisebox 16.5pt; Figure 3 at `.57/.40\textwidth`), and it
already had `\usepackage[table]{xcolor}` and `\definecolor{scmblue}{...}` in the
preamble (added in Section 31), so the `\rowcolor` error could not reproduce here.
Appears to be hand-editing a separate copy of the file (not the one synced
to this device's Downloads folder) -- flagged directly rather than silently
assumed reconciled.

**Figure 2 raisebox recalibration for `.50/.45\textwidth`.** The previous 16.5pt
raisebox was calibrated for the old `.45/.45` split and goes stale whenever either
panel's width changes, because the two images have different aspect ratios (so
widening one panel changes its rendered height, and LaTeX's `[t]`-minipage
side-by-side placement effectively pins the two boxes' *bottoms* together when a
box's reference point is a single non-text graphic, not their tops -- confirmed by
a diagnostic fbox-boundary test showing both panels' fbox bottoms coincide while
their tops differ by the height difference). Rebuilt the same fbox-boundary +
pixel-scan test used in Section 31, first re-validating the method reproduces
16.5pt at the old `.45/.45` split (got 16.98pt, matching to within measurement
noise), then measured the new `.50/.45` split: 42.6pt. Applied, then rendered the
actual two-panel figure (no fbox) and visually confirmed "A"/"B" baked-in titles
sit on the same line.

**Figure 3 error fix, font, and width ratio.** Added nothing new to the real
`main2.tex` preamble (xcolor/helvet/scmblue were already present from Section 31)
-- the fix that mattered was pasting back the complete preamble block again
in case my own separate working copy is missing it. Shrunk the table font
`\footnotesize` -> `\scriptsize` and tightened `\arraystretch` from 1.15 to 1.12
to match; rebalanced the two-panel widths from `.57/.40\textwidth` to
`.60/.35\textwidth` so the SCM diagram is visibly larger than the coefficient
table, as requested. Verified in a cloud-sandbox test build with the real
`fig5_scm_hierarchical.pdf` that this combination produces zero overfull/underfull
warnings (the original `.60/.36` split had produced an 8.23pt overfull hbox before
the font was shrunk in Section 31; shrinking the font this round freed enough
width that `.60/.35` is now clean).

**Full re-verification compile** (fresh `pdflatex -> bibtex -> pdflatex ->
pdflatex`, `references.bib` copied from `~/Downloads/` specifically): 39 pages,
zero fatal errors, same 8 known missing BibTeX keys, same 2 known `??` citations
(`pearl2009causality`/`peters2017elements`), zero new overfull/underfull warnings
(the pre-existing ones at lines 927-931 and 1320-1323, and a cluster of
underfull hboxes in the Supplement's item-coding table around lines 454-587, are
unrelated to this edit and were already present). Rendered pages 8 and 10 directly
and visually confirmed both figures look correct at real document scale.

Backup on disk (device, `~/Downloads/`): `main2.tex.bak_before_fig2fig3_relayout`.


## Section 33: Figure 3 spacing/header polish, old-Figure-S4 removed, FCI+PC stability matrices combined into new S5, dashed-edge count double-checked

**Figure 3 further polish** (three quick follow-ups to Section 32, each cloud-sandbox-tested before applying): removed the `\hfill` between the SCM diagram and the coefficient table minipages so they sit immediately adjacent (my own note: "no space between scm and table") -- some visual whitespace remains on the graph's right edge, but that is margin baked into `fig5_scm_hierarchical.pdf` itself, not closeable from LaTeX spacing; would need the source PDF cropped if I want it gone entirely. Shrunk the table font further, `\scriptsize` -> `\tiny`, and rebalanced the two panels from `.60/.35\textwidth` to `.63/.32\textwidth` so the graph reads more dominant -- verified zero overfull/underfull at this combination. Bumped just the header row ("Edge"/"$\beta$") to `\small` while data rows stay `\tiny`, so the header is visually distinct; confirmed this doesn't reintroduce any width warnings. Also confirmed, by inspecting the render, that the table's font is Helvetica (`\usepackage{helvet}` + `\sffamily`), matching `FIG_FONT <- "Arial"` in `03_figure_style.R` -- not Palatino.

**SCM dashed-edge count double-checked, no redraw needed.** Wondered whether the SCM diagram should now show 4 dashed ("substantive completion") edges instead of 3, worried something had drifted out of sync during the layout edits. Traced this to source: read `pipeline_outputs/scm_edges_finalized.csv` directly and counted `final_tier == "substantive"` -- exactly 3 rows (politics->belief_concern, policy_support->social_norms, social_norms->climate_behavior), matching the main-text sentence enumerating "the remaining three edges" verbatim. Also zoomed into the actual `fig5_scm_hierarchical.pdf` at 300dpi and visually confirmed exactly 3 dashed lines. The "4" I was recalling is a different count entirely: Section 3.5's structural-orientation sensitivity analysis varies 4 edges (the same 3, plus future harm->present harm), because that edge's *pooled* bootstrap asymmetry is positive (so it renders solid in Figure 3) even though the sign reverses when the two alpha thresholds are checked separately -- this distinction (pooled sign for the diagram vs. per-threshold consistency for the sensitivity analysis) is already spelled out explicitly in both the main-text paragraph introducing "four remained directionally unresolved" and in Figure 3's own caption. Confirmed nothing was changed inconsistently; no data or figure edit was needed, just an explanation.

**Old Figure S4 (Wave-5 behavior-item correlation heatmap) removed.** Flagged `figS_behavior_corr.pdf` as visually out of place (plain base-R heatmap + dendrogram, generic `cc_behavior_*` item codes, no relation to the paper's styled-figure system) and considered whether it could just be dropped. Confirmed it has no in-text cross-reference beyond its own label, and that its substantive content (two-cluster item structure, all-positive correlations, supporting the one-composite decision) is already stated in the surrounding S2 prose, so removed the figure environment outright. All later supplementary figures renumber automatically via `\ref`/`\label` (no hardcoded "Figure S_n_" numerals anywhere, confirmed by grep as in every previous reorder this project). Document is now 38 pages (was 39).

**New Figure S5: FCI's detailed stability matrix combined with PC-stable's, side by side.** Noticed the main-text dot-matrix (Figure 2A) doesn't show exact bootstrap proportions or edge-type (directed/bidirected/undirected) the way the existing PC-stable supplement figure does, and wanted an FCI version in that same style, shown alongside PC-stable's. Before building anything, checked whether a matching, current FCI figure already existed rather than reusing a possibly-stale one: `r_patches/10_cci_addition_and_figure4_pc_appendix.R` turned out to already generate exactly this -- a shared `plot_edge_matrix()` function renders both `figures/figS_stability_matrix_detailed.pdf` (FCI) and `figures/figS_stability_pc.pdf` (PC-stable) in the same restyled "v4" visual system (tint fill + colored stroke/text, real ggplot legend, matches Figure 1's grammar), from the same current bootstrap data (`fci_props_combined`/`pc_props_combined`), same script run (same file mtimes) -- it was already built and intended for exactly this pairing (the script's own comment: "moves to the Supplement instead of being discarded, alongside the PC-stable matrix below"), just never actually wired into `main2.tex` before now. Ruled out two other candidates first: the current main-text dot-matrix (different visual encoding, extended 9-node scope vs. PC-stable's 7-node scope -- would look and scope-mismatched) and an older archived discrete-cell FCI matrix, `fig4_stability_matrix.pdf` (matching style and scope, but superseded/archived as of 2026-09-05 per its own script's location in `_archive_2026-09-05/`, and its "ugly as f, not pub-ready" v3 styling per that script's own header notes -- explicitly what the v4 rewrite fixed). Replaced the old standalone `fig:pc_stability_supp` figure (PC-stable only) with a two-panel minipage version (FCI left, PC-stable right, `.48/.48\textwidth`), rewrote its caption to describe both panels and explicitly connect it to the main-text Figure 2A dot-matrix ("the exact per-endpoint-mark proportions underlying the continuous stability summary"), kept the same `\label{fig:pc_stability_supp}` so the existing main-text Figure 2 cross-reference still resolves correctly. Now Figure S5 (post-S4-removal renumbering).

**Full re-verification compile** after all of the above (fresh `pdflatex -> bibtex -> pdflatex -> pdflatex`, `references.bib` from `~/Downloads/` specifically): 38 pages, zero fatal errors, same 8 known missing BibTeX keys, same known overfull/underfull list (unchanged from before this batch). Rendered the actual Figure 3 page and the new Figure S5 page directly and visually confirmed both.

Backups on disk (device, `~/Downloads/`): `main2.tex.bak_before_fig3_tiny`, `main2.tex.bak_before_fig3_nogap`, `main2.tex.bak_before_fig3_header_size`, `main2.tex.bak_before_remove_s4`, `main2.tex.bak_before_s5_combine`.


## Section 34: orientation-uncertainty criterion redefined, two SCM edges reversed to follow bootstrap direction, sensitivity flip-set changed

Went back to the directional-sensitivity setup after discussing it directly: the
old approach was mixing two genuinely different situations under one "directionally
unresolved" label. `politics->belief_concern` and `policy_support->social_norms`
have LARGE pooled orientation asymmetry (-.336 / -.325, consistently negative at
both alpha=.05 and alpha=.01, not a staleness artifact) -- the bootstrap isn't
undecided about these, it's actively saying the other direction. That's a
theory-vs-data conflict, not directional uncertainty, and lumping them into the
same flip set as edges with genuinely weak evidence made the sensitivity analysis
harder to defend than it needed to be.

**New rule, applied uniformly to all 16 edges:** an edge counts as directionally
uncertain (dashed in Figure 3/5, and varied in the sensitivity enumeration) if its
pooled orientation asymmetry falls under a .10 magnitude band, or if its sign
disagrees between alpha=.05 and alpha=.01. Implemented as a single threshold change
-- `ORIENTATION_ASYMMETRY_EPS` in `clean_pipeline/00_config.R`, .01 -> .10 -- since
`classify_orientation()` in `04_scm_finalize.R` already rejects negative asymmetry
regardless of magnitude, so the sign-disagreement case (harm_future->harm_present,
pooled +.076 but -.108/+.261 per-alpha, per `r_patches/21_orientation_crossalpha_
table.R` and `27_diagnose_hf_hp_edge.R`) is already caught by the magnitude band
alone for the current bootstrap array; didn't wire in a second, separate per-alpha
check since it wouldn't change the outcome and would be one more thing to keep in
sync.

Checked this threshold change against all 16 edges' live asymmetry (Sep-10 verified
bootstrap array) before applying it, specifically to make sure nothing else crossed
the new band unintentionally -- confirmed only the intended edges move.

**Two edges reversed in the working SCM itself**, not just the sensitivity set:
`politics->belief_concern` becomes `belief_concern->politics`, and
`policy_support->social_norms` becomes `social_norms->policy_support` -- following
the bootstrap's actual preferred direction instead of asserting the old
theory-driven one against it. Updated `CURRENT_SCM_EDGES` in `04_scm_finalize.R`
and the canonical `base_edges` in `05_scm_intervention_helpers.R` (kept identical
sets -- 05's own provenance check would fail loudly against `scm_edges_finalized.csv`
otherwise). Confirmed the reversed 16-edge set is still acyclic. These two edges
now render solid (data_aligned) rather than dashed, since their asymmetry, correctly
signed, clears the new .10 band easily.

**Sensitivity flip-set changed accordingly.** Out: `politics->belief_concern`,
`policy_support->social_norms` (no longer uncertain -- now correctly oriented, not
flip candidates). In: `belief_concern->weather_risk_prep` (+.065) and
`harm_present->weather_risk_prep` (+.078), both genuinely near zero. Unchanged:
`social_norms->climate_behavior` (0) and `harm_future->harm_present` (sign-flip).
Updated `flip_candidates` in `05_scm_intervention_helpers.R` (canonical, used by
clean_pipeline 06-14). The r_patches copies of this same table (`02_full_
orientation_enumeration_v4.R`, `30`/`31`/`32_deterministic_*_ates.R`) still have
the OLD flip set and have NOT been updated -- they predate 05's consolidation and
aren't sourced by `run_all.R`, but they'll give stale/wrong answers if run as-is.
Worth archiving or updating before anyone reaches for them again.

**Still needed before any of this is trustworthy:** refit the lavaan SCM with the
two reversed edges and paste the new `belief_concern->politics` /
`social_norms->policy_support` coefficients into `07_figure5_scm_hierarchical_v3.R`'s
`BETA_TR` (currently `NA_real_` placeholders -- the script will refuse to plot
until they're filled in, by design) and check global fit hasn't degraded; rerun
`04_scm_finalize.R` (new final_tier for all 16 edges), then the reference pipeline,
then the 16-specification sensitivity enumeration with the new flip set. Redraw
Figure 3/5 only after the refit, not before -- the dashed/solid pattern depends on
the real coefficients existing, not just the edge list being right.

Not yet touched: `main2.tex`'s Table 2 (still shows the old politics->belief_concern
/ policy_support->social_norms direction and asymmetry) and its Methods/Section 3.5
text (still describes the old 4-edge flip set). Both need updating once the refit
is in, and are being left alone until then per standing instruction not to edit
`main2.tex` without being asked directly.

## Section 35: refit complete, betas filled in, direction fix propagated to 15/17/29; 11_interaction_moderation.R flagged as a real methodological question, not a mechanical fix

Kyuri reran `01_data_prep.R` + the refit (N=870, MLR, `fixed.x=FALSE`): converged,
CFI=.978, TLI=.959, RMSEA=.089, SRMR=.046 (vs. old CFI=.979, TLI=.962, RMSEA=.086,
SRMR=.039 -- essentially unchanged, no fit degradation). `politics ~ belief_concern`
=.530 (matches the old reverse-direction coefficient exactly, as expected for a
bivariate/single-predictor standardized regression -- symmetric with the correlation,
not an error); `policy_support ~ belief_concern + trust_science + politics +
social_norms`, social_norms coefficient=.135. Both pasted into
`r_patches/07_figure5_scm_hierarchical_v3.R`'s `BETA_TR` (replacing the `NA_real_`
placeholders), global-fit comment updated with both old and new numbers.

**Propagated the two-edge reversal to the other places that had their own hardcoded
copies**, beyond 04/05/07 (already done in Section 34):
- `15_confound_sensitivity_diagnostic.R`: `confound_candidates`'
  `policy_support->social_norms` row reversed to `social_norms->policy_support` --
  otherwise `drop_edges()`'s exact-direction match would silently fail to remove
  this edge from `base_edges` once the direction changed, and the "all-8-confounds"
  specification would silently include an edge it was supposed to drop.
- `17_edge_confounding_classification.R`: its own `confound_candidates_in_15`
  cross-check copy updated the same way, so it doesn't fire a false "does NOT
  match 15" warning.
- `r_patches/29_joint_pag_edgetype_audit.R`: `retained_edges` (the 16-edge
  from/to labeling convention for the joint PAG audit table) and `current_flip_set`
  (descriptive cross-check column) both updated to the new direction/flip set.
  Purely a reporting-convention change -- the audit itself still estimates all
  p*(p-1)/2 pairs regardless of this list, so nothing about what's computed
  changes, only which physical direction gets labeled "from_to" vs "to_from" for
  these two edges and which edges get flagged `in_current_flip_set` in the printout.

Verified programmatically afterward that `04`'s `CURRENT_SCM_EDGES`, `05`'s
`base_edges`, and `29`'s `retained_edges` now have byte-for-byte identical
(from,to) key sets (16/16 match) -- no accidental third copy left inconsistent.

**`clean_pipeline/11_interaction_moderation.R` was deliberately NOT fixed**, and
this is a bigger deal than the label/tribble mismatches above. Its `edge_abbrev`
tribble has the same stale-direction problem as 15/17/29's lists (still lists
`"politics","belief_concern","pol_bc"` and `"policy_support","social_norms","ps_sn"`),
so its own guard (`stopifnot(identical(key(base_edges), key(edge_abbrev[...])))`)
will fail loudly the next time this script runs. But swapping just those two
tribble rows would silently produce a wrong result rather than fix anything,
because this script's whole Part D analysis (the politics-conditional SHIFT
intervention feeding Table S9) hardcodes a causal story that the edge reversal
has now inverted:

- `scm_resid_sd` is built with the comment "politics is exogenous (SD=1)" and
  `endogenous_nodes <- setdiff(all_nodes, "politics")` -- true under the OLD SCM
  (politics had no incoming edges), false under the NEW one: `belief_concern` is
  now the node with zero incoming edges in `base_edges` (politics's only equation
  is now `politics ~ belief_concern`), so `belief_concern`, not `politics`, is the
  actual exogenous root.
- `simulate_scm_shift()`/`exact_shift_mean()` treat `politics_level` as a fixed
  external input and compute `belief_concern` FROM it (`bc <- p$pol_bc * pol + ...`).
  This is backwards from the new working SCM regardless of what `p$pol_bc`'s
  numeric value is -- the .530 coefficient is correct in magnitude (bivariate
  symmetry, see above) but using it to generate belief_concern from a fixed
  politics level re-asserts the OLD causal direction inside the simulation itself.
- Part D's entire stated rationale -- "a shift, not a common absolute target,
  because politics -> belief_concern is a real fitted path and a common
  do(belief_concern=0.5) would give the two political groups very different-sized
  treatments" -- no longer holds. Under the new SCM politics doesn't confound
  belief_concern (it's a descendant of it, not a common cause), so the original
  reason for using a shift-conditional-on-politics design instead of a plain
  absolute do(belief_concern=X) may no longer apply, or may need to be reframed
  entirely (e.g. conditioning on politics as a downstream consequence rather than
  an upstream confound changes what the S9 comparison is actually testing).

This needs a real decision from Kyuri/the supervisor about what Part D and Table S9
should now measure, not a code fix -- flagged to her directly rather than guessed at,
the same way the earlier "belief_concern ~ ... + politics" wording contradiction was
flagged instead of silently implemented. Parts A/B/C of the same script (the
politics x belief_concern -> policy_support interaction significance test itself)
are NOT affected by this -- that test only needs both variables as joint predictors
of policy_support, which is unchanged, and doesn't depend on which of the two is
causally upstream of the other.

Committed as `d75a542`: the 07/15/17/29 fixes plus the regenerated
`pipeline_outputs/scm_edges_finalized.csv`. `11_interaction_moderation.R` left
unmodified and unrun.

Not yet touched, per standing instruction: `main2.tex`.

## Section 36: cyclic orientation combinations stay excluded (12 of 16 valid); propagated the "15 alternatives" staleness to Figure 4's script; flagged run_all.R's locked baseline numbers as due for an update

Decision (Kyuri, in chat): keep the 4 structurally-cyclic combinations (combo_1/5/9/13,
see Section 35 addendum below / chat) excluded from every downstream table and figure,
exactly as `is_acyclic()` already does. This is standard practice for this kind of
orientation-sensitivity enumeration, not a selective exclusion needing justification
beyond "these aren't valid SCMs" -- a graph with a directed cycle can't be fit as a
recursive path model in the first place, so there's no ATE to report for it. All
reported ranges/figures are based on the 12 valid (acyclic, converged) combinations
out of 16 possible.

**Propagated to the one other place that had a matching hardcoded assumption**:
`r_patches/10_figure7_uncertainty_pub.R` (the live Figure 4 script -- confirmed current
via file date and its own header referencing the deterministic-8node CSV that
`clean_pipeline/07_intervention_ates_8node.R` writes; the sibling `_v2`/`_hump`/
`_jitter`/`_violin` files in the same directory are earlier exploratory variants, not
the one in use). Its plotting logic was already robust (filters `scenario !=
"combo_0"` generically, no hardcoded row count), so no functional change was needed --
but its header comments claimed "15 alternative orientations" and named the OLD
4-edge flip set (politics->belief_concern etc.), which would have quietly misdescribed
the figure to anyone reading the script later. Updated to the correct count (11
alternatives, 12 total valid combos) and the current flip set, with a pointer to
`clean_pipeline/07_intervention_ates_8node.R`'s header for the exact cyclic-combo
mechanism instead of repeating it.

**Flagged, not yet fixed** (needs a fresh run, not a guess): `clean_pipeline/run_all.R`'s
own `expected_baseline` vector (a locked snapshot of the 8 single-node baseline ATEs
from 2026-09-08, used as a reference-run sanity check) predates the edge reversal and
will not match the new model -- Kyuri's own already-completed `06` run confirms this
for the 6 nodes it covers (e.g. policy_support locked=.0162, new=0 exactly; trust_science
locked=.0342, new=.0437; social_norms locked=.0571, new=.0634). This means the next
`run_all.R` execution will correctly print "REFERENCE RUN MISMATCH" -- that's the guard
doing its job, not a bug, and not something to silence. Once Kyuri runs `run_all.R`
fresh (restart R first, per its own header instructions) and gets a clean end-to-end
pass, `expected_baseline` should be updated to that run's real combo_0 values for all
8 nodes (need harm_future/politics from a successful `07` run, which requires this
session's `07_intervention_ates_8node.R` fix -- not yet confirmed run by Kyuri as of
this entry).

Not yet touched, per standing instruction: `main2.tex`.

## Section 37: cyclic combinations stay a clean exclusion, not a separate analysis -- decision finalized, adopted Methods/Supplement framing recorded

Resolved (Kyuri, relaying a colleague's review of Section 36's cyclic-combo question):
we do NOT fold the 4 structurally-cyclic orientation combinations into the same
sensitivity analysis as the 12 acyclic ones, and we do NOT pursue a feedback/
equilibrium-SEM treatment of them for this paper. Two reasons, both already true of
the existing implementation, now made explicit for Methods: (1) the working SCM's
intervention calculation (`scm_mean_propagate()`, topological forward propagation)
assumes a recursive/acyclic structure -- a cyclic completion isn't a variant of the
same calculation, it's a calculation that isn't defined; (2) FCI itself represents
causal structure under an acyclic-graph framework, so a completion implying a directed
feedback cycle is outside the causal model class the discovery stage is answering
questions about, not just computationally inconvenient.

**Adopted framing, to go into Methods** (colleague's wording, verbatim): "We considered
all 2^4=16 combinations of the four directionally uncertain relationships. Twelve
combinations produced acyclic graphs and were retained for the directional-sensitivity
analysis. Four produced directed cycles and were excluded because both the working SCM
and the causal structures represented by FCI assume an acyclic causal system. Each of
the 12 retained specifications was refitted and the intervention analysis was repeated."

**Adopted framing, Supplement**: "The four excluded specifications were not treated as
failed models. They implied directed feedback cycles and therefore fell outside the
recursive SCM used for intervention analysis and the acyclic causal-graph framework
assumed by FCI." Plus a small supplementary table listing which four combinations
generated cycles, for transparency -- **no new script needed for this**:
`clean_pipeline/10_orientation_enumeration_fit.R` already writes one row per combo,
including the 4 excluded ones, tagged `status="cyclic_skipped"` with their exact
flipped-edge labels (`orientation_enumeration_fit.csv`) -- just filter to those 4 rows.
The copy of that file currently on disk (2026-09-09) predates the reversal and still
shows all 16 as "ok" with the old flip labels; needs a rerun of `10` post-reversal
before those 4 rows reflect the current model.

**Why all 4 excluded combos share one root cause, not four independent ones** (worth a
sentence in the Supplement rather than leaving it looking like 4 unrelated failures):
all four have `belief_concern -> weather_risk_prep` flipped while `harm_present ->
weather_risk_prep` stays unflipped -- combined with the fixed (non-flippable)
`belief_concern -> harm_present` edge, this always closes the 3-cycle `belief_concern
-> harm_present -> weather_risk_prep -> belief_concern`. The other two flip candidates
(`social_norms <-> climate_behavior`, `harm_future <-> harm_present`) vary freely
across the four excluded combos without affecting cyclicity -- confirmed against
Kyuri's real run (combo_1/5/9/13 exactly, matching this structural account bit-for-bit).

**Explicitly NOT pursued for this paper, logged for the record as a future-paper
direction**: a linear feedback treatment (X = BX + eps, equilibrium X = (I-B)^{-1}eps
if stable) was discussed as technically possible but requires, at minimum: the cyclic
equations be statistically identified (an instrument or covariance constraint per loop
variable), a unique and dynamically stable equilibrium (checked via the loop's path-
coefficient product having |product| < 1 -- worked through informally in chat using
raw correlations as a quick proxy, not yet done as a rigorous check), and a
substantive argument that an equilibrium interpretation is appropriate for these
attitude constructs. That's a second causal model class, not four more rows on the
current table -- explicitly out of scope here per Kyuri's colleague's review, noted as
a natural fit for future work given other feedback-systems interests.

Not yet touched, per standing instruction: `main2.tex` (the Methods/Supplement text
above is ready to paste in whenever Kyuri does that pass).

## Section 38 (2026-09-11): Interaction analysis/S9 removed entirely; theory-overridden/uncertain mixture retired; manuscript update order locked

**Superseded**: the "rebuild Part D using observed politics as a real subgroup" plan
from earlier the same day (Section 37 follow-up, chosen via the 4-option question
about Part D's redesign) is superseded by this decision. The two small uncommitted
edits made toward that rebuild (`edge_abbrev`'s two rows renamed to the post-reversal
direction/abbreviation; the dead `scm_resid_sd` block removed) are moot now that the
whole file is gone, but are harmless -- the file's entire content was replaced anyway.

**Decision, in Kyuri's own words (2026-09-11), locking the revised analysis**:
- Follow the data-supported direction whenever it meets the orientation criterion:
  keep `belief_concern -> politics` and `social_norms -> policy_support` as currently
  implemented, rather than overriding either back to the old theory-asserted direction.
- Do not call theory-overridden edges "uncertain." That whole three-way mixture
  (bootstrap-resolved / theory-asserted-but-now-also-bootstrap-supported /
  theory-asserted-and-contradicted-by-bootstrap) is retired for framing purposes.
- Directional sensitivity varies only the four genuinely weak/inconsistent directions:
  `belief_concern<->weather_risk_prep`, `harm_present<->weather_risk_prep`,
  `social_norms<->climate_behavior`, `harm_future<->harm_present`.
- Remove the belief x politics moderation analysis and Supplementary Section S9
  entirely -- with belief/concern now upstream of politics, it no longer contributes
  cleanly to the paper (politics is not exogenous to belief_concern anymore, and
  policy_support is a sink node, so the interaction has no downstream route to
  climate_behavior regardless).
- Figure 3 (`fig:scm`, `figures/fig5_scm_hierarchical.pdf`) and Table 2
  (`tab:orientation`) get a matching conceptual simplification: solid arrows /
  "Data-supported direction" = the 12 directions retained from the working SCM;
  dashed arrows / "Directionally weak/inconsistent" = the four varied in sensitivity
  analysis. No third "theory-overridden" category.
- Update order once the clean baseline refit is done: baseline coefficients/effects ->
  Figure 3/Table 2 -> directional sensitivity/Figure 4 -> confounding sensitivity (if
  affected by the new baseline) -> cyclic feedback extension (proposed Supp. S15,
  Section 37) -> manuscript prose/Supplement -> one full clean-pipeline run.

**What's done as of this entry**:
- `clean_pipeline/11_interaction_moderation.R` replaced in full with a short
  deprecation stub explaining the removal and pointing here; left in place (not
  `git rm`'d) so the rationale stays attached to a real file at its old path.
  Confirmed via grep across `clean_pipeline/` and `r_patches/` that nothing else reads
  its outputs (`interaction_coefficient.csv`, `interaction_shift_results_montecarlo.csv`,
  `interaction_shift_results_exact.csv`) or sources the file itself -- safe to remove
  outright rather than just gut Part D.
- `clean_pipeline/run_all.R`'s exclusion of `"11_interaction_moderation.R"` from
  `scripts_in_order` changed from "temporarily excluded pending a decision" to
  "permanently removed," comment updated accordingly.
- `r_patches/07_figure5_scm_hierarchical_v3.R`: checked `pipeline_outputs/
  scm_edges_finalized.csv`'s `final_tier` column directly (already regenerated
  post-reversal) -- it already classifies exactly the 4 flip-candidate edges as
  `substantive` (dashed) and all other 12 as `data_aligned` (solid), so the
  line-style *logic* already matched the new rule with no code change needed. Only
  the legend text was stale: "Data-supported orientation" / "Data-uncertain /
  theory-completed" -> "Retained from working SCM" / "Varied in orientation-
  sensitivity analysis." Same header-comment cleanup (two spots) removing the
  "theory-completed"/"picked on substantive grounds" framing.

**Not yet done, and deliberately deferred per the update order above**:
- Figure 3's actual redraw (rerun `r_patches/07_figure5_scm_hierarchical_v3.R` in R --
  the legend/label logic is fixed now, but Kyuri needs to run it to regenerate the PDF)
  and Table 2's rebuild with the two-category "Basis" scheme are next in line, but
  Table 2 itself is hand-written LaTeX (`main2.tex` lines ~264-289, `tab:orientation`),
  not script-generated -- rewriting it is a manuscript-prose edit, held for the
  manuscript step per the order above, not done now.
- Confounding-sensitivity rerun check (`17_edge_confounding_classification.R`,
  `29_joint_pag_edgetype_audit.R`) -- both already carry the post-reversal edge
  directions (Section 35's fixes), so likely nothing further needed here, but not
  re-verified against a fresh end-to-end run this entry.
- The cyclic feedback extension (proposed Supp. S15) -- not built yet.
- Manuscript prose/Supplement pass -- **not started, per standing instruction not to
  touch `main2.tex` without explicit sign-off, and because it's explicitly last in
  Kyuri's own stated order.** Located for when that pass happens:
  - Abstract/intro framing that references "theory-augmented"/"theory-completed"
    edges as a class: `main2.tex` lines ~30-31, ~70-72.
  - Figure 3: `main2.tex` lines ~239-268 (`fig:scm`, sources
    `figures/fig5_scm_hierarchical.pdf`); caption currently describes the old
    bootstrap-vs-theory two-line-style scheme and needs the retained/varied wording
    to match the figure once it's redrawn.
  - Table 2: `main2.tex` lines ~264-289 (`tab:orientation`) -- five hand-written
    evidential-basis groups (bootstrap-aligned, bootstrap+theory-aligned,
    theory-then-bootstrap-supported, theory-only, threshold-sensitive, weak/symmetric)
    collapse to two: "Data-supported direction" (12 rows) / "Directionally
    weak/inconsistent" (4 rows, the flip candidates) -- footnotes referencing the old
    categories (`$^{\dagger}$`, `$^{\ddagger}$`, etc.) will need rewording or dropping.
  - Results prose immediately after Table 2 (`main2.tex` line ~271) restates the old
    grouping in words and needs the same simplification; also still has pre-reversal
    baseline numbers (`.206`/`.195`/`.093`/`.034`/`.057`/`.016` at lines ~322, ~362) --
    these are stale versus the confirmed post-reversal values already locked into
    `run_all.R`'s `expected_baseline` (Section 35): `.201`/`.187`/`.093`/`.044`/`.063`/
    `0`. Same section (`\subsection{Intervention implications}` onward) also still
    reports old pair/triple bootstrap percentages (`38.3%`/`37.8%`/`23.9%`) which
    should be re-checked against the current `intervention_bootstrap_pair_rank1_freq.csv`
    (`45.5%`/`38.8%`/`15.7%`, confirmed this session) before pasting in.
  - Supplementary Section S9 itself: `main2.tex` line 1168
    (`\subsection*{S9. Belief/concern $\times$ political orientation interaction}`)
    through line 1264 (next `\subsection*{S10...}` at line 1265) -- delete this whole
    span, then renumber S10-S14 down by one (or repoint labels/refs if numbering is
    manual rather than automatic; not yet checked which).
  - Also line 351 (main-text sentence reporting the interaction's point estimate and
    pointing to `\nameref{sec:interactionappendix}`) and any other inline
    cross-references to S9/`sec:interactionappendix`/`supp:interaction` -- not yet
    fully enumerated beyond this one hit; a fresh grep for `interactionappendix\|
    supp:interaction\|S9\b` right before the prose pass will catch stragglers,
    including any in `main9.tex` if that's a separate/parallel file.
  These are recorded here as a checklist for that pass, not acted on.

**Uncommitted at the time of this entry**: this file plus
`clean_pipeline/11_interaction_moderation.R`, `clean_pipeline/run_all.R`, and
`r_patches/07_figure5_scm_hierarchical_v3.R` -- committing together right after this
entry is written.
