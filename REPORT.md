# Geographic origin of KLMA Yellow Rail feathers: stable isotope assignment

## Question
Do body feathers and flight feathers from Yellow Rails (*Coturnicops
noveboracensis*) captured at Klamath Marsh NWR, Oregon in July grow in
different places, consistent with Cornell's account of the species' molt
strategy (flight feathers replaced in a complete molt primarily on/near the
breeding grounds; body feathers replaced in a molt that *can* continue into
the wintering period, though this is not guaranteed)? And does habitat
(wetland/wet agriculture only, per Cornell's habitat account) help narrow
down where that is?

## Major revision (this version)

This version follows an external code/results review (conducted via
GitHub Copilot, on the previous version of this pipeline and its outputs)
that identified several analytical issues materially affecting the
biological conclusions. Two changes were made as a direct result, both
requested by the study's lead author after reviewing that feedback:

1. **The two non-KLMA individuals ("San Diego ref bird" and `1272-31742`)
   are removed from the analysis entirely**, not merely excluded from
   population-level averages as in the previous version. They were
   previously framed as "reference individuals," but that framing doesn't
   hold up: neither was captured as part of the KLMA sampling effort,
   their capture dates relative to the KLMA birds are unknown, and -- most
   importantly -- there is no independent knowledge of where either bird
   actually bred or molted for their isotope-based assignment to be
   checked against. (The external review flagged the San Diego bird's
   assignment to the Pacific Northwest rather than the California coast as
   a "red flag" for exactly this reason: without knowing where that
   feather was actually grown, an assignment result can't be judged right
   or wrong, so it wasn't functioning as a reference at all.) The
   consolidated lab spreadsheet now used as this pipeline's data source
   (see "Data" below) is already scoped to just the 12 KLMA individuals,
   so there is no longer a separate excluded-records file to maintain --
   these two individuals simply never enter the pipeline.
2. **d2H alone (not the dual d2H+d18O model) is now the primary
   assignment model.** The previous version used the dual-isotope model as
   primary and reported the d2H-only model only as a sensitivity check.
   That was backwards: d18O's calibration here is weak (R^2 = 0.08) and
   its raw isoscape is r^2 = 0.99 correlated with d2H's across this study
   area -- a condition `assignR::pdRaster()` itself warns about
   ("specificity of assignments may be inflated"). Section 1b below
   quantifies the consequence directly: switching between d2H-only and
   dual-isotope changes the population-level body-feather result from
   Pacific-Northwest-majority to Mexico-majority. A result that depends on
   including a weak, highly collinear second marker should not be the
   default model. The dual-isotope model is still run and reported, now
   explicitly as one of several sensitivity checks rather than as the
   primary result.

This version also adds, directly from the external review's remaining
recommendations: a habitat-prior-strength sensitivity check (weak/current/
strong), a calibration-reference-group sensitivity check (Passerine vs.
Water bird), and a check of whether the habitat prior itself is balanced
between the Mexico and Pacific Northwest zones (independent of any isotope
data). See "Results, section 1b" for all of these.

## Data
- `data/SIA_Results_KLMA.xlsx` (H/O isotopes, SIRFER lab): the consolidated
  lab spreadsheet, read directly by `scripts/00_clean_data.R`. This
  replaced an earlier two-file version of this analysis (a separately
  exported, decimal-rounded "clean" results file cross-referenced against
  a raw QC text export) after the lead author pointed out the two-file
  split was error-prone -- it was how a non-KLMA individual (`1272-31742`)
  ended up silently included in this pipeline for a while, recovered only
  by text-regex from the raw QC file rather than being a deliberate
  inclusion. The consolidated file is already scoped to the 12 KLMA
  individuals only (no San Diego bird, no `1272-31742`), and gives full
  floating-point isotope precision instead of the previous file's values
  rounded to ~1 decimal place. `scripts/00_clean_data.R` includes a sanity
  check that stops with an error if any non-KLMA envelope ID is ever found
  in this file, rather than silently filtering it out.
- Switching to full precision changed population-level results
  negligibly (all changes <0.2 percentage points; see change log) but did
  flip one individual bird's body-feather result across the 50% Mexico
  threshold (`1422-02214`: 49.7% -> 50.5%, now a marginal rather than
  clear non-majority) -- a reminder that individual-level results this
  close to a round-number threshold were never meaningfully different in
  the first place.
- `data/SIBS_CN_230621_Kellerman.csv` (C/N isotopes, separate SIBS lab,
  same envelope-ID numbering), bias-corrected using the lab's own QA
  standards in `data/SIBS_CN_230621_Kellerman_QC.xlsx` (see "Carbon/nitrogen
  isotopes" in Results). Individual `1272-31742`'s C/N samples are
  excluded here too, for the same reason.
- Cleaned/merged tables: `data/klma_clean.csv` (`scripts/00_clean_data.R`,
  43 feather samples from 12 KLMA individuals), `data/cn_clean.csv` /
  `data/isotopes_combined_summary.csv` (`scripts/04_clean_cn_data.R`).
- Three samples were QC-flagged during cleaning: a Nape sample from bird
  1422-02212 and a Nape sample from bird 1422-02214 (both implausibly high
  Wt% H/O, possible double-loading), and one Rectrix from 1422-02220 (O:H
  combustion ratio 8.5, far outside the normal 3.6-4.5 range). Separately,
  the C/N lab's own QA standards showed a consistent ~+2.05 permil
  enrichment bias for d13C specifically (not d15N); this is now corrected
  in the reported values -- see "Carbon/nitrogen isotopes" in Results.
- Growing-season precipitation isoscapes `d2h_GS.tif` / `d18o_GS.tif` (mean
  surfaces only, no prediction-SE grid -- see Limitations).
- `data/habitat_prior.tif`: wetland/wet-cropland habitat suitability, built
  from ESA WorldCover 10m land cover (`scripts/05_build_habitat_prior.R`).

## Methods
All analysis in R (`scripts/00`-`06`) using **assignR** (installed from
GitHub source, `SPATIAL-Lab/assignR`) plus `terra`/`sf`. True run order is
00 -> 05 -> 01 -> 02 -> 03 -> 04 -> 06 (see `METHODS_WALKTHROUGH.md`).

1. **Primary model: d2H only.** assignR's built-in known-origin database
   contains **no Rallidae** and essentially no North American
   non-passerine isotope data, so a Yellow-Rail-specific calibration isn't
   possible. The North-America-masked **Passerine** group (n=699 sites) is
   used as the best available proxy for the primary model (R^2 = 0.76,
   slope = 0.71, p < 2e-16). d18O is deliberately **not** included in the
   primary model -- see "Major revision" above and Results section 1b for
   why, and for the dual-isotope, Water-bird-calibration, and
   habitat-prior-strength variants that are still run as named sensitivity
   checks rather than folded into one default result.
2. **Isoscape uncertainty.** The supplied raster has no SE band (the
   authoritative source, wateriso.utah.edu, wasn't reachable from this
   sandboxed environment). A spatially uniform SE was assumed (5 permil
   for d2H), sized against assignR's own low-res example SE raster.
3. **Study area.** CA/NV/AZ/NM south through Baja California and the
   western Mexican coast (the original frame of reference), extended
   north to include OR/WA/ID/UT so Klamath Marsh itself (~43 N) remains a
   candidate flight-feather origin.
4. **Habitat suitability prior.** Yellow Rails occupy wet sedge meadow /
   marsh / wet agriculture only (Cornell's Birds of the World account) --
   a suitability layer was built from ESA WorldCover 2021 10m land cover
   (herbaceous wetland = 1.0, cropland = 0.3, permanent water = 0.1,
   everything else = 0), aggregated to the isoscape's ~9 km grid, Gaussian-
   smoothed, and supplied to `assignR::pdRaster()` as a multiplicative
   `prior`. 85 of 108 candidate 10m tiles across the study area had land
   data; after smoothing, 101,813 of 1,209,600 cells (~8.4%) across the
   full North American grid are nonzero (most of the continent stays
   exactly zero); within the narrower study area the smoothing kernel's
   finite radius means almost every cell has *some* nonzero value, but the
   distribution stays heavily skewed toward near-zero (median ~0.003; only
   about a third of study-area cells exceed 0.01).
   **This required a workaround for a real bug in assignR itself**:
   `pdRaster.isoStack()` crops the isoscape to `mask` internally but never
   crops `prior` to match, so the elementwise `assign * prior` step
   silently recycled two differently-sized vectors. The fix was to
   pre-crop/mask both the isoscape and the prior to the study area
   ourselves and pass `mask = NULL` to `pdRaster()`.
   Sanity check: Mt. Whitney (4,421 m alpine peak, no habitat) samples at
   ~8.5x10^-15 probability in the population-body surface; Sacramento
   National Wildlife Refuge (a real Central Valley wetland complex)
   samples at ~6.2x10^-4 -- more than ten orders of magnitude apart,
   confirming the prior is doing its job. **New this version**: the
   prior's own balance between the two zones under comparison was also
   checked directly, independent of any isotope data (Results 1b) --
   39.3% of total habitat-suitability mass falls in the Mexico zone versus
   49.6% in the Pacific Northwest zone, so the mask is not dramatically
   skewed toward either region.
5. **Assignment.** `assignR::pdRaster()` on all 43 feather samples (12
   KLMA individuals only), normalized to the study area.
6. **Combining feathers per bird.** Where a bird had >1 feather in the same
   category, posteriors were multiplied and renormalized (shared-origin
   assumption). Where that produced no shared support at all, the mean of
   the individual posteriors was used instead and flagged (`combine_method`
   column in the summary table).
7. **Population summary.** Equal-weighted mean of KLMA individuals'
   posteriors per category.
8. **Credible regions.** `assignR::qtlRaster(thresholdType = "prob")` --
   the smallest area containing 50/75/90% of the cumulative posterior
   probability.
9. **C/N isotopes** (`scripts/04_clean_cn_data.R`) were cleaned and
   cross-referenced against the H/O results, interpreted qualitatively
   (there's no established continental d13C/d15N isoscape product at
   usable resolution the way there is for d2H/d18O).
10. **Sensitivity analysis** (`scripts/06_sensitivity_check.R`, new this
    version): reruns the assignment under five alternative configurations
    -- dual-isotope, no habitat prior, weak habitat prior (prior^0.5),
    strong habitat prior (prior^2), and Water-bird calibration -- and
    compares each against the primary model. See Results 1b.

## Results

### 1. Flight feathers robustly point to the Pacific Northwest; body feathers are genuinely unresolved between the Pacific Northwest and Mexico

| | Body feathers | Flight feathers |
|---|---|---|
| Population probability, Mexico (Baja + mainland), primary model | 34.7% | 0.8% |
| Population probability, Pacific NW (OR+WA+CA+ID+UT), primary model | **59.5%** | **98.2%** |
| Range across all 6 model variants (Results 1b), % Mexico | 27.9-53.1% | 0.01-4.0% |
| Range across all 6 model variants, % Pacific NW | 36.6-64.1% | 74.1-99.4% |
| 50%-credible-region area, population summary, primary model | 112,733 km^2 | 27,858 km^2 |

Under the primary model (d2H only, Passerine calibration, current habitat
prior), **both** body and flight feathers have a population-level majority
probability in the Pacific Northwest -- but the two are still clearly
different: flight feathers are far more concentrated there (98.2%) than
body feathers (59.5%), and this gap (flight consistently far more
Pacific-Northwest-concentrated than body) holds across **every one of the
six model variants tested** (see Results 1b), making it the most robust
finding in this dataset. **Individual-level credible areas are large and
mostly not trustworthy at face value under this primary model -- see
"Match quality" below.**

- **Flight feathers** are robustly assigned to the Pacific Northwest
  interior/coast under every model variant tested (74.1-99.4% across all
  six), consistent with the complete prebasic molt of remiges happening
  primarily on/near breeding habitat (Cornell's account).
- **Body feathers** are not robustly assigned to any single zone. Under
  the primary model they lean Pacific Northwest (59.5% vs. 34.7% Mexico);
  under the dual-isotope and Water-bird-calibration sensitivity variants
  they lean Mexico (53.1% and 50.5% respectively); under the
  habitat-prior-strength variants (weak/strong) the result stays close to
  the primary model's (33.5-34.7% Mexico) regardless of how strongly the
  habitat mask is weighted. In other words: **the habitat prior's
  weighting is not what's driving the uncertainty in the body-feather
  result -- the choice of isotope marker(s) and calibration reference
  group is.** This is consistent with Cornell's account that body
  (preformative) molt *can* continue into the wintering period but isn't
  guaranteed to -- body feathers may represent breeding grounds,
  migration stopover sites, or wintering areas depending on age and molt
  completion timing, so an unresolved, mixed population-level signal is
  biologically plausible, not just a modeling artifact.
- **Unlike the previous (dual-isotope) version, no feather-combination
  cases are flagged as isotopically inconsistent under the primary model**:
  all 20 multi-feather combinations (two birds had only a single feather
  in one category and needed no combination) succeeded via the standard
  joint-product method, with no `inconsistent_feathers_mean_fallback`
  cases. This makes sense mechanically -- a single, less spatially
  discriminating marker (d2H alone) produces broader individual-feather
  posteriors that are more likely to numerically overlap somewhere in the
  study area, even for feathers that would have looked "inconsistent"
  under the sharper dual-isotope posteriors. This is not evidence that the
  shared-origin assumption is better supported under the primary model --
  only that it's less often detectably violated by this particular check.

### 1b. Sensitivity analysis: six model variants, compared directly

Six configurations were run and compared, individual-by-individual and at
the population level (`scripts/06_sensitivity_check.R`,
`outputs/tables/sensitivity_comparison.csv`):

| Variant | What changed from primary | Body: % Mexico | Body: % Pacific NW | Body: 50%-CA (km^2) |
|---|---|---|---|---|
| **Primary** | -- | 34.7 | 59.5 | 112,733 |
| Dual-isotope | + d18O | 53.1 | 36.6 | 68,257 |
| No habitat prior | prior removed entirely | 27.9 | 40.8 | 1,176,595 |
| Weak habitat prior | prior^0.5 (compressed) | 34.3 | 55.0 | 261,107 |
| Strong habitat prior | prior^2 (exaggerated) | 33.5 | 64.1 | 46,927 |
| Water-bird calibration | Passerine -> Water bird ref. group | 50.5 | 42.2 | 111,772 |

(Flight-feather numbers across the same six variants range only 74.1-99.4%
Pacific NW and 0.01-4.0% Mexico -- stable regardless of variant, and
omitted from the table above for space; see the CSV for full figures.)

Three conclusions follow directly from this table:

1. **Habitat-prior strength has comparatively little effect on the
   body-feather zone split**: a four-fold range in prior strength (weak to
   strong) moves % Mexico by only about 1 percentage point (33.5-34.7%),
   all close to the primary model's 34.7%. It has a much larger effect on
   credible-region *size* (46,927 to 261,107 km^2) -- the prior mainly
   controls precision, not which zone wins.
2. **Removing the habitat prior entirely drops Mexico probability further
   (27.9%) and explodes the credible area to 1.17 million km^2** --
   confirming the prior is necessary for any usable precision, even though
   it isn't the main driver of the Mexico-vs-Pacific-NW split itself.
3. **The isotope marker and calibration-group choices are what actually
   move the body-feather conclusion**: switching to dual-isotope or to a
   Water-bird calibration reference group both push the body-feather
   result to a Mexico-majority (53.1%, 50.5%), while the Passerine-based,
   habitat-prior-strength variants all agree on a Pacific-Northwest-leaning
   result (55.0-64.1%). Four of six variants lean Pacific NW for body
   feathers; two (both involving either d18O or a different calibration
   reference group) lean Mexico. This is genuinely unresolved with the
   present data -- not a matter of picking the "right" prior weighting.

**Habitat-prior zone balance** (independent of any isotope data,
addressing whether the mask itself could be biasing results toward one
region): 39.3% of the habitat prior's total suitability mass falls within
the Mexico zone, versus 49.6% within the Pacific Northwest zone (the
remainder is in Nevada/Arizona/New Mexico); suitable area (prior > 0.01)
is 459,286 km^2 in Mexico versus 512,430 km^2 in the Pacific Northwest.
The mask is not dramatically skewed toward either region, so the
difference in results between variants is not simply an artifact of
Mexico habitat being underrepresented in ESA WorldCover.

### 2. Carbon/nitrogen isotopes: an independent habitat cross-check

**d13C values reported here are bias-corrected** (new this version). The
SIBS lab's own QA standards (`data/SIBS_CN_230621_Kellerman_QC.xlsx`) show
d13C reading a consistent ~+2.05 permil too enriched relative to their
known consensus values (SIBS-pCo: +2.10 permil; SIBS-pEc: +2.01 permil --
consistent to within 0.09 permil despite the two standards sitting ~27
permil apart, so a simple additive correction is justified). An earlier
version of this report *noted* this bias in passing but never actually
corrected the data for it; `scripts/04_clean_cn_data.R` now applies the
correction (subtracting 2.053 permil, the mean of the two standards'
observed bias) to every sample, and both the raw and corrected values are
kept in `data/cn_clean.csv` for transparency. d15N showed no comparable
bias (+0.31, +0.15 permil on the two standards -- within normal
measurement noise) and was not corrected.

- **d13C is solidly C3 across every KLMA bird** (-27.6 to -24.3 permil,
  corrected) -- nowhere near the C4 range (roughly -12 to -17 permil) that
  would flag *Spartina*/*Distichlis* salt-marsh vegetation. The correction
  moves every value *further* from the C4 range, not closer, so this
  conclusion is if anything more solidly supported than before. These
  feathers reflect freshwater sedge-marsh/wet-meadow-type diets, not
  brackish salt marsh, for whichever molt event each feather represents.
- **d15N does not reliably track the H/O-based Mexico probability.** This
  was checked directly against the data in an earlier version of this
  report and did not hold up even before the primary-model change (one
  bird, 1422-02216, directly contradicted the pattern: low d15N, high
  Mexico probability). With the primary model now generally assigning
  lower Mexico probabilities to begin with, and only 3 KLMA individuals
  having both C/N and geographic-assignment data for body feathers, this
  comparison remains too small to support a trend claim either way.
- One bird (1422-02220) had four body-feather C/N subsamples split sharply
  between freshwater-typical (d15N = 8.1) and elevated (d15N = 14.9-15.3)
  -- another individual whose different body feathers evidently didn't
  grow together, on top of the H/O-flagged cases above.

## Outputs
- `outputs/maps/population_body_vs_flight.png` -- the main figure (primary
  model: d2H only, current habitat prior).
- `outputs/maps/habitat_prior.png` -- the wetland/wet-cropland suitability
  layer itself.
- `outputs/maps/individuals_body.png`, `individuals_flight.png` -- main
  per-bird figures (unflagged birds only; under the primary model, this is
  a small minority of birds -- see Match quality).
- `outputs/maps/individuals_body_supplementary_flagged.png`,
  `individuals_flight_supplementary_flagged.png` -- per-bird figures for
  flagged birds.
- `outputs/maps/calibration_regression_*.png` -- calibration diagnostics.
- `outputs/tables/assignment_summary.csv` -- per-individual and population
  peak location, 50/75/90% credible areas, probability mass by zone, and
  match-quality columns (`peak_share`, `low_precision_flag`,
  `single_cell_dominant_flag`, `shared_peak_flag`).
- `outputs/tables/sensitivity_comparison.csv` -- the full 6-variant
  comparison behind Results section 1b, every individual and population
  surface.
- `outputs/tables/habitat_prior_by_zone.csv` -- the habitat-prior
  zone-balance check behind Results section 1b.
- `data/isotopes_combined_summary.csv`, `data/cn_vs_geography_comparison.csv`
  -- C/N results cross-referenced against H/O values and the geographic
  assignment, respectively (KLMA individuals only).
- `outputs/rasters/*.tif` -- full posterior surfaces, for further GIS work.
- `data/habitat_prior.tif` -- the habitat suitability layer, for reuse.

## Match quality: not every tight credible region means a precise result

This issue, first found in an earlier (dual-isotope-primary) version of
this analysis, is **more severe under the current (d2H-only) primary
model**, not less: 17 of 22 individual-category surfaces (77%) now carry
at least one match-quality flag, up from 14 of 26 (54%) in the
dual-isotope version. This is an expected and honestly-reported trade-off,
not a new bug: d2H alone is less spatially discriminating than d2H+d18O
(that's exactly why the two isotopes gave different population-level
results in section 1b), so more individuals' weakly-informative likelihood
surfaces default to whichever single cell the habitat prior favors most
locally -- and because the habitat mask is sparse, many unrelated
individuals converge on the *same* cell. The practical upshot is unchanged
from before: **population-level zone probabilities (section 1 and 1b) are
the trustworthy output of this analysis; individual-level peak locations
and credible areas, under the primary model, mostly are not**, and should
not be cited as precise results without checking the flag columns in
`assignment_summary.csv` first.

The three flag columns (unchanged in definition from the previous
version): `peak_share` (fraction of a surface's probability mass in its
single top cell), `single_cell_dominant_flag` (peak_share > 0.3 -- 0 of 22
surfaces under the primary model, since d2H alone never concentrates that
much probability in one cell), `shared_peak_flag` (peak cell identical to
another individual's -- all 17 flagged surfaces are flagged for this
reason alone), `low_precision_flag` (50%-credible area < 1,000 km^2 -- 0
of 22 surfaces; the precision problem under d2H-only looks different,
showing up as
large shared credible regions rather than small ones).

## Key limitations (read before over-interpreting)
1. **No Rallidae calibration data exist anywhere** in the standard
   known-origin isotope databases used by assignR. The Passerine-based
   calibration is the best available proxy, not a validated one for rails.
   Section 1b shows this choice matters: switching to a Water-bird
   reference group changes the body-feather population result from
   Pacific-Northwest-majority to Mexico-majority.
2. **d18O is weakly calibrated (R^2 = 0.08) and nearly collinear with d2H
   (r^2 = 0.99 in the raw isoscape)** -- this is why it was removed from
   the primary model (see "Major revision"). Section 1b shows its
   inclusion alone is enough to flip the body-feather population result.
3. **Isoscape SE was assumed, not measured.** This is a real limitation,
   but it is not an analytical-precision problem: the SIRFER lab's own
   per-batch secondary-standard (POW) replicate measurements show d2H
   standard uncertainty of 0.76-1.20 permil across the two job batches
   used here -- far tighter than the calibration regression's own
   residual scatter (~18 permil) or the assumed 5 permil isoscape SD.
   Measurement noise is not what's limiting this analysis.
4. **The habitat prior is approximate**: a single year (2021) of 10m land
   cover, a hand-picked (not statistically fit) weighting scheme, and it
   cannot distinguish flooded/irrigated cropland (real Yellow Rail
   habitat) from dry cropland (not habitat). However, section 1b's
   prior-strength sensitivity check found this weighting scheme is *not*
   the main driver of the body-feather result's instability -- a
   four-fold range in prior strength moved the result by only ~1
   percentage point.
5. **No age data.** Cornell's account makes age class (HY/juvenile vs.
   AHY/adult) important, and separately, the preformative (body) molt's
   timing relative to migration/wintering is itself variable (see Results
   section 1). Without age class, neither the body-vs-flight contrast nor
   the specific body-feather result can be tied to a single molt event
   with confidence.
6. **C/N isotopes are a qualitative cross-check, not folded into the
   spatial model**, and the sample available for cross-referencing against
   geography is very small (3 KLMA individuals with both measurements for
   body feathers).
7. **Sample size** is small (12 KLMA individuals only, as of this
   version). Several feathers per bird were pooled under a shared-origin
   assumption (20 multi-feather combinations); under the primary model
   this assumption was not detectably violated for any individual (see
   "Combining feathers" above), though this check is weaker under d2H
   alone than it was under the dual-isotope model, where 4 combinations
   failed outright.
8. **Precision implied by reported percentages.** Given limitations 1-3
   above, percentages in this report (e.g. "59.5% Pacific Northwest") are
   reported to one decimal place for traceability against the output
   tables, not as a claim of that level of real precision. Readers should
   treat the ranges in section 1b's table, not any single percentage, as
   the honest summary of what the data support for body feathers
   specifically.

## Suggested next steps
- **Highest priority, unresolved**: the body-feather Mexico-vs-Pacific-
  Northwest question cannot be resolved with the present data and
  reference materials -- it depends on which of two defensible but
  different modeling choices (isotope marker set, calibration reference
  group) is used. Resolving it would need either a validated Rallidae- or
  wetland-bird-specific calibration, or an isoscape uncertainty product
  that properly models the d2H/d18O covariance instead of treating the
  two as independent evidence.
- Obtain the official prediction-SE isoscape grids from waterisotopes.org
  (blocked from this sandboxed environment) and rerun with real,
  spatially-varying SE.
- If plumage-based ages are assignable, rerun stratified by HY vs AHY.
- If additional truly independent reference individuals (known capture
  date and context matching the KLMA birds, ideally with some independent
  knowledge of origin) become available, re-introduce them properly rather
  than treating unvalidated specimens as references.
- Consider re-running the analysis excluding/re-measuring 1422-02212 given
  its QC flag, and treating 1422-02219/1422-02251 as reflecting genuinely
  discordant feathers rather than combining them.

## Change log
- Initial version: assignR pipeline, no habitat mask.
- Correction 1: fixed `qtlRaster` credible-region math and the
  Baja/mainland Mexico zone split; fixed map color scaling and missing
  legends.
- Correction 2: added the ESA WorldCover habitat suitability prior
  (fixing a real assignR bug in the process), added two individuals as
  "reference" birds (since removed -- see below), added the C/N
  cross-check.
- Correction 3: map colors classified directly by credible region instead
  of by quantile-of-pixel-value.
- Correction 4: found and partially mitigated the single-cell convergence
  issue described in "Match quality" above.
- Correction 5: fixed several numbers that had drifted out of sync with
  the pipeline (QC-flag description, inconsistent-feather count, habitat-
  prior percentages, Mt. Whitney/Sacramento NWR sanity-check values).
- Correction 6: full pipeline rerun confirmed byte-for-byte determinism;
  found and quantified the d2H/d18O collinearity sensitivity issue (see
  former Results 1b, now superseded by this version's expanded section
  1b).
- **This version (major revision)**: following external code/results
  review, (a) removed the two non-KLMA "reference" individuals entirely,
  since they were never true reference individuals; (b) made d2H alone
  the primary model instead of dual-isotope, moving dual-isotope to a
  named sensitivity check; (c) added habitat-prior-strength and
  calibration-reference-group sensitivity checks; (d) added a
  habitat-prior zone-balance check; (e) rewrote the Results and
  Limitations sections around the finding that the body-feather
  Mexico-vs-Pacific-Northwest question is genuinely unresolved by the
  present data, while the body-vs-flight contrast itself remains robust
  across every model variant tested.
- **This version (data-source upgrade)**: switched from the original
  two-file input (a decimal-rounded "clean" results spreadsheet plus a
  separately matched raw QC text export) to a single consolidated lab
  spreadsheet (`data/SIA_Results_KLMA.xlsx`) with full floating-point
  isotope precision and already scoped to the 12 KLMA individuals.
  Re-ran the full pipeline: every population-level result changed by
  less than 0.2 percentage points, confirming the earlier rounding
  wasn't masking anything. One individual (`1422-02214`) crossed the 50%
  Mexico-probability threshold (49.7% -> 50.5%), changing the
  individual-majority count for body feathers from 1 of 10 to 2 of 10 --
  flagged explicitly as a marginal case, not read as a new finding. Also
  confirmed via this file's per-batch QA/QC that analytical (instrument)
  precision is not a limiting factor in this analysis (POW standard
  uncertainty 0.76-1.20 permil for d2H, far tighter than the ~18 permil
  calibration-regression residual scatter that dominates the overall
  uncertainty budget).
- **This version (C/N bias correction)**: the C/N lab's own QA-standard
  file (`data/SIBS_CN_230621_Kellerman_QC.xlsx`, reattached by the lead
  author alongside the original C/N data when asked whether the C/N
  results could be strengthened) was used to implement an actual
  correction for the ~+2.05 permil d13C enrichment bias that an earlier
  version of this report had only noted in passing. `scripts/04_clean_cn_data.R`
  now parses the QA-standard values directly and applies the correction
  to every sample (raw values kept alongside corrected ones in
  `data/cn_clean.csv` for transparency). Corrected d13C values (-27.6 to
  -24.3 permil) are further from the C4 salt-marsh signature than the
  uncorrected values were, so the existing "solidly C3" conclusion is
  strengthened, not changed in direction. This does not add any new
  individuals to the C/N-vs-geography cross-check (still only 3 KLMA
  individuals with both measurement types for body feathers) and does
  not change the Mexico-vs-Pacific-Northwest determination -- there is
  still no continental d13C/d15N isoscape to convert a diet signal into
  a geographic probability the way there is for d2H.
