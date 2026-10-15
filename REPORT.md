# Geographic origin of KLMA Yellow Rail feathers: stable isotope assignment

## Question
Do body feathers and flight feathers from Yellow Rails (*Coturnicops noveboracensis*)
captured at Klamath Marsh NWR, Oregon in July grow in different places, consistent
with Cornell's account of the species' molt strategy (flight feathers replaced in a
complete molt primarily on/near the breeding grounds; body feathers replaced in a
molt that can complete on the wintering grounds)? And does habitat (wetland/wet
agriculture only, per Cornell's habitat account) and diet (carbon/nitrogen isotopes)
help narrow down and cross-check where that is?

## Data
- `data/SIA_Results_KLMA_2026.csv` (H/O isotopes, SIRFER lab): 46 feather samples,
  12 KLMA-banded individuals + 1 "San Diego" reference bird, cross-referenced
  against `data/SIA_Results_QAQC.csv` by SIRFER sample number.
- Bird **`1272-31742`**: present in the original SIRFER QC file (job 24-215.2)
  but never included in the KLMA study CSV -- its band-number prefix doesn't
  match the KLMA birds' (1422- vs 1272-). Per user confirmation, this is a
  reference specimen, not part of the core KLMA population; it's included
  throughout for comparison (same treatment as the San Diego bird) but
  excluded from the KLMA population-level summaries.
- `data/SIBS_CN_230621_Kellerman.csv` (C/N isotopes, separate SIBS lab, same
  envelope-ID numbering) cross-referenced against `data/SIBS_CN_QC_230621_Kellerman.csv`.
- Cleaned/merged tables: `data/klma_clean.csv` (`scripts/00_clean_data.R`),
  `data/cn_clean.csv` / `data/isotopes_combined_summary.csv` (`scripts/04_clean_cn_data.R`).
- Three H/O samples were QC-flagged during cleaning: a Nape sample from bird
  1422-02212 and a Nape sample from bird 1422-02214 (both implausibly high
  Wt% H/O, possible double-loading), and one Rectrix from 1422-02220 (O:H
  combustion ratio 8.5, far outside the normal 3.6-4.5 range). Separately,
  the C/N lab's own QC standards showed a consistent ~+2 permil enrichment
  bias for d13C specifically (both QA standards read ~2 permil less
  negative than their known consensus value).
- Growing-season precipitation isoscapes `d2h_GS.tif` / `d18o_GS.tif` (mean
  surfaces only, no prediction-SE grid -- see Limitations).
- `data/habitat_prior.tif`: wetland/wet-cropland habitat suitability, built
  from ESA WorldCover 10m land cover (`scripts/05_build_habitat_prior.R`).

## Methods
All analysis in R (`scripts/00`-`05`) using **assignR** (installed from GitHub
source, `SPATIAL-Lab/assignR`, since CRAN was not reachable from this
environment) plus `terra`/`sf`.

1. **Calibration (precipitation -> feather keratin).** assignR's built-in
   known-origin database contains **no Rallidae** and essentially no North
   American non-passerine oxygen data, so a Yellow-Rail-specific calibration
   isn't possible. We used the North-America-masked **Passerine** group
   (n=699 sites for d2H, n=167 sites for d18O) as the best available proxy.
   The d2H regression is solid (R^2=0.76, slope=0.71, p<2e-16); the **d18O
   regression is weak** (R^2=0.08, though still significant, p=0.0002). A
   d2H-only "Water bird" (mallard/scaup) calibration was run as a
   sensitivity check (R^2=0.55, similar slope, different intercept).
   **A more consequential sensitivity issue, found on this review:**
   `assignR::pdRaster()` itself warns "two or more isoscapes have shared
   variance > 0.7, added information will be limited, and specificity of
   assignments may be inflated." Checking this directly: the raw d2H and
   d18O growing-season isoscapes are correlated at r^2=0.99 across this
   study area (both trace the same precipitation fractionation process),
   so despite its weak calibration fit, d18O is not statistically
   independent of d2H the way `pdRaster`'s likelihood treats it. This
   turns out to matter a lot in practice -- see "Results, 1b" below, where
   dropping d18O and reassigning with d2H alone is shown to change one of
   the headline results.
2. **Isoscape uncertainty.** The supplied rasters have no SE band (the
   authoritative source, wateriso.utah.edu, wasn't reachable from this
   sandboxed environment). A spatially uniform SE was assumed (5 permil for
   d2H, 0.6 permil for d18O), sized against assignR's own low-res example SE
   raster.
3. **Study area.** CA/NV/AZ/NM south through Baja California and the western
   Mexican coast (the original frame of reference), extended north to
   include OR/WA/ID/UT so Klamath Marsh itself (~43 N) remains a candidate
   flight-feather origin.
4. **Habitat suitability prior.** Yellow Rails occupy wet sedge meadow /
   marsh / wet agriculture only (Cornell's Birds of the World account) --
   they are not found in the mountain, forest, or desert terrain that the
   H/O-only maps were placing high probability on, purely because those
   pixels happen to match the target d2H/d18O values by elevation. A
   suitability layer was built from ESA WorldCover 2021 10m land cover
   (herbaceous wetland = 1.0, cropland = 0.3 -- Cornell explicitly notes
   "wet agricultural areas" as habitat, but WorldCover doesn't distinguish
   flooded/irrigated cropland from dry cropland so this is a partial-credit
   weight, permanent water = 0.1 for marsh-edge effects, everything else =
   0), aggregated to the isoscape's ~9 km grid, and supplied to
   `assignR::pdRaster()` as a multiplicative `prior`. 85 of 108 candidate
   10m tiles across the study area had land data; after Gaussian smoothing,
   101,813 of 1,209,600 cells (~8.4%) across the full North American grid
   are nonzero -- most of the continent stays exactly zero. Within the
   narrower study area itself the smoothing kernel's finite radius means
   almost every cell picks up some nonzero value (it's not a hard binary
   mask there), but the distribution stays heavily skewed toward near-zero
   (median ~0.002; only about a third of study-area cells exceed 0.01).
   **This required a workaround for a real bug in assignR itself**:
   `pdRaster.isoStack()` crops the isoscape to `mask` internally but never
   crops `prior` to match, so the elementwise `assign * prior` step
   silently recycled two differently-sized vectors -- every single sample's
   posterior came out as exactly zero when first tested this way. The fix
   was to pre-crop/mask both the isoscape and the prior to the study area
   ourselves and pass `mask = NULL` to `pdRaster()`.
   Sanity check: Mt. Whitney (4,421 m alpine peak, no habitat) samples at
   ~8.5x10^-15 probability in the population-body surface; Sacramento
   National Wildlife Refuge (a real Central Valley wetland complex) samples
   at ~6.2x10^-4 -- more than ten orders of magnitude apart, confirming the
   prior is doing its job.
5. **Assignment.** `assignR::pdRaster()` (dual d2H+d18O joint-normal
   likelihood x habitat prior) on all 53 feather samples (43 from the 12
   KLMA-banded individuals, 3 from the "San Diego" reference bird, and 7
   from the `1272-31742` reference bird), normalized to the study area.
6. **Combining feathers per bird.** Where a bird had >1 feather in the same
   category, posteriors were multiplied and renormalized (shared-origin
   assumption). Where that produced no shared support at all (checked
   numerically), the mean of the individual posteriors was used instead and
   flagged (`combine_method` column in the summary table).
7. **Population summary.** Equal-weighted mean of KLMA individuals'
   posteriors per category (San Diego and 1272-31742 reference birds
   excluded).
8. **Credible regions.** `assignR::qtlRaster(thresholdType = "prob")` --
   the smallest area containing 50/75/90% of the cumulative posterior
   probability.
9. **C/N isotopes** (`scripts/04_clean_cn_data.R`) were cleaned and
   cross-referenced against the H/O results by envelope ID + feather
   category, and interpreted qualitatively (there's no established
   continental d13C/d15N isoscape product at usable resolution the way
   there is for d2H/d18O, so these aren't folded into the spatial model
   directly -- they're a habitat/diet cross-check).

## Results

### 1. Body vs. flight feathers point to different places, and it's now a much sharper picture

| | Body feathers | Flight feathers |
|---|---|---|
| Population probability, Mexico (Baja + mainland) | **58%** | 6% |
| Population probability, OR+WA+CA+ID+UT | 31% | **86%** |
| KLMA individuals with Mexico > 50% of their own probability | 6 of 10 | -- |
| KLMA individuals with OR+WA+CA > 50% | -- | 11 of 12 |
| 50%-credible-region area, population summary | **57,000 km^2** | **1,900 km^2** |

Adding the habitat prior shrank the credible regions substantially
(population body 50%-area: 389,000 km^2 pre-habitat-mask -> 57,000 km^2;
flight: 3,900 -> 1,900 km^2) without changing the basic body-vs-flight story
-- the broad zone-level split was already coming from the isotope gradient
itself, and the habitat mask sharpens *where within* that zone, excluding
the mountain/forest/desert terrain that isotope value alone couldn't rule
out. **Individual-level credible areas vary hugely and many are not
trustworthy at face value -- see "Match quality" below before citing any
single bird's area as a precise result.**

- **Body feathers** lean toward Mexico (mainland + Baja combined) in the
  main dual-isotope model, consistent with body (Formative/Definitive
  Basic) molt completing on the wintering grounds -- but see **1b below**,
  where this specific result (not the body-vs-flight contrast itself) is
  shown to be sensitive to the d18O collinearity issue and should be read
  with real caution. Several individuals (1422-02211, -02214, -02217, and
  especially the reference bird 1272-31742 at 95%) show most of that
  signal specifically on wintering-appropriate wetland patches along the
  Mexican coast rather than diffusely across the whole region.
- **Flight feathers** remain concentrated in the Pacific Northwest
  interior/coast, consistent with the complete prebasic molt of remiges
  happening primarily on/near breeding habitat. Unlike the body-feather
  result, this one is *not* sensitive to dropping d18O (see 1b) -- it
  actually gets stronger (86% -> 98% Pacific Northwest probability) under
  d2H alone.
- The same birds flagged before remain outliers under the habitat-masked
  model too -- the QC/consistency flags and the geographic outlier status
  still line up. This now includes a fourth case found on review:
  1422-02212's own flight feathers are isotopically inconsistent with each
  other (`combine_method = inconsistent_feathers_mean_fallback`), plausibly
  linked to its QC-flagged Nape sample above, in addition to 1422-02219's
  body feathers and 1422-02251's flight feathers.

### 1b. Sensitivity check: how much of the body-feather Mexico result is really coming from d18O?

Prompted by `pdRaster()`'s own collinearity warning (see Methods, item 1),
the full assignment was rerun using **d2H alone** (same habitat prior, same
AOI, same calibration approach) to see how much the dual-isotope model's
results actually depend on adding d18O, given that d18O's own spatial
pattern is 99%-correlated with d2H's in this study area. The two runs were
compared directly, bird by bird and population-level, not just by
eyeballing maps.

**Result: the flight-feather conclusion is robust; the body-feather
Mexico-majority conclusion is not.**

| | Population Body | Population Flight |
|---|---|---|
| P(Mexico), dual-isotope (main model) | 57.7% | 6.1% |
| P(Mexico), d2H-only | 34.6% | 0.8% |
| P(Pacific NW: OR/WA/ID/UT/CA), dual-isotope | 31.4% | 85.8% |
| P(Pacific NW: OR/WA/ID/UT/CA), d2H-only | 59.6% | 98.2% |

For **flight feathers**, dropping d18O barely changes anything -- if
anything the Pacific Northwest signal gets even stronger (86% -> 98%),
meaning that result does not depend on the collinearity issue at all.

For **body feathers**, dropping d18O flips the population-level majority
zone: the dual-isotope model's 58%-Mexico / 31%-Pacific-NW split becomes a
35%-Mexico / 60%-Pacific-NW split under d2H alone -- the *opposite*
majority zone. This isn't just a population-averaging artifact: checking
individual KLMA birds' body-feather posteriors the same way, most
(1422-02211, -02214, -02215, -02216, -02217) show their own Mexico
probability roughly cut in half when d18O is dropped (e.g. 1422-02211:
84% -> 40%; 1422-02217: 83% -> 35%), while a few (1422-02212, -02218,
-02221) stay Pacific-NW-majority under both models, and only 1422-02213
stays robustly Mexico-majority under both (93% -> 89%).

**One important exception: the reference bird 1272-31742's body-feather
result is robust.** Checked the same way, its 94.9% dual-isotope Mexico
probability barely moves under d2H alone (95.6%) -- unlike the more
marginal ~80-90%-Mexico KLMA individuals, its isotope values are extreme
enough that either model points the same way with high confidence. (Its
flight-feather result, already flagged elsewhere as isotopically
inconsistent between its own feathers, is *not* robust -- it swings from
29% to 0% Mexico probability -- but that result was already being
excluded from firm interpretation for the unrelated inconsistency reason.)
So this collinearity issue specifically undermines the *marginal* cases
(most of the KLMA population's body feathers), not every Mexico-leaning
result in the dataset.

**What this does and doesn't mean:** it does not mean the dual-isotope
result is simply "wrong" -- d18O's calibration, while weak (R^2=0.08), is
still statistically significant, and `calRaster()` does correctly inflate
its uncertainty to reflect that weak fit (rescaled sd ~2.6 permil vs the
raw assumed 0.6 permil). d18O is carrying *some* real information. But
`pdRaster`'s likelihood treats d2H and d18O as independent evidence when
computing the joint posterior, and at r^2=0.99 they very much are not --
so d18O's contribution is, per assignR's own warning, likely inflating the
model's apparent confidence in the Mexico direction for body feathers more
than is statistically justified. **The body-vs-flight contrast itself
survives this check and remains the headline result** (flight feathers are
far more Pacific-Northwest-concentrated than body feathers under both
models, 98% vs 60% even in the more conservative d2H-only version) -- but
the specific claim "body feathers show a Mexico-wintering majority" should
be presented as suggestive, not confirmed, pending either an isoscape SE
product that properly accounts for the d2H/d18O covariance, or treating
d2H-only as the primary model with dual-isotope as a secondary sensitivity
check (rather than the reverse, as this analysis currently does).

### 2. Reference bird 1272-31742: a clean, confident result
This bird (not part of the KLMA population, included for comparison) has
the tightest, most confident result in the whole dataset: its 4 combined
body feathers give a 94.9% probability of a mainland-western-Mexico origin,
in a 50%-credible area of just 715 km^2. Unlike most of the KLMA
population's more marginal body-feather results, this one is **robust** to
the d2H/d18O collinearity issue in section 1b below -- it barely moves
(94.9% -> 95.6%) when reassigned with d2H alone, because its isotope
values are extreme enough that either model points the same way with high
confidence. Its 3 flight feathers are more ambiguous (isotopically
inconsistent with each other -- flagged) and split across Mexico,
California, and the Pacific Northwest interior; that particular result is
*not* robust to the same check, but it was already excluded from firm
interpretation for the unrelated feather-inconsistency reason.

### 3. Carbon/nitrogen isotopes: an independent habitat cross-check
- **d13C is solidly C3 across every bird** (-18.4 to -25.5 permil) --
  nowhere near the C4 range (roughly -12 to -17 permil) that would flag
  Spartina/Distichlis salt-marsh vegetation (the classic Gulf Coast
  wintering habitat Cornell describes for some populations). These feathers
  reflect freshwater sedge-marsh/wet-meadow-type diets, not brackish salt
  marsh, for whichever molt event each feather represents.
- **d15N does *not* reliably track the H/O-based Mexico probability --
  this claim from an earlier version of this report was checked directly
  against the data and does not hold up.** The full set of body-feather
  d15N values alongside their dual-isotope Mexico probability
  (`data/cn_vs_geography_comparison.csv`):

  | Bird | d15N (permil) | P(Mexico), dual-isotope |
  |---|---|---|
  | 1272-31742 | 11.55 | 94.9% |
  | 1422-02221 | 8.84 | 0.4% |
  | 1422-02216 | 8.28 | **85.4%** |
  | 1422-02219 | 11.61 | 43.9% (flagged, inconsistent feathers) |

  1272-31742 and 1422-02221 do fit a "higher d15N -> higher Mexico
  probability" pattern, but 1422-02216 directly contradicts it (low d15N,
  85.4% Mexico), and 1422-02219 is a weak/ambiguous case already flagged
  for an unrelated reason. With only 4 usable data points and one clear
  counterexample, this is not a real trend -- it should not have been
  reported as one. We're leaving this correction here rather than quietly
  deleting the original claim, since it's a concrete example of the kind
  of small-n pattern-matching this analysis needs to guard against.
- One bird (1422-02220) had four body-feather C/N subsamples split sharply
  between freshwater-typical (d15N = 8.1) and elevated (d15N = 14.9-15.3)
  -- another case of an individual whose different body feathers evidently
  didn't grow together, on top of the H/O-flagged cases above.

## Outputs
- `outputs/maps/population_body_vs_flight.png` -- the main figure (now
  habitat-masked).
- `outputs/maps/habitat_prior.png` -- the wetland/wet-cropland suitability
  layer itself.
- `outputs/maps/individuals_body.png`, `individuals_flight.png` -- per-bird
  maps, including the reference birds.
- `outputs/maps/calibration_regression_*.png` -- calibration diagnostics.
- `outputs/tables/assignment_summary.csv` -- per-individual and population
  peak location, 50/75/90% credible areas, probability mass by zone, and
  three match-quality columns (see "Match quality" below):
  `peak_share`, `low_precision_flag`, `single_cell_dominant_flag`,
  `shared_peak_flag`.
- `data/isotopes_combined_summary.csv` -- C/N results cross-referenced
  against raw H/O isotope values.
- `data/cn_vs_geography_comparison.csv` -- C/N results cross-referenced
  against the H/O-based geographic assignment (zone probabilities, peak
  location, match-quality flags) -- the table behind Results sections 1b
  and 3.
- `outputs/rasters/*.tif` -- full posterior surfaces, for further GIS work.
- `data/habitat_prior.tif` -- the habitat suitability layer, for reuse.

## Match quality: not every tight credible region means a precise result

Manual review after posting the last version found that several individuals'
strikingly small credible regions (near-single-pixel, ~61 km^2) weren't
independent precise results -- **unrelated birds with different measured
isotope values were landing on the exact same grid cell.** Traced to its
source: parts of the study area (eastern Idaho/Wyoming in particular) have
an unusually *flat* isoscape (d2H standard deviation ~4.8 permil over a
5x4 degree box there, vs. ~17 permil over the same size area in California),
so the isotope data barely discriminates locations within it. Combined with
a habitat-suitability layer that's sparse to begin with, several birds'
broad, weakly-informative likelihood surfaces were all getting tipped onto
whichever single cell had the locally highest habitat weight in that flat
zone -- not a real match to that bird's own isotope signature. This is a
real property of the data and calibration (see limitations below, esp. the
assumed rather than measured isoscape SE) interacting with a sparse habitat
mask, not a coding bug, but it means **a small credible area is not on its
own evidence of a confident result.**

Two mitigations are in the pipeline now:
1. The habitat prior is Gaussian-smoothed (~18 km) before use, so nearby
   cells get partial credit instead of a hard cliff to zero. This helped
   partially (e.g. reference bird 1272-31742's flight credible area grew
   from 7,062 to 14,907 km^2) but did not eliminate the effect for the most
   affected individuals -- the underlying flat-isoscape problem is still
   there.
2. Three flag columns in `assignment_summary.csv`:
   - `peak_share`: fraction of a surface's entire probability mass held by
     its single top cell. High values (>0.3, flagged as
     `single_cell_dominant_flag`) mean the "credible region" is really
     one dominant pixel.
   - `shared_peak_flag`: TRUE if this surface's peak cell is identical to
     another bird's peak cell -- the direct signature of the failure mode
     above.
   - `low_precision_flag`: TRUE if area50_km2 < 1,000 km^2 (a coarser,
     easier-to-eyeball version of the same concern).

At last count, 14 of 26 individual surfaces (roughly half) carry at least
one of these flags. **For a manuscript, individual-level "this bird's
origin was pinpointed to N km^2" claims should not be made for flagged
birds** -- lean on the population-level zone probabilities instead, which
are far more robust (they're area-weighted sums over the whole zone, not
dependent on one argmax pixel).

## Key limitations (read before over-interpreting)
1. **No Rallidae calibration data exist anywhere** in the standard
   known-origin isotope databases used by assignR. The Passerine-based
   calibration is the best available proxy, not a validated one for rails.
2. **d18O is weakly calibrated (R^2=0.08) AND nearly collinear with d2H
   (r^2=0.99 in the raw isoscape)**, which `pdRaster()` itself flags as a
   condition where "specificity of assignments may be inflated." This is
   not theoretical: Results section 1b above shows that dropping d18O and
   reassigning with d2H alone flips the population body-feather result
   from Mexico-majority (58%) to Pacific-Northwest-majority (60%), while
   leaving the flight-feather result intact (in fact stronger). **The
   body-feather Mexico-wintering claim specifically should be treated as
   suggestive, not confirmed, until this is resolved** -- e.g. with an
   isoscape SE product that models the d2H/d18O covariance rather than
   treating the two as independent.
3. **Isoscape SE was assumed, not measured.**
4. **The habitat prior is approximate**: a single year (2021) of 10m land
   cover, a coarse hand-picked weighting scheme (wetland=1.0/cropland=0.3/
   water=0.1) rather than a validated habitat model, and it can't
   distinguish flooded/irrigated cropland (real Yellow Rail habitat per
   Cornell) from dry cropland (not habitat) -- both get the same 0.3
   weight. It also can't see small/narrow marsh patches below its ~9 km
   effective resolution particularly well (the 100m intermediate
   decimation step used for computational tractability could miss the
   smallest patches, though the largest, most likely-relevant marshes
   should survive).
5. **No age data.** Cornell's account makes age class (HY/juvenile vs.
   AHY/adult) important: a July-caught juvenile's feathers were grown on
   the natal territory, while an adult's most recent flight feathers were
   grown during the *previous* molt cycle, wherever that bird was at the
   time. Without age class, the body-vs-flight contrast found here is
   still consistent with the hypothesis, but can't distinguish "grown here
   this year" from "grown somewhere similar last year."
6. **C/N isotopes are a qualitative cross-check, not folded into the
   spatial model** -- there's no established fine-resolution continental
   d13C/d15N isoscape the way there is for d2H/d18O.
7. **Sample size** is small (12 KLMA individuals + 2 reference birds) and
   several feathers per bird were pooled under a shared-origin assumption
   that failed outright for 4 of the 14 birds with more than one feather in
   some category (1422-02212, 1422-02219, and 1422-02251 among the KLMA
   birds, plus the reference bird 1272-31742's flight feathers).

## Suggested next steps
- **Highest priority**: resolve the d2H/d18O collinearity issue in
  "Match quality... 1b" and "Key limitations" item 2 before treating the
  body-feather Mexico-majority result as a finding for publication. At
  minimum, report both the dual-isotope and d2H-only versions side by
  side rather than only the dual-isotope one; ideally, obtain or derive an
  isoscape product whose uncertainty model accounts for the known d2H/d18O
  covariance (the Global/Local Meteoric Water Line) rather than treating
  the two as independent evidence, which is what `pdRaster()`'s own
  warning is flagging.
- Obtain the official prediction-SE isoscape grids from waterisotopes.org
  (blocked from this sandboxed environment) and rerun with real,
  spatially-varying SE.
- Refine the habitat prior: a validated wetland classification (e.g.
  National Wetlands Inventory for the US portion) and a way to distinguish
  flooded/irrigated cropland from dryland cropland would both sharpen the
  Mexico results specifically, where the wet-agriculture signal matters
  most.
- If plumage-based ages are assignable, rerun stratified by HY vs AHY.
- Consider re-running the analysis excluding/re-measuring 1422-02212 given
  its QC flag, and treating 1422-02219/1422-02251/1272-31742-Flight as
  reflecting genuinely discordant feathers rather than combining them.

## Change log
- Initial version: assignR pipeline, no habitat mask.
- Correction 1: fixed `qtlRaster` credible-region math (was using
  `thresholdType = "area"`, which is not a real credible region) and the
  Baja/mainland Mexico zone split (was silently returning 0 for Baja on
  every individual); fixed map color scaling and missing legends.
- Correction 2: added the ESA WorldCover habitat suitability prior
  (fixing a real assignR bug in the process), added bird 1272-31742 as a
  reference individual, added the carbon/nitrogen isotope cross-check.
- Correction 3 (this version): map colors are now classified directly by
  credible region (50% / 50-75% / 75-90% / outside 90%) instead of by
  quantile-of-pixel-value. The quantile approach made broad,
  actually-low-probability areas look "hot" by relative rank, so the
  map's coloring and the credible-region area stats told visually
  different stories (a "28,000 km^2 credible area" stat next to a map
  that looked like most of the Southwest was lit up). Now the colored
  area on the map matches the credible-region area in the stats table by
  construction. Also dropped `per_feather_grid.png` -- an uncurated
  byproduct of assignR's own internal plotting (terra caps it at 16 of
  the 53 samples, with a color scheme unrelated to this project's maps)
  that was never a real deliverable; `individuals_body.png`/
  `individuals_flight.png` (one panel per bird, feathers already
  combined) are the per-bird figures to use.
- Correction 4: found and partially mitigated the single-cell convergence
  issue described in "Match quality" above -- Gaussian-smoothed the habitat
  prior and added `peak_share`, `low_precision_flag`,
  `single_cell_dominant_flag`, and `shared_peak_flag` columns to the
  summary table.
- Correction 5 (this version): a line-by-line check of this report against
  the actual data found several numbers had drifted out of sync with the
  current pipeline output (mostly stale from before the habitat-prior
  smoothing step was added). Fixed: the QC-flagged-sample count and
  description (it's one Nape sample each from 1422-02212 and 1422-02214,
  not two from 1422-02212); the isotopically-inconsistent-feather count
  (4 cases, not 3 -- 1422-02212's flight feathers were missing from the
  original list); the habitat-prior nonzero-cell percentages (the old
  "~5.8%/~58%" figures were stale pre-smoothing numbers); and the Mt.
  Whitney/Sacramento NWR sanity-check values (updated to the current
  raster's actual values, which if anything show an even larger
  mountain-vs-wetland contrast than originally reported).
- Correction 6 (this version): reran the full pipeline end to end from the
  raw data and confirmed it is fully deterministic (every output file --
  rasters, tables, and derived CSVs -- came back byte-for-byte identical
  to the previously committed version; `git status` showed zero diff
  after the full rerun). While investigating an `assignR::pdRaster()`
  warning surfaced during that rerun ("two or more isoscapes have shared
  variance > 0.7, added information will be limited, and specificity of
  assignments may be inflated"), found and quantified a substantive,
  previously-undocumented sensitivity issue: the raw d2H/d18O isoscapes
  are r^2=0.99 correlated in this study area, and a d2H-only reassignment
  (same habitat prior, same calibration approach) flips the population
  body-feather result from Mexico-majority to Pacific-Northwest-majority,
  while leaving the flight-feather result intact. Added as "Results,
  section 1b" and promoted to the #2 item in "Key limitations" and the
  top item in "Suggested next steps," since it bears directly on how
  confidently the body-feather Mexico-wintering claim can be stated.
