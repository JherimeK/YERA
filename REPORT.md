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
- Four samples were QC-flagged during cleaning: two Nape samples from bird
  1422-02212 (implausibly high Wt% H/O, possible double-loading), one
  Rectrix from 1422-02220 (O:H combustion ratio 8.5, far outside the normal
  3.6-4.5 range), and a consistent ~+2 permil enrichment bias on the C/N
  lab's own QC standards for d13C specifically (both QA standards read ~2
  permil less negative than their known consensus value).
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
   regression is weak** (R^2=0.08, though still significant, p=0.0002) --
   oxygen is doing comparatively little work in these maps. A d2H-only
   "Water bird" (mallard/scaup) calibration was run as a sensitivity check
   (R^2=0.55, similar slope, different intercept).
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
   10m tiles across the study area had land data; ~5.8% of the full domain
   and ~58% of the study-area cells have nonzero suitability.
   **This required a workaround for a real bug in assignR itself**:
   `pdRaster.isoStack()` crops the isoscape to `mask` internally but never
   crops `prior` to match, so the elementwise `assign * prior` step
   silently recycled two differently-sized vectors -- every single sample's
   posterior came out as exactly zero when first tested this way. The fix
   was to pre-crop/mask both the isoscape and the prior to the study area
   ourselves and pass `mask = NULL` to `pdRaster()`.
   Sanity check: Mt. Whitney (4,421 m alpine peak, no habitat) samples at
   ~4x10^-14 probability in the population-body surface; Sacramento
   National Wildlife Refuge (a real Central Valley wetland complex) samples
   at ~1.4x10^-3 -- a >10,000x difference, confirming the prior is doing
   its job.
5. **Assignment.** `assignR::pdRaster()` (dual d2H+d18O joint-normal
   likelihood x habitat prior) on all 53 feather samples (46 KLMA + 7 from
   the reference bird), normalized to the study area.
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

- **Body feathers** remain concentrated in Mexico (mainland + Baja
  combined), consistent with body (Formative/Definitive Basic) molt
  completing on the wintering grounds. Several individuals (1422-02211,
  -02214, -02217, and especially the reference bird 1272-31742 at 99%)
  show most of that signal specifically on wintering-appropriate wetland
  patches along the Mexican coast rather than diffusely across the whole
  region.
- **Flight feathers** remain concentrated in the Pacific Northwest
  interior/coast, consistent with the complete prebasic molt of remiges
  happening primarily on/near breeding habitat.
- The same three birds flagged before (1422-02212, whose Nape sample was
  QC-flagged; 1422-02219 and 1422-02251, whose feathers are isotopically
  inconsistent with each other) remain outliers under the habitat-masked
  model too -- the QC/consistency flags and the geographic outlier status
  still line up.

### 2. Reference bird 1272-31742: a clean, confident result
This bird (not part of the KLMA population, included for comparison) has
the tightest, most confident result in the whole dataset: its 4 combined
body feathers give a 99.3% probability of a mainland-western-Mexico origin,
in a 50%-credible area of just 238 km^2. Its 3 flight feathers are more
ambiguous (isotopically inconsistent with each other -- flagged) and split
across Mexico, California, and the Pacific Northwest interior.

### 3. Carbon/nitrogen isotopes: an independent habitat cross-check
- **d13C is solidly C3 across every bird** (-18.4 to -25.5 permil) --
  nowhere near the C4 range (roughly -12 to -17 permil) that would flag
  Spartina/Distichlis salt-marsh vegetation (the classic Gulf Coast
  wintering habitat Cornell describes for some populations). These feathers
  reflect freshwater sedge-marsh/wet-meadow-type diets, not brackish salt
  marsh, for whichever molt event each feather represents.
- **d15N varies by bird in a way that lines up with the H/O-based
  geography**: low d15N (~7-9 permil, freshwater-marsh-typical) tracks with
  lower Mexico probability (1422-02216, -02221); elevated d15N (~11.5-15.4
  permil, plausibly marking brackish/estuarine influence or fertilized wet
  agriculture -- both consistent with Cornell's account of wet agricultural
  habitat, and common in coastal Sinaloa/Nayarit shrimp-aquaculture and
  irrigated-cropland wetlands) tracks with higher Mexico probability
  (1272-31742 at 99.3% Mexico has d15N = 11.55, the second-highest
  reference-bird-adjacent value in the dataset).
- One bird (1422-02220) had four body-feather C/N subsamples split sharply
  between freshwater-typical (d15N = 8.1) and elevated (d15N = 13.1-15.3)
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
  against H/O geographic assignment.
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
2. **d18O calibration is weak (R^2=0.08)**; d2H is carrying most of the
   signal.
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
7. **Sample size** is small (12 KLMA individuals + 1 reference) and several
   feathers per bird were pooled under a shared-origin assumption that
   failed outright for 3 of 13 birds (including the reference bird's
   flight feathers).

## Suggested next steps
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
- Correction 4 (this version): found and partially mitigated the
  single-cell convergence issue described in "Match quality" above --
  Gaussian-smoothed the habitat prior and added `peak_share`,
  `low_precision_flag`, `single_cell_dominant_flag`, and
  `shared_peak_flag` columns to the summary table.
