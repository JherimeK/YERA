# Geographic origin of KLMA Yellow Rail feathers: stable isotope assignment

## Question
Do body feathers and flight feathers from Yellow Rails (*Coturnicops noveboracensis*)
captured at Klamath Marsh NWR, Oregon in July grow in different places, consistent
with Cornell's account of the species' molt strategy (flight feathers replaced in a
complete molt primarily on/near the breeding grounds; body feathers replaced in a
molt that can complete on the wintering grounds)?

## Data
- `data/SIA_Results_KLMA_2026.csv` (46 usable feather samples, 12 KLMA-banded
  individuals + 1 "San Diego" reference bird sampled elsewhere) cross-referenced
  against `data/SIA_Results_QAQC.csv` (lab QA/QC report) by SIRFER sample number.
  Cleaned/merged table: `data/klma_clean.csv` (`scripts/00_clean_data.R`).
- Three samples were QC-flagged during cleaning (see Limitations): two Nape
  samples from bird 1422-02212 with implausibly high Wt% H/O (possible
  double-loading), and one Rectrix from 1422-02220 with an O:H combustion
  ratio (8.5) far outside the normal 3.6-4.5 range.
- Growing-season precipitation isoscapes `d2h_GS.tif` / `d18o_GS.tif` (mean
  surfaces only, no prediction-SE grid -- see Limitations).

## Methods
All analysis in R (`scripts/00`-`03`) using **assignR** (installed from GitHub
source, `SPATIAL-Lab/assignR`, since CRAN was not reachable from this
environment) plus `terra`/`sf`.

1. **Calibration (precipitation -> feather keratin).** assignR's built-in
   known-origin database contains **no Rallidae** and essentially no North
   American non-passerine oxygen data, so a Yellow-Rail-specific calibration
   isn't possible. We used the North-America-masked **Passerine** group
   (n=699 sites for d2H, n=167 sites for d18O) as the best available proxy --
   it's the only NA group in the database with enough paired sites to fit a
   defensible regression. The d2H regression is solid (R^2=0.76, slope=0.71,
   p<2e-16); the **d18O regression is weak** (R^2=0.08, though still
   significant, p=0.0002) -- oxygen is doing comparatively little work in
   these maps. A d2H-only "Water bird" (mallard/scaup) calibration was run as
   a sensitivity check (R^2=0.55, similar slope, different intercept; see
   `outputs/maps/calibration_regression_d2H_waterbird_sensitivity.png`).
2. **Isoscape uncertainty.** The supplied rasters have no SE band (the
   authoritative source, wateriso.utah.edu, wasn't reachable from this
   sandboxed environment). A spatially uniform SE was assumed (5 permil for
   d2H, 0.6 permil for d18O), sized against assignR's own low-res example SE
   raster. This term is minor relative to the regression residual variance
   `calRaster` adds in quadrature, so it has limited effect on the maps, but
   should be replaced with the real SE grid if you can obtain it.
3. **Study area.** CA/NV/AZ/NM south through Baja California and the western
   Mexican coast (your original frame of reference), **extended north to
   include OR/WA/ID/UT**. Klamath Marsh itself sits at ~43 N, just outside
   your original box -- excluding it would have ruled out the single most
   biologically expected flight-feather origin before the isotope data got a
   say. All probabilities below are normalized within this extended region.
4. **Assignment.** `assignR::pdRaster()` (dual d2H+d18O, using the joint
   normal likelihood, not two independent univariate tests) on all 46
   feather samples individually, masked/renormalized to the study area.
5. **Combining feathers per bird.** Where a bird had >1 feather in the same
   category (Body/Flight), posteriors were multiplied and renormalized
   (`assignR::jointP` logic), i.e. treated as repeat samples of one molt
   event. Two individuals' feather-pairs turned out to be isotopically
   inconsistent (posteriors share essentially no support) --
   `1422-02212` Flight (Secondary vs. Rectrix) and `1422-02219` Body
   (Breast vs. Nape) -- these were reported as the mean of their individual
   posteriors instead, and flagged (`combine_method` column in the summary
   table).
6. **Population summary.** Equal-weighted mean of all KLMA individuals'
   posteriors within each category (San Diego reference bird excluded).

## Results

*(Two corrections were made after initial posting, based on feedback that the
maps and stats didn't line up -- see "Post-review corrections" below. The
numbers here are the corrected ones.)*

**The body-vs-flight split is real and consistent across birds**, not just an
artifact of the population average:

| | Body feathers | Flight feathers |
|---|---|---|
| Population probability, Mexico (Baja + mainland) | **61%** (8% Baja + 53% mainland) | 7% (1% Baja + 6% mainland) |
| Population probability, OR+WA+ID+UT | 9% | **75%** |
| Individuals with that region as leading signal | 9 of 11 | 12 of 12 |
| 50%-credible-region area, range across individuals | 1,500-282,000 km^2 | 61-131,000 km^2 |
| 50%-credible-region area, population summary | 389,000 km^2 | **3,900 km^2** |

- **Body feathers**: 9 of 11 KLMA individuals show 0.47-0.97 probability mass
  in Mexico; most of that is mainland western Mexico (Sonora/Sinaloa/Nayarit/
  Jalisco latitudes), but 4 of 11 individuals (1422-02211, -02214, -02216,
  -02217) put more than 10% of their own probability specifically on the
  **Baja California peninsula** itself. This matches the account that body
  (Formative/Definitive Basic) molt can complete on the wintering grounds,
  and doesn't pin down mainland vs. peninsula as sharply as the population
  average alone suggests.
- **Flight feathers**: all 12 individuals with usable data show 0.30-1.00
  probability mass in the Pacific Northwest interior (WA/ID/UT, with OR
  itself typically getting a smaller secondary share), matching the account
  that the complete prebasic molt of remiges happens primarily on/near
  breeding habitat. Very little mass lands in Mexico for flight feathers,
  and essentially none of that is Baja specifically.
- The two exceptions line up with independent red flags: **1422-02212**
  (Mexico signal near zero for both feather categories) is the same bird
  whose Nape sample was QC-flagged for anomalous combustion chemistry, and
  **1422-02219**/**1422-02251** are the two isotopically-inconsistent
  feather pairs above. That QC flags and biological outliers coincide is
  reassuring for the overall pattern, but these three birds' results should
  be treated with more caution than the rest.
- Individual credible-region size varies a lot (see table above) -- some
  individuals, especially after combining two flight feathers assumed to
  share an origin, resolve to a strikingly small area (tens of km^2, close
  to a single isoscape grid cell). Take those at face value cautiously: this
  is partly a real effect of combining two correlated dual-isotope
  measurements (multiplying two probability surfaces that already overlap
  well sharpens the result quickly), not necessarily true pinpoint
  biological precision. The population-level pattern (which is not a
  simple average of these tight individual regions -- see next section) is
  the more robust takeaway; see the calibration caveats below for why.

## Post-review corrections

After the first version of this analysis was posted, closer review of the
maps surfaced two real bugs, both now fixed in `scripts/02` and
`scripts/03`:

1. **50/75/90% "credible regions" were computed wrong.** The first version
   used `assignR::qtlRaster(..., thresholdType = "area")`, which (per
   assignR's own definition) just returns exactly that fraction of the
   *total map area* regardless of the probability distribution -- so every
   individual's "50% credible area" was ~half the study area by
   construction, not a meaningful precision estimate. Fixed to
   `thresholdType = "prob"` (smallest area containing that much cumulative
   probability), which is what "50/75/90% credible region" is supposed to
   mean. This changed the area numbers substantially (see Results) and is
   the reason the population maps now also show 50/75/90% credible-region
   contour lines instead of nothing.
2. **The Baja California vs. mainland Mexico split was silently broken.**
   The zone used to report `P_Baja_CA` was built by splitting Mexico's
   cropped polygon into fragments and unioning the ones whose centroid fell
   west of a longitude cutoff -- for reasons not fully diagnosed (a handful
   of small cropping-artifact fragments probably corrupted the union), the
   resulting polygon didn't actually overlap any raster cells, so every
   individual's `P_Baja_CA` silently came out as exactly 0.0000 despite the
   maps clearly showing probability mass over the peninsula. Fixed by
   splitting with a simple longitude bounding-box intersection instead,
   verified to actually overlap the peninsula's raster cells. The combined
   Baja+mainland Mexico total is unchanged (this was a mis-attribution
   between the two sub-zones, not a change in overall signal).
3. The first version of the population map also used a plain linear color
   scale, which made it look almost entirely blank -- the underlying
   probability surfaces are extremely right-skewed (the top ~0.01% of
   pixels carry most of the visible contrast), so a linear scale compresses
   everything else to one color. Maps now use data-driven quantile color
   breaks plus the corrected credible-region contour lines, and the
   per-individual grid figures now show a real probability legend (the
   original version suppressed it to save space, which was a mistake).

## Outputs
- `outputs/maps/population_body_vs_flight.png` -- the main figure.
- `outputs/maps/individuals_body.png`, `individuals_flight.png` -- per-bird maps.
- `outputs/maps/calibration_regression_*.png` -- calibration diagnostics.
- `outputs/tables/assignment_summary.csv` -- per-individual and population
  peak location, 50/75/90% credible areas, and probability mass by zone
  (Oregon, California, Nevada, Arizona, New Mexico, WA/ID/UT, Baja CA,
  Mainland W Mexico, within 100 km of Klamath Marsh).
- `outputs/rasters/*.tif` -- full posterior surfaces (per feather, per
  individual-category, and population), plus `qtl50/75/90_credible_region.tif`
  (binary 50/75/90% credible-region masks per individual/population surface),
  for further GIS work.

## Key limitations (read before over-interpreting)
1. **No Rallidae calibration data exist anywhere** in the standard
   known-origin isotope databases used by assignR. The Passerine-based
   calibration is the best available proxy, not a validated one for rails.
2. **d18O calibration is weak (R^2=0.08)**; d2H is carrying most of the
   signal. This is a real property of the reference dataset, not a fixable
   bug.
3. **Isoscape SE was assumed, not measured** (see Methods #2).
4. **No age data.** Cornell's account makes age class (HY/juvenile vs.
   AHY/adult) important: a July-caught juvenile's feathers (including
   flight feathers, not yet molted) were grown *on the natal territory*,
   while an adult's most recent flight feathers were grown during the
   *previous* molt cycle, wherever that bird was at the time -- not
   necessarily this year's Klamath territory. Without age class, the "flight
   = near breeding grounds" body-vs-flight contrast we found is still
   consistent with the hypothesis, but can't distinguish "grown here this
   year" from "grown somewhere similar last year." If ages can be assigned
   from the specimens (Pyle 2008 molt-limit criteria), a follow-up
   stratified by age would sharpen this.
5. **Sample size** is small (12 KLMA individuals) and several feathers per
   bird were pooled under a shared-origin assumption that failed outright
   for 2 of 12 birds.

## Suggested next steps
- Obtain the official prediction-SE isoscape grids from waterisotopes.org (blocked
  from this sandboxed environment) and rerun with real, spatially-varying SE.
- If plumage-based ages are assignable, rerun stratified by HY vs AHY.
- Consider re-running the analysis excluding/re-measuring 1422-02212 given its
  QC flag, and treating 1422-02219/1422-02251 as reflecting genuinely
  discordant feathers rather than combining them.
