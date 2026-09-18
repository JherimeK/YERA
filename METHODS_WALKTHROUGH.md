# Step-by-step methods walkthrough: KLMA Yellow Rail geographic origin analysis

This document walks through the full pipeline, in the order things actually
need to run, from the two raw lab data sets to the final maps and tables.
It's meant to be read alongside the R code in `scripts/`, so a second reviewer
can follow (and re-run) exactly what was done. For the short publication-style
version, see `MANUSCRIPT_METHODS.md`.

## The R files, and the order to run them in

The scripts are numbered `00`-`05`, but the numbering is *not* the run order
for a full pipeline, because script `04` has a final step that depends on
script `02`'s output. **Run them in this order:**

| Order | File | What it needs | What it produces |
|---|---|---|---|
| 1 | `scripts/00_clean_data.R` | `data/SIA_Results_KLMA_2026.csv`, `data/SIA_Results_QAQC.csv` | `data/klma_clean.csv` |
| 2 | `scripts/05_build_habitat_prior.R` | internet access (downloads ESA WorldCover tiles), `d2h_GS.tif` | `data/habitat_prior.tif` |
| 3 | `scripts/01_run_assignment.R` | `data/klma_clean.csv`, `data/habitat_prior.tif`, `d2h_GS.tif`, `d18o_GS.tif` | `outputs/rasters/pd_per_feather.tif`, calibration diagnostic plots |
| 4 | `scripts/02_summarize_assignment.R` | `data/klma_clean.csv`, `outputs/rasters/pd_per_feather.tif` | `outputs/tables/assignment_summary.csv`, per-individual/population posterior rasters, credible-region rasters |
| 5 | `scripts/03_make_maps.R` | everything from step 4 | all files in `outputs/maps/`, `outputs/tables/Table1_manuscript.csv` |
| 6 | `scripts/04_clean_cn_data.R` | `data/klma_clean.csv`, `data/SIBS_CN_230621_Kellerman.csv`, and (for its last step only) `outputs/tables/assignment_summary.csv` | `data/cn_clean.csv`, `data/cn_summary_by_individual.csv`, `data/isotopes_combined_summary.csv`, `data/cn_vs_geography_comparison.csv` |

Script `05` only needs to be run once (its output, `data/habitat_prior.tif`,
is checked into the repo already, so step 2 can be skipped on a re-run unless
you want to rebuild the habitat layer from scratch — it re-downloads ~85
ESA WorldCover tiles and takes several minutes).

Script `04` is written so it's still useful to run early (right after step 1)
if you just want to look at the raw C/N data — its first few outputs only
need `klma_clean.csv`. But its *last* block, which cross-references C/N
against the geographic assignment, checks whether
`outputs/tables/assignment_summary.csv` exists yet, and simply prints a
message and skips that one output if it doesn't. So running `04` twice
(once early, once at the end) is harmless and even expected.

## Step-by-step: what each script actually does

### Step 1 — `00_clean_data.R`: build the master feather-sample table

- Reads the KLMA study spreadsheet (`SIA_Results_KLMA_2026.csv`: one row per
  feather sample, with δ2H, δ18O, %H, %O, the H:O combustion ratio, bird
  envelope ID, and feather type) and the SIRFER lab's raw QA/QC report
  (`SIA_Results_QAQC.csv`, a stacked text export covering 4 analytical
  batches/"jobs").
- Cross-references every sample to its batch by SIRFER sample number
  (`Original_ID`), pulling in the run date and that batch's own QC-standard
  precision (from the POW secondary reference material) for traceability.
- Independently re-derives Body vs. Flight category from the feather-type
  text (matching against rectrix/rectrice/secondary/primary/tertial/covert
  keywords) and cross-checks it against the supplied category column, so any
  mislabeling gets flagged rather than silently trusted.
- Pulls in one extra individual, band `1272-31742`, directly out of the raw
  QC text. This bird's band-number prefix (`1272-`) doesn't match the KLMA
  study's `1422-` prefix — it was analyzed in the same SIRFER batch but never
  entered into the KLMA study spreadsheet. Per your confirmation, it's a
  reference specimen, not a KLMA-population bird, so it's carried through the
  whole pipeline labeled `Group = "Reference_1272-31742"` and excluded from
  population-level averages, but included everywhere individual results are
  shown for comparison (the same treatment given to the existing "San Diego"
  reference bird already in the KLMA file).
- Flags (but does not drop) samples with an implausible O:H combustion ratio
  (outside ~3.0–5.5) or implausible %H/%O (possible double-loading). Four
  samples were flagged this way: two Nape samples from bird `1422-02212`,
  one Rectrix from `1422-02220`.
- Writes `data/klma_clean.csv`: one row per feather sample, every bird
  (KLMA + both reference birds) in one table.

### Step 2 — `05_build_habitat_prior.R`: build the wetland/wet-cropland mask

**Why this step exists:** Yellow Rails are a wetland/wet-agriculture obligate
(Cornell's Birds of the World account: wet sedge meadow, marsh, wet
agricultural areas — not forest, not desert, not high mountains). The very
first version of this analysis used the isotope data alone, and it produced
"credible" origins on mountain peaks and desert terrain, purely because those
places happen to match a bird's δ2H/δ18O values by elevation. This step
builds a raster that down-weights everywhere that isn't plausible rail
habitat, so that elevation-driven false matches get suppressed.

- Downloads ESA WorldCover 2021 10 m global land-cover tiles covering the
  whole study area (108 candidate tiles; 85 actually contain land, the rest
  are open ocean and are skipped).
- Each 10 m tile is first decimated to 100 m with a fast nearest-neighbor
  `gdalwarp` resample (full 10 m resolution across the whole domain would be
  ~3.2 billion pixels — far too slow to reclassify directly in R), then
  reclassified in R to a suitability weight: **herbaceous wetland = 1.0**
  (primary habitat), **cropland = 0.3** (Cornell explicitly lists "wet
  agricultural areas" as habitat, but WorldCover's cropland class doesn't
  distinguish flooded/irrigated cropland from dry cropland, so this is a
  partial-credit weight, not a confirmed-habitat weight), **permanent water
  = 0.1** (marsh-edge effect; rails don't use open water, but wetlands are
  often fringed by or misclassified as water), **everything else = 0**.
- Each tile's suitability raster is then averaged down onto the isoscape's
  native grid (~0.083°, ~9 km), and all tiles are mosaicked into one
  domain-wide raster.
- The raw mask is then **Gaussian-smoothed** (~2 grid cells, ~18 km). This
  matters and is explained more in Step 3 below — without it, isolated
  wetland pixels surrounded by hard zeros caused a real problem where several
  unrelated birds' posteriors all collapsed onto the same single pixel.
- The result is force-aligned (same CRS, same exact grid, via
  `resample(..., method="near")`) to the isoscape template, to avoid a subtle
  floating-point/CRS-string mismatch that otherwise breaks a downstream
  `assignR` sanity check (see Step 3).
- Writes `data/habitat_prior.tif`. About 5.8% of the full raster domain (and
  ~58% of the narrower study-area cells) end up with nonzero suitability
  after smoothing.

### Step 3 — `01_run_assignment.R`: the core isotope assignment

This is the main `assignR` step: converting each feather's measured δ2H/δ18O
into a map of "how likely is it this feather grew at pixel X," for every
pixel in the study area.

- **Study area.** California, Nevada, Arizona, New Mexico, Oregon,
  Washington, Idaho, Utah, plus Baja California and the western Mexican
  coast (bounded roughly 14–33° N, 101–118° W). Your original frame of
  reference was CA/NV/AZ/NM south through Baja/western Mexico; Oregon (plus
  WA/ID/UT for geographic continuity) was added because Klamath Marsh NWR
  itself — the capture/breeding site — sits at ~43° N, just north of the
  California line, and Cornell's account says flight feathers grow primarily
  on or near the breeding grounds. Excluding Oregon would have ruled out the
  single most biologically expected flight-feather origin before the isotope
  data got a say.
- **Isoscapes.** Uses the supplied growing-season precipitation δ2H and δ18O
  rasters (`d2h_GS.tif`, `d18o_GS.tif`) as the base surfaces. These are
  mean-only rasters with no prediction-uncertainty band (the official SE
  grids live at wateriso.utah.edu, which wasn't reachable from this
  environment). A spatially uniform standard deviation is assumed instead (5‰
  for δ2H, 0.6‰ for δ18O), sized to be roughly consistent with `assignR`'s
  own bundled example uncertainty raster. This is a real approximation —
  see limitations.
- **Calibration (precipitation isotope value → feather keratin isotope
  value).** `assignR`'s reference database of species with known-origin
  isotope measurements contains no Rallidae (rails) at all, and essentially
  no North American non-passerine δ18O data — so a rail-specific calibration
  isn't possible with the tools/data available. The North-America-masked
  **"Passerine"** group (n = 699 known-origin sites for δ2H, a subset with
  δ18O) is used as the best available stand-in. The δ2H regression is solid
  (R² = 0.76). The δ18O regression is weak (R² = 0.08, though still
  statistically significant, p = 0.0002) — meaning δ18O contributes
  comparatively little discriminating power to the final maps; δ2H is doing
  most of the work. As a sensitivity check, a second, δ2H-only calibration
  using the "Water bird" reference group (dabbling/diving waterfowl —
  ecologically closer to a marsh-dwelling rail than a generic passerine) is
  also run and saved (R² = 0.55, similar slope, different intercept) but not
  used in the main results.
- **Habitat prior integration — and a real `assignR` bug found and worked
  around.** The habitat suitability raster from Step 2 is supplied to
  `assignR::pdRaster()` as a multiplicative prior, so the final posterior at
  each pixel is (isotope-match likelihood) × (habitat suitability). Testing
  this the straightforward way (passing both `mask` and `prior` to
  `pdRaster()`) produced every single sample's posterior as exactly zero
  everywhere. Diagnosis: `pdRaster.isoStack()` crops the isoscape to `mask`
  internally, but never applies that same crop to `prior` — so the
  elementwise `likelihood * prior` multiplication silently multiplies two
  differently-sized/misaligned vectors (R's default behavior for this is
  silent recycling, not an error, which is why it failed quietly rather than
  throwing an error message). The fix used here: pre-crop and pre-mask both
  the calibrated isoscape and the habitat prior to the study area ourselves,
  then call `pdRaster(..., mask = NULL, prior = habitat_prior)` so nothing
  gets cropped a second time inside the function.
- **Sanity check performed:** in the resulting population body-feather
  surface, Mt. Whitney (a 4,421 m alpine peak with zero rail habitat) samples
  at a posterior probability of ~4×10⁻¹⁴, while Sacramento National Wildlife
  Refuge (a real Central Valley wetland complex) samples at ~1.4×10⁻³ — a
  >10,000-fold difference, confirming the habitat prior is meaningfully
  reshaping the result rather than being a no-op.
- Runs `pdRaster()` once per feather sample (53 total: 46 KLMA feathers + 7
  from the two reference birds), producing one posterior-probability raster
  per feather. Writes the full stack to
  `outputs/rasters/pd_per_feather.tif`.

### Step 4 — `02_summarize_assignment.R`: combine feathers, summarize, flag

- **Combine feathers within a bird.** Where a bird contributed more than one
  feather of the same category (e.g., two Nape samples), those feathers'
  posteriors are multiplied together and renormalized — the standard
  assumption that feathers of the same category share one molt origin. If
  that product comes out effectively zero everywhere (i.e., the two
  feathers' independent posteriors don't actually overlap anywhere in the
  study area — meaning the feathers are isotopically inconsistent with each
  other), the script falls back to the simple mean of the two feathers'
  posteriors instead, and records `combine_method =
  "inconsistent_feathers_mean_fallback"` in the output table so this can be
  reviewed rather than silently masked. This triggered for 3 of the
  multi-feather cases (bird `1422-02219`, bird `1422-02251`, and the
  reference bird `1272-31742`'s flight feathers).
- **Population summaries.** For KLMA birds only (both reference birds
  excluded), an unweighted mean posterior is computed separately for Body and
  Flight categories, each renormalized to sum to 1 over the study area.
- **Credible regions.** For every individual×category surface and both
  population surfaces, computes the 50%, 75%, and 90% credible regions using
  `assignR::qtlRaster(thresholdType = "prob")` — the smallest-area region
  containing that much cumulative posterior probability. (An earlier version
  of this analysis used `thresholdType = "area"`, which is **not** a real
  credible region — it just returns exactly that fraction of the map's total
  area regardless of where probability actually sits. That bug was caught
  and fixed before this version.)
- **Zone probability table.** For 8 named zones (Oregon, California, Nevada,
  Arizona, New Mexico, a combined Washington/Idaho/Utah zone, Baja
  California, and mainland western Mexico), sums each surface's posterior
  probability mass falling inside that zone's polygon.
- **Peak location and match-quality flags** (this is the important part for
  interpreting individual results correctly — see below): for every surface,
  records the single highest-probability pixel's coordinates and the
  fraction of that surface's *entire* probability mass held by that one
  pixel (`peak_share`). Three flag columns are derived from this:
  - `low_precision_flag`: 50%-credible area < 1,000 km².
  - `single_cell_dominant_flag`: `peak_share` > 0.3 (i.e., one grid cell
    alone holds more than 30% of a surface's total probability).
  - `shared_peak_flag`: this surface's peak pixel is the *exact same*
    pixel as another, unrelated bird's peak pixel.
- Writes `outputs/tables/assignment_summary.csv` (the master results table)
  and the underlying rasters (`pd_per_individual_category.tif`,
  `pd_population_by_category.tif`, `qtl50/75/90_credible_region.tif`).

**Why the flag columns matter — please read before trusting any single
bird's "precise" result.** Manual review found that several individuals'
strikingly small credible regions (some as small as a single ~61 km² pixel)
were not, in fact, independent precise results — different birds with
different measured isotope values were landing on the *identical* grid cell.
The cause: parts of the study area (eastern Idaho/Wyoming, in particular)
have an unusually flat isoscape (δ2H standard deviation ~4.8‰ over a 5°×4°
box there, versus ~17‰ over a similarly sized box in California) — meaning
the isotope data barely distinguishes locations within that region. Combined
with a habitat mask that's inherently sparse (wetlands are a small fraction
of any landscape), several birds whose true isotope values fell outside what
was achievable anywhere within the masked wetland cells ended up defaulting
to whichever single wetland pixel was the "least bad" available match — and
because wetland pixels are sparse, unrelated birds could converge on the
exact same "least bad" pixel. This is a genuine property of the data and
calibration (see limitations) interacting with the sparse habitat mask, not
a coding bug — but it means a small credible area or a confident-looking
peak is not, on its own, evidence of a precise result. The Gaussian smoothing
in Step 2 and the flag columns here are mitigations, not a full fix: as of
this version, roughly half of the 26 individual×category surfaces (14 of 26)
carry at least one of the three flags.

### Step 5 — `03_make_maps.R`: final maps and the publication table

- Classifies every surface into 4 discrete color classes **directly from the
  credible-region rasters** from Step 4 (top 50% / 50–75% / 75–90% / outside
  the 90% region), rather than by relative pixel-value quantile. (An earlier
  version colored by quantile-of-pixel-value, which — because these
  posterior surfaces are extremely right-skewed — made broad, actually
  low-probability areas look visually "hot," so the map's coloring and the
  credible-region area statistics told two different stories. Coloring
  directly by credible region means the colored area on the map *is* the
  credible-region area reported in the stats table, by construction.)
- Produces the main population-level Body-vs-Flight comparison figure
  (`population_body_vs_flight.png`).
- Produces individual-bird map grids, **split into a main figure (unflagged
  birds only) and a supplementary figure (flagged birds, each panel labeled
  "[low match quality]")** for both Body and Flight — this split was a
  deliberate choice discussed with you, so the primary figures for
  publication aren't cluttered with results that carry the match-quality
  caveat above, while nothing is hidden (the supplementary figures are still
  produced and are just as complete).
- Produces a habitat-prior diagnostic map (`habitat_prior.png`), so the
  suitability layer itself can be visually sanity-checked against known
  marsh locations.
- Produces `outputs/tables/Table1_manuscript.csv`, a cleaned publication
  table: one row per bird × feather category, with peak location, 50%/90%
  credible area, % probability mass in Mexico vs. the Pacific Northwest, and
  a plain-text quality note for flagged rows.

### Step 6 — `04_clean_cn_data.R`: carbon/nitrogen cross-check

- Cleans the SIBS lab's carbon/nitrogen results (a separate lab and separate
  analytical run from the SIRFER δ2H/δ18O data, but keyed to the same
  envelope IDs), parsing sample labels like `"1422-02216 Breast 1"` into bird
  ID + feather category.
- Computes per-bird, per-feather-category mean δ13C and δ15N.
- Cross-references those means against the δ2H/δ18O values for the matching
  bird/feather-category (`data/isotopes_combined_summary.csv`).
- **Final step (only runs once `outputs/tables/assignment_summary.csv`
  exists, i.e., after Step 4 above):** merges the C/N means with each
  bird/feather-category's Mexico probability, Pacific-Northwest probability,
  peak location, and match-quality flags from the geographic assignment,
  writing `data/cn_vs_geography_comparison.csv`. This is the table used to
  ask "does the independent diet/habitat isotope signal line up with the
  H/O-based geographic assignment?" If `assignment_summary.csv` isn't present
  yet, this step is skipped with a printed message rather than erroring, so
  the script still works as a standalone early look at the C/N data.

**Important: the carbon/nitrogen data are not folded into the spatial model
itself.** There is no established, fine-resolution continental δ13C/δ15N
isoscape the way there is for δ2H/δ18O, so C/N values can't be converted into
a geographic probability surface the way H/O can. Instead, they're used
purely as an independent qualitative cross-check: δ13C indicates C3 vs. C4
vegetation base (all 6 birds here are solidly C3, i.e., freshwater
sedge-marsh, not C4 salt-marsh), and δ15N tracks trophic level / marine vs.
freshwater / agricultural nitrogen enrichment, which is compared against
each bird's H/O-based Mexico probability as a plausibility check, not
combined mathematically with it.

## Full output inventory

| File | What it is |
|---|---|
| `data/klma_clean.csv` | Cleaned master feather-sample table (H/O + metadata) |
| `data/habitat_prior.tif` | Wetland/wet-cropland habitat suitability raster |
| `data/cn_clean.csv` | Cleaned C/N sample-level data |
| `data/cn_summary_by_individual.csv` | C/N means per bird × feather category |
| `data/isotopes_combined_summary.csv` | C/N means joined to H/O means |
| `data/cn_vs_geography_comparison.csv` | C/N means joined to geographic assignment results |
| `outputs/rasters/pd_per_feather.tif` | Raw per-feather posterior probability surfaces |
| `outputs/rasters/pd_per_individual_category.tif` | Combined per-bird × category posteriors |
| `outputs/rasters/pd_population_by_category.tif` | Population-mean Body/Flight posteriors |
| `outputs/rasters/qtl50/75/90_credible_region.tif` | Credible-region rasters |
| `outputs/tables/assignment_summary.csv` | Master results table: peak location, credible areas, zone probabilities, quality flags |
| `outputs/tables/Table1_manuscript.csv` | Cleaned publication-ready summary table |
| `outputs/maps/population_body_vs_flight.png` | Main population-level figure |
| `outputs/maps/individuals_body.png`, `individuals_flight.png` | Main individual-bird figures (unflagged birds) |
| `outputs/maps/individuals_body_supplementary_flagged.png`, `individuals_flight_supplementary_flagged.png` | Supplementary individual-bird figures (flagged birds) |
| `outputs/maps/habitat_prior.png` | Habitat suitability layer diagnostic map |
| `outputs/maps/calibration_regression_*.png` | Calibration regression diagnostic plots |

## Known limitations (worth your colleague's attention)

1. No Rallidae calibration data exist in any standard known-origin isotope
   database — the Passerine-based calibration is the best available proxy,
   not a validated one for rails.
2. The δ18O calibration is weak (R² = 0.08); δ2H is carrying most of the
   discriminating power.
3. Isoscape prediction uncertainty was assumed (uniform 5‰/0.6‰), not
   obtained from the official spatially varying SE grids.
4. The habitat prior is approximate: a single year of 10 m land cover, a
   hand-picked (not statistically fit) weighting scheme, and it cannot
   distinguish flooded/irrigated cropland (real habitat) from dry cropland
   (not habitat) — both get the same partial-credit weight.
5. No age-class data — a July-caught juvenile's feathers reflect the natal
   territory, while an adult's most recent flight feathers were grown during
   the *previous* molt cycle wherever that bird was at the time.
6. Small sample size (12 KLMA individuals + 2 reference birds), and the
   shared-molt-origin assumption for combining same-category feathers failed
   outright for 3 of 13 birds.
7. As detailed above, roughly half of individual-level results carry a
   match-quality flag and should not be cited as precise on their own —
   population-level zone probabilities are the more robust result.

See `REPORT.md` in the repository for the full results narrative and
`MANUSCRIPT_METHODS.md` for the condensed version intended for submission.
