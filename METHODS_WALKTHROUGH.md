# Step-by-step methods walkthrough: KLMA Yellow Rail geographic origin analysis

This document walks through the full pipeline, in the order things actually
need to run, from the two raw lab data sets to the final maps and tables.
It's meant to be read alongside the R code in `scripts/`, so a second reviewer
can follow (and re-run) exactly what was done. For the short publication-style
version, see `MANUSCRIPT_METHODS.md`.

**Major revision note:** this version follows an external code/results
review (via GitHub Copilot) of the previous version of this pipeline. Two
changes were made as a direct result: (1) the two non-KLMA individuals
("San Diego ref bird" and `1272-31742`) are now excluded from the analysis
entirely, since they were never true reference individuals (their actual
breeding/molt origin is unknown, so there was nothing for their assignment
to be validated against); (2) **δ2H alone, not the dual δ2H+δ18O model, is
now the primary assignment model**, because δ18O's calibration here is
weak and nearly collinear with δ2H, and including it was shown to change
one of the two headline results. A new script,
`scripts/06_sensitivity_check.R`, runs five additional model
configurations (dual-isotope, three habitat-prior strengths, and an
alternate calibration reference group) and compares every one against the
primary model — this is now a required part of interpreting the results,
not an optional extra. See Step 3 and Step 7 below, and `REPORT.md`'s
Results section 1b, for what it found.

## The R files, and the order to run them in

The scripts are numbered `00`-`06`, but the numbering is *not* the run
order for a full pipeline, because script `04` has a final step that
depends on script `02`'s output, and script `06` depends on both `01` and
`02`. **Run them in this order:**

| Order | File | What it needs | What it produces |
|---|---|---|---|
| 1 | `scripts/00_clean_data.R` | `data/SIA_Results_KLMA_2026.csv`, `data/SIA_Results_QAQC.csv` | `data/klma_clean.csv`, `data/excluded_individuals.csv` |
| 2 | `scripts/05_build_habitat_prior.R` | internet access (downloads ESA WorldCover tiles), `d2h_GS.tif` | `data/habitat_prior.tif` |
| 3 | `scripts/01_run_assignment.R` | `data/klma_clean.csv`, `data/habitat_prior.tif`, `d2h_GS.tif`, `d18o_GS.tif` | `outputs/rasters/pd_per_feather.tif` (primary model), calibration diagnostic plots |
| 4 | `scripts/02_summarize_assignment.R` | `data/klma_clean.csv`, `outputs/rasters/pd_per_feather.tif` | `outputs/tables/assignment_summary.csv`, per-individual/population posterior rasters, credible-region rasters |
| 5 | `scripts/03_make_maps.R` | everything from step 4 | all files in `outputs/maps/`, `outputs/tables/Table1_manuscript.csv` |
| 6 | `scripts/04_clean_cn_data.R` | `data/klma_clean.csv`, `data/SIBS_CN_230621_Kellerman.csv`, and (for its last step only) `outputs/tables/assignment_summary.csv` | `data/cn_clean.csv`, `data/cn_summary_by_individual.csv`, `data/isotopes_combined_summary.csv`, `data/cn_vs_geography_comparison.csv` |
| 7 | `scripts/06_sensitivity_check.R` | `outputs/tables/assignment_summary.csv`, `data/klma_clean.csv`, `data/habitat_prior.tif` | `outputs/tables/sensitivity_comparison.csv`, `outputs/tables/habitat_prior_by_zone.csv` |

Step 7 reruns the whole assignment under five alternative configurations
(dual-isotope; no habitat prior; weak and strong habitat prior; an
alternate calibration reference group) and compares each against the
primary model from steps 3-4. **This step is required reading before
citing the body-feather result specifically** — see Step 7 and `REPORT.md`
section 1b. It is not required to produce the main maps/tables (those come
from steps 1-6), but the body-feather result from those steps should not be
cited on its own without the context step 7 provides.

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

`scripts/full_pipeline_combined.R` is also available: it's all seven
scripts above concatenated verbatim into one file, already in this exact
run order, for sharing or reviewing the whole pipeline top to bottom in
one place. The seven numbered scripts remain the source of truth for
running or editing any one step on its own.

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
- Also identifies two individuals that are **not part of the KLMA sampling
  effort**: the "San Diego ref bird" rows already present in the KLMA study
  spreadsheet (a different capture/collection context), and band
  `1272-31742`, pulled directly out of the raw SIRFER QC text (present in
  that lab's job batch 24-215.2 but never entered into the KLMA study
  spreadsheet — its band-number prefix, `1272-`, doesn't match the KLMA
  birds' `1422-`). **Both are excluded from the analysis entirely as of
  this version** — see the note at the top of this document for why. Their
  10 combined feather records are written to
  `data/excluded_individuals.csv` for transparency, but nothing downstream
  reads that file.
- Flags (but does not drop) samples with an implausible O:H combustion ratio
  (outside ~3.0–5.5) or implausible %H/%O (possible double-loading). Three
  samples were flagged this way: a Nape sample from bird `1422-02212`
  (WtH/WtO 11.4/39.2) and a Nape sample from bird `1422-02214` (11.6/44.9),
  both flagged for implausibly high %H/%O (possible double-loading); and a
  Rectrix from `1422-02220` (O:H ratio 8.5, far outside the normal ~3.6–4.5
  range).
- Writes `data/klma_clean.csv`: 43 feather samples from the 12 KLMA
  individuals only.

### Step 2 — `05_build_habitat_prior.R`: build the wetland/wet-cropland mask

**Why this step exists:** Yellow Rails are a wetland/wet-agriculture obligate
(Cornell's Birds of the World account: wet sedge meadow, marsh, wet
agricultural areas — not forest, not desert, not high mountains). The very
first version of this analysis used the isotope data alone, and it produced
"credible" origins on mountain peaks and desert terrain, purely because those
places happen to match a bird's δ2H value by elevation. This step builds a
raster that down-weights everywhere that isn't plausible rail habitat, so
that elevation-driven false matches get suppressed.

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
  matters and is explained more in Step 4 below — without it, isolated
  wetland pixels surrounded by hard zeros caused a real problem where several
  unrelated birds' posteriors all collapsed onto the same single pixel.
- The result is force-aligned (same CRS, same exact grid, via
  `resample(..., method="near")`) to the isoscape template, to avoid a subtle
  floating-point/CRS-string mismatch that otherwise breaks a downstream
  `assignR` sanity check.
- Writes `data/habitat_prior.tif`. Of the 1,209,600 cells across the full
  North American grid, 101,813 (~8.4%) are nonzero after smoothing — most
  of the continent stays exactly zero. Within the narrower 8-state-plus-
  Mexico study area itself, the Gaussian kernel's finite radius means
  almost every cell picks up *some* nonzero value, so it isn't a hard
  binary mask there — but the distribution stays heavily skewed toward
  near-zero (median suitability ~0.003; only about a third of study-area
  cells exceed 0.01), so it still meaningfully down-weights most of the
  study area in practice.
- **New this version** (used by Step 7, not here): two additional variants
  of this same prior are built by a simple power transform —
  `habitat_prior^0.5` ("weak," compressed toward a uniform/no-op prior) and
  `habitat_prior^2` ("strong," exaggerated toward a hard mask) — to test how
  much the prior's *strength*, independent of its underlying
  classification, drives the results.

### Step 3 — `01_run_assignment.R`: the core isotope assignment (primary model)

This is the main `assignR` step: converting each feather's measured δ2H
into a map of "how likely is it this feather grew at pixel X," for every
pixel in the study area.

- **Primary model: δ2H alone.** An earlier version of this analysis used
  the dual-isotope (δ2H+δ18O) model as primary. That changed for two
  reasons, found on review: (1) δ18O's calibration here is weak (R² =
  0.08), and (2) the raw δ2H and δ18O growing-season isoscapes are r² =
  0.99 correlated across this study area — `assignR::pdRaster()` itself
  warns about exactly this condition ("two or more isoscapes have shared
  variance > 0.7 ... specificity of assignments may be inflated"). Treating
  a weak, highly collinear second marker as independent evidence can
  narrow posteriors and shift zone-level probabilities in a way that isn't
  statistically justified. **Quantified directly in Step 7**: switching
  between δ2H-only and dual-isotope changes the population-level
  body-feather result from Pacific-Northwest-majority (primary model) to
  Mexico-majority (dual-isotope) — a large enough swing that the
  weaker-calibrated, collinear model should not be the default. The
  dual-isotope model is still run, now explicitly as one of five named
  sensitivity checks in Step 7, not folded into the headline result.
- **Study area.** California, Nevada, Arizona, New Mexico, Oregon,
  Washington, Idaho, Utah, plus Baja California and the western Mexican
  coast (bounded roughly 14–33° N, 101–118° W for the Mexican portion).
  Your original frame of reference was CA/NV/AZ/NM south through
  Baja/western Mexico; Oregon (plus WA/ID/UT for geographic continuity)
  was added because Klamath Marsh NWR itself — the capture/breeding site —
  sits at ~43° N, just north of the California line, and Cornell's account
  says flight feathers grow primarily on or near the breeding grounds.
  Excluding Oregon would have ruled out the single most biologically
  expected flight-feather origin before the isotope data got a say.
- **Isoscape.** Uses the supplied growing-season precipitation δ2H raster
  (`d2h_GS.tif`) as the base surface. This is a mean-only raster with no
  prediction-uncertainty band (the official SE grid lives at
  wateriso.utah.edu, which wasn't reachable from this environment). A
  spatially uniform standard deviation is assumed instead (5‰), sized to
  be roughly consistent with `assignR`'s own bundled example uncertainty
  raster. This is a real approximation — see limitations. (The δ18O
  isoscape, `d18o_GS.tif`, is still loaded and calibrated here too, purely
  so Step 7's dual-isotope sensitivity variant has it available.)
- **Calibration (precipitation isotope value → feather keratin isotope
  value).** `assignR`'s reference database of species with known-origin
  isotope measurements contains no Rallidae (rails) at all — so a
  rail-specific calibration isn't possible with the tools/data available.
  The North-America-masked **"Passerine"** group (n = 699 known-origin
  sites) is used as the primary model's calibration (R² = 0.76). A second
  calibration using the **"Water bird"** reference group (dabbling/diving
  waterfowl — ecologically closer to a marsh-dwelling rail than a generic
  passerine; R² = 0.55, similar slope, different intercept) is also built
  here and used in Step 7's calibration-sensitivity check — it is *not*
  used in the primary model, but Step 7 found that switching to it changes
  the body-feather result about as much as adding δ18O does, which is an
  important finding in its own right (see Step 7).
- **Habitat prior integration — and a real `assignR` bug found and worked
  around.** The habitat suitability raster from Step 2 is supplied to
  `assignR::pdRaster()` as a multiplicative prior, so the final posterior at
  each pixel is (isotope-match likelihood) × (habitat suitability). Testing
  this the straightforward way (passing both `mask` and `prior` to
  `pdRaster()`) produced every single sample's posterior as exactly zero
  everywhere. Diagnosis: `pdRaster.isoStack()` crops the isoscape to `mask`
  internally, but never applies that same crop to `prior` — so the
  elementwise `likelihood * prior` multiplication silently multiplies two
  differently-sized/misaligned vectors. The fix used here: pre-crop and
  pre-mask both the calibrated isoscape and the habitat prior to the study
  area ourselves, then call `pdRaster(..., mask = NULL, prior =
  habitat_prior)` so nothing gets cropped a second time inside the
  function.
- **Sanity check performed:** in the resulting population body-feather
  surface, Mt. Whitney (a 4,421 m alpine peak with zero rail habitat) samples
  at a posterior probability of ~8.5×10⁻¹⁵, while Sacramento National
  Wildlife Refuge (a real Central Valley wetland complex) samples at
  ~6.2×10⁻⁴ — more than ten orders of magnitude apart, confirming the
  habitat prior is meaningfully reshaping the result rather than being a
  no-op.
- Runs `pdRaster()` once per feather sample (43 total, all from the 12
  KLMA individuals — the two non-KLMA individuals from Step 1 are gone from
  the dataset by this point, so nothing further is needed to exclude
  them), producing one posterior-probability raster per feather under the
  primary (δ2H-only) model. Writes the full stack to
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
  reviewed rather than silently masked. **Under the primary (δ2H-only)
  model, this did not trigger for any of the 20 multi-feather
  combinations** — all succeeded via the standard joint-product method.
  (Under the earlier dual-isotope model, 4 combinations had failed this
  way, including bird `1422-02212`'s flight feathers, plausibly linked to
  its QC-flagged Nape sample from Step 1.) This is expected, not
  reassuring on its own: a single, less spatially discriminating marker
  produces broader per-feather posteriors that are more likely to overlap
  somewhere in the study area even when the underlying feathers would have
  looked inconsistent under the sharper dual-isotope posteriors, so this
  check is simply weaker under the primary model, not passed more
  convincingly.
- **Population summaries.** An unweighted mean posterior is computed
  separately for Body and Flight categories across all 12 KLMA individuals,
  each renormalized to sum to 1 over the study area.
- **Credible regions.** For every individual×category surface and both
  population surfaces, computes the 50%, 75%, and 90% credible regions using
  `assignR::qtlRaster(thresholdType = "prob")` — the smallest-area region
  containing that much cumulative posterior probability.
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
  - `single_cell_dominant_flag`: `peak_share` > 0.3.
  - `shared_peak_flag`: this surface's peak pixel is the *exact same*
    pixel as another, unrelated bird's peak pixel.
- Writes `outputs/tables/assignment_summary.csv` (the master results table)
  and the underlying rasters (`pd_per_individual_category.tif`,
  `pd_population_by_category.tif`, `qtl50/75/90_credible_region.tif`).

**Why the flag columns matter — please read before trusting any single
bird's "precise" result, and note this is now worse, not better, than in
the previous (dual-isotope-primary) version.** Under the primary δ2H-only
model, **17 of 22 individual-category surfaces (77%) carry at least one
match-quality flag**, up from 14 of 26 (54%) under the old dual-isotope
model. This is an expected trade-off of moving to the more statistically
defensible single-marker model, not a new bug: δ2H alone is less spatially
discriminating than δ2H+δ18O, so more individuals' weakly-informative
likelihood surfaces default to whichever single cell the habitat prior
favors most locally — and because the habitat mask is sparse, many
unrelated individuals converge on the exact same cell. The practical
upshot is the same as before, just more so: **population-level zone
probabilities are the trustworthy output of this analysis; individual-level
peak locations and credible areas, under the primary model, mostly are
not.**

### Step 5 — `03_make_maps.R`: final maps and the publication table

- Classifies every surface into 4 discrete color classes directly from the
  credible-region rasters from Step 4 (top 50% / 50–75% / 75–90% / outside
  the 90% region), rather than by relative pixel-value quantile.
- Produces the main population-level Body-vs-Flight comparison figure
  (`population_body_vs_flight.png`), now titled to note it's the δ2H-only
  primary model.
- Produces individual-bird map grids, split into a main figure (unflagged
  birds only) and a supplementary figure (flagged birds, each panel labeled
  "[low match quality]") for both Body and Flight. Under the primary model
  this split is lopsided — most birds land in the supplementary/flagged
  figure (see Step 4) — which is itself informative about how much weight
  individual-level maps should be given.
- Produces a habitat-prior diagnostic map (`habitat_prior.png`).
- Produces `outputs/tables/Table1_manuscript.csv`, a cleaned publication
  table: one row per bird × feather category (all KLMA, no reference-bird
  column any more), with peak location, 50%/90% credible area, %
  probability mass in Mexico vs. the Pacific Northwest, and a plain-text
  quality note for flagged rows.

### Step 6 — `04_clean_cn_data.R`: carbon/nitrogen cross-check

- Cleans the SIBS lab's carbon/nitrogen results (a separate lab and separate
  analytical run from the SIRFER δ2H/δ18O data, but keyed to the same
  envelope IDs), parsing sample labels like `"1422-02216 Breast 1"` into bird
  ID + feather category.
- **Excludes individual `1272-31742`** from the C/N data too (it's a
  separate raw file from the H/O data, so this needed its own explicit
  filter — see the script's comments), for the same reason it's excluded
  everywhere else in this version.
- Computes per-bird, per-feather-category mean δ13C and δ15N (5 KLMA
  individuals, 14 samples).
- Cross-references those means against the δ2H values for the matching
  bird/feather-category (`data/isotopes_combined_summary.csv`).
- **Final step (only runs once `outputs/tables/assignment_summary.csv`
  exists):** merges the C/N means with each bird/feather-category's Mexico
  probability, Pacific-Northwest probability, peak location, and
  match-quality flags from the geographic assignment, writing
  `data/cn_vs_geography_comparison.csv` (4 rows under the current KLMA-only
  dataset: 3 individuals with body-feather data, 1 with flight-feather
  data).

**Important: the carbon/nitrogen data are not folded into the spatial model
itself.** There is no established, fine-resolution continental δ13C/δ15N
isoscape the way there is for δ2H, so C/N values can't be converted into a
geographic probability surface. They're used purely as an independent
qualitative cross-check: δ13C indicates C3 vs. C4 vegetation base (all 6
KLMA birds here are solidly C3), and δ15N tracks trophic level / marine
vs. freshwater / agricultural nitrogen enrichment — but with only 3-4 data
points for the geography comparison, this doesn't support a trend claim
either way (see `REPORT.md`).

### Step 7 — `06_sensitivity_check.R`: does the body-feather result survive reasonable alternative choices?

This is the step added in response to external review, and it's required
reading before citing the body-feather result. It reruns the full
assignment under five alternative configurations and compares each to the
primary model (δ2H, Passerine calibration, current habitat prior),
individual-by-individual and at the population level:

| Variant | What changes | Why it's tested |
|---|---|---|
| `dual_isotope` | adds δ18O | the model this pipeline used to treat as primary; see Step 3 |
| `no_prior` | habitat prior removed entirely | tests whether the prior is doing any work at all |
| `weak_prior` | habitat_prior^0.5 (compressed) | tests sensitivity to prior *strength* |
| `strong_prior` | habitat_prior^2 (exaggerated) | tests sensitivity to prior *strength*, other direction |
| `waterbird_calibration` | Water bird reference group instead of Passerine | tests sensitivity to calibration-group choice |

**What it found** (population-level body feathers; flight feathers are
stable — Pacific-Northwest-majority — across every variant, so the table
below is about the unresolved result specifically):

| Variant | % Mexico | % Pacific NW | 50%-credible area (km²) |
|---|---|---|---|
| Primary (δ2H, Passerine, current prior) | 34.6 | 59.6 | 112,407 |
| Dual-isotope | 53.3 | 36.3 | 70,477 |
| No habitat prior | 27.8 | 40.8 | 1,172,234 |
| Weak habitat prior | 34.2 | 55.0 | 260,612 |
| Strong habitat prior | 33.4 | 64.1 | 46,672 |
| Water-bird calibration | 50.4 | 42.2 | 111,795 |

Two things stand out:

1. **Habitat-prior strength (weak vs. strong, a four-fold range) barely
   moves the body-feather split** (33.4–34.2% Mexico, all close to the
   primary model's 34.6%) — it mainly controls credible-region *size*
   (46,672 to 260,612 km²), not which zone wins. This directly answers a
   concern raised in review: the subjective habitat-weighting scheme is
   *not* the main source of instability in this result.
2. **The isotope-marker and calibration-group choices are what actually
   move it**: both the dual-isotope and Water-bird-calibration variants
   push the body-feather result to a Mexico-majority (53.3%, 50.4%),
   while the Passerine-based, prior-strength variants all agree on a
   Pacific-Northwest-leaning result (55.0–64.1%). Four of six variants
   lean Pacific Northwest; two lean Mexico. **This is genuinely unresolved
   with the present data** — see `REPORT.md` Results section 1b for the
   full discussion and individual-level numbers.

The script also checks (independent of any isotope data) whether the
habitat prior itself is balanced between the two zones under comparison:
39.3% of its total suitability mass falls in the Mexico zone versus 49.6%
in the Pacific Northwest zone (`outputs/tables/habitat_prior_by_zone.csv`)
— not dramatically skewed toward either, so the prior itself isn't
manufacturing the Mexico-vs-Pacific-Northwest disagreement among the other
variants.

## Full output inventory

| File | What it is |
|---|---|
| `data/klma_clean.csv` | Cleaned master feather-sample table, KLMA individuals only |
| `data/excluded_individuals.csv` | Raw records for the 2 excluded non-KLMA individuals (not used downstream) |
| `data/habitat_prior.tif` | Wetland/wet-cropland habitat suitability raster |
| `data/cn_clean.csv` | Cleaned C/N sample-level data, KLMA individuals only |
| `data/cn_summary_by_individual.csv` | C/N means per bird × feather category |
| `data/isotopes_combined_summary.csv` | C/N means joined to H/O means |
| `data/cn_vs_geography_comparison.csv` | C/N means joined to geographic assignment results |
| `outputs/rasters/pd_per_feather.tif` | Raw per-feather posterior probability surfaces, primary model |
| `outputs/rasters/pd_per_individual_category.tif` | Combined per-bird × category posteriors, primary model |
| `outputs/rasters/pd_population_by_category.tif` | Population-mean Body/Flight posteriors, primary model |
| `outputs/rasters/qtl50/75/90_credible_region.tif` | Credible-region rasters, primary model |
| `outputs/tables/assignment_summary.csv` | Master results table: peak location, credible areas, zone probabilities, quality flags |
| `outputs/tables/Table1_manuscript.csv` | Cleaned publication-ready summary table |
| `outputs/tables/sensitivity_comparison.csv` | All 6 model variants compared, every individual and population surface |
| `outputs/tables/habitat_prior_by_zone.csv` | Habitat-prior zone-balance check |
| `outputs/maps/population_body_vs_flight.png` | Main population-level figure |
| `outputs/maps/individuals_body.png`, `individuals_flight.png` | Main individual-bird figures (unflagged birds) |
| `outputs/maps/individuals_body_supplementary_flagged.png`, `individuals_flight_supplementary_flagged.png` | Supplementary individual-bird figures (flagged birds) |
| `outputs/maps/habitat_prior.png` | Habitat suitability layer diagnostic map |
| `outputs/maps/calibration_regression_*.png` | Calibration regression diagnostic plots |

## Known limitations (worth your colleague's attention)

1. No Rallidae calibration data exist in any standard known-origin isotope
   database — the Passerine-based calibration is the best available proxy,
   not a validated one for rails. Step 7 shows this choice matters:
   switching to a Water-bird reference group changes the body-feather
   population result from Pacific-Northwest-majority to Mexico-majority.
2. δ18O is weakly calibrated (R² = 0.08) **and** nearly collinear with δ2H
   (r² = 0.99 in the raw isoscape) — this is why it was removed from the
   primary model. Step 7 shows its inclusion alone is enough to flip the
   body-feather population result.
3. Isoscape prediction uncertainty was assumed (uniform 5‰), not obtained
   from the official spatially varying SE grid.
4. The habitat prior is approximate: a single year of 10 m land cover, a
   hand-picked (not statistically fit) weighting scheme, and it cannot
   distinguish flooded/irrigated cropland (real habitat) from dry cropland
   (not habitat). Step 7's prior-strength sensitivity check found this is
   *not* the main driver of the body-feather result's instability, though.
5. No age-class data, and separately, the preformative (body) molt's
   timing relative to migration/wintering is itself variable per Cornell's
   account — body feathers may represent breeding grounds, migration
   stopovers, or wintering areas depending on age and molt-completion
   timing. An unresolved population-level body-feather signal is
   biologically plausible for this reason, not just a modeling artifact.
6. Small sample size: 12 KLMA individuals only (as of this version). The
   shared-molt-origin assumption for combining same-category feathers
   was not detectably violated for any individual under the primary
   (d2H-only) model (0 of 20 multi-feather combinations failed) — but
   this check is weaker under a single, less discriminating marker than
   it was under the dual-isotope model (where 4 combinations failed
   outright), so this should not be read as stronger support for the
   shared-origin assumption, just a less sensitive test of it.
7. Under the primary model, 77% of individual-level results (17 of 22)
   carry a match-quality flag and should not be cited as precise on their
   own — population-level zone probabilities are the robust result, and
   even those are only robust for flight feathers; the body-feather result
   is genuinely unresolved (Step 7).

See `REPORT.md` in the repository for the full results narrative and
`MANUSCRIPT_METHODS.md`/`MANUSCRIPT_RESULTS.md` for the condensed,
publication-style versions.
