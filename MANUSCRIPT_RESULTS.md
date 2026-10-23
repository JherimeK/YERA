# Results (draft, for manuscript — target journal: *Ornithology* or *The Wilson Journal of Ornithology*)

This is a full Results section, written to follow directly from
`MANUSCRIPT_METHODS.md`, with every statistic pulled from the current
pipeline output. This version reflects a substantial revision: δ2H alone
(not the dual δ2H+δ18O model) is now the primary assignment model, and two
individuals previously treated as "reference" individuals have been
removed from the dataset entirely (see `MANUSCRIPT_METHODS.md`'s revision
note). Table and figure numbers are placeholders (T1, T2, F1, F2...) for
you to renumber once the manuscript's full figure/table order is set.
Bracketed items are still open and need your input before submission.

---

## Results

### Sample composition and data quality

Stable hydrogen and oxygen isotope ratios (δ2H, δ18O) were obtained from
43 feather samples from 12 individually banded Yellow Rails captured at
Klamath Marsh National Wildlife Refuge (18 body feathers, 25 flight
feathers). Three feather samples (7.0% of the total) were flagged during
quality control: two Nape samples (individuals 1422-02212 and 1422-02214)
for implausibly high %H/%O yield consistent with possible double-loading,
and one Rectrix (1422-02220) for an O:H combustion ratio (8.5) well
outside the expected keratin range (~3.6–4.5). Flagged samples were
retained in the dataset but are noted where relevant below.

Two additional individuals sampled in the same laboratory batches were
excluded from the analysis entirely (see Methods): neither had an
independently known breeding or molt origin, so there was no basis to
treat them as reference individuals.

Where more than one feather of the same category (body or flight) was
available for an individual, their posterior probability surfaces were
combined by multiplication and renormalization under a shared-origin
assumption (20 such combinations; 2 individuals had only a single feather
in one category and needed no combination). Under the primary (δ2H-only)
model, this produced a numerically consistent joint posterior for all 20
combinations (0 failures). This is a weaker check than it would be under
a more spatially discriminating marker set: under the dual-isotope
sensitivity model (below), 4 of the same combinations failed this check
(no numerical overlap between the individual feathers' posteriors), so
the absence of failures under the primary model should not be read as
stronger support for the shared-origin assumption — only that δ2H alone
is less likely to detect a violation of it.

### Isotope calibration

Because no Rallidae, and essentially no North American non-passerine,
known-origin data exist in the reference database used for calibration,
feather values were calibrated against the precipitation isoscape using a
North America-masked Passerine reference group. The δ2H calibration (n =
699 known-origin sites) explained approximately 76% of variance in tissue
isotope values (R² ≈ 0.76, F1,697 ≈ 2,190, *P* < 0.001; slope ≈ 0.71,
intercept ≈ −14‰). This is the calibration used in the primary model. The
δ18O calibration (n = 167 sites) was comparatively weak (R² ≈ 0.08, F1,165
≈ 14.4, *P* < 0.001; slope ≈ 0.22, intercept ≈ 15.4‰) and is used only in
the dual-isotope sensitivity variant (below), not the primary model. A
further sensitivity calibration using a Water bird reference group (δ2H
only, n = 75 sites) produced a qualitatively similar fit (R² ≈ 0.55,
F1,73 ≈ 88.7, *P* < 0.001; slope ≈ 0.71, intercept ≈ −23‰) and is used in
one of the sensitivity variants below, not the primary model.
[Calibration-regression coefficients showed small (≤0.3 R²-point)
run-to-run variation in this computing environment, consistent with
internal resampling in assignR's calibration routine; values above are
rounded accordingly and should not be over-interpreted to more than 2
significant figures.]

### Habitat suitability prior

After reclassification of ESA WorldCover 10 m land cover (85 of 108
candidate tiles contained land within the study area) and Gaussian
smoothing, 101,813 of 1,209,600 grid cells (8.4%) retained nonzero
wetland/wet-cropland suitability across the full North American domain;
within the narrower study area, suitability values remained strongly
right-skewed (median 0.003; only 32.1% of study-area cells exceeded
0.01). The prior's effect on the posterior was validated against two
locations of known, opposite habitat status: in the population-level
body-feather surface, Mt. Whitney (4,421 m alpine peak; no rail habitat)
received a posterior probability of 8.5 × 10⁻¹⁵, compared with 6.2 × 10⁻⁴
at Sacramento National Wildlife Refuge (a documented Central Valley
wetland complex) — a difference exceeding ten orders of magnitude.

As an additional check (prompted by the possibility that the prior's
suitability classification could be systematically biased toward one
geographic region, independent of any isotope signal), the prior's own
distribution of mass and area was tabulated between the two zones under
comparison: 39.3% of its total suitability mass fell within the Mexico
zone (Baja California + mainland western Mexico) versus 49.6% within the
Pacific Northwest zone (Oregon, Washington, Idaho, Utah, California);
suitable area (prior > 0.01) was comparable between the two (459,286 km²
in Mexico vs. 512,430 km² in the Pacific Northwest). The prior is
therefore not dramatically skewed toward either region.

### Geographic assignment: body versus flight feathers

Population-level posterior surfaces (equal-weighted means across the 12
individuals, primary model: δ2H only, Passerine calibration,
habitat-constrained) assigned different probable origins to body and
flight feathers, but less starkly than an earlier dual-isotope analysis
had suggested (Table T1; Figure F1). Flight feathers were assigned a
large majority of posterior probability (98.2%) to the Pacific Northwest
(Oregon, Washington, Idaho, Utah, California combined), versus 0.8% to
Mexico (Baja California + mainland western Mexico); the 50%-credible
region covered 27,858 km². Body feathers leaned toward the Pacific
Northwest as well, but far less decisively (59.5% Pacific Northwest vs.
34.7% Mexico; 50%-credible region 112,733 km²). At the individual level,
12 of 12 individuals with flight-feather data had a majority (>50%) of
their own posterior probability in the Pacific Northwest (range
89.0–99.5%, median 99.1%), while only 2 of 10 individuals with
body-feather data had a majority in Mexico (range 0.04–89.3%, median
37.6%; one of the two, 1422-02214, is a marginal case at 50.5%, not a
clear majority) — the body-feather population-level lean toward the
Pacific Northwest is driven by most individuals sitting in an ambiguous
30–60% range rather than a consistent majority either way (Table T1; Figures F2,
F3).

### Sensitivity analysis: the body-feather result is not robust to isotope-marker or calibration-group choice

The primary model's result for flight feathers was essentially unchanged
across five alternative configurations tested (Table T2): a dual-isotope
(δ2H+δ18O) model, a model with no habitat prior, "weak" and "strong"
habitat-prior-strength variants (prior raised to the 0.5 and 2 power,
respectively), and a Water-bird-calibration variant. Across all six
configurations (including the primary model), the flight-feather
population-level Pacific Northwest probability ranged only 74.1–99.4% (and
Mexico probability 0.01–4.0%).

The body-feather result was not stable (Table T2). Two of the five
sensitivity variants reversed the primary model's Pacific-Northwest lean
into a Mexico majority: the dual-isotope model (53.1% Mexico vs. 36.6%
Pacific Northwest) and the Water-bird-calibration model (50.5% vs. 42.2%).
The habitat-prior-strength variants, by contrast, stayed close to the
primary model regardless of direction: weak prior gave 34.3% Mexico /
55.0% Pacific Northwest, and strong prior gave 33.5% / 64.1% — a
four-fold range in prior strength moved the Mexico probability by only
about 1 percentage point. Removing the habitat prior entirely dropped
Mexico probability further (27.9%) and increased the 50%-credible area
from 112,733 km² to 1,176,595 km², confirming the prior's necessity for
any usable precision even though it was not the primary driver of the
Mexico-vs-Pacific-Northwest disagreement among variants. In sum: **the
choice of isotope marker(s) and calibration reference group, not the
habitat-prior weighting scheme, determines whether the body-feather result
reads as Pacific-Northwest-leaning or Mexico-leaning**, and the present
data do not resolve which of these choices is correct.

### Reliability screening of individual-level assignments

Individual-category posterior surfaces were screened for a match-quality
artifact in which spatially flat regions of the isoscape, combined with
the sparseness of the habitat mask, produce coincident credible regions
across unrelated individuals. Under the primary model, 17 of 22
individual-category surfaces (77%) shared an identical peak-probability
grid cell with at least one other individual's surface
(`shared_peak_flag`); no surface had a single cell accounting for more
than 30% of its posterior probability (`single_cell_dominant_flag`: 0 of
22), and no surface had a 50%-credible area under 1,000 km²
(`low_precision_flag`: 0 of 22). This is a larger fraction of flagged
surfaces than under the dual-isotope model evaluated as a sensitivity
check (14 of 26, 54%), because δ2H alone is less spatially discriminating,
making individual-level peak locations and credible areas under the
primary model a weaker basis for location claims than population-level
zone probabilities.

### Carbon and nitrogen isotopes

Carbon and nitrogen isotope values were obtained from 14 feather samples
across 5 of the 12 individuals. δ13C values were corrected for a
systematic +2.05‰ enrichment bias identified from the analytical lab's
in-house QA standards (two independent plant standards, both reading
δ13C enriched relative to known consensus values by a consistent ~2‰
despite spanning a ~27‰ range: +2.10‰ and +2.01‰); δ15N showed no
comparable bias and was not corrected. Corrected δ13C values across all
individual-by-category means ranged from −27.6‰ to −24.3‰, consistent
with a C3-plant-based diet throughout and showing no individual
approaching the C4 signature (approximately −12‰ to −17‰) associated
with *Spartina*/*Distichlis* salt-marsh vegetation — the correction moves
every value further from the C4 range, not closer, so this conclusion is
strengthened rather than altered by it. δ15N values ranged from
7.05‰ to 15.35‰ across individual means; one individual (1422-02220)
showed marked within-individual heterogeneity, with four body-feather
subsamples splitting between a freshwater-typical value (8.1‰) and three
elevated values (14.9–15.3‰). A comparison of δ15N against the
geographic assignment, restricted to the 3 individuals with both
measurements available for body feathers, did not show a consistent
relationship and is too small a sample to support a trend claim in either
direction.

---

## Table T1. Per-individual and population-level summary (primary model)

*(from `outputs/tables/Table1_manuscript.csv`)*

| Bird | Feathers | n | Peak (lon, lat) | 50%-CA (km²) | 90%-CA (km²) | % Mexico | % Pacific NW | Note |
|---|---|---|---|---|---|---|---|---|
| 1422-02211 | Body | 2 | −121.96, 39.38 | 63,840 | 293,712 | 39.5 | 54.9 | Low match quality |
| 1422-02211 | Flight | 2 | −117.12, 47.29 | 19,420 | 94,692 | 0.0 | 98.8 | Low match quality |
| 1422-02212 | Body | 2 | −118.96, 46.96 | 43,149 | 220,674 | 5.7 | 92.4 | Low match quality |
| 1422-02212 | Flight | 2 | −118.96, 46.96 | 49,321 | 248,034 | 8.6 | 89.0 | Low match quality |
| 1422-02213 | Body | 2 | −108.38, 25.54 | 16,464 | 105,902 | 89.3 | 0.7 | |
| 1422-02213 | Flight | 2 | −117.29, 47.12 | 23,827 | 99,277 | 0.0 | 99.3 | Low match quality |
| 1422-02214 | Body | 2 | −122.04, 39.29 | 61,898 | 287,742 | 50.5 | 43.0 | |
| 1422-02214 | Flight | 2 | −117.29, 47.12 | 27,710 | 100,757 | 0.0 | 99.4 | Low match quality |
| 1422-02215 | Body | 2 | −122.04, 39.38 | 61,773 | 290,955 | 44.4 | 49.6 | Low match quality |
| 1422-02215 | Flight | 2 | −117.12, 47.29 | 20,129 | 95,716 | 0.0 | 98.9 | Low match quality |
| 1422-02216 | Body | 1 | −122.04, 39.38 | 198,287 | 886,452 | 46.8 | 42.0 | Low match quality |
| 1422-02216 | Flight | 2 | −117.29, 47.12 | 22,903 | 98,881 | 0.0 | 99.2 | Low match quality |
| 1422-02217 | Body | 2 | −121.96, 39.38 | 66,802 | 296,251 | 35.6 | 59.1 | Low match quality |
| 1422-02217 | Flight | 2 | −117.29, 47.12 | 24,161 | 99,613 | 0.0 | 99.3 | Low match quality |
| 1422-02218 | Body | 2 | −118.96, 46.96 | 43,685 | 224,169 | 6.0 | 92.0 | Low match quality |
| 1422-02218 | Flight | 2 | −117.29, 46.96 | 31,637 | 129,138 | 0.8 | 98.4 | |
| 1422-02219 | Body | 2 | −117.29, 47.04 | 28,848 | 102,057 | 0.0 | 99.4 | |
| 1422-02219 | Flight | 2 | −117.29, 47.12 | 20,488 | 96,296 | 0.0 | 99.0 | Low match quality |
| 1422-02220 | Flight | 2 | −117.12, 47.29 | 18,931 | 94,002 | 0.0 | 98.7 | Low match quality |
| 1422-02221 | Body | 1 | −121.96, 39.38 | 197,771 | 911,324 | 29.0 | 61.9 | Low match quality |
| 1422-02221 | Flight | 2 | −117.29, 47.12 | 23,134 | 98,987 | 0.0 | 99.2 | Low match quality |
| 1422-02251 | Flight | 3 | −118.88, 46.96 | 19,192 | 71,713 | 0.2 | 99.5 | |
| **All KLMA (population mean)** | Body | 10 | −121.96, 39.38 | 112,733 | 483,474 | **34.7** | 59.5 | |
| **All KLMA (population mean)** | Flight | 12 | −117.29, 47.12 | 27,858 | 125,774 | 0.8 | **98.2** | |

"% Mexico" = Baja California + mainland western Mexico zone probabilities
combined. "% Pacific NW" = Oregon + Washington/Idaho/Utah + California
zone probabilities combined. 50%-CA / 90%-CA = area of the 50%/90%
credible region. "Low match quality" = `shared_peak_flag` (see Methods;
`single_cell_dominant_flag` and `low_precision_flag` did not trigger for
any surface under this model).

## Table T2. Population-level results across all six model configurations

| | Body: % Mexico | Body: % Pacific NW | Body: 50%-CA (km²) | Flight: % Mexico | Flight: % Pacific NW |
|---|---|---|---|---|---|
| **Primary (δ2H, Passerine, current prior)** | **34.7** | **59.5** | 112,733 | 0.8 | **98.2** |
| Dual-isotope | 53.1 | 36.6 | 68,257 | 0.01 | 90.0 |
| No habitat prior | 27.9 | 40.8 | 1,176,595 | 0.7 | 74.1 |
| Weak habitat prior | 34.3 | 55.0 | 261,107 | 0.9 | 93.0 |
| Strong habitat prior | 33.5 | 64.1 | 46,927 | 0.6 | 99.4 |
| Water-bird calibration | 50.5 | 42.2 | 111,772 | 4.0 | 94.6 |

Source: `outputs/tables/sensitivity_comparison.csv` (full per-individual
version also in that file).

---

## Notes for Discussion (not drafted here)

- The flight-feather result (Pacific Northwest origin, consistent with
  breeding-ground molt) is robust across all six model configurations and
  should be reported as a confirmed finding.
- The body-feather result is **not** robust: it depends on whether δ18O
  or a Water-bird calibration reference group is used, and the present
  data do not resolve which choice is more defensible. We recommend
  reporting the full range (Table T2) rather than a single percentage,
  and explicitly stating that the body-feather geographic origin remains
  an open question.
- The body-vs-flight *contrast* itself (flight consistently and
  substantially more Pacific-Northwest-concentrated than body, across
  every one of the six configurations) is the most defensible finding in
  this dataset and is consistent with Cornell's account that flight
  feathers molt on/near the breeding grounds while body (preformative)
  molt timing relative to migration/wintering is more variable.
- No age-class data were available; this limits how confidently either
  result can be tied to a single molt event.
- The two individuals removed from this version (formerly "reference"
  individuals) should not be reintroduced without independent knowledge of
  their breeding/molt origin, or some other basis for treating them as a
  genuine comparison point.
- Carbon/nitrogen results support a freshwater (not salt-marsh) diet
  across all individuals, but the small cross-reference sample (n = 3 for
  body feathers) does not support a quantitative diet-geography
  relationship claim.

## Figure/Table inventory referenced above

| Placeholder | Source file |
|---|---|
| Figure F1 | `outputs/maps/population_body_vs_flight.png` |
| Figure F2 | `outputs/maps/individuals_body.png` |
| Figure F3 | `outputs/maps/individuals_flight.png` |
| Figure S1 (supp.) | `outputs/maps/habitat_prior.png` |
| Figure S2 (supp.) | `outputs/maps/individuals_body_supplementary_flagged.png`, `individuals_flight_supplementary_flagged.png` |
| Figure S3 (supp.) | `outputs/maps/calibration_regression_d2H.png`, `_d18O.png`, `_d2H_waterbird_sensitivity.png` |
| Table T1 | `outputs/tables/Table1_manuscript.csv` |
| Table T2 | `outputs/tables/sensitivity_comparison.csv` (population rows) |
| Table S1 (supp.) | `outputs/tables/sensitivity_comparison.csv` (full per-individual version) |
| Table S2 (supp.) | `outputs/tables/habitat_prior_by_zone.csv` |
| Table S3 (supp.) | `data/cn_vs_geography_comparison.csv` |
