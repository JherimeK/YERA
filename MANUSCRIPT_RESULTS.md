# Results (draft, for manuscript — target journal: *Ornithology* or *The Wilson Journal of Ornithology*)

This is a full Results section, written to follow directly from
`MANUSCRIPT_METHODS.md`, with every statistic pulled from the current
pipeline output rather than estimated or recalled. Numbers were
regenerated directly from `outputs/tables/assignment_summary.csv`,
`data/cn_summary_by_individual.csv`, `data/cn_vs_geography_comparison.csv`,
and two additional sensitivity-analysis reruns described below, as part of
the same end-to-end pipeline check described in `REPORT.md`. Table and
figure numbers are placeholders (T1, T2, F1, F2...) for you to renumber
once the manuscript's full figure/table order is set. Bracketed items are
still open ([YEAR], [capture dates]) and need your input before
submission.

---

## Results

### Sample composition and data quality

Stable hydrogen and oxygen isotope ratios (δ2H, δ18O) were obtained from
53 feather samples: 43 from 12 individually banded Yellow Rails captured
at Klamath Marsh National Wildlife Refuge ("KLMA" individuals; 18 body
feathers, 25 flight feathers), plus two non-KLMA reference individuals
included for comparison — a "San Diego" bird (3 feathers) and an
additional individual (band 1272-31742, not part of the KLMA capture
series; 7 feathers). Three feather samples (5.7% of the total) were
flagged during quality control: two Nape samples (individuals 1422-02212
and 1422-02214) for implausibly high %H/%O yield consistent with possible
double-loading, and one Rectrix (1422-02220) for an O:H combustion ratio
(8.5) well outside the expected keratin range (~3.6–4.5). Flagged samples
were retained in the dataset but are noted where relevant below.

Where more than one feather of the same category (body or flight) was
available for an individual, their posterior probability surfaces were
combined by multiplication and renormalization under a shared-origin
assumption (n = 23 such combinations). This produced a numerically
consistent joint posterior in 19 of 23 cases (82.6%); in the remaining 4
cases (17.4%: individuals 1422-02212's flight feathers, 1422-02219's body
feathers, 1422-02251's flight feathers, and reference individual
1272-31742's flight feathers) the two feathers' posteriors shared no
numerical overlap anywhere in the study area, indicating isotopically
inconsistent values between feathers of the same nominal category; these
were instead represented by the simple mean of the two surfaces and are
flagged in Table T1.

### Isotope calibration

Because no Rallidae, and essentially no North American non-passerine,
known-origin δ2H/δ18O data exist in the reference database used for
calibration, feather values were calibrated against the precipitation
isoscape using a North America-masked Passerine reference group. The δ2H
calibration (n = 699 known-origin sites) explained 75.9% of variance in
tissue isotope values (R² = 0.759, F1,697 = 2195, *P* < 0.001; slope =
0.708 ± 0.015, intercept = −14.16 ± 1.18‰). The δ18O calibration (n = 167
sites) was comparatively weak (R² = 0.080, F1,165 = 14.4, *P* < 0.001;
slope = 0.225 ± 0.059, intercept = 15.40 ± 0.69‰). A sensitivity
calibration using a Water bird reference group (δ2H only, n = 75 sites;
dabbling/diving waterfowl, ecologically closer to a marsh-dwelling rail
than a generic passerine) produced a qualitatively similar fit (R² =
0.549, F1,73 = 88.7, *P* < 0.001; slope = 0.705 ± 0.075, intercept =
−23.07 ± 9.54‰), supporting the general robustness of the δ2H calibration
slope despite the taxonomic mismatch inherent to both reference groups.

### Habitat suitability prior

After reclassification of ESA WorldCover 10 m land cover (85 of 108
candidate tiles contained land within the study area) and Gaussian
smoothing, 101,813 of 1,209,600 grid cells (8.4%) retained nonzero
wetland/wet-cropland suitability across the full North American domain;
within the narrower study area the smoothing kernel's finite radius left
few cells at exactly zero, but suitability values remained strongly
right-skewed (median 0.003; only 32.1% of study-area cells exceeded
0.01), so the prior functioned as a graded weighting rather than a hard
mask. The prior's effect on the posterior was validated against two
locations of known, opposite habitat status: in the population-level
body-feather surface, Mt. Whitney (4,421 m alpine peak; no rail habitat)
received a posterior probability of 8.5 × 10⁻¹⁵, compared with 6.2 × 10⁻⁴
at Sacramento National Wildlife Refuge (a documented Central Valley
wetland complex) — a difference exceeding ten orders of magnitude,
confirming the prior substantially reshapes the assignment rather than
acting as a negligible weighting term. Incorporating the habitat prior
reduced the 50%-credible-region area of the population-level body-feather
surface from 389,188 km² (isotope data alone) to 57,055 km² (habitat-
constrained), and the flight-feather surface from 3,877 km² to 1,906
km², without altering the underlying zone-level geographic pattern
described below (Figure S1).

### Geographic assignment: body versus flight feathers

Population-level posterior surfaces (equal-weighted means across the 12
KLMA individuals, habitat-constrained, dual-isotope model) assigned
substantially different probable origins to body and flight feathers
(Table T1; Figure F1). Body feathers were assigned a majority of
posterior probability (57.7%) to Mexico (Baja California + mainland
western Mexico combined), versus 31.4% to the Pacific Northwest region
(Oregon, Washington, Idaho, Utah, and California combined); the
50%-credible region covered 57,055 km², the 90%-credible region 393,051
km². Flight feathers showed the reverse pattern: 85.8% of posterior
probability in the Pacific Northwest versus 6.1% in Mexico, with a much
smaller 50%-credible region (1,906 km²; 90%-credible region 153,110 km²).
Neither population-level surface assigned appreciable probability (<0.2%)
to a 100-km radius around the Klamath Marsh capture site itself.

At the individual level, 6 of 10 KLMA individuals with usable body-feather
data had a majority (>50%) of their own posterior probability in Mexico
(range 0.2–92.8%, median 82.5%), and 11 of 12 individuals with
flight-feather data had a majority of their posterior probability in the
Pacific Northwest (range 50.0–99.99%, median 90.0%) (Table T1; Figures
F2, F3).

### Sensitivity of the body-feather result to δ18O

`pdRaster()` issues a built-in warning when input isoscapes share more
than 70% of their spatial variance, flagging that joint-likelihood
specificity may be inflated; this warning was triggered here. Direct
calculation confirmed the raw δ2H and δ18O growing-season isoscapes are
correlated at r² = 0.99 across the study area, meaning the two isotope
dimensions are far from the statistically independent evidence that
`pdRaster()`'s joint-normal likelihood assumes. To quantify the practical
consequence, the full assignment was repeated using δ2H alone (identical
calibration approach, habitat prior, and study area) and compared
directly against the dual-isotope model (Table T2).

The flight-feather result was unaffected by this substitution: the
Pacific Northwest majority strengthened slightly under δ2H alone (85.8% →
98.2% population-level; 11 of 12 → 12 of 12 individuals with a Pacific
Northwest majority). The body-feather result was not robust: the
population-level majority zone reversed, from 57.7% Mexico / 31.4% Pacific
Northwest under the dual-isotope model to 34.6% Mexico / 59.6% Pacific
Northwest under δ2H alone. This reversal held at the individual level for
most KLMA birds — e.g., individual 1422-02211's body-feather Mexico
probability fell from 84.4% to 39.5%, and 1422-02217's from 82.6% to
35.1% — while a minority of individuals (1422-02212, 1422-02218,
1422-02221) remained Pacific-Northwest-majority under both models, and
only one (1422-02213) remained robustly Mexico-majority under both (92.8%
→ 89.3%). The non-KLMA reference individual 1272-31742 was a further
exception: its body-feather Mexico probability was essentially unchanged
by the substitution (94.9% → 95.6%), consistent with its isotope values
being extreme enough that either model converges on the same result.
50%-credible-region areas were systematically smaller under the
dual-isotope model than under δ2H alone across nearly all individual
surfaces (median ratio 0.12, i.e., dual-isotope credible areas averaged
roughly one-eighth the size of the corresponding δ2H-only areas; Table
S2), consistent with `pdRaster()`'s own warning that joint specificity may
be overstated under isoscape collinearity of this magnitude.

### Reliability screening of individual-level assignments

Individual-category posterior surfaces were screened for a
match-quality artifact in which spatially flat regions of the isoscape
(e.g., portions of the northern Rocky Mountain interior, where δ2H
standard deviation is roughly a third of its value in California over a
comparable area) combine with the sparseness of the habitat mask to
produce spuriously narrow or coincident credible regions. Of 26
individual-category surfaces, 6 (23.1%) had a single grid cell accounting
for more than 30% of total posterior probability, and 14 (53.8%) shared
an identical peak-probability grid cell with at least one other,
unrelated individual's surface — a signature of isotope values that fall
outside what is achievable anywhere within the habitat-masked domain,
causing the model to default to the same "least-bad" available cell. In
total, 14 of 26 surfaces (53.8%) met at least one of these two criteria,
and 11 (42.3%) had a 50%-credible-region area under 1,000 km². Flagged
surfaces are noted in Table T1 and are interpreted at the population level
only; individual-level credible-area and peak-location statistics for
these surfaces should not be read as precise point estimates.

### Carbon and nitrogen isotopes

Carbon and nitrogen isotope values were obtained from 17 feather samples
across 6 individuals (2 KLMA-banded birds with multiple feather
categories, 2 with body feathers only, and reference individual
1272-31742). δ13C values across all individual-by-category means ranged
from −25.5‰ to −18.4‰, consistent with a C3-plant-based diet throughout
and showing no individual approaching the C4 signature (approximately
−12‰ to −17‰) associated with *Spartina*/*Distichlis* salt-marsh
vegetation. δ15N values ranged from 7.05‰ to 15.35‰ across individual
means; one individual (1422-02220) showed marked within-individual
heterogeneity, with four body-feather subsamples splitting between a
freshwater-typical value (8.1‰) and three elevated values (14.9–15.3‰),
indicating that its sampled body feathers did not grow under a single
isotopic regime. A comparison of δ15N against the dual-isotope
geographic assignment, restricted to the 4 individuals with both
measurements available for body feathers, did not show a consistent
relationship (Table S3): two individuals fit a pattern of higher δ15N
co-occurring with higher assigned Mexico probability, but one
individual (1422-02216: δ15N = 8.3‰, 85.4% Mexico probability) directly
contradicted it; with only 4 data points this is not interpreted as a
reliable trend.

---

## Table T1. Per-individual and population-level summary

*(from `outputs/tables/Table1_manuscript.csv`; reproduced here at full
precision available — percentages rounded to 1 decimal in the source
file)*

| Bird | Type | Feathers | n | Peak (lon, lat) | 50%-CA (km²) | 90%-CA (km²) | % Mexico | % Pacific NW | Note |
|---|---|---|---|---|---|---|---|---|---|
| 1422-02211 | KLMA | Body | 2 | −103.29, 25.79 | 14,124 | 65,648 | 84.4 | 0.0 | |
| 1422-02211 | KLMA | Flight | 2 | −119.96, 39.13 | 929 | 9,033 | 0.0 | 68.3 | |
| 1422-02212 | KLMA | Body | 2 | −121.79, 39.21 | 8,011 | 39,101 | 0.2 | 96.6 | |
| 1422-02212 | KLMA | Flight | 2 | −113.63, 45.13 | 427 | 130,792 | 44.4 | 50.0 | Low match quality |
| 1422-02213 | KLMA | Body | 2 | −105.29, 21.96 | 1,033 | 3,100 | 92.8 | 0.0 | Low match quality |
| 1422-02213 | KLMA | Flight | 2 | −120.88, 45.63 | 2,887 | 17,229 | 0.0 | 90.0 | |
| 1422-02214 | KLMA | Body | 2 | −103.29, 25.71 | 9,727 | 69,664 | 82.4 | 0.0 | Low match quality |
| 1422-02214 | KLMA | Flight | 2 | −119.54, 45.96 | 6,967 | 25,673 | 0.0 | 99.8 | |
| 1422-02215 | KLMA | Body | 2 | −103.38, 24.04 | 22,425 | 105,510 | 91.4 | 0.5 | |
| 1422-02215 | KLMA | Flight | 2 | −113.63, 45.13 | 61 | 183 | 0.0 | 99.2 | Low match quality |
| 1422-02216 | KLMA | Body | 1 | −102.54, 20.29 | 72,254 | 322,700 | 85.4 | 1.1 | |
| 1422-02216 | KLMA | Flight | 2 | −113.63, 45.13 | 61 | 61 | 0.0 | 100.0 | Low match quality |
| 1422-02217 | KLMA | Body | 2 | −103.29, 25.71 | 9,889 | 70,617 | 82.6 | 0.0 | Low match quality |
| 1422-02217 | KLMA | Flight | 2 | −113.63, 45.13 | 437 | 2,638 | 0.0 | 78.4 | Low match quality |
| 1422-02218 | KLMA | Body | 2 | −121.96, 38.13 | 12,451 | 63,155 | 13.4 | 82.8 | |
| 1422-02218 | KLMA | Flight | 2 | −122.96, 44.46 | 5,639 | 28,298 | 0.1 | 90.1 | |
| 1422-02219 | KLMA | Body | 2 | −113.63, 45.13 | 973 | 97,821 | 43.9 | 50.0 | Low match quality; inconsistent feathers |
| 1422-02219 | KLMA | Flight | 2 | −112.63, 44.13 | 5,608 | 44,627 | 0.0 | 89.7 | |
| 1422-02220 | KLMA | Flight | 2 | −113.63, 45.13 | 61 | 183 | 0.0 | 99.7 | Low match quality |
| 1422-02221 | KLMA | Body | 1 | −122.71, 44.96 | 20,852 | 98,761 | 0.4 | 82.6 | |
| 1422-02221 | KLMA | Flight | 2 | −113.63, 45.13 | 61 | 183 | 0.0 | 99.8 | Low match quality |
| 1422-02251 | KLMA | Flight | 3 | −105.29, 21.96 | 28,755 | 183,486 | 29.3 | 64.5 | Low match quality; inconsistent feathers |
| **All KLMA (population mean)** | KLMA | Body | 10 | −113.63, 45.13 | 57,055 | 393,051 | **57.7** | 31.4 | |
| **All KLMA (population mean)** | KLMA | Flight | 12 | −113.63, 45.13 | 1,906 | 153,110 | 6.2 | **85.8** | |
| 1272-31742 | Reference | Body | 4 | −105.29, 21.88 | 715 | 2,623 | 94.9 | 0.0 | |
| 1272-31742 | Reference | Flight | 3 | −113.63, 45.13 | 14,907 | 104,157 | 29.1 | 59.1 | Low match quality; inconsistent feathers |
| San Diego ref. | Reference | Body | 1 | −113.63, 45.13 | 61 | 427 | 0.0 | 99.4 | Low match quality |
| San Diego ref. | Reference | Flight | 2 | −113.63, 45.13 | 61 | 61 | 0.0 | 100.0 | Low match quality |

"% Mexico" = Baja California + mainland western Mexico zone probabilities
combined. "% Pacific NW" = Oregon + Washington/Idaho/Utah + California
zone probabilities combined. 50%-CA / 90%-CA = area of the 50%/90%
credible region. "Low match quality" = `single_cell_dominant_flag` and/or
`shared_peak_flag` (see Methods); "inconsistent feathers" =
`combine_method = inconsistent_feathers_mean_fallback`.

## Table T2. Dual-isotope versus δ2H-only population-level results

| | Dual-isotope (primary model) | δ2H-only (sensitivity check) |
|---|---|---|
| Body: % Mexico | 57.7 | 34.6 |
| Body: % Pacific Northwest | 31.4 | 59.6 |
| Body: 50%-credible area (km²) | 57,055 | 112,407 |
| Flight: % Mexico | 6.1 | 0.8 |
| Flight: % Pacific Northwest | 85.8 | 98.2 |
| Flight: 50%-credible area (km²) | 1,906 | 27,858 |

---

## Notes for Discussion (not drafted here)

- The body-vs-flight contrast (the paper's central hypothesis test) holds
  under both the dual-isotope and δ2H-only models and should be reported
  as the robust finding.
- The specific framing "body feathers show a majority Mexico-wintering
  signal" should be presented as sensitivity-dependent (Table T2), not as
  a confirmed result, pending either an isoscape uncertainty product that
  accounts for δ2H/δ18O covariance, or adopting δ2H-only as the primary
  model with the dual-isotope result reported as a secondary check.
- Reference individual 1272-31742's body-feather result is an exception
  worth highlighting on its own: unlike the broader KLMA sample, it is
  robust to the δ18O sensitivity check, likely because its isotope values
  are far enough from the regional "flat zone" described in Reliability
  screening above.
- No age-class data were available to distinguish feathers grown during
  the current versus a previous molt cycle; this should be flagged as a
  limitation on interpreting the body-vs-flight contrast as evidence of
  this year's movements specifically.
- Carbon/nitrogen results support a freshwater (not salt-marsh) diet
  across all individuals regardless of assigned region, but the small
  sample (n = 4 for the geography cross-check) does not support a
  quantitative diet-geography relationship claim.

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
| Table T2 | population rows of `outputs/tables/sensitivity_d2Honly_vs_dual.csv` |
| Table S2 (supp.) | `outputs/tables/sensitivity_d2Honly_vs_dual.csv` — full per-individual dual-vs-δ2H-only comparison, all 26 surfaces plus the 2 population rows |
| Table S3 (supp.) | `data/cn_vs_geography_comparison.csv` |
