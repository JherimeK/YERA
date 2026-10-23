# Methods (draft, for manuscript — target journal: *Ornithology*)

This is a condensed methods draft, written for a manuscript submission. It
summarizes the same pipeline documented in full in `METHODS_WALKTHROUGH.md`
and `REPORT.md`. Numbers should be spot-checked against the current
`outputs/tables/assignment_summary.csv` and
`outputs/tables/sensitivity_comparison.csv` before submission. Citations
are given in full below since they aren't yet in a reference manager;
adjust formatting to the journal's style on submission.

**Revision note:** this draft reflects a substantial revision made after
an external methodological review. The dual-isotope (δ2H+δ18O) model used
as primary in an earlier draft is now reported only as a sensitivity
variant; δ2H alone is the primary model. Two individuals previously
described as "reference" individuals have been removed from the analysis
entirely, since neither had an independently known breeding or molt origin
for their assignment to be checked against.

---

### Sample collection and stable isotope analysis

Feather samples were collected from 12 Yellow Rails (*Coturnicops
noveboracensis*) captured at Klamath Marsh National Wildlife Refuge,
Oregon, in July [YEAR — fill in]. Body feathers (grown during the
preformative/definitive basic molt, which can continue into the wintering
period but does not necessarily do so) and flight feathers (grown during
the complete prebasic molt, primarily on or near the breeding grounds)
were sampled separately from each bird to test for a difference in
molt-origin location between these two feather generations [insert your
capture/handling/IACUC/permit details]. Hydrogen and oxygen stable isotope
ratios (δ2H, δ18O) were measured at the University of Utah SIRFER
facility on 43 feathers [insert instrumentation/method details from the
lab report if required by the journal]. Carbon and nitrogen stable
isotope ratios (δ13C, δ15N) were measured separately at [SIBS lab — fill
in full name/location] on 14 feathers from 5 of the 12 individuals.
Samples with combustion H:O ratios or %H/%O yields outside the expected
range for keratin were flagged during QA/QC review (3 of 43 feathers);
results are reported with these individuals noted rather than excluded.

Two additional individuals sampled in the same laboratory batches (a bird
collected near San Diego, California, and an individual with a
non-matching band-number prefix) were excluded from this analysis. Both
were initially treated as reference individuals in early analyses, but
neither has an independently known breeding or molt origin, so there was
no basis against which to validate their isotope-based assignments; they
are not part of the Klamath Marsh capture series and their capture dates
relative to the study individuals are unknown. [If capture/collection
details for either individual become available, consider whether
reintroducing them as genuine comparison points is warranted.]

### Geographic assignment

Geographic origin was estimated using the Bayesian isotope-assignment
framework implemented in the R package **assignR** (Ma et al. 2020),
applied to a growing-season precipitation δ2H isoscape. Because no
Rallidae, and essentially no North American non-passerine, known-origin
isotope data exist in assignR's reference database, feather keratin
values were calibrated against precipitation isoscape values using the
package's built-in North America-masked Passerine reference group (n =
699 known-origin sites) as the primary model's calibration (R² = 0.76).
Isoscape prediction uncertainty was not available for the region at the
time of analysis and was approximated as spatially uniform (5‰); we
recommend this be revisited with region-specific prediction standard
errors before final publication (see Discussion/Limitations).

δ18O was measured for all samples and an initial version of this analysis
used both δ2H and δ18O jointly, following common practice for this
assignment framework. On review, this was changed: the δ18O calibration
against the Passerine reference group was weak (R² = 0.08), and the raw
δ2H and δ18O growing-season isoscapes are correlated at r² = 0.99 across
this study area (both trace the same precipitation fractionation
process). `assignR::pdRaster()` issues a built-in warning under exactly
this condition (isoscapes sharing more than 70% of their spatial
variance), flagging that joint-likelihood specificity may be inflated
when the two markers are this collinear. We therefore adopted δ2H alone
as the primary model and report the dual-isotope model as one of several
named sensitivity analyses (below), rather than treating it as the
default.

Because Yellow Rails are a wetland and wet-agriculture habitat specialist
(Bookhout 1995; [Birds of the World account citation]) and unconstrained
isotope-only assignment produced biologically implausible high-probability
regions in montane and desert terrain, a habitat suitability layer was
incorporated as a multiplicative prior in the assignment following the
approach of Vander Zanden et al. (2014), implemented via assignR's `prior`
argument to `pdRaster()`. The suitability layer was derived from ESA
WorldCover 10 m land cover (v200, 2021; Zanaga et al. 2022), reclassified
to weight herbaceous wetland highest, cropland at partial weight
(reflecting this species' use of wet agricultural habitat, which
WorldCover's cropland class cannot further resolve into flooded/dry),
permanent open water at low weight, and all other land cover at zero,
then aggregated to the isoscape's native resolution (~9 km) and lightly
spatially smoothed to avoid discretization artifacts at wetland-patch
edges.

Feathers of the same category (body or flight) sampled from the same
individual were combined by multiplying and renormalizing their individual
posterior probability surfaces, under the assumption of a shared molt
origin; where two feathers from the same individual and category produced
no numerically overlapping posterior support, the simple mean of their
posteriors was used instead and the individual was flagged for that
category. Under the primary (δ2H-only) model, this did not occur for any
of the 20 multi-feather combinations (all combined successfully via the
joint-product method); under the dual-isotope sensitivity model, by
comparison, 4 combinations failed this check, since a single, less
spatially discriminating marker produces broader per-feather posteriors
that are more likely to numerically overlap regardless of whether the
underlying feathers are truly consistent. Population-level surfaces were
computed as the equal-weighted mean of
individual posteriors within each feather category across all 12
individuals. Fifty, 75, and 90% credible regions were computed for every
individual and population-level surface as the smallest-area region
containing that fraction of cumulative posterior probability.

### Sensitivity analysis

Because the primary model's calibration and habitat-prior specification
both involve subjective choices, the full assignment was rerun under five
alternative configurations and compared against the primary model,
individual-by-individual and at the population level
(`scripts/06_sensitivity_check.R`): (1) the dual-isotope (δ2H+δ18O) model
described above; (2) no habitat prior; (3) a "weak" habitat prior (the
primary prior raised to the 0.5 power, compressing it toward a uniform,
non-informative weighting); (4) a "strong" habitat prior (the primary
prior squared, exaggerating it toward a near-binary mask); and (5) a
calibration using the "Water bird" reference group (dabbling/diving
waterfowl, ecologically closer to a marsh-dwelling rail than a generic
passerine; R² = 0.55) in place of the Passerine group, with δ2H only. The
flight-feather result was stable across all five variants. The
body-feather result was not: it depended on which isotope marker(s) and
calibration reference group were used (ranging from 27.9% to 53.1%
population-level probability assigned to Mexico across the six total
configurations including the primary model), while habitat-prior strength
specifically had comparatively little effect (33.5-34.7% across a
four-fold range in prior strength). This is reported as a result in its
own right (see Results) rather than used to select a single "best"
configuration.

As a check on whether the habitat prior itself could be mechanically
biasing results toward one region (e.g., if ESA WorldCover underrepresents
ephemeral or seasonally flooded wetlands in Mexico relative to the United
States), the prior's own distribution of suitability mass and area between
the Mexico and Pacific Northwest zones was tabulated independent of any
isotope data (`outputs/tables/habitat_prior_by_zone.csv`).

Individual-level results were additionally screened for a match-quality
artifact arising from the interaction of spatially flat regions of the
isoscape with the sparseness of the habitat mask, in which the assignment
can default to a single "best available" pixel independent of a
genuinely close isotopic match, producing spuriously large or coincident
credible regions shared across unrelated individuals. Under the primary
(δ2H-only) model this affected the majority of individual-level surfaces
(17 of 22, 77%) — a larger fraction than under the dual-isotope model
previously evaluated (14 of 26, 54%), since δ2H alone is less spatially
discriminating. Flagged surfaces are interpreted at the population level
only.

### Carbon and nitrogen isotopes as an independent habitat cross-check

Because no fine-resolution continental δ13C/δ15N isoscape suitable for
spatial assignment exists, carbon and nitrogen isotope values were not
incorporated into the spatial assignment model directly. δ13C values were
corrected for a systematic +2.05‰ enrichment bias identified from the
analytical lab's in-house QA plant standards (two standards, each run in
triplicate, both reading δ13C enriched relative to their known consensus
values by a consistent amount — +2.10‰ and +2.01‰ respectively — despite
spanning a ~27‰ range, supporting a simple additive correction); δ15N
showed no comparable bias. δ13C was used to screen for a C4-plant (e.g.,
*Spartina*/*Distichlis* salt marsh) dietary signature, and δ15N was
interpreted qualitatively as an indicator of trophic position and
marine/agricultural nitrogen enrichment, compared against the independent
δ2H-based geographic assignment for the same individual and feather
category. Only 3-4 individuals had both measurements available for a
given feather category, limiting this comparison to a qualitative check
rather than a quantitative trend test.

### Software

All analyses were conducted in R [version — fill in] using the packages
**assignR** (Ma et al. 2020), **terra**, and **sf** for spatial data
handling, and **dplyr** for data management.

---

## References (verify formatting/completeness against journal style)

- Ma, C., Vander Zanden, H. B., Wunder, M. B., & Bowen, G. J. (2020). assignR:
  An R package for isotope-based geographic assignment. *Methods in Ecology
  and Evolution*, 11(9), 996–1001. https://doi.org/10.1111/2041-210X.13426
- Wunder, M. B. (2010). Using isoscapes to model probability surfaces for
  determining geographic origins. In *Isoscapes: Understanding movement,
  pattern, and process on Earth through isotope mapping* (pp. 251–270).
  Springer. ISBN 9789048133536.
- Vander Zanden, H. B., Wunder, M. B., Hobson, K. A., Van Wilgenburg, S. L.,
  Wassenaar, L. I., Welker, J. M., & Bowen, G. J. (2014). Contrasting
  assignment of migratory organisms to geographic origins using long-term
  versus year-specific precipitation isotope maps. *Methods in Ecology and
  Evolution*, 5(9), 891–900. https://doi.org/10.1111/2041-210X.12229
- Zanaga, D., Van De Kerchove, R., De Keersmaecker, W., Souverijns, N.,
  Brockmann, C., Quast, R., Wevers, J., Grosu, A., Paccini, A., Vergnaud, S.,
  Cartus, O., Santoro, M., Fritz, S., Georgieva, I., Lesiv, M., Carter, S.,
  Herold, M., Li, Linlin, Tsendbazar, N. E., Ramoino, F., & Arino, O. (2022).
  ESA WorldCover 10 m 2021 v200. [Dataset DOI to verify before submission —
  the v100/2020 product's DOI is 10.5281/zenodo.5571936; confirm the correct
  v200/2021 DOI directly from esa-worldcover.org before citing.]
- [Cornell Lab of Ornithology / Birds of the World Yellow Rail account —
  insert full citation per the account's requested format.]
- Bookhout, T. A. (1995). Yellow Rail (*Coturnicops noveboracensis*), version
  2.0. In *The Birds of North America* (A. F. Poole & F. B. Gill, Eds.).
  Cornell Lab of Ornithology. [Confirm current Birds of the World DOI/version.]

---

**Note on what's still placeholder:** capture dates/permits, instrumentation
details for both labs, R package version numbers, and the exact WorldCover
DOI all need to be filled in/verified before submission — these are flagged
in brackets above rather than guessed. Everything else describes methods
actually implemented in the current pipeline.

**Note for Results/Discussion (not drafted in full here, but important for
whoever writes them):** the sensitivity analysis is not a minor technical
aside — it determines which population-level zone has the majority of
posterior probability for body feathers (Pacific Northwest under the
primary model; Mexico under two of the five sensitivity variants), while
leaving the flight-feather result (Pacific Northwest, stable across all
variants) unchanged. We recommend the Results section report the full
range across variants for body feathers rather than only the primary
model's figure, and the Discussion/Limitations explicitly state that the
body-feather geographic origin is unresolved with the present data, while
the body-vs-flight contrast itself (the paper's core hypothesis test) is
robust across every variant tested. See `REPORT.md` and
`MANUSCRIPT_RESULTS.md` for the full numeric breakdown.
