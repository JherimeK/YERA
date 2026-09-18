# Methods (draft, for manuscript — target journal: *Ornithology*)

This is a condensed methods draft, written for a manuscript submission. It
summarizes the same pipeline documented in full in `METHODS_WALKTHROUGH.md`
and `REPORT.md`. Numbers should be spot-checked against the current
`outputs/tables/assignment_summary.csv` before submission, since the pipeline
may still be refined. Citations are given in full below since they aren't
yet in a reference manager; adjust formatting to the journal's style on
submission.

---

### Sample collection and stable isotope analysis

Feather samples were collected from Yellow Rails (*Coturnicops
noveboracensis*) captured at Klamath Marsh National Wildlife Refuge, Oregon,
in July [YEAR — fill in]. Body feathers (assumed to be replaced during a molt
that can complete on the wintering grounds) and flight feathers (replaced
during the complete prebasic molt, primarily on or near the breeding
grounds) were sampled separately from each bird to test for a difference in
molt-origin location between these two feather generations [insert your
capture/handling/IACUC/permit details]. Hydrogen and oxygen stable isotope
ratios (δ2H, δ18O) were measured at the University of Utah SIRFER facility
[insert instrumentation/method details from the lab report if required by
the journal]. Carbon and nitrogen stable isotope ratios (δ13C, δ15N) were
measured separately at [SIBS lab — fill in full name/location] on a subset
of individuals. Samples with combustion H:O ratios or %H/%O yields outside
the expected range for keratin were flagged during QA/QC review; results are
reported with and without these individuals noted (see Results).

### Geographic assignment

Geographic origin was estimated using the dual-isotope (δ2H, δ18O) Bayesian
assignment framework implemented in the R package **assignR** (Ma et al.
2020), applied to a growing-season precipitation isoscape. Because no
Rallidae, and essentially no North American non-passerine, known-origin
δ18O keratin data exist in assignR's reference database, feather keratin
values were calibrated against precipitation isoscape values using the
package's built-in North America-masked Passerine reference group (n = 699
known-origin sites for δ2H) as the best available proxy calibration (δ2H:
R² = 0.76; δ18O: R² = 0.08). The weak δ18O calibration means δ2H contributed
the majority of the geographic discriminating power in the resulting
assignments. A sensitivity calibration using a Water bird reference group
(δ2H only) produced a qualitatively similar result (R² = 0.55). Isoscape
prediction uncertainty was not available for the region at the time of
analysis and was approximated as spatially uniform (5‰ for δ2H, 0.6‰ for
δ18O); we recommend this be revisited with region-specific prediction
standard errors before final publication (see Discussion/Limitations).

The assignment domain covered California, Nevada, Arizona, New Mexico,
Oregon, Washington, Idaho, and Utah, and extended south through Baja
California and mainland western Mexico to approximately 14° N, encompassing
both the capture location and the broader western North American and Mexican
range plausible for this species' migratory connectivity.

Because Yellow Rails are a wetland and wet-agriculture habitat specialist
(Bookhout 1995; [Birds of the World account citation]) and unconstrained
isotope-only assignment produced biologically implausible high-probability
regions in montane and desert terrain, a habitat suitability layer was
incorporated as a multiplicative prior in the assignment following the
approach of Vander Zanden et al. (2014), implemented via assignR's `prior`
argument to `pdRaster()`. The suitability layer was derived from ESA
WorldCover 10 m land cover (v200, 2021; Zanaga et al. 2022), reclassified to
weight herbaceous wetland highest, cropland at partial weight (reflecting
this species' use of wet agricultural habitat, which WorldCover's cropland
class cannot further resolve into flooded/dry), permanent open water at low
weight, and all other land cover at zero, then aggregated to the isoscape's
native resolution (~9 km) and lightly spatially smoothed to avoid
discretization artifacts at wetland-patch edges. Incorporating this prior
substantially reduced the spatial extent of the resulting credible regions
(50%-credible-region area for the population-level body-feather surface:
389,000 km² without the habitat prior vs. 57,000 km² with it) without
changing the underlying isotope-driven geographic pattern.

Feathers of the same category (body or flight) sampled from the same
individual were combined by multiplying and renormalizing their individual
posterior probability surfaces, under the assumption of a shared molt
origin; where two feathers from the same individual and category produced
no numerically overlapping posterior support (i.e., were isotopically
inconsistent with one another), the simple mean of their posteriors was used
instead and the individual was flagged for that category. Population-level
surfaces were computed as the equal-weighted mean of individual posteriors
within each feather category. Fifty, 75, and 90% credible regions were
computed for every individual and population-level surface as the
smallest-area region containing that fraction of cumulative posterior
probability (assignR's `qtlRaster()`, probability-based threshold).

Individual-level results were additionally screened for a match-quality
artifact arising from the interaction of spatially flat regions of the
isoscape (e.g., portions of the northern Rocky Mountain interior have markedly
lower δ2H spatial variance than California) with the sparseness of the
wetland habitat mask: where a bird's measured isotope value fell outside the
range achievable within any wetland-masked pixel, the assignment could
default to a single "best available" pixel independent of a genuinely close
isotopic match, occasionally producing spuriously narrow, and occasionally
coincident, individual credible regions. Individuals for which a single grid
cell accounted for a disproportionate share of the total posterior
probability, or whose peak-probability cell coincided exactly with another
individual's, were flagged and are interpreted at the population level only
(see Results/Discussion); we report this transparently as a limitation of
applying a fine-scale habitat mask to a coarse-resolution isoscape with
uniform assumed uncertainty, rather than treating narrow individual credible
regions as precise results in all cases.

### Carbon and nitrogen isotopes as an independent habitat cross-check

Because no fine-resolution continental δ13C/δ15N isoscape suitable for
spatial assignment exists, carbon and nitrogen isotope values were not
incorporated into the spatial assignment model directly. Instead, δ13C was
used to screen for a C4-plant (e.g., *Spartina*/*Distichlis* salt marsh)
dietary signature, and δ15N was interpreted qualitatively as an indicator of
trophic position and marine/agricultural nitrogen enrichment, each compared
against the independent δ2H/δ18O-based geographic assignment for the same
individual and feather category as a cross-validation of habitat and origin
inference.

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
  insert full citation per the account's requested format; BOW citations are
  typically author, year, species account title, in *Birds of the World*
  (Cornell Lab of Ornithology, Ithaca, NY), with a DOI.]
- Bookhout, T. A. (1995). Yellow Rail (*Coturnicops noveboracensis*), version
  2.0. In *The Birds of North America* (A. F. Poole & F. B. Gill, Eds.).
  Cornell Lab of Ornithology. [Confirm current Birds of the World DOI/version
  — this account has likely been revised since the original BNA number.]

---

**Note on what's still placeholder:** capture dates/permits, instrumentation
details for both labs, R package version numbers, and the exact WorldCover
DOI all need to be filled in/verified before submission — these are flagged
in brackets above rather than guessed. Everything else describes methods
actually implemented in the current pipeline.
