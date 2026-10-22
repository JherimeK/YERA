# 01_run_assignment.R
# assignR-based geographic origin assignment for KLMA Yellow Rail feathers.
#
# PRIMARY MODEL (this script's main output, outputs/rasters/pd_per_feather.tif):
# d2H ALONE, Passerine calibration, current (Gaussian-smoothed) habitat
# prior, KLMA individuals only (12 birds, 43 feathers).
#
# This changed from an earlier version of this analysis, which used the
# dual-isotope (d2H+d18O) model as primary. On review (prompted by an
# external code/results review and confirmed with our own sensitivity
# check -- see scripts/06_sensitivity_check.R and REPORT.md), that was the
# wrong default for two independent reasons:
#   1. d18O calibration against this North America-masked Passerine
#      reference group is weak (R^2 = 0.08), and the raw d2H/d18O growing-
#      season isoscapes are r^2 = 0.99 correlated across this study area --
#      assignR::pdRaster() itself warns about exactly this condition
#      ("shared variance > 0.7 ... specificity of assignments may be
#      inflated"). Treating d18O as a second independent evidence
#      dimension when it is this collinear with d2H can narrow posteriors
#      (and shift zone-level probabilities) in a way that isn't
#      statistically justified.
#   2. Dropping d18O and reassigning with d2H alone reversed the
#      population-level body-feather result (Mexico-majority ->
#      Pacific-Northwest-majority) while leaving the flight-feather result
#      intact. A result that reverses when a weakly-calibrated, highly
#      collinear second marker is removed should not be the default/primary
#      analysis -- it belongs in a sensitivity check, not the headline
#      model. The dual-isotope model is still run and reported (as a
#      sensitivity variant, in scripts/06_sensitivity_check.R), not
#      dropped.
#
# Calibration: assignR ships no Rallidae (or any non-passerine) known-origin
# keratin data for North America, so a species-specific calibration is not
# possible. The North-America-masked "Passerine" group (n=699 d2H sites) is
# used as the best-available proxy for the primary model. This is a real
# limitation (Rallidae keratin-water fractionation may differ from
# passerines) and is reported as such. A "Water bird" d2H-only calibration
# (dabbling/diving waterfowl, ecologically closer to a marsh rail) is run
# as a further sensitivity check, also in script 06.
#
# Isoscape uncertainty: the supplied d2h_GS.tif / d18o_GS.tif are mean-only
# growing-season precipitation isoscapes (no SE band could be downloaded in
# this environment -- wateriso.utah.edu is not reachable). A spatially
# uniform isoscape SE is assumed (5 permil H, 0.6 permil O), based on the
# magnitude of assignR's own example prediction-SE raster. This is a real
# approximation and should be replaced with the official SE grid if it
# becomes available.
#
# Reference individuals: the two non-KLMA individuals formerly carried
# through this pipeline ("San Diego ref bird", "1272-31742") have been
# removed entirely as of this version -- see the note in
# scripts/00_clean_data.R. They do not appear in data/klma_clean.csv, so
# nothing further is needed here to exclude them.
#
# Habitat prior sensitivity: this script also builds, and saves for script
# 06 to use, "weak" and "strong" variants of the habitat suitability prior
# (prior^0.5 and prior^2, respectively -- a one-parameter family that
# preserves the underlying ESA WorldCover-derived suitability pattern while
# varying how strongly it discounts non-habitat terrain), so the effect of
# the prior's strength on the results can be checked directly rather than
# asserted.

suppressMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(assignR)
  library(maps)
})

set.seed(42)
sf_use_s2(FALSE)

dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/maps", showWarnings = FALSE)
dir.create("outputs/rasters", showWarnings = FALSE)

ASSUMED_SD_D2H <- 5.0

# ---- 1. Study area (AOI) ----------------------------------------------------
# NOTE: the user's original frame of reference (CA/NV/AZ/NM south through
# Baja and western Mexico) omits Oregon -- but Klamath Marsh NWR (the
# capture/breeding site) is itself in Oregon at ~43 N, just north of the
# California line. Cornell's molt account says flight feathers are grown
# "primarily on or near breeding grounds," so excluding Oregon would a
# priori rule out the single most biologically expected flight-feather
# origin before the isotope data get a say. Oregon (plus Washington/Idaho/
# Utah for geographic continuity of the western breeding range) is added
# to the AOI so that possibility remains on the table.
us <- st_as_sf(maps::map("state", regions = c("california", "nevada", "arizona", "new mexico",
                                               "oregon", "washington", "idaho", "utah"),
                          fill = TRUE, plot = FALSE))
mx <- st_as_sf(maps::map("world", "Mexico", fill = TRUE, plot = FALSE))
st_crs(us) <- 4326; st_crs(mx) <- 4326
us <- st_make_valid(us); mx <- st_make_valid(mx)
bbox_mx <- st_bbox(c(xmin = -118, xmax = -101, ymin = 14, ymax = 33), crs = 4326)
mx_clip <- suppressWarnings(st_crop(mx, bbox_mx))

aoi_sf <- st_as_sf(st_union(st_union(st_geometry(us)), st_union(st_geometry(mx_clip))))
aoi <- vect(aoi_sf)
crs(aoi) <- "EPSG:4326"

# split Mexico clip into Baja peninsula vs mainland west-coast parts (by
# centroid longitude, since the Gulf of California naturally separates the
# two polygon pieces in the source land polygon)
mx_parts <- st_cast(st_make_valid(st_union(mx_clip)), "POLYGON")
mx_cent <- st_coordinates(st_centroid(mx_parts))
baja_sf <- st_as_sf(st_union(mx_parts[mx_cent[, 1] < -111.5]))
mainland_mx_sf <- st_as_sf(st_union(mx_parts[mx_cent[, 1] >= -111.5]))

zones <- list(
  California   = vect(st_as_sf(us[us$ID == "california", ])),
  Nevada       = vect(st_as_sf(us[us$ID == "nevada", ])),
  Arizona      = vect(st_as_sf(us[us$ID == "arizona", ])),
  New_Mexico   = vect(st_as_sf(us[us$ID == "new mexico", ])),
  Baja_CA      = vect(baja_sf),
  Mainland_W_Mexico = vect(mainland_mx_sf)
)
for (z in names(zones)) crs(zones[[z]]) <- "EPSG:4326"

focal_points <- data.frame(
  site = c("Klamath_Marsh", "SF_Bay", "San_Diego"),
  lon = c(-121.6008333, -122.0262222, -117.2288611),
  lat = c(42.9766389, 37.4533611, 32.7931111)
)
klamath_buf <- st_as_sf(st_buffer(st_sfc(st_point(c(-121.6008333, 42.9766389)), crs = 4326),
                                   dist = 100000))
klamath_buf_v <- vect(klamath_buf)

# ---- 2. Isoscape (d2H only -- the primary model's only marker) ------------
# d18O and the Water-bird calibration are NOT loaded here any more: the
# primary model doesn't use them, and scripts/06_sensitivity_check.R
# recomputes everything it needs from scratch anyway (SpatRaster objects
# don't survive saveRDS()/readRDS() across the separate R processes each
# script runs as -- see that script's header comment), so building them
# here was pure wasted computation. Their calibration diagnostic plots
# (calibration_regression_d18O.png, _d2H_waterbird_sensitivity.png) are
# produced by script 06 instead, right where those models are actually
# used.
na_ext <- ext(-170, -50, 5, 75)
d2h_mean <- crop(rast("d2h_GS.tif"), na_ext)
d2h_iso <- c(d2h_mean, setValues(d2h_mean, ASSUMED_SD_D2H))
names(d2h_iso) <- c("mean", "sd")

# ---- 3. Calibration (tissue <- precipitation), Passerine group, d2H -------
data(naMap)
known_H <- subOrigData(marker = "d2H", group = "Passerine", mask = naMap,
                        niter = 2000, genplot = FALSE)
png("outputs/maps/calibration_regression_d2H.png", width = 900, height = 700)
cal_H <- calRaster(known_H, d2h_iso, genplot = TRUE, verboseLM = TRUE)
dev.off()

# ---- 3b. Habitat prior (wetland / wet-cropland suitability) ----------------
# Built by scripts/05_build_habitat_prior.R from ESA WorldCover 10m land
# cover (herbaceous wetland = 1.0, cropland = 0.3, permanent water = 0.1,
# everything else = 0), on the exact same grid as the isoscape above. This
# is what keeps the assignment out of mountain/forest/desert terrain that
# happens to match the target d2H/d18O values by elevation alone but isn't
# habitat Yellow Rails would ever occupy (wet sedge meadow / marsh / wet
# agriculture only, per Cornell's Birds of the World account).
#
# assignR::pdRaster.isoStack() has a real bug when mask and prior are both
# supplied: it crops the isoscape to mask internally, but never crops prior
# to match, so the elementwise assign*prior multiplication silently
# recycles two differently-sized vectors (confirmed: every single sample's
# posterior came out exactly zero when tested this way). The fix is to
# pre-crop/mask BOTH the isoscape and the prior to the AOI ourselves and
# pass mask = NULL to pdRaster, so nothing gets cropped a second time
# inside the function.
aoi_crop_mask <- function(r) {
  r <- crop(r, aoi)
  mask(r, aoi)
}
cal_H$isoscape.rescale <- aoi_crop_mask(cal_H$isoscape.rescale)

habitat_prior <- rast("data/habitat_prior.tif")
habitat_prior <- resample(aoi_crop_mask(habitat_prior), cal_H$isoscape.rescale, method = "near")
crs(habitat_prior) <- crs(cal_H$isoscape.rescale)

# "Weak"/"strong" prior-strength variants for script 06's habitat-prior
# sensitivity check: a power transform preserves the underlying suitability
# pattern (same relative ranking of cells) while compressing (weak,
# exponent < 1, closer to a uniform/no-op prior) or exaggerating (strong,
# exponent > 1, closer to a hard mask) how much it discounts non-habitat.
habitat_prior_weak <- habitat_prior^0.5
habitat_prior_strong <- habitat_prior^2

# ---- 4. Load cleaned sample data --------------------------------------------
samp <- read.csv("data/klma_clean.csv", stringsAsFactors = FALSE)
samp$sample_id <- paste(samp$Envelope_ID, samp$Feather_type_clean, samp$Original_ID, sep = "__")

unknown_H <- data.frame(ID = samp$sample_id, d2H = samp$d2H)

message("Running primary per-feather assignment (d2H only, Passerine calibration, ",
        "current habitat prior) for ", nrow(unknown_H), " samples...")
# genplot = FALSE: pdRaster's own internal plot is a raw, uncurated
# byproduct (terra caps it at 16 panels, with assignR's default color
# scheme unrelated to this project's maps) -- not a real deliverable.
# scripts/03_make_maps.R builds the actual per-bird figures.
pd_feather <- pdRaster(cal_H, unknown = unknown_H, mask = NULL, prior = habitat_prior, genplot = FALSE)
names(pd_feather) <- unknown_H$ID
writeRaster(pd_feather, "outputs/rasters/pd_per_feather.tif", overwrite = TRUE)

# Saved for same-session inspection only (e.g. interactively in RStudio).
# terra::SpatRaster objects wrap external C++ pointers that do not survive
# saveRDS()/readRDS() across separate R processes -- they come back with
# dead pointers and error on first use -- so scripts/06_sensitivity_check.R
# recomputes these objects itself rather than reading this file (it runs as
# its own `Rscript` process, like every script in this pipeline).
saveRDS(list(cal_H = cal_H,
             habitat_prior = habitat_prior,
             habitat_prior_weak = habitat_prior_weak,
             habitat_prior_strong = habitat_prior_strong,
             aoi = aoi_sf, zones_names = names(zones)),
        "outputs/calibration_objects.rds")

message("Step 1-4 complete. Saved primary per-feather posterior stack (d2H-only).")
