# 01_run_assignment.R
# assignR-based geographic origin assignment for KLMA Yellow Rail feathers.
#
# Calibration: assignR ships no Rallidae (or any non-passerine) known-origin
# d18O keratin data for North America, so a species-specific calibration is
# not possible. We use the built-in North-America-masked "Passerine" group
# (n=699 d2H sites, n=699/187 usable d18O sites) as the best-available proxy
# calibration -- it is the only NA group in the database with enough paired
# d2H/d18O known-origin sites to fit a robust regression. This is a real
# limitation (Rallidae keratin-water fractionation may differ from
# passerines) and is reported as such. A "Water bird" d2H-only calibration
# (dabbling/diving waterfowl, ecologically closer to a marsh rail) is run
# as a sensitivity check.
#
# Isoscape uncertainty: the supplied d2h_GS.tif / d18o_GS.tif are mean-only
# growing-season precipitation isoscapes (no SE band could be downloaded in
# this environment -- wateriso.utah.edu is not reachable). A spatially
# uniform isoscape SE is assumed (5 permil H, 0.6 permil O), based on the
# magnitude of assignR's own example prediction-SE raster. This SE term is
# minor relative to the tissue-calibration regression variance that
# calRaster adds in quadrature, so its effect on final maps is small, but
# it should be replaced with the official SE grid if the user obtains it.

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
ASSUMED_SD_D18O <- 0.6

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

# ---- 2. Isoscapes -----------------------------------------------------------
na_ext <- ext(-170, -50, 5, 75)
d2h_mean <- crop(rast("d2h_GS.tif"), na_ext)
d18o_mean <- crop(rast("d18o_GS.tif"), na_ext)
d18o_mean <- resample(d18o_mean, d2h_mean, method = "bilinear")

d2h_iso <- c(d2h_mean, setValues(d2h_mean, ASSUMED_SD_D2H))
names(d2h_iso) <- c("mean", "sd")
d18o_iso <- c(d18o_mean, setValues(d18o_mean, ASSUMED_SD_D18O))
names(d18o_iso) <- c("mean", "sd")

# ---- 3. Calibration (tissue <- precipitation) -------------------------------
data(naMap)

known_H <- subOrigData(marker = "d2H", group = "Passerine", mask = naMap,
                        niter = 2000, genplot = FALSE)
known_O <- subOrigData(marker = "d18O", group = "Passerine", ref_scale = "VSMOW_O",
                        mask = naMap, niter = 2000, genplot = FALSE)

png("outputs/maps/calibration_regression_d2H.png", width = 900, height = 700)
cal_H <- calRaster(known_H, d2h_iso, genplot = TRUE, verboseLM = TRUE)
dev.off()

png("outputs/maps/calibration_regression_d18O.png", width = 900, height = 700)
cal_O <- calRaster(known_O, d18o_iso, genplot = TRUE, verboseLM = TRUE)
dev.off()

# sensitivity calibration: Water bird group, d2H only
known_H_wb <- subOrigData(marker = "d2H", group = "Water bird", mask = naMap,
                           niter = 2000, genplot = FALSE)
png("outputs/maps/calibration_regression_d2H_waterbird_sensitivity.png", width = 900, height = 700)
cal_H_wb <- calRaster(known_H_wb, d2h_iso, genplot = TRUE, verboseLM = TRUE)
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
cal_O$isoscape.rescale <- aoi_crop_mask(cal_O$isoscape.rescale)
iso_stack <- isoStack(cal_H, cal_O)

habitat_prior <- rast("data/habitat_prior.tif")
habitat_prior <- resample(aoi_crop_mask(habitat_prior), iso_stack[[1]], method = "near")
crs(habitat_prior) <- crs(iso_stack[[1]])

# ---- 4. Load cleaned sample data --------------------------------------------
samp <- read.csv("data/klma_clean.csv", stringsAsFactors = FALSE)
samp$sample_id <- paste(samp$Envelope_ID, samp$Feather_type_clean, samp$Original_ID, sep = "__")

unknown_df <- data.frame(ID = samp$sample_id, d2H = samp$d2H, d18O = samp$d18O)

message("Running per-feather dual-isotope assignment for ", nrow(unknown_df), " samples...")
png("outputs/maps/per_feather_grid.png", width = 1400, height = 1400)
pd_feather <- pdRaster(iso_stack, unknown = unknown_df, mask = NULL, prior = habitat_prior, genplot = TRUE)
dev.off()
names(pd_feather) <- unknown_df$ID
writeRaster(pd_feather, "outputs/rasters/pd_per_feather.tif", overwrite = TRUE)

saveRDS(list(cal_H = cal_H, cal_O = cal_O, cal_H_wb = cal_H_wb, aoi = aoi_sf,
             zones_names = names(zones)),
        "outputs/calibration_objects.rds")

message("Step 1-4 complete. Saved per-feather posterior stack.")
