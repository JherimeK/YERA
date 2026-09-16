# 05_build_habitat_prior.R
# Build a wetland / wet-cropland habitat suitability raster from ESA
# WorldCover 10m land cover, matching the isoscape grid, for use as a prior
# in assignR::pdRaster(). This is what excludes the biologically implausible
# mountain/desert hotspots from the H/O-only maps: those pixels can match
# the target d2H/d18O values by elevation alone, but Yellow Rails don't
# occupy that terrain (wet sedge meadow / marsh / wet agriculture only, per
# Cornell's Birds of the World account).
#
# ESA WorldCover classes used (others get weight 0):
#   90 Herbaceous wetland -> weight 1.0 (primary habitat)
#   40 Cropland           -> weight 0.3 (Cornell notes "wet agricultural
#                             areas"; WorldCover doesn't distinguish
#                             flooded/irrigated cropland from dry cropland,
#                             so this is a partial-credit weight, not a
#                             confirmed-habitat weight)
#   80 Permanent water     -> weight 0.1 (marsh edge effect only; rails
#                             don't use open water, but wetlands are often
#                             fringed by/misclassified as water)
#
# Each 10m tile (3.24 billion pixels' worth across the full study area) is
# far too slow/large to reclassify directly in R, so each tile is first
# decimated to 100m with a fast nearest-neighbor gdalwarp (a few hundred ms
# per tile, since it's simple block-strided reading, not a full-resolution
# computation), THEN reclassified and averaged down to the isoscape's
# ~0.083 degree grid in R. Raw tiles and the 100m intermediate are deleted
# immediately after each tile is processed to keep disk usage low.

suppressMessages({
  library(terra)
})

dir.create("data/worldcover_tmp", showWarnings = FALSE)

# ---- 1. isoscape template grid (same grid pdRaster's prior must match) -----
na_ext <- ext(-170, -50, 5, 75)
template <- crop(rast("d2h_GS.tif"), na_ext)

# ---- 2. enumerate candidate tiles covering the study area -------------------
lon_starts <- seq(-126, -102, by = 3)
lat_starts <- seq(15, 48, by = 3)
tiles <- expand.grid(lon = lon_starts, lat = lat_starts)

tile_name <- function(lon, lat) {
  ns <- if (lat >= 0) sprintf("N%02d", lat) else sprintf("S%02d", -lat)
  ew <- if (lon >= 0) sprintf("E%03d", lon) else sprintf("W%03d", -lon)
  paste0(ns, ew)
}
tiles$name <- mapply(tile_name, tiles$lon, tiles$lat)
base_url <- "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_"

reclass_weights <- matrix(c(
  10, 0,   # tree cover
  20, 0,   # shrubland
  30, 0,   # grassland
  40, 0.3, # cropland
  50, 0,   # built-up
  60, 0,   # bare/sparse
  70, 0,   # snow/ice
  80, 0.1, # permanent water
  90, 1.0, # herbaceous wetland
  95, 0,   # mangroves (outside study area, but just in case)
  100, 0   # moss/lichen
), ncol = 2, byrow = TRUE)

pieces <- list()
n_ok <- 0; n_missing <- 0

for (i in seq_len(nrow(tiles))) {
  tn <- tiles$name[i]
  url <- paste0(base_url, tn, "_Map.tif")
  raw_f <- file.path("data/worldcover_tmp", paste0(tn, "_raw.tif"))
  dec_f <- file.path("data/worldcover_tmp", paste0(tn, "_100m.tif"))

  dl <- tryCatch(
    download.file(url, raw_f, quiet = TRUE, mode = "wb", method = "libcurl"),
    error = function(e) 1L
  )
  if (!isTRUE(dl == 0) || !file.exists(raw_f) || file.info(raw_f)$size < 1000) {
    n_missing <- n_missing + 1
    unlink(raw_f)
    next
  }
  n_ok <- n_ok + 1

  ok <- system2("gdalwarp",
                c("-r", "near", "-tr", "0.0008333333", "0.0008333333",
                  "-overwrite", "-co", "COMPRESS=LZW", "-q", raw_f, dec_f))
  unlink(raw_f)
  if (ok != 0 || !file.exists(dec_f)) { unlink(dec_f); next }

  r <- tryCatch(rast(dec_f), error = function(e) NULL)
  if (is.null(r)) { unlink(dec_f); next }

  suit <- classify(r, reclass_weights, others = 0)
  tmpl_sub <- crop(template, ext(suit))
  suit_lr <- resample(suit, tmpl_sub, method = "average")
  pieces[[length(pieces) + 1]] <- suit_lr

  unlink(dec_f)
  message(sprintf("[%d/%d] %s: ok (n_ok=%d, n_missing=%d)", i, nrow(tiles), tn, n_ok, n_missing))
}

message("Total tiles with land data: ", n_ok, "; empty/ocean-only: ", n_missing)

# ---- 3. mosaic all tiles onto the full template grid ------------------------
habitat_prior <- rast(template)
values(habitat_prior) <- 0
names(habitat_prior) <- "habitat_suitability"

mos <- do.call(merge, pieces)
habitat_prior <- merge(mos, habitat_prior)  # fills any gaps outside tile coverage with 0

# Force an exact grid + CRS match to the isoscape template: repeated
# merge()/resample() through this pipeline leaves tiny (~1e-4 degree)
# floating point drift in the extent, and terra's CRS WKT for d2h_GS.tif
# ("unknown" datum based on WGS84) is textually different from (though
# equivalent to) the WorldCover tiles' full EPSG:4326 WKT. Both are enough
# to fail assignR::pdRaster()'s internal compareGeom() check on the prior,
# so pin both explicitly rather than relying on approximate equality.
habitat_prior <- resample(habitat_prior, template, method = "near")
crs(habitat_prior) <- crs(template)
names(habitat_prior) <- "habitat_suitability"

writeRaster(habitat_prior, "data/habitat_prior.tif", overwrite = TRUE)
unlink("data/worldcover_tmp", recursive = TRUE)

message("Wrote data/habitat_prior.tif")
v <- values(habitat_prior, mat = FALSE)
message("Nonzero suitability cells: ", sum(v > 0, na.rm = TRUE), " of ", length(v))
