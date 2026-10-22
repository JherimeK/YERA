# 06_sensitivity_check.R
# Comprehensive sensitivity analysis for the primary model (d2H only,
# Passerine calibration, current habitat prior -- scripts/01-02). Reruns
# the assignment under five alternative configurations and compares each
# against the primary model, individual-by-individual and at the
# population level:
#
#   A. dual_isotope      -- d2H + d18O (the model this pipeline used to
#                            treat as primary; see scripts/01_run_assignment.R
#                            for why that changed)
#   B. no_prior           -- d2H only, no habitat suitability prior at all
#   C. weak_prior         -- d2H only, habitat_prior^0.5 (compressed toward
#                            a uniform/no-op prior)
#   D. strong_prior       -- d2H only, habitat_prior^2 (exaggerated toward
#                            a hard mask)
#   E. waterbird_calib    -- d2H only, Water bird reference-group
#                            calibration instead of Passerine, current prior
#
# This directly answers three questions raised in external review of this
# pipeline: (1) how much does including d18O change the result despite its
# weak, highly-collinear-with-d2H calibration; (2) how much does the
# subjective habitat-prior weighting scheme drive the result; (3) how much
# does the passerine-vs-waterbird calibration choice matter. If results are
# stable across B-E, the specific choice of prior strength/calibration group
# is not doing the work; if results swing the way they did for A in the
# previous (dual-isotope-primary) version of this analysis, that choice
# matters and should be reported as such.
#
# RUN ORDER: after 01 and 02 (needs outputs/rasters/pd_per_individual_category.tif,
# outputs/rasters/pd_population_by_category.tif, outputs/calibration_objects.rds,
# and outputs/tables/assignment_summary.csv for the primary-model numbers).
# Writes outputs/tables/sensitivity_comparison.csv (per-individual and
# population rows) and outputs/tables/habitat_prior_by_zone.csv (the
# Major-Issue-4 check: how much of the habitat prior's suitable area/mass
# falls in each zone, independent of any isotope assignment).

suppressMessages({
  library(terra); library(sf); library(dplyr); library(assignR); library(maps)
})
set.seed(42)
sf_use_s2(FALSE)

# NOTE: this script recomputes calibration/isoscape objects from scratch
# (mirroring scripts/01_run_assignment.R) rather than reading
# outputs/calibration_objects.rds. terra::SpatRaster objects wrap external
# C++ pointers that do not survive saveRDS()/readRDS() across separate R
# sessions/processes (they come back with dead/NULL pointers and error on
# first use) -- the .rds file is kept for inspection within a single live
# R session, but every script that runs as its own `Rscript` process (as
# all of these do) needs to rebuild its own raster objects.
samp <- read.csv("data/klma_clean.csv", stringsAsFactors = FALSE)
samp$sample_id <- paste(samp$Envelope_ID, samp$Feather_type_clean, samp$Original_ID, sep = "__")
unknown_H  <- data.frame(ID = samp$sample_id, d2H = samp$d2H)
unknown_HO <- data.frame(ID = samp$sample_id, d2H = samp$d2H, d18O = samp$d18O)

# ---- rebuild AOI + zones (same logic as scripts 01-02) ---------------------
us <- st_as_sf(maps::map("state", regions = c("california", "nevada", "arizona", "new mexico",
                                               "oregon", "washington", "idaho", "utah"),
                          fill = TRUE, plot = FALSE))
mx <- st_as_sf(maps::map("world", "Mexico", fill = TRUE, plot = FALSE))
st_crs(us) <- 4326; st_crs(mx) <- 4326
us <- st_make_valid(us); mx <- st_make_valid(mx)
bbox_mx <- st_bbox(c(xmin = -118, xmax = -101, ymin = 14, ymax = 33), crs = 4326)
mx_clip <- suppressWarnings(st_crop(mx, bbox_mx))
aoi_sf <- st_as_sf(st_union(st_union(st_geometry(us)), st_union(st_geometry(mx_clip))))
aoi <- vect(aoi_sf); crs(aoi) <- "EPSG:4326"
mx_clip_u <- st_union(mx_clip)
baja_bbox <- st_as_sfc(st_bbox(c(xmin = -118, xmax = -111.5, ymin = 14, ymax = 33), crs = 4326))
mainland_bbox <- st_as_sfc(st_bbox(c(xmin = -111.5, xmax = -101, ymin = 14, ymax = 33), crs = 4326))
baja_sf <- st_as_sf(suppressWarnings(st_intersection(mx_clip_u, baja_bbox)))
mainland_mx_sf <- st_as_sf(suppressWarnings(st_intersection(mx_clip_u, mainland_bbox)))
mexico_v <- vect(st_as_sf(st_union(st_geometry(baja_sf), st_geometry(mainland_mx_sf))))
pnw_v <- vect(st_as_sf(us[us$ID %in% c("oregon", "washington", "idaho", "utah", "california"), ]))
crs(mexico_v) <- "EPSG:4326"; crs(pnw_v) <- "EPSG:4326"

# ---- rebuild isoscapes + calibration (identical to script 01) -------------
ASSUMED_SD_D2H <- 5.0
ASSUMED_SD_D18O <- 0.6
na_ext <- ext(-170, -50, 5, 75)
d2h_mean <- crop(rast("d2h_GS.tif"), na_ext)
d18o_mean <- crop(rast("d18o_GS.tif"), na_ext)
d18o_mean <- resample(d18o_mean, d2h_mean, method = "bilinear")
d2h_iso <- c(d2h_mean, setValues(d2h_mean, ASSUMED_SD_D2H)); names(d2h_iso) <- c("mean", "sd")
d18o_iso <- c(d18o_mean, setValues(d18o_mean, ASSUMED_SD_D18O)); names(d18o_iso) <- c("mean", "sd")

data(naMap)
known_H <- subOrigData(marker = "d2H", group = "Passerine", mask = naMap, niter = 2000, genplot = FALSE)
known_O <- subOrigData(marker = "d18O", group = "Passerine", ref_scale = "VSMOW_O", mask = naMap, niter = 2000, genplot = FALSE)
known_H_wb <- subOrigData(marker = "d2H", group = "Water bird", mask = naMap, niter = 2000, genplot = FALSE)
cal_H <- calRaster(known_H, d2h_iso, genplot = FALSE, verboseLM = FALSE)
# These two calibrations are only used by the sensitivity variants below
# (dual-isotope and waterbird_calibration respectively), not the primary
# model -- their diagnostic plots are produced here rather than in
# scripts/01_run_assignment.R for exactly that reason.
dir.create("outputs/maps", showWarnings = FALSE)
png("outputs/maps/calibration_regression_d18O.png", width = 900, height = 700)
cal_O <- calRaster(known_O, d18o_iso, genplot = TRUE, verboseLM = TRUE)
dev.off()
png("outputs/maps/calibration_regression_d2H_waterbird_sensitivity.png", width = 900, height = 700)
cal_H_wb <- calRaster(known_H_wb, d2h_iso, genplot = TRUE, verboseLM = TRUE)
dev.off()

aoi_crop_mask <- function(r) mask(crop(r, aoi), aoi)
cal_H$isoscape.rescale <- aoi_crop_mask(cal_H$isoscape.rescale)
cal_O$isoscape.rescale <- aoi_crop_mask(cal_O$isoscape.rescale)
cal_H_wb$isoscape.rescale <- aoi_crop_mask(cal_H_wb$isoscape.rescale)

habitat_prior <- rast("data/habitat_prior.tif")
habitat_prior <- resample(aoi_crop_mask(habitat_prior), cal_H$isoscape.rescale, method = "near")
crs(habitat_prior) <- crs(cal_H$isoscape.rescale)
habitat_prior_weak <- habitat_prior^0.5
habitat_prior_strong <- habitat_prior^2

cal <- list(cal_H = cal_H, cal_O = cal_O, cal_H_wb = cal_H_wb,
            habitat_prior = habitat_prior, habitat_prior_weak = habitat_prior_weak,
            habitat_prior_strong = habitat_prior_strong)

mass_in_zone <- function(layer, zonevect) {
  e <- terra::extract(layer, zonevect, fun = sum, na.rm = TRUE, ID = FALSE)
  s <- sum(unlist(e), na.rm = TRUE)
  if (length(s) == 0 || is.na(s)) 0 else s
}
combine_feathers <- function(sub) {
  if (nlyr(sub) == 1) return(sub[[1]])
  prod <- sub[[1]]
  for (i in 2:nlyr(sub)) prod <- prod * sub[[i]]
  s <- global(prod, sum, na.rm = TRUE)[1, 1]
  if (is.na(s) || s <= 0) {
    m <- app(sub, fun = mean, na.rm = TRUE)
    m / global(m, sum, na.rm = TRUE)[1, 1]
  } else prod / s
}
area_km2 <- function(qtl_layer, cellkm2) {
  vals <- values(qtl_layer, mat = FALSE); ok <- which(vals == 1)
  sum(values(cellkm2, mat = FALSE)[ok], na.rm = TRUE)
}
groups <- samp %>% distinct(Envelope_ID, Feather_cat)
build_combo <- function(pd_feather) {
  layers <- list()
  for (i in seq_len(nrow(groups))) {
    eid <- groups$Envelope_ID[i]; fcat <- groups$Feather_cat[i]
    ids <- samp$sample_id[samp$Envelope_ID == eid & samp$Feather_cat == fcat]
    sub <- pd_feather[[match(ids, names(pd_feather))]]
    outname <- paste(eid, fcat, sep = "__")
    combo <- combine_feathers(sub); names(combo) <- outname
    layers[[outname]] <- combo
  }
  rast(layers)
}
population_surface <- function(combo_stack, fcat) {
  sel <- grep(paste0("__", fcat, "$"), names(combo_stack))
  m <- app(combo_stack[[sel]], fun = mean, na.rm = TRUE)
  m / global(m, sum, na.rm = TRUE)[1, 1]
}
summarize_variant <- function(pd_feather, label) {
  combo <- build_combo(pd_feather)
  cellkm2 <- cellSize(combo[[1]], unit = "km")
  qtl <- qtlRaster(combo, threshold = 0.5, thresholdType = "prob", genplot = FALSE)
  rows <- lapply(names(combo), function(n) {
    data.frame(surface = n, variant = label,
               P_Mexico = round(mass_in_zone(combo[[n]], mexico_v), 4),
               P_PacNW = round(mass_in_zone(combo[[n]], pnw_v), 4),
               area50_km2 = round(area_km2(qtl[[n]], cellkm2)))
  })
  pop_rows <- lapply(c("Body", "Flight"), function(fcat) {
    pb <- population_surface(combo, fcat)
    qp <- qtlRaster(pb, threshold = 0.5, thresholdType = "prob", genplot = FALSE)
    data.frame(surface = paste0("Population_", fcat), variant = label,
               P_Mexico = round(mass_in_zone(pb, mexico_v), 4),
               P_PacNW = round(mass_in_zone(pb, pnw_v), 4),
               area50_km2 = round(area_km2(qp, cellkm2)))
  })
  rbind(do.call(rbind, rows), do.call(rbind, pop_rows))
}

all_rows <- list()

# ---- PRIMARY (reuse scripts/01-02 output; no rerun needed) ----------------
message("Summarizing PRIMARY model (d2H only, Passerine, current prior) from existing output...")
primary_tbl <- read.csv("outputs/tables/assignment_summary.csv", stringsAsFactors = FALSE)
primary_tbl$P_Mexico <- primary_tbl$P_Baja_CA + primary_tbl$P_Mainland_W_Mexico
primary_tbl$P_PacNW <- primary_tbl$P_Oregon + primary_tbl$P_Washington_Idaho_Utah + primary_tbl$P_California
all_rows$primary <- data.frame(surface = primary_tbl$surface, variant = "primary_d2H_passerine_currentprior",
                                P_Mexico = round(primary_tbl$P_Mexico, 4), P_PacNW = round(primary_tbl$P_PacNW, 4),
                                area50_km2 = primary_tbl$area50_km2)

# ---- A. dual-isotope (current prior) ---------------------------------------
message("A. Running dual-isotope (d2H+d18O) sensitivity variant...")
iso_stack <- isoStack(cal$cal_H, cal$cal_O)
pd_dual <- pdRaster(iso_stack, unknown = unknown_HO, mask = NULL, prior = cal$habitat_prior, genplot = FALSE)
names(pd_dual) <- unknown_HO$ID
all_rows$dual <- summarize_variant(pd_dual, "dual_isotope")

# ---- B. no habitat prior (d2H only) ----------------------------------------
message("B. Running no-habitat-prior sensitivity variant...")
pd_noprior <- pdRaster(cal$cal_H, unknown = unknown_H, mask = NULL, prior = NULL, genplot = FALSE)
names(pd_noprior) <- unknown_H$ID
all_rows$noprior <- summarize_variant(pd_noprior, "no_prior")

# ---- C. weak habitat prior (d2H only, prior^0.5) ---------------------------
message("C. Running weak-habitat-prior sensitivity variant...")
pd_weak <- pdRaster(cal$cal_H, unknown = unknown_H, mask = NULL, prior = cal$habitat_prior_weak, genplot = FALSE)
names(pd_weak) <- unknown_H$ID
all_rows$weak <- summarize_variant(pd_weak, "weak_prior")

# ---- D. strong habitat prior (d2H only, prior^2) ---------------------------
message("D. Running strong-habitat-prior sensitivity variant...")
pd_strong <- pdRaster(cal$cal_H, unknown = unknown_H, mask = NULL, prior = cal$habitat_prior_strong, genplot = FALSE)
names(pd_strong) <- unknown_H$ID
all_rows$strong <- summarize_variant(pd_strong, "strong_prior")

# ---- E. Water bird calibration (d2H only, current prior) -------------------
message("E. Running Water bird calibration sensitivity variant...")
pd_wb <- pdRaster(cal$cal_H_wb, unknown = unknown_H, mask = NULL, prior = cal$habitat_prior, genplot = FALSE)
names(pd_wb) <- unknown_H$ID
all_rows$waterbird <- summarize_variant(pd_wb, "waterbird_calibration")

# ---- combine into one wide comparison table --------------------------------
long <- do.call(rbind, all_rows)
wide <- reshape(long, idvar = "surface", timevar = "variant",
                 v.names = c("P_Mexico", "P_PacNW", "area50_km2"), direction = "wide")
wide <- wide[order(wide$surface), ]
dir.create("outputs/tables", showWarnings = FALSE)
write.csv(wide, "outputs/tables/sensitivity_comparison.csv", row.names = FALSE)
message("\nWrote outputs/tables/sensitivity_comparison.csv (", nrow(wide), " rows x ", length(all_rows), " model variants)")

pop_wide <- wide[grepl("^Population", wide$surface), ]
print(pop_wide, row.names = FALSE)

# ---- Major Issue 4: habitat prior's own distribution across zones ---------
# Independent of any isotope data: how much of the habitat suitability
# prior's mass (and nonzero area) falls in each zone? If Mexico habitat is
# systematically under-represented relative to the US side (e.g. because
# ESA WorldCover classifies ephemeral/seasonal Mexican wetlands less
# completely than US wetlands), that alone -- with no isotope signal
# needed -- would bias assignments northward, and vice versa.
message("\nComputing habitat-prior distribution by zone...")
hp <- cal$habitat_prior
zone_list <- list(Mexico = mexico_v, Pacific_NW = pnw_v)
cell_km2_hp <- cellSize(hp, unit = "km")
hp_zone_rows <- lapply(names(zone_list), function(zn) {
  zv <- zone_list[[zn]]
  mass <- mass_in_zone(hp, zv)
  total_mass <- global(hp, sum, na.rm = TRUE)[1, 1]
  area_cells <- terra::extract(hp > 0.01, zv, fun = sum, na.rm = TRUE, ID = FALSE)
  area_km2_suitable <- sum(unlist(area_cells), na.rm = TRUE) *
    mean(values(cell_km2_hp, mat = FALSE), na.rm = TRUE)
  data.frame(zone = zn, pct_of_total_habitat_mass = round(100 * mass / total_mass, 2),
             suitable_area_km2_gt_0.01 = round(area_km2_suitable))
})
hp_zone_tbl <- do.call(rbind, hp_zone_rows)
write.csv(hp_zone_tbl, "outputs/tables/habitat_prior_by_zone.csv", row.names = FALSE)
print(hp_zone_tbl, row.names = FALSE)
message("Wrote outputs/tables/habitat_prior_by_zone.csv")
