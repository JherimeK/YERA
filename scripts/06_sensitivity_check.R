# 06_sensitivity_check.R
# Quantifies how much the dual-isotope (d2H+d18O) assignment depends on
# d18O specifically, given that assignR::pdRaster() warns when input
# isoscapes share >70% of their spatial variance ("two or more isoscapes
# have shared variance > 0.7, added information will be limited, and
# specificity of assignments may be inflated") -- a warning this project's
# raw d2h_GS.tif/d18o_GS.tif isoscapes trigger (r^2 = 0.99 correlated
# across the study area). Reruns the full assignment using d2H alone
# (same calibration approach, habitat prior, and study area as script 01)
# and compares both population-level and individual-level results against
# the primary dual-isotope model from scripts 01-02.
#
# RUN ORDER: after 01 and 02 (needs outputs/rasters/pd_per_feather.tif and
# data/habitat_prior.tif). Writes outputs/tables/sensitivity_d2Honly_vs_dual.csv.
# See REPORT.md "Results, section 1b" and MANUSCRIPT_RESULTS.md for the
# narrative writeup of what this found.

suppressMessages({
  library(terra); library(sf); library(dplyr); library(assignR); library(maps)
})
set.seed(42)
sf_use_s2(FALSE)

ASSUMED_SD_D2H <- 5.0

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

mexico_v <- vect(st_as_sf(st_union(st_geometry(mx_clip))))
pnw_v <- vect(st_as_sf(us[us$ID %in% c("oregon", "washington", "idaho", "utah", "california"), ]))

# ---- 1. check the collinearity directly, as pdRaster's own warning flags ---
na_ext <- ext(-170, -50, 5, 75)
d2h_mean <- crop(rast("d2h_GS.tif"), na_ext)
d18o_mean <- crop(rast("d18o_GS.tif"), na_ext)
d18o_mean <- resample(d18o_mean, d2h_mean, method = "bilinear")
v1 <- values(d2h_mean, mat = FALSE); v2 <- values(d18o_mean, mat = FALSE)
ok <- !is.na(v1) & !is.na(v2)
r2_iso <- cor(v1[ok], v2[ok])^2
message(sprintf("Raw d2H/d18O isoscape correlation: r^2 = %.4f", r2_iso))

# ---- 2. rerun the assignment with d2H alone (mirrors script 01) -----------
d2h_iso <- c(d2h_mean, setValues(d2h_mean, ASSUMED_SD_D2H))
names(d2h_iso) <- c("mean", "sd")

data(naMap)
known_H <- subOrigData(marker = "d2H", group = "Passerine", mask = naMap, niter = 2000, genplot = FALSE)
cal_H <- calRaster(known_H, d2h_iso, genplot = FALSE, verboseLM = FALSE)

aoi_crop_mask <- function(r) mask(crop(r, aoi), aoi)
cal_H$isoscape.rescale <- aoi_crop_mask(cal_H$isoscape.rescale)

habitat_prior <- rast("data/habitat_prior.tif")
habitat_prior <- resample(aoi_crop_mask(habitat_prior), cal_H$isoscape.rescale, method = "near")
crs(habitat_prior) <- crs(cal_H$isoscape.rescale)

samp <- read.csv("data/klma_clean.csv", stringsAsFactors = FALSE)
samp$sample_id <- paste(samp$Envelope_ID, samp$Feather_type_clean, samp$Original_ID, sep = "__")
unknown_df <- data.frame(ID = samp$sample_id, d2H = samp$d2H)

message("Running d2H-only per-feather assignment for ", nrow(unknown_df), " samples (sensitivity check)...")
pd_h <- pdRaster(cal_H, unknown = unknown_df, mask = NULL, prior = habitat_prior, genplot = FALSE)
names(pd_h) <- unknown_df$ID

# ---- 3. combine feathers per individual x category, both models -----------
pd_dual <- rast("outputs/rasters/pd_per_feather.tif")

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
mass_in_zone <- function(layer, zonevect) {
  e <- terra::extract(layer, zonevect, fun = sum, na.rm = TRUE, ID = FALSE)
  s <- sum(unlist(e), na.rm = TRUE)
  if (length(s) == 0 || is.na(s)) 0 else s
}
area_km2 <- function(qtl_layer, cellkm2) {
  vals <- values(qtl_layer, mat = FALSE); ok <- which(vals == 1)
  sum(values(cellkm2, mat = FALSE)[ok], na.rm = TRUE)
}

groups <- samp %>% distinct(Envelope_ID, Feather_cat, Group)
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
combo_dual <- build_combo(pd_dual)
combo_h <- build_combo(pd_h)
cellkm2 <- cellSize(combo_dual[[1]], unit = "km")
qtl_dual <- qtlRaster(combo_dual, threshold = 0.5, thresholdType = "prob", genplot = FALSE)
qtl_h <- qtlRaster(combo_h, threshold = 0.5, thresholdType = "prob", genplot = FALSE)

res <- data.frame(
  surface = names(combo_dual), Group = groups$Group,
  P_Mexico_dual = sapply(names(combo_dual), function(n) mass_in_zone(combo_dual[[n]], mexico_v)),
  P_PacNW_dual  = sapply(names(combo_dual), function(n) mass_in_zone(combo_dual[[n]], pnw_v)),
  P_Mexico_d2Honly = sapply(names(combo_h), function(n) mass_in_zone(combo_h[[n]], mexico_v)),
  P_PacNW_d2Honly  = sapply(names(combo_h), function(n) mass_in_zone(combo_h[[n]], pnw_v)),
  area50_dual_km2    = sapply(names(combo_dual), function(n) area_km2(qtl_dual[[n]], cellkm2)),
  area50_d2Honly_km2 = sapply(names(combo_h), function(n) area_km2(qtl_h[[n]], cellkm2))
)
res[3:6] <- round(res[3:6], 4)
res$area50_dual_km2 <- round(res$area50_dual_km2)
res$area50_d2Honly_km2 <- round(res$area50_d2Honly_km2)
res$ratio_area50_dual_over_d2Honly <- round(res$area50_dual_km2 / res$area50_d2Honly_km2, 3)
res <- res[order(res$surface), ]

# ---- 4. population-level comparison ----------------------------------------
klma_ids <- unique(samp$Envelope_ID[samp$Group == "KLMA"])
pop_for <- function(pd_feather, fcat) {
  layers <- list()
  for (eid in klma_ids) {
    ids <- samp$sample_id[samp$Envelope_ID == eid & samp$Feather_cat == fcat]
    if (length(ids) == 0) next
    layers[[eid]] <- combine_feathers(pd_feather[[match(ids, names(pd_feather))]])
  }
  st <- rast(layers); m <- app(st, fun = mean, na.rm = TRUE)
  m / global(m, sum, na.rm = TRUE)[1, 1]
}
pop_rows <- list()
for (fcat in c("Body", "Flight")) {
  pb_dual <- pop_for(pd_dual, fcat); pb_h <- pop_for(pd_h, fcat)
  qd <- qtlRaster(pb_dual, threshold = 0.5, thresholdType = "prob", genplot = FALSE)
  qh <- qtlRaster(pb_h, threshold = 0.5, thresholdType = "prob", genplot = FALSE)
  pop_rows[[fcat]] <- data.frame(
    surface = paste0("Population_", fcat), Group = "KLMA_summary",
    P_Mexico_dual = round(mass_in_zone(pb_dual, mexico_v), 4),
    P_PacNW_dual = round(mass_in_zone(pb_dual, pnw_v), 4),
    P_Mexico_d2Honly = round(mass_in_zone(pb_h, mexico_v), 4),
    P_PacNW_d2Honly = round(mass_in_zone(pb_h, pnw_v), 4),
    area50_dual_km2 = round(area_km2(qd, cellkm2)),
    area50_d2Honly_km2 = round(area_km2(qh, cellkm2)),
    ratio_area50_dual_over_d2Honly = NA
  )
}
res_full <- rbind(do.call(rbind, pop_rows), res)

dir.create("outputs/tables", showWarnings = FALSE)
write.csv(res_full, "outputs/tables/sensitivity_d2Honly_vs_dual.csv", row.names = FALSE)
message("Wrote outputs/tables/sensitivity_d2Honly_vs_dual.csv (", nrow(res_full), " rows)")
message(sprintf("Population Body:   dual P(Mexico)=%.1f%%  d2H-only P(Mexico)=%.1f%%",
                 100 * pop_rows$Body$P_Mexico_dual, 100 * pop_rows$Body$P_Mexico_d2Honly))
message(sprintf("Population Flight: dual P(Mexico)=%.1f%%  d2H-only P(Mexico)=%.1f%%",
                 100 * pop_rows$Flight$P_Mexico_dual, 100 * pop_rows$Flight$P_Mexico_d2Honly))
message(sprintf("Median individual area50 ratio (dual/d2H-only): %.3f", median(res$ratio_area50_dual_over_d2Honly, na.rm = TRUE)))
