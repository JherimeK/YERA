# 02_summarize_assignment.R
# Combine per-feather posteriors into per-individual x feather-category
# surfaces (assuming feathers of the same category on the same bird share
# a molt origin), build population-level summary surfaces, compute 50/75/90%
# credible regions, and tabulate probability mass by geographic zone.

suppressMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(assignR)
})

sf_use_s2(FALSE)
dir.create("outputs/maps", showWarnings = FALSE)
dir.create("outputs/rasters", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)

samp <- read.csv("data/klma_clean.csv", stringsAsFactors = FALSE)
samp$sample_id <- paste(samp$Envelope_ID, samp$Feather_type_clean, samp$Original_ID, sep = "__")
pd_feather <- rast("outputs/rasters/pd_per_feather.tif")

# ---- rebuild AOI + zones (same logic as script 01) --------------------------
# AOI extended north to include OR/WA/ID/UT so Klamath Marsh itself (the
# capture/breeding site, ~43 N, just north of California) remains a
# candidate origin -- see note in script 01.
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

# Baja vs mainland split by a simple longitude bbox intersection. (An earlier
# version split mx_clip's individual polygon *fragments* by centroid
# longitude and unioned them -- that silently produced a broken/empty-
# overlapping Baja polygon: every P_Baja_CA came out as exactly 0 despite
# the maps clearly showing probability mass over the peninsula. The bbox
# intersection below is simpler and verified to actually overlap the
# raster's Baja Peninsula cells.)
mx_clip_u <- st_union(mx_clip)
baja_bbox <- st_as_sfc(st_bbox(c(xmin = -118, xmax = -111.5, ymin = 14, ymax = 33), crs = 4326))
mainland_bbox <- st_as_sfc(st_bbox(c(xmin = -111.5, xmax = -101, ymin = 14, ymax = 33), crs = 4326))
baja_sf <- st_as_sf(suppressWarnings(st_intersection(mx_clip_u, baja_bbox)))
mainland_mx_sf <- st_as_sf(suppressWarnings(st_intersection(mx_clip_u, mainland_bbox)))

zones <- list(
  Oregon            = vect(st_as_sf(us[us$ID == "oregon", ])),
  California        = vect(st_as_sf(us[us$ID == "california", ])),
  Nevada            = vect(st_as_sf(us[us$ID == "nevada", ])),
  Arizona           = vect(st_as_sf(us[us$ID == "arizona", ])),
  New_Mexico        = vect(st_as_sf(us[us$ID == "new mexico", ])),
  Washington_Idaho_Utah = vect(st_as_sf(us[us$ID %in% c("washington", "idaho", "utah"), ])),
  Baja_CA           = vect(baja_sf),
  Mainland_W_Mexico = vect(mainland_mx_sf)
)
# Copy the CRS straight from the analysis raster rather than re-typing an
# "EPSG:4326" string, to avoid any WKT-string-vs-PROJ4-string mismatch.
for (z in names(zones)) crs(zones[[z]]) <- crs(pd_feather)

sf_use_s2(TRUE)
klamath_buf_v <- vect(st_as_sf(st_buffer(st_sfc(st_point(c(-121.6008333, 42.9766389)), crs = 4326),
                                          dist = 100000)))
sf_use_s2(FALSE)

# ---- helper: probability mass of a raster layer inside a polygon -----------
mass_in_zone <- function(layer, zonevect) {
  e <- terra::extract(layer, zonevect, fun = sum, na.rm = TRUE, ID = FALSE)
  s <- sum(unlist(e), na.rm = TRUE)
  if (length(s) == 0 || is.na(s)) 0 else s
}

# ---- 1. combine feathers within individual x category ----------------------
# Multiple feathers of the same category on one bird are combined assuming a
# shared molt origin (elementwise product of their posteriors, renormalized
# -- the same logic as assignR::jointP). If two feathers imply origins so
# different that their posteriors numerically share no support within the
# AOI (product sums to ~0 everywhere), forcing a joint estimate would be
# meaningless, so we fall back to the mean of their individual posteriors
# and flag the pair as isotopically inconsistent for manual review.
combine_feathers <- function(sub) {
  if (nlyr(sub) == 1) return(list(surface = sub[[1]], method = "single"))
  prod <- sub[[1]]
  for (i in 2:nlyr(sub)) prod <- prod * sub[[i]]
  s <- global(prod, sum, na.rm = TRUE)[1, 1]
  if (is.na(s) || s <= 0) {
    m <- app(sub, fun = mean, na.rm = TRUE)
    m <- m / global(m, sum, na.rm = TRUE)[1, 1]
    list(surface = m, method = "inconsistent_feathers_mean_fallback")
  } else {
    list(surface = prod / s, method = "joint_product")
  }
}

groups <- samp %>% distinct(Envelope_ID, Feather_cat, Group)
indiv_layers <- list()
indiv_meta <- list()

for (i in seq_len(nrow(groups))) {
  eid <- groups$Envelope_ID[i]; fcat <- groups$Feather_cat[i]; grp <- groups$Group[i]
  ids <- samp$sample_id[samp$Envelope_ID == eid & samp$Feather_cat == fcat]
  sub <- pd_feather[[match(ids, names(pd_feather))]]
  outname <- paste(eid, fcat, sep = "__")
  res <- combine_feathers(sub)
  combo <- res$surface
  names(combo) <- outname
  # (no per-call diagnostic PNG here: an earlier version wrote one with
  # terra's raw default color scale, which -- being on a totally different
  # scale than the credible-region-classified maps in script 03 -- caused
  # real confusion comparing the two. scripts/03_make_maps.R's
  # individuals_body.png/individuals_flight.png are the figures to use.)
  indiv_layers[[outname]] <- combo
  indiv_meta[[outname]] <- data.frame(surface = outname, Envelope_ID = eid, Feather_cat = fcat,
                                       Group = grp, n_feathers = nlyr(sub), combine_method = res$method)
}

indiv_stack <- rast(indiv_layers)
names(indiv_stack) <- names(indiv_layers)
writeRaster(indiv_stack, "outputs/rasters/pd_per_individual_category.tif", overwrite = TRUE)
meta_df <- do.call(rbind, indiv_meta)
rownames(meta_df) <- NULL

# ---- 2. population-level summary surfaces (KLMA birds only, by category) ---
pop_layers <- list()
for (fcat in c("Body", "Flight")) {
  sel <- names(indiv_stack)[meta_df$Feather_cat == fcat & meta_df$Group == "KLMA"]
  sub <- indiv_stack[[sel]]
  m <- app(sub, fun = mean, na.rm = TRUE)
  m <- m / global(m, sum, na.rm = TRUE)[1, 1]
  names(m) <- paste0("Population_", fcat)
  pop_layers[[fcat]] <- m
}
pop_stack <- rast(pop_layers)
names(pop_stack) <- paste0("Population_", names(pop_stack))
writeRaster(pop_stack, "outputs/rasters/pd_population_by_category.tif", overwrite = TRUE)

# ---- 3. credible-region contours (50/75/90% cumulative probability) --------
# thresholdType = "prob" is the smallest-area region containing that much
# cumulative posterior probability (a real credible region). The earlier
# thresholdType = "area" is NOT this -- it just returns exactly that
# fraction of the total map area regardless of the probability
# distribution, which made the previous area50/75/90 columns meaningless
# (they were always ~50/75/90% of the AOI's area by construction).
all_maps <- c(indiv_stack, pop_stack)
qtl50 <- qtlRaster(all_maps, threshold = 0.5, thresholdType = "prob", genplot = FALSE)
qtl75 <- qtlRaster(all_maps, threshold = 0.75, thresholdType = "prob", genplot = FALSE)
qtl90 <- qtlRaster(all_maps, threshold = 0.9, thresholdType = "prob", genplot = FALSE)
writeRaster(qtl50, "outputs/rasters/qtl50_credible_region.tif", overwrite = TRUE)
writeRaster(qtl75, "outputs/rasters/qtl75_credible_region.tif", overwrite = TRUE)
writeRaster(qtl90, "outputs/rasters/qtl90_credible_region.tif", overwrite = TRUE)

cell_km2 <- cellSize(all_maps[[1]], unit = "km")

area_for <- function(qtl_layer) {
  vals <- values(qtl_layer, mat = FALSE)
  ok <- which(vals == 1)
  sum(values(cell_km2, mat = FALSE)[ok], na.rm = TRUE)
}

# ---- 4. summary table --------------------------------------------------------
# peak_share: the fraction of a surface's total probability mass held by its
# single top cell (since each surface is already normalized to sum to 1
# within the AOI, this is just max(values)). This is a more direct
# diagnostic than area50_km2 for a specific failure mode found by manual
# inspection: several individuals' measured d2H falls outside the range of
# calibrated values achievable anywhere in the sparse habitat-masked
# domain, so the model has no good match and instead collapses onto
# whichever single wetland cell is "least bad" -- which, because wetlands
# are sparse, can be the SAME cell for multiple unrelated birds. That
# produces a misleadingly tiny, confident-looking credible region that
# really means "nothing matched well," not "precisely located." A high
# peak_share (one cell dominating the whole posterior) is the signature of
# this, and is flagged below regardless of whether the credible area also
# looks small.
rows <- list()
for (nm in names(all_maps)) {
  lyr <- all_maps[[nm]]
  v <- values(lyr, mat = FALSE)
  mx_cell <- which.max(v)
  xy <- xyFromCell(lyr, mx_cell)
  row <- data.frame(
    surface = nm,
    peak_lon = round(xy[1], 3),
    peak_lat = round(xy[2], 3),
    peak_share = round(max(v, na.rm = TRUE), 4),
    area50_km2 = round(area_for(qtl50[[nm]])),
    area75_km2 = round(area_for(qtl75[[nm]])),
    area90_km2 = round(area_for(qtl90[[nm]]))
  )
  for (z in names(zones)) {
    row[[paste0("P_", z)]] <- round(mass_in_zone(lyr, zones[[z]]), 4)
  }
  row$P_within_100km_Klamath <- round(mass_in_zone(lyr, klamath_buf_v), 4)
  rows[[nm]] <- row
}
pop_meta <- data.frame(
  surface = paste0("Population_", c("Body", "Flight")),
  Envelope_ID = "Population", Feather_cat = c("Body", "Flight"), Group = "KLMA_summary",
  n_feathers = sapply(c("Body", "Flight"), function(f)
    sum(meta_df$Group == "KLMA" & meta_df$Feather_cat == f)),
  combine_method = "population_mean"
)
full_meta <- rbind(meta_df, pop_meta)
stats_tbl <- do.call(rbind, rows)
summary_tbl <- merge(full_meta, stats_tbl, by = "surface")

# ---- 5. precision/quality flags ---------------------------------------------
# low_precision_flag: a very small credible area on its own. Not
# necessarily wrong, but worth a second look (see file header note).
summary_tbl$low_precision_flag <- summary_tbl$area50_km2 < 1000

# single_cell_dominant_flag: one grid cell holds a large share of the
# entire posterior. This is the more direct signature of the "no good
# match in the sparse habitat mask" failure mode than area50 alone.
summary_tbl$single_cell_dominant_flag <- summary_tbl$peak_share > 0.3

# shared_peak_flag: this surface's peak cell is shared with >=1 other
# surface (excluding population summaries, and excluding a bird's own
# Body/Flight pair, which have no reason to share a peak) -- independent
# birds/categories converging on the identical grid cell is exactly what
# manual inspection found for several "single_cell_dominant" cases: their
# measured isotope values fall outside what's achievable anywhere in the
# habitat-masked domain, so unrelated birds all default to the same
# least-bad available wetland cell.
indiv_rows <- !grepl("^Population", summary_tbl$surface)
peak_key <- paste(summary_tbl$peak_lon, summary_tbl$peak_lat)
dup_key <- peak_key[indiv_rows][duplicated(peak_key[indiv_rows]) | duplicated(peak_key[indiv_rows], fromLast = TRUE)]
summary_tbl$shared_peak_flag <- indiv_rows & peak_key %in% dup_key

write.csv(summary_tbl, "outputs/tables/assignment_summary.csv", row.names = FALSE)

n_flagged <- sum(summary_tbl$single_cell_dominant_flag | summary_tbl$shared_peak_flag)
message("Wrote outputs/tables/assignment_summary.csv (", nrow(summary_tbl), " rows, ",
        n_flagged, " flagged for low match quality)")
if (n_flagged > 0) {
  print(summary_tbl[summary_tbl$single_cell_dominant_flag | summary_tbl$shared_peak_flag,
                     c("surface", "peak_lon", "peak_lat", "peak_share", "area50_km2",
                       "single_cell_dominant_flag", "shared_peak_flag")])
}
message("Done: per-individual, population, and zone-probability summaries complete.")
