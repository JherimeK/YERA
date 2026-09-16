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

mx_parts <- st_cast(st_make_valid(st_union(mx_clip)), "POLYGON")
mx_cent <- st_coordinates(st_centroid(mx_parts))
baja_sf <- st_as_sf(st_union(mx_parts[mx_cent[, 1] < -111.5]))
mainland_mx_sf <- st_as_sf(st_union(mx_parts[mx_cent[, 1] >= -111.5]))

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
for (z in names(zones)) crs(zones[[z]]) <- "EPSG:4326"

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
  if (nlyr(sub) > 1) {
    png(paste0("outputs/maps/combo_", gsub("[^A-Za-z0-9_]", "_", outname), ".png"),
        width = 800, height = 600)
    plot(combo, main = paste(outname, "(", res$method, ")"))
    dev.off()
  }
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

# ---- 3. credible-region contours (50/75/90% by area) ------------------------
all_maps <- c(indiv_stack, pop_stack)
qtl50 <- qtlRaster(all_maps, threshold = 0.5, thresholdType = "area", genplot = FALSE)
qtl75 <- qtlRaster(all_maps, threshold = 0.25, thresholdType = "area", genplot = FALSE)
qtl90 <- qtlRaster(all_maps, threshold = 0.10, thresholdType = "area", genplot = FALSE)

cell_km2 <- cellSize(all_maps[[1]], unit = "km")

area_for <- function(qtl_layer) {
  vals <- values(qtl_layer, mat = FALSE)
  ok <- which(vals == 1)
  sum(values(cell_km2, mat = FALSE)[ok], na.rm = TRUE)
}

# ---- 4. summary table --------------------------------------------------------
rows <- list()
for (nm in names(all_maps)) {
  lyr <- all_maps[[nm]]
  mx_cell <- which.max(values(lyr, mat = FALSE))
  xy <- xyFromCell(lyr, mx_cell)
  row <- data.frame(
    surface = nm,
    peak_lon = round(xy[1], 3),
    peak_lat = round(xy[2], 3),
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
write.csv(summary_tbl, "outputs/tables/assignment_summary.csv", row.names = FALSE)

message("Wrote outputs/tables/assignment_summary.csv (", nrow(summary_tbl), " rows)")
message("Done: per-individual, population, and zone-probability summaries complete.")
