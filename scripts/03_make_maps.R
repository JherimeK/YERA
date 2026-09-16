# 03_make_maps.R
# Final presentation maps: population-level Body vs Flight comparison, and a
# faceted grid of all individual-category posteriors.
#
# Colors are classified directly by credible region (50% / 50-75% / 75-90% /
# outside 90%, from script 02's thresholdType = "prob" rasters), not by
# quantile-of-pixel-value. An earlier version used quantile color breaks,
# which -- because the underlying probability surfaces are extremely
# right-skewed -- made broad, actually-low-probability areas look "hot" by
# relative rank even though they hold very little of the true cumulative
# probability. That made the map's coloring and the credible-region area
# stats tell visually different stories. Classifying by credible region
# directly means the colored area on the map IS the credible-region area
# in the stats table, by construction.

suppressMessages({
  library(terra)
  library(sf)
  library(maps)
})
sf_use_s2(FALSE)

pop <- rast("outputs/rasters/pd_population_by_category.tif")
indiv <- rast("outputs/rasters/pd_per_individual_category.tif")
qtl50 <- rast("outputs/rasters/qtl50_credible_region.tif")
qtl75 <- rast("outputs/rasters/qtl75_credible_region.tif")
qtl90 <- rast("outputs/rasters/qtl90_credible_region.tif")

us_sf <- st_as_sf(maps::map("state", regions = c("california", "nevada", "arizona",
                                                   "new mexico", "oregon", "utah", "idaho",
                                                   "washington", "colorado", "wyoming", "texas"),
                             fill = TRUE, plot = FALSE))
mx_sf <- st_as_sf(maps::map("world", "Mexico", fill = TRUE, plot = FALSE))
st_crs(us_sf) <- 4326; st_crs(mx_sf) <- 4326

focal_points <- data.frame(
  site = c("Klamath Marsh\n(capture site)", "SF Bay", "San Diego"),
  lon = c(-121.6008333, -122.0262222, -117.2288611),
  lat = c(42.9766389, 37.4533611, 32.7931111)
)

# Credible-region class: 4 = top 50% (most likely), 3 = 50-75%, 2 = 75-90%,
# 1 = outside the 90% credible region (essentially negligible probability).
cr_class <- function(qname) {
  c50 <- qtl50[[qname]]; c75 <- qtl75[[qname]]; c90 <- qtl90[[qname]]
  cls <- terra::ifel(c50 == 1, 4, terra::ifel(c75 == 1, 3, terra::ifel(c90 == 1, 2, 1)))
  names(cls) <- qname
  cls
}

cr_colors <- c("grey94", "#FED976", "#FD8D3C", "#B10026")
cr_labels <- c("outside 90% region", "75-90%", "50-75%", "top 50%")

plot_one <- function(qname, main) {
  cls <- cr_class(qname)
  plot(cls, breaks = c(0.5, 1.5, 2.5, 3.5, 4.5), col = cr_colors, type = "interval",
       main = main, mar = c(2, 2, 3, 9), axes = TRUE,
       plg = list(cex = 0.65, legend = cr_labels, title = "credible region"))
  plot(st_geometry(us_sf), add = TRUE, border = "grey30", lwd = 0.6)
  plot(st_geometry(mx_sf), add = TRUE, border = "grey30", lwd = 0.6)
  points(focal_points$lon, focal_points$lat, pch = 17, col = "blue", cex = 1.1)
  text(focal_points$lon, focal_points$lat, labels = focal_points$site,
       pos = 4, cex = 0.6, col = "blue4", offset = 0.4)
}

png("outputs/maps/population_body_vs_flight.png", width = 1800, height = 950, res = 130)
par(mfrow = c(1, 2))
plot_one("Population_Body", "KLMA Yellow Rail -- Body feathers\n(population mean origin probability)")
plot_one("Population_Flight", "KLMA Yellow Rail -- Flight feathers\n(population mean origin probability)")
dev.off()

# ---- individual grids: same credible-region classification per panel -------
# type = "interval" (not "classes") for the multi-layer case: terra's
# "classes" plotting auto-detects factor levels per layer and can assign
# inconsistent colors/legends across panels when different birds don't
# have all 4 classes present. "interval" with explicit breaks avoids that
# entirely and gives every panel the same fixed color mapping.
make_grid <- function(idx, outfile, w = 1900, h = 1500) {
  cls_list <- lapply(names(indiv)[idx], cr_class)
  cls_stack <- rast(cls_list)
  names(cls_stack) <- names(indiv)[idx]
  png(outfile, width = w, height = h, res = 140)
  plot(cls_stack, breaks = c(0.5, 1.5, 2.5, 3.5, 4.5), col = cr_colors, type = "interval",
       nc = 4, plg = list(cex = 0.5, legend = cr_labels, title = "credible region"),
       mar = c(1.5, 1.5, 2, 7))
  dev.off()
}

nm <- names(indiv)
make_grid(grep("Body$", nm), "outputs/maps/individuals_body.png")
make_grid(grep("Flight$", nm), "outputs/maps/individuals_flight.png")

# ---- habitat prior diagnostic map -------------------------------------------
# Shows the wetland/wet-cropland suitability layer itself (before it's
# combined with the isotope likelihood), so it's possible to sanity-check
# that it actually captures Klamath Marsh and known coastal marsh areas.
if (file.exists("data/habitat_prior.tif")) {
  hp <- rast("data/habitat_prior.tif")
  hp_aoi <- crop(hp, ext(-125, -101, 17, 49))
  png("outputs/maps/habitat_prior.png", width = 1100, height = 950, res = 130)
  plot(hp_aoi, col = hcl.colors(50, "Greens 3", rev = TRUE),
       main = "Wetland / wet-cropland habitat suitability prior\n(ESA WorldCover 2021)",
       mar = c(2, 2, 3, 5), axes = TRUE, plg = list(cex = 0.7, title = "suitability"))
  plot(st_geometry(us_sf), add = TRUE, border = "grey40", lwd = 0.6)
  plot(st_geometry(mx_sf), add = TRUE, border = "grey40", lwd = 0.6)
  points(focal_points$lon, focal_points$lat, pch = 17, col = "blue", cex = 1.1)
  text(focal_points$lon, focal_points$lat, labels = focal_points$site,
       pos = 4, cex = 0.65, col = "blue4", offset = 0.4)
  dev.off()
}

message("Maps written to outputs/maps/")
