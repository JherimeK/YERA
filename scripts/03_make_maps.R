# 03_make_maps.R
# Final presentation maps: population-level Body vs Flight comparison, and a
# faceted grid of all individual-category posteriors.
#
# The raw posterior surfaces are extremely right-skewed (the 99.9th
# percentile pixel is often 100-1000x the median pixel), so a plain linear
# color scale makes almost the entire map read as a single flat color and
# hides the pattern. Quantile-based color breaks (computed from the data
# itself) are used instead so the full color ramp is actually used, and the
# 50/75/90% credible regions (from script 02, thresholdType = "prob") are
# drawn as contour outlines for a quantitatively meaningful "legend".

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

quantile_breaks <- function(v, n = 8) {
  v <- v[!is.na(v) & v > 0]
  b <- unique(quantile(v, probs = seq(0, 1, length.out = n), na.rm = TRUE))
  c(0, b)
}

contour_outline <- function(qtl_layer, name, ...) {
  bin <- qtl_layer[[name]]
  if (isTRUE(global(bin, sum, na.rm = TRUE)[1, 1] > 0)) {
    poly <- as.polygons(bin, dissolve = TRUE)
    poly <- poly[poly[[1]][[1]] == 1, ]
    if (nrow(poly) > 0) plot(poly, add = TRUE, ...)
  }
}

plot_one <- function(layer, qname, main) {
  breaks <- quantile_breaks(values(layer, mat = FALSE))
  pal <- hcl.colors(length(breaks) - 1, "YlOrRd", rev = TRUE)
  plot(layer, breaks = breaks, col = pal, type = "interval", main = main,
       mar = c(2, 2, 3, 7), axes = TRUE, plg = list(cex = 0.55, title = "prob./cell"))
  plot(st_geometry(us_sf), add = TRUE, border = "grey30", lwd = 0.6)
  plot(st_geometry(mx_sf), add = TRUE, border = "grey30", lwd = 0.6)
  contour_outline(qtl90, qname, border = "grey20", lwd = 0.8, lty = 3)
  contour_outline(qtl75, qname, border = "grey10", lwd = 1.1, lty = 2)
  contour_outline(qtl50, qname, border = "black", lwd = 1.6, lty = 1)
  points(focal_points$lon, focal_points$lat, pch = 17, col = "blue", cex = 1.1)
  text(focal_points$lon, focal_points$lat, labels = focal_points$site,
       pos = 4, cex = 0.6, col = "blue4", offset = 0.4)
  legend("bottomleft", legend = c("50% credible region", "75%", "90%"),
         lty = c(1, 2, 3), lwd = c(1.6, 1.1, 0.8), col = c("black", "grey10", "grey20"),
         bg = "white", box.col = "grey50", cex = 0.6, inset = 0.02)
}

png("outputs/maps/population_body_vs_flight.png", width = 1700, height = 950, res = 130)
par(mfrow = c(1, 2))
plot_one(pop[["Population_Body"]], "Population_Body",
         "KLMA Yellow Rail -- Body feathers\n(population mean origin probability)")
plot_one(pop[["Population_Flight"]], "Population_Flight",
         "KLMA Yellow Rail -- Flight feathers\n(population mean origin probability)")
dev.off()

# ---- individual grids: shared quantile breaks + terra's native per-panel legend
make_grid <- function(idx, outfile, w = 1900, h = 1500) {
  sub <- indiv[[idx]]
  breaks <- quantile_breaks(values(sub, mat = FALSE))
  pal <- hcl.colors(length(breaks) - 1, "YlOrRd", rev = TRUE)
  png(outfile, width = w, height = h, res = 140)
  plot(sub, breaks = breaks, col = pal, type = "interval", nc = 4,
       plg = list(cex = 0.5, title = "prob./cell"), mar = c(1.5, 1.5, 2, 5))
  dev.off()
}

nm <- names(indiv)
make_grid(grep("Body$", nm), "outputs/maps/individuals_body.png")
make_grid(grep("Flight$", nm), "outputs/maps/individuals_flight.png")

message("Maps written to outputs/maps/")
