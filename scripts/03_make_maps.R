# 03_make_maps.R
# Final presentation maps: population-level Body vs Flight comparison, and a
# faceted grid of all individual-category posteriors.

suppressMessages({
  library(terra)
  library(sf)
  library(maps)
})
sf_use_s2(FALSE)

pop <- rast("outputs/rasters/pd_population_by_category.tif")
indiv <- rast("outputs/rasters/pd_per_individual_category.tif")

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

plot_one <- function(layer, main, zlim = NULL) {
  plot(layer, main = main, col = hcl.colors(100, "YlOrRd", rev = TRUE),
       mar = c(2, 2, 3, 4), axes = TRUE)
  plot(st_geometry(us_sf), add = TRUE, border = "grey30", lwd = 0.6)
  plot(st_geometry(mx_sf), add = TRUE, border = "grey30", lwd = 0.6)
  points(focal_points$lon, focal_points$lat, pch = 17, col = "blue", cex = 1.1)
  text(focal_points$lon, focal_points$lat, labels = focal_points$site,
       pos = 4, cex = 0.65, col = "blue4", offset = 0.4)
}

png("outputs/maps/population_body_vs_flight.png", width = 1600, height = 900, res = 130)
par(mfrow = c(1, 2))
plot_one(pop[["Population_Body"]], "KLMA Yellow Rail -- Body feathers\n(population mean origin probability)")
plot_one(pop[["Population_Flight"]], "KLMA Yellow Rail -- Flight feathers\n(population mean origin probability)")
dev.off()

# grid of individuals, split by category
nm <- names(indiv)
body_idx <- grep("Body$", nm)
flight_idx <- grep("Flight$", nm)

png("outputs/maps/individuals_body.png", width = 1800, height = 1400, res = 130)
n <- length(body_idx); nc <- 4; nr <- ceiling(n / nc)
par(mfrow = c(nr, nc), mar = c(1.5, 1.5, 2, 1))
for (i in body_idx) {
  plot(indiv[[i]], main = names(indiv)[i], legend = FALSE,
       col = hcl.colors(100, "YlOrRd", rev = TRUE), axes = FALSE)
  plot(st_geometry(us_sf), add = TRUE, border = "grey40", lwd = 0.4)
  plot(st_geometry(mx_sf), add = TRUE, border = "grey40", lwd = 0.4)
}
dev.off()

png("outputs/maps/individuals_flight.png", width = 1800, height = 1400, res = 130)
n <- length(flight_idx); nc <- 4; nr <- ceiling(n / nc)
par(mfrow = c(nr, nc), mar = c(1.5, 1.5, 2, 1))
for (i in flight_idx) {
  plot(indiv[[i]], main = names(indiv)[i], legend = FALSE,
       col = hcl.colors(100, "YlOrRd", rev = TRUE), axes = FALSE)
  plot(st_geometry(us_sf), add = TRUE, border = "grey40", lwd = 0.4)
  plot(st_geometry(mx_sf), add = TRUE, border = "grey40", lwd = 0.4)
}
dev.off()

message("Maps written to outputs/maps/")
