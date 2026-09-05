# run_molt_assignments_by_feathertype.R
# Per-individual, per-feather-type isotope assignment maps (δ2H + δ18O)
# Uses isoscape rasters d2h_GS.tif and d18o_GS.tif in repo root and SIA_Results.csv
# Outputs written to molt_assignments_outputs/ and committed back to the repo by the Action

library(raster)
library(sf)
library(dplyr)
library(assignR)
library(sp)

# ---- User inputs ----
input_csv <- "SIA_Results.csv"
out_dir <- "molt_assignments_outputs"
dir.create(out_dir, showWarnings = FALSE)

# focal points (change as needed)
focal_points <- data.frame(
  site = c("Klamath", "SF_Bay", "San_Diego"),
  lon = c(-121.6008333, -122.4194, -117.1611),
  lat = c(42.9766389, 37.7749, 32.7157),
  stringsAsFactors = FALSE
)
radius_km <- 50

# calibrations
sd_d2H_precip_calib <- 18.3 / 0.83   # ≈ 22.05 ‰
sd_d18O_precip_calib <- 1.37

# ---- Load data ----
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = TRUE)

# detect columns
feather_col <- names(df)[grepl("^Feather$", names(df), ignore.case=TRUE)][1]
envelope_col <- names(df)[grepl("ID.*envelop|Envelope|ID_on_envelope|ID.on.envelope|ID on envelope", names(df), ignore.case=TRUE)][1]
if(is.na(feather_col)) stop("No 'Feather' column found in CSV (case-insensitive).")
if(is.na(envelope_col)) {
  # fallback
  envelope_col <- names(df)[1]
  warning(paste("Could not detect 'ID on envelope' column; using", envelope_col))
}

d2H_col <- names(df)[grepl("d2h|d2H|delta2H", names(df), ignore.case=TRUE)][1]
d18O_col <- names(df)[grepl("d18o|d18O|delta18O", names(df), ignore.case=TRUE)][1]
if(is.na(d2H_col) || is.na(d18O_col)) stop("Could not find d2H and/or d18O columns in CSV.")

# load isoscapes
d2H_file <- list.files(".", pattern="(?i)d2h.*tif$", full.names=TRUE, perl=TRUE)[1]
d18O_file <- list.files(".", pattern="(?i)d18o.*tif$", full.names=TRUE, perl=TRUE)[1]
if(is.na(d2H_file) || is.na(d18O_file)) stop("d2h_GS.tif or d18o_GS.tif not found in repo root.")

iso_d2H <- raster(d2H_file)
iso_d18O <- raster(d18O_file)
study_ext <- extent(-130, -110, 20, 55)
iso_d2H <- crop(iso_d2H, study_ext)
iso_d18O <- crop(iso_d18O, study_ext)
iso_d18O <- resample(iso_d18O, iso_d2H, method="bilinear")

# compute predicted precipitation values
suppressWarnings(
  df <- df %>% mutate(
    d2H_meas = as.numeric(get(d2H_col)),
    d18O_meas = as.numeric(get(d18O_col)),
    d2H_precip_pred = (d2H_meas + 18.0) / 0.83,
    d18O_precip_pred = 0.35 * d18O_meas + 4.55,
    sd_d2H_precip = sd_d2H_precip_calib,
    sd_d18O_precip = sd_d18O_precip_calib
  )
)

# names for outputs
if(!("Original.ID" %in% names(df))) df$Original.ID <- seq_len(nrow(df))
df$sample_name <- paste0(gsub("[^A-Za-z0-9_\\-]", "_", df[[envelope_col]]),
                         "_ID", df$Original.ID, "_", gsub("[^A-Za-z0-9_\\-]", "_", df[[feather_col]]))

# classify feather types
df$feathertype <- ifelse(grepl("flight|rectrix|secondary|primary|tertial|primary|rectrice|rectrice|rectric", df[[feather_col]], ignore.case=TRUE),
                         "Flight", "Body")

# per-sample pd rasters
for(i in seq_len(nrow(df))) {
  s <- df[i, ]
  if(is.na(s$d2H_precip_pred) || is.na(s$d18O_precip_pred)) next
  pH <- pdRaster(iso_d2H, s$d2H_precip_pred, total.sd = s$sd_d2H_precip)
  pO <- pdRaster(iso_d18O, s$d18O_precip_pred, total.sd = s$sd_d18O_precip)
  pH[pH[]==0] <- 1e-12; pO[pO[]==0] <- 1e-12
  ll <- log(pH) + log(pO)
  comb <- exp(ll - max(values(ll), na.rm=TRUE))
  comb <- comb / cellStats(comb, stat='sum', na.rm=TRUE)
  sname <- s$sample_name
  writeRaster(pH, filename=file.path(out_dir, paste0("pd_d2H_", sname, ".tif")), overwrite=TRUE)
  writeRaster(pO, filename=file.path(out_dir, paste0("pd_d18O_", sname, ".tif")), overwrite=TRUE)
  writeRaster(comb, filename=file.path(out_dir, paste0("pd_dual_", sname, ".tif")), overwrite=TRUE)
  png(file.path(out_dir, paste0("map_dual_", sname, ".png")), width=900, height=600)
  plot(comb, main=paste("Dual p(δ2H+δ18O) ", sname))
  dev.off()
}

# combine per-individual × feather-type
ids <- unique(df[[envelope_col]])
summary_rows <- list()
for(id in ids) {
  for(ft in c("Flight","Body")) {
    rows <- df %>% filter(get(envelope_col) == id & feathertype == ft)
    if(nrow(rows)==0) next
    rasters <- list()
    for(r in rows$sample_name) {
      fpath <- file.path(out_dir, paste0("pd_dual_", r, ".tif"))
      if(file.exists(fpath)) rasters[[r]] <- raster(fpath)
    }
    if(length(rasters)==0) next
    stack_r <- stack(rasters)
    logsum <- calc(stack_r, fun=function(...) {
      vals <- c(...)
      if(all(is.na(vals))) return(NA)
      return(sum(log(vals + 1e-12)))
    })
    joint <- exp(logsum - max(values(logsum), na.rm=TRUE))
    joint <- joint / cellStats(joint, stat='sum', na.rm=TRUE)
    outname <- paste0(gsub("[^A-Za-z0-9_\\-]","_", id), "_", ft)
    writeRaster(joint, filename=file.path(out_dir, paste0("joint_pd_", outname, ".tif")), overwrite=TRUE)
    png(file.path(out_dir, paste0("map_joint_", outname, ".png")), width=1000, height=700)
    plot(joint, main=paste("Joint dual p (", id, " / ", ft, ")"))
    dev.off()
    # contours
    for(p in c(0.5, 0.75, 0.9)) {
      lvl <- getPDRlevel(joint, p = p)
      bin <- joint >= lvl
      writeRaster(bin, filename=file.path(out_dir, paste0("joint_", outname, "_", p*100, "pc.tif")), overwrite=TRUE)
    }
    maxcell <- which.max(values(joint))
    coords_mat <- xyFromCell(joint, maxcell)
    areas <- sapply(c(0.5,0.75,0.9), function(p) {
      lvl <- getPDRlevel(joint, p=p)
      bin <- joint[] >= lvl
      if(all(is.na(bin))) return(NA)
      cell_area_km2 <- area(joint) / 1e6
      sum(cell_area_km2[which(bin)], na.rm=TRUE)
    })
    # probability mass within focal buffers
    pts <- rasterToPoints(joint, spatial=TRUE)
    pts_sf <- st_as_sf(pts)
    probs <- pts_sf$layer
    pmass <- c()
    pts_ll <- st_transform(pts_sf, 4326)
    for(fp in seq_len(nrow(focal_points))) {
      pt_sf <- st_sfc(st_point(c(focal_points$lon[fp], focal_points$lat[fp])), crs = 4326)
      d <- as.numeric(st_distance(pts_ll, pt_sf))
      sel <- d <= (radius_km * 1000)
      pmass[focal_points$site[fp]] <- sum(probs[which(sel)], na.rm=TRUE)
    }
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      envelope_id = id, feather_type = ft,
      n_feathers = length(rasters),
      centroid_lon = coords_mat[1], centroid_lat = coords_mat[2],
      area_50_km2 = areas[1], area_75_km2 = areas[2], area_90_km2 = areas[3],
      prob_in_Klamath_50km = pmass["Klamath"],
      prob_in_SF_Bay_50km = pmass["SF_Bay"],
      prob_in_San_Diego_50km = pmass["San_Diego"],
      stringsAsFactors = FALSE
    )
  }
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file=file.path(out_dir, "molt_assignment_summary_by_individual_and_feathertype.csv"), row.names=FALSE)

# end
