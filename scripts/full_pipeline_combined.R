# full_pipeline_combined.R
#
# The complete KLMA Yellow Rail geographic-origin pipeline, as a single file,
# for sharing. This is scripts/00_clean_data.R, 05_build_habitat_prior.R,
# 01_run_assignment.R, 02_summarize_assignment.R, 03_make_maps.R,
# 04_clean_cn_data.R, and 06_sensitivity_check.R concatenated VERBATIM, in
# the order they actually need to run (00 -> 05 -> 01 -> 02 -> 03 -> 04 -> 06
# -- NOT numeric file order, since script 04's last step depends on script
# 02's output, and script 06 depends on both 01 and 02; see the run-order
# note in METHODS_WALKTHROUGH.md). Nothing has been changed from the
# individual files beyond adding the section banners below, consolidating
# each file's own library()-loading block into one at the top, and removing
# a short list of now-redundant duplicate setup calls (sf_use_s2(FALSE),
# dir.create(), a repeated extent/AOI/calibration object assignment) that
# are already established earlier in this combined run.
#
# PRIMARY MODEL: d2H alone, Passerine calibration, current habitat prior,
# the 12 KLMA-banded individuals only (43 feathers). Step 7 reruns the
# assignment under five alternative configurations (dual-isotope, no/weak/
# strong habitat prior, Water-bird calibration) and compares each against
# the primary model -- required reading before citing the body-feather
# result specifically. See REPORT.md Results section 1b.
#
# The individual numbered files in scripts/ remain the source of truth for
# running or modifying any one step on its own. This combined file is
# provided purely so the whole pipeline can be reviewed or re-run top to
# bottom from one place. See METHODS_WALKTHROUGH.md for the full narrative
# explanation of each step.
#
# Run from the repository root (working directory must be the repo root,
# since every path below is relative to it, e.g. "data/...", "outputs/...").

suppressMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(stringr)
  library(maps)
  library(assignR)
})


# ==============================================================================
# STEP 1 (scripts/00_clean_data.R)
# ==============================================================================


klma <- read.csv("data/SIA_Results_KLMA_2026.csv", stringsAsFactors = FALSE)
qaqc <- read.csv("data/SIA_Results_QAQC.csv", stringsAsFactors = FALSE, header = TRUE)

# ---- 1. Standardize KLMA study data ----------------------------------------
names(klma) <- c("Original_ID", "d2H", "d18O", "WtH", "WtO", "OHratio",
                  "Envelope_ID", "Feather_type_raw", "Feather_cat")

klma <- klma %>%
  mutate(
    Envelope_ID = str_trim(Envelope_ID),
    Envelope_ID = ifelse(str_detect(Envelope_ID, regex("san\\s*de[ai]go", ignore_case = TRUE)),
                          "San Diego ref bird", Envelope_ID),
    Feather_type_raw = str_trim(Feather_type_raw),
    # strip trailing replicate index ("Breast 1" -> "Breast", "Wing covert 1"/"2" kept distinct
    # only for within-individual bookkeeping; category classification ignores the index)
    Feather_type_clean = str_trim(str_remove(Feather_type_raw, "\\s*\\d+$")),
    Feather_cat = str_trim(Feather_cat),
    Group = ifelse(Envelope_ID == "San Diego ref bird", "San_Diego_ref", "KLMA")
  )

# Independent re-derivation of feather category from feather type, to check
# the supplied Feather_cat column rather than trust it blindly.
flight_kw <- "rectrix|rectrice|secondary|primary|tertial|covert"
klma <- klma %>%
  mutate(
    Feather_cat_check = ifelse(str_detect(Feather_type_clean, regex(flight_kw, ignore_case = TRUE)),
                                "Flight", "Body")
  )
mismatch <- klma %>% filter(Feather_cat != Feather_cat_check)
if (nrow(mismatch) > 0) {
  message("Feather category mismatches (using supplied Feather_cat): ")
  print(mismatch[, c("Original_ID", "Envelope_ID", "Feather_type_raw", "Feather_cat", "Feather_cat_check")])
}

# ---- 2. Pull QA/QC context (job batch, run date, per-batch precision) ------
# The QAQC file is a stacked lab report (4 job batches). Extract each batch's
# header (job id, run date) and its "standard uncertainty" rows (DS, ORX, POW
# reference materials), then map each Original_ID to its batch by row position.

qraw <- readLines("data/SIA_Results_QAQC.csv")
job_starts <- grep("^SIRFER TCEA-IRMS-HO", qraw)
job_ends <- c(job_starts[-1] - 1, length(qraw))

batch_info <- list()
sample_batch <- list()

for (i in seq_along(job_starts)) {
  chunk <- qraw[job_starts[i]:job_ends[i]]
  job_line <- chunk[grep("^Job #", chunk)]
  job_id <- str_trim(str_split_fixed(job_line, ",", 2)[1])
  job_id <- str_remove(job_id, "^Job #\\s*")
  date_line <- chunk[grep("Date analyzed:", chunk, fixed = TRUE)][1]
  parts_dl <- str_split(date_line, ",")[[1]]
  da_pos <- which(str_detect(parts_dl, "Date analyzed:"))
  run_date <- str_trim(parts_dl[da_pos + 2])

  # standard-uncertainty rows: "standard uncertainty,,<d2H>,<d18O>,..."
  su_idx <- grep("^standard uncertainty", chunk)
  su_vals <- lapply(su_idx, function(j) {
    parts <- str_split(chunk[j], ",")[[1]]
    c(d2H_su = suppressWarnings(as.numeric(parts[3])),
      d18O_su = suppressWarnings(as.numeric(parts[4])))
  })
  # order in file: DS (primary ref), ORX (primary ref 2), POW (secondary/QC)
  su_labels <- c("DS", "ORX", "POW")[seq_along(su_vals)]
  names(su_vals) <- su_labels

  batch_info[[job_id]] <- list(run_date = run_date, su = su_vals)

  # sample rows in this chunk: lines starting with ",<number>,"
  samp_lines <- chunk[str_detect(chunk, "^,[0-9]+,")]
  for (ln in samp_lines) {
    parts <- str_split(ln, ",")[[1]]
    sample_batch[[parts[2]]] <- job_id
  }
}

klma$Job_ID <- sapply(as.character(klma$Original_ID), function(x) {
  v <- sample_batch[[x]]
  if (is.null(v)) NA_character_ else v
})
klma$Run_date <- sapply(klma$Job_ID, function(j) if (is.na(j)) NA_character_ else batch_info[[j]]$run_date)
klma$POW_d2H_su <- sapply(klma$Job_ID, function(j) if (is.na(j)) NA_real_ else batch_info[[j]]$su[["POW"]]["d2H_su"])
klma$POW_d18O_su <- sapply(klma$Job_ID, function(j) if (is.na(j)) NA_real_ else batch_info[[j]]$su[["POW"]]["d18O_su"])

# ---- 2b. Pull in the "1272-31742" reference bird ----------------------------
# This individual is in the original SIRFER QC file (job 24-215.2) but was
# never included in the KLMA study CSV -- its band-number prefix (1272- vs.
# the KLMA birds' 1422-) doesn't match, and per user confirmation it's a
# reference specimen, not part of the core KLMA population. It's kept for
# comparison (same treatment as the "San Diego ref bird" rows).
ref_lines <- qraw[str_detect(qraw, ",1272-31742,")]
ref_df <- do.call(rbind, lapply(ref_lines, function(ln) {
  parts <- str_split(ln, ",")[[1]]
  data.frame(Original_ID = as.integer(parts[2]), d2H = as.numeric(parts[3]), d18O = as.numeric(parts[4]),
             WtH = as.numeric(parts[5]), WtO = as.numeric(parts[6]), OHratio = as.numeric(parts[7]),
             Envelope_ID = parts[8], Feather_type_raw = parts[9], Feather_cat = NA_character_,
             stringsAsFactors = FALSE)
}))
ref_df <- ref_df %>%
  mutate(
    Feather_type_raw = str_trim(Feather_type_raw),
    Feather_type_clean = str_trim(str_remove(Feather_type_raw, "\\s*\\d+$")),
    Feather_cat = ifelse(str_detect(Feather_type_clean, regex(flight_kw, ignore_case = TRUE)),
                          "Flight", "Body"),
    Feather_cat_check = Feather_cat,
    Group = "Reference_1272-31742"
  )
klma <- bind_rows(klma, ref_df)

# job batch / precision lookup for the reference bird's samples too
klma$Job_ID <- sapply(as.character(klma$Original_ID), function(x) {
  v <- sample_batch[[x]]
  if (is.null(v)) NA_character_ else v
})
klma$Run_date <- sapply(klma$Job_ID, function(j) if (is.na(j)) NA_character_ else batch_info[[j]]$run_date)
klma$POW_d2H_su <- sapply(klma$Job_ID, function(j) if (is.na(j)) NA_real_ else batch_info[[j]]$su[["POW"]]["d2H_su"])
klma$POW_d18O_su <- sapply(klma$Job_ID, function(j) if (is.na(j)) NA_real_ else batch_info[[j]]$su[["POW"]]["d18O_su"])

# ---- 3. QA flags -------------------------------------------------------------
# %O:%H ratio for keratin combustion is typically ~3.6-4.5. Flag rows well
# outside that range (possible incomplete combustion / contamination) and
# rows with implausibly high Wt% H or Wt% O (possible double-loading).
klma <- klma %>%
  mutate(
    QC_flag = case_when(
      OHratio > 5.5 | OHratio < 3.0 ~ "flag_OHratio_outlier",
      WtH > 8 | WtO > 30            ~ "flag_wt_pct_outlier",
      TRUE ~ "ok"
    )
  )

if (any(klma$QC_flag != "ok")) {
  message("QC-flagged rows:")
  print(klma %>% filter(QC_flag != "ok") %>%
          select(Original_ID, Envelope_ID, Feather_type_raw, d2H, d18O, WtH, WtO, OHratio, QC_flag))
}

# ---- 4. Exclude the two non-KLMA individuals from the analysis -------------
# "San Diego ref bird" and "1272-31742" were originally carried through this
# pipeline as "reference individuals" for comparison. On review (see
# REPORT.md change log), that framing doesn't hold up: neither was captured
# as part of the KLMA sampling effort, their capture/collection dates are
# unknown relative to the KLMA birds', and -- critically -- we have no
# independent knowledge of where either one actually bred or molted, so
# there is nothing for their isotope-based assignment to be "referenced"
# against. Treating them as reference individuals implied a validation they
# don't provide. Their raw records are preserved in
# data/excluded_individuals.csv for transparency, but they are dropped here
# and do not appear anywhere downstream (assignment, maps, tables, C/N
# cross-check).
excluded <- klma %>% filter(Group != "KLMA")
if (nrow(excluded) > 0) {
  write.csv(excluded %>%
              select(Original_ID, Envelope_ID, Group, Feather_type_raw, Feather_type_clean,
                     Feather_cat, d2H, d18O, WtH, WtO, OHratio, Job_ID, Run_date,
                     POW_d2H_su, POW_d18O_su, QC_flag),
            "data/excluded_individuals.csv", row.names = FALSE)
  message("\nExcluded ", nrow(excluded), " feather samples from ",
          length(unique(excluded$Envelope_ID)),
          " non-KLMA individuals (not true reference individuals -- breeding/molt",
          " origin unknown, different capture context): ",
          paste(unique(excluded$Envelope_ID), collapse = ", "),
          ". Preserved in data/excluded_individuals.csv.")
}
klma <- klma %>% filter(Group == "KLMA")

# ---- 5. Write cleaned dataset -----------------------------------------------
out <- klma %>%
  select(Original_ID, Envelope_ID, Group, Feather_type_raw, Feather_type_clean,
         Feather_cat, d2H, d18O, WtH, WtO, OHratio, Job_ID, Run_date,
         POW_d2H_su, POW_d18O_su, QC_flag)

write.csv(out, "data/klma_clean.csv", row.names = FALSE)
message("\nWrote data/klma_clean.csv with ", nrow(out), " feather samples from ",
        length(unique(out$Envelope_ID)), " KLMA individuals.")
message("Individuals: ", paste(unique(out$Envelope_ID), collapse = ", "))



# ==============================================================================
# STEP 2 (scripts/05_build_habitat_prior.R)
# ==============================================================================


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

# ---- 3b. smooth to remove hard mask edges -----------------------------------
# Without this, the mask is a sparse patchwork of isolated wetland cells
# surrounded by hard zeros. Manual inspection after the first version of
# this prior found that many DIFFERENT birds' independent posteriors were
# landing on the exact same single grid cell despite having different
# measured isotope values -- because whenever a bird's true value wasn't
# achievable at any wetland cell (common, since wetlands are such a small
# fraction of the domain), pdRaster would default to whichever isolated
# wetland cell was the "least bad" available match, and unrelated birds
# with different-but-similarly-off-target values were all funneled onto
# the same isolated cell. A Gaussian-smoothed version (sigma ~2 grid
# cells, ~18 km) gives cells near a wetland partial credit instead of a
# hard cliff to zero, so nearby non-wetland cells with a better isotopic
# match can compete instead of the whole posterior collapsing onto one
# isolated pixel.
gauss_w <- focalMat(habitat_prior, 0.15, type = "Gauss")
habitat_prior <- focal(habitat_prior, w = gauss_w, fun = "sum", na.rm = TRUE)

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



# ==============================================================================
# STEP 3 (scripts/01_run_assignment.R)
# ==============================================================================


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

# Water bird calibration: d2H only, saved here for script 06's calibration
# sensitivity check (Passerine vs. Water bird), not used in the primary model.
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
cal_H_wb$isoscape.rescale <- aoi_crop_mask(cal_H_wb$isoscape.rescale)

habitat_prior <- rast("data/habitat_prior.tif")
habitat_prior <- resample(aoi_crop_mask(habitat_prior), cal_H$isoscape.rescale, method = "near")
crs(habitat_prior) <- crs(cal_H$isoscape.rescale)

# "Weak"/"strong" prior-strength variants for script 06's habitat-prior
# sensitivity check: a power transform preserves the underlying suitability
# pattern (same relative ranking of cells) while compressing (weak,
# exponent < 1, closer to a uniform/no-op prior) or exaggerating (strong,
# exponent > 1, closer to a hard mask) how much it discounts non-habitat.
habitat_prior_weak <- habitat_prior^0.5
habitat_prior_strong <- habitat_prior^2

# ---- 4. Load cleaned sample data --------------------------------------------
samp <- read.csv("data/klma_clean.csv", stringsAsFactors = FALSE)
samp$sample_id <- paste(samp$Envelope_ID, samp$Feather_type_clean, samp$Original_ID, sep = "__")

unknown_H <- data.frame(ID = samp$sample_id, d2H = samp$d2H)

message("Running primary per-feather assignment (d2H only, Passerine calibration, ",
        "current habitat prior) for ", nrow(unknown_H), " samples...")
# genplot = FALSE: pdRaster's own internal plot is a raw, uncurated
# byproduct (terra caps it at 16 panels, with assignR's default color
# scheme unrelated to this project's maps) -- not a real deliverable.
# scripts/03_make_maps.R builds the actual per-bird figures.
pd_feather <- pdRaster(cal_H, unknown = unknown_H, mask = NULL, prior = habitat_prior, genplot = FALSE)
names(pd_feather) <- unknown_H$ID
writeRaster(pd_feather, "outputs/rasters/pd_per_feather.tif", overwrite = TRUE)

# Saved for same-session inspection only (e.g. interactively in RStudio).
# terra::SpatRaster objects wrap external C++ pointers that do not survive
# saveRDS()/readRDS() across separate R processes -- they come back with
# dead pointers and error on first use -- so scripts/06_sensitivity_check.R
# recomputes these objects itself rather than reading this file (it runs as
# its own `Rscript` process, like every script in this pipeline).
saveRDS(list(cal_H = cal_H, cal_O = cal_O, cal_H_wb = cal_H_wb,
             habitat_prior = habitat_prior,
             habitat_prior_weak = habitat_prior_weak,
             habitat_prior_strong = habitat_prior_strong,
             aoi = aoi_sf, zones_names = names(zones)),
        "outputs/calibration_objects.rds")

message("Step 1-4 complete. Saved primary per-feather posterior stack (d2H-only).")



# ==============================================================================
# STEP 4 (scripts/02_summarize_assignment.R)
# ==============================================================================


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



# ==============================================================================
# STEP 5 (scripts/03_make_maps.R)
# ==============================================================================

sf_use_s2(FALSE)

pop <- rast("outputs/rasters/pd_population_by_category.tif")
indiv <- rast("outputs/rasters/pd_per_individual_category.tif")
qtl50 <- rast("outputs/rasters/qtl50_credible_region.tif")
qtl75 <- rast("outputs/rasters/qtl75_credible_region.tif")
qtl90 <- rast("outputs/rasters/qtl90_credible_region.tif")
qual <- read.csv("outputs/tables/assignment_summary.csv", stringsAsFactors = FALSE)
flagged <- with(qual, surface[single_cell_dominant_flag | shared_peak_flag])

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
plot_one("Population_Body", "KLMA Yellow Rail -- Body feathers\n(population mean origin probability, d2H-only primary model)")
plot_one("Population_Flight", "KLMA Yellow Rail -- Flight feathers\n(population mean origin probability, d2H-only primary model)")
dev.off()

# ---- individual grids: split into main (clean) figure + supplementary -----
# (flagged) figure. type = "interval" (not "classes") for the multi-layer
# case: terra's "classes" plotting auto-detects factor levels per layer and
# can assign inconsistent colors/legends across panels when different birds
# don't have all 4 classes present. "interval" with explicit breaks avoids
# that entirely and gives every panel the same fixed color mapping.
make_grid <- function(idx, outfile, w = 1900, h = 1500, mark_flagged = FALSE) {
  nms <- names(indiv)[idx]
  if (length(nms) == 0) return(invisible(NULL))
  cls_list <- lapply(nms, cr_class)
  cls_stack <- rast(cls_list)
  titles <- if (mark_flagged) ifelse(nms %in% flagged, paste0(nms, "  [low match quality]"), nms) else nms
  names(cls_stack) <- titles
  nc <- min(4, length(nms))
  nr <- ceiling(length(nms) / nc)
  png(outfile, width = w, height = h * nr / 3, res = 140)
  plot(cls_stack, breaks = c(0.5, 1.5, 2.5, 3.5, 4.5), col = cr_colors, type = "interval",
       nc = nc, plg = list(cex = 0.5, legend = cr_labels, title = "credible region"),
       mar = c(1.5, 1.5, 2, 7))
  dev.off()
}

nm <- names(indiv)
body_idx <- grep("Body$", nm); flight_idx <- grep("Flight$", nm)
body_main <- body_idx[!nm[body_idx] %in% flagged]; body_flag <- body_idx[nm[body_idx] %in% flagged]
flight_main <- flight_idx[!nm[flight_idx] %in% flagged]; flight_flag <- flight_idx[nm[flight_idx] %in% flagged]

make_grid(body_main, "outputs/maps/individuals_body.png")
make_grid(flight_main, "outputs/maps/individuals_flight.png")
make_grid(body_flag, "outputs/maps/individuals_body_supplementary_flagged.png")
make_grid(flight_flag, "outputs/maps/individuals_flight_supplementary_flagged.png")

message("Body: ", length(body_main), " main + ", length(body_flag), " flagged; ",
        "Flight: ", length(flight_main), " main + ", length(flight_flag), " flagged")

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

# ---- publication-ready summary table ----------------------------------------
# All rows are KLMA individuals -- the two non-KLMA individuals formerly
# shown here as "reference" birds were removed from the analysis entirely
# (see scripts/00_clean_data.R), so no Origin_type/reference column is
# needed any more.
tbl <- qual[order(match(qual$Group, c("KLMA", "KLMA_summary")), qual$Envelope_ID, qual$Feather_cat), ]
tbl$Bird <- ifelse(tbl$Group == "KLMA_summary", "All KLMA (population mean)", tbl$Envelope_ID)
tbl$Quality_note <- ifelse(tbl$single_cell_dominant_flag | tbl$shared_peak_flag,
                            "Low match quality -- see Methods", "")
pub_tbl <- data.frame(
  Bird = tbl$Bird,
  Feathers = tbl$Feather_cat,
  n_feathers = tbl$n_feathers,
  `Peak_lon` = tbl$peak_lon,
  `Peak_lat` = tbl$peak_lat,
  `50pct_credible_km2` = tbl$area50_km2,
  `90pct_credible_km2` = tbl$area90_km2,
  `P_Mexico_pct` = round(100 * (tbl$P_Baja_CA + tbl$P_Mainland_W_Mexico), 1),
  `P_PacificNW_pct` = round(100 * (tbl$P_Oregon + tbl$P_Washington_Idaho_Utah + tbl$P_California), 1),
  Quality_note = tbl$Quality_note,
  check.names = FALSE
)
write.csv(pub_tbl, "outputs/tables/Table1_manuscript.csv", row.names = FALSE)

message("Maps written to outputs/maps/")
message("Publication table written to outputs/tables/Table1_manuscript.csv")



# ==============================================================================
# STEP 6 (scripts/04_clean_cn_data.R)
# ==============================================================================


cn <- read.csv("data/SIBS_CN_230621_Kellerman.csv", stringsAsFactors = FALSE)
names(cn)[1:7] <- c("SIBS_ID", "Sample_ID_raw", "Feather_code", "d15N", "N_pct", "d13C", "C_pct")
cn <- cn %>% filter(!is.na(SIBS_ID) & SIBS_ID != "" & SIBS_ID != "#REF!")
cn <- cn[, 1:7]

# parse "1422-02216 Breast 1" -> Envelope_ID "1422-02216", Feather_type_raw "Breast 1"
cn <- cn %>%
  mutate(
    Sample_ID_raw = str_trim(Sample_ID_raw),
    Envelope_ID = str_extract(Sample_ID_raw, "^\\S+"),
    Feather_type_raw = str_trim(str_remove(Sample_ID_raw, "^\\S+\\s*")),
    Feather_type_clean = str_trim(str_remove(Feather_type_raw, "\\s*(partial|\\d+)$")),
    Feather_cat = case_when(
      Feather_code == "B" ~ "Body",
      Feather_code == "N" ~ "Body",
      Feather_code == "R" ~ "Flight",
      TRUE ~ NA_character_
    ),
    # band-number prefix 1272- does not match the KLMA 1422- prefix pattern
    # used elsewhere in this study
    Group = ifelse(str_starts(Envelope_ID, "1422-"), "KLMA", "Non_KLMA_excluded")
  )

# Exclude individual 1272-31742, same as scripts/00_clean_data.R: not a
# true reference individual (unknown breeding/molt origin, different
# capture context from the KLMA birds), so dropped from the analysis
# entirely rather than carried through as a comparison point.
n_excluded <- sum(cn$Group == "Non_KLMA_excluded")
if (n_excluded > 0) {
  message("Excluding ", n_excluded, " C/N sample(s) from non-KLMA individual(s): ",
          paste(unique(cn$Envelope_ID[cn$Group == "Non_KLMA_excluded"]), collapse = ", "))
}
cn <- cn %>% filter(Group == "KLMA")

write.csv(cn, "data/cn_clean.csv", row.names = FALSE)
message("Wrote data/cn_clean.csv: ", nrow(cn), " C/N samples from ",
        length(unique(cn$Envelope_ID)), " individuals: ",
        paste(unique(cn$Envelope_ID), collapse = ", "))

# ---- per-individual x feather-category means -------------------------------
cn_summary <- cn %>%
  group_by(Envelope_ID, Group, Feather_cat) %>%
  summarise(
    n = n(),
    d13C_mean = mean(d13C), d13C_sd = sd(d13C),
    d15N_mean = mean(d15N), d15N_sd = sd(d15N),
    .groups = "drop"
  )
write.csv(cn_summary, "data/cn_summary_by_individual.csv", row.names = FALSE)
print(cn_summary)

# ---- cross-reference against H/O results (data/klma_clean.csv) -------------
ho <- read.csv("data/klma_clean.csv", stringsAsFactors = FALSE)
ho_summary <- ho %>%
  group_by(Envelope_ID, Feather_cat) %>%
  summarise(d2H_mean = mean(d2H), d18O_mean = mean(d18O), .groups = "drop")

combined <- cn_summary %>%
  left_join(ho_summary, by = c("Envelope_ID", "Feather_cat"))
write.csv(combined, "data/isotopes_combined_summary.csv", row.names = FALSE)

message("\n--- Combined C/N + H/O summary ---")
print(as.data.frame(combined))

message("\nIndividuals in C/N data with NO matching H/O record for that feather category:")
print(combined %>% filter(is.na(d2H_mean)) %>% select(Envelope_ID, Feather_cat, n))

# ---- cross-reference against the geographic assignment (script 02 output) --
# Requires outputs/tables/assignment_summary.csv, which only exists after
# 02_summarize_assignment.R has run (see run-order note at top of file).
# This is the table used to sanity-check the H/O-based spatial assignment
# against the independent C/N diet/habitat signal: e.g. a bird assigned
# high P_Mainland_W_Mexico/P_Baja_CA probability should plausibly show a
# more C4/marine-enriched d13C and/or higher d15N than a bird assigned to
# the Pacific Northwest, if the geographic assignment is picking up a real
# signal rather than isoscape noise.
assign_path <- "outputs/tables/assignment_summary.csv"
if (file.exists(assign_path)) {
  assign_tbl <- read.csv(assign_path, stringsAsFactors = FALSE)

  cn_vs_geo <- cn_summary %>%
    inner_join(assign_tbl, by = c("Envelope_ID", "Feather_cat")) %>%
    transmute(
      Envelope_ID, Group = Group.x, Feather_cat,
      d13C_mean, d13C_sd, d15N_mean, d15N_sd,
      P_Mexico_total = round(P_Baja_CA + P_Mainland_W_Mexico, 4),
      P_PacificNW = round(P_Oregon + P_Washington_Idaho_Utah + P_California, 4),
      P_within_100km_Klamath,
      peak_lon, peak_lat, area50_km2, peak_share,
      single_cell_dominant_flag, shared_peak_flag
    )

  write.csv(cn_vs_geo, "data/cn_vs_geography_comparison.csv", row.names = FALSE)
  message("\nWrote data/cn_vs_geography_comparison.csv: ", nrow(cn_vs_geo),
          " C/N samples matched to a geographic assignment record")
  print(cn_vs_geo)
} else {
  message("\n", assign_path, " not found -- skipping C/N-vs-geography comparison.",
          " Run 01_run_assignment.R and 02_summarize_assignment.R first, then ",
          "re-run this script to produce data/cn_vs_geography_comparison.csv.")
}



# ==============================================================================
# STEP 7 (scripts/06_sensitivity_check.R)
# ==============================================================================

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
cal_O <- calRaster(known_O, d18o_iso, genplot = FALSE, verboseLM = FALSE)
cal_H_wb <- calRaster(known_H_wb, d2h_iso, genplot = FALSE, verboseLM = FALSE)

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


