# 00_clean_data.R
# Clean and merge the KLMA Yellow Rail feather isotope results with the
# SIRFER lab QA/QC file, cross-referenced by the shared "Original ID"
# (SIRFER sample number). Output is a single tidy CSV used by the
# assignR pipeline.

suppressMessages({
  library(dplyr)
  library(stringr)
})

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

# ---- 4. Write cleaned dataset -----------------------------------------------
out <- klma %>%
  select(Original_ID, Envelope_ID, Group, Feather_type_raw, Feather_type_clean,
         Feather_cat, d2H, d18O, WtH, WtO, OHratio, Job_ID, Run_date,
         POW_d2H_su, POW_d18O_su, QC_flag)

write.csv(out, "data/klma_clean.csv", row.names = FALSE)
message("\nWrote data/klma_clean.csv with ", nrow(out), " feather samples from ",
        length(unique(out$Envelope_ID)), " individuals.")
message("Individuals: ", paste(unique(out$Envelope_ID), collapse = ", "))
