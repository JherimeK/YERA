# 00_clean_data.R
# Clean the KLMA Yellow Rail feather isotope results, read directly from
# the single consolidated lab spreadsheet (data/SIA_Results_KLMA.xlsx),
# which includes both the per-feather results and the SIRFER QA/QC context
# (job batches, reference-standard runs, run dates, per-batch precision)
# in one file. Output is a single tidy CSV used by the assignR pipeline.
#
# This replaced an earlier two-file version of this script (a separately
# exported, decimal-rounded "clean" results spreadsheet cross-referenced
# against a raw QC text export) after the lead author pointed out the
# original two-file split was error-prone -- it was how a non-KLMA
# individual (1272-31742) ended up silently included in this pipeline for
# a while, recovered only by text-regex from the raw QC file rather than
# being a deliberate inclusion. The consolidated file removes that
# failure mode entirely (there's only one file, and it's already scoped
# to just the 12 KLMA individuals), and as a bonus gives full floating-
# point isotope precision instead of the previous file's values rounded
# to ~1 decimal place.

suppressMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})

path <- "data/SIA_Results_KLMA.xlsx"
raw <- suppressMessages(read_excel(path, sheet = 1, col_names = FALSE))
raw <- as.data.frame(lapply(raw, as.character), stringsAsFactors = FALSE)
names(raw) <- paste0("c", seq_len(ncol(raw)))

# ---- 1. Split into per-job blocks -------------------------------------------
# The sheet is a stacked lab report: one block per SIRFER job batch, each
# with a header (job ID, analyst, run date), a table of per-feather sample
# rows, then a "Quality Assurance"/"Quality Control" section with replicate
# runs of the DS/ORX (primary) and POW (secondary) reference standards,
# their averages, and standard uncertainties.
job_starts <- which(raw$c1 == "SIRFER TCEA-IRMS-HO")
job_ends <- c(job_starts[-1] - 1, nrow(raw))

sample_rows <- list()
batch_info <- list()

for (i in seq_along(job_starts)) {
  chunk <- raw[job_starts[i]:job_ends[i], ]

  job_line <- chunk$c1[str_detect(chunk$c1, "^Job #")]
  job_id <- str_trim(str_remove(job_line[1], "^Job #\\s*"))

  date_row <- chunk[str_detect(chunk$c4, "Date analyzed:") & !is.na(chunk$c4), ]
  run_date <- str_trim(date_row$c6[1])

  # standard-uncertainty rows (one per reference material: DS, ORX, POW, in
  # that order) have c1 == "standard uncertainty" and the d2H/d18O SU in
  # columns 3/4.
  su_rows <- chunk[!is.na(chunk$c1) & chunk$c1 == "standard uncertainty", ]
  su_labels <- c("DS", "ORX", "POW")[seq_len(nrow(su_rows))]
  su <- setNames(
    lapply(seq_len(nrow(su_rows)), function(j) {
      c(d2H_su = as.numeric(su_rows$c3[j]), d18O_su = as.numeric(su_rows$c4[j]))
    }),
    su_labels
  )
  batch_info[[job_id]] <- list(run_date = run_date, su = su)

  # sample rows: c1 (SIRFER #) is blank and c2 (Original ID) is numeric
  samp <- chunk[is.na(chunk$c1) & !is.na(chunk$c2) &
                  suppressWarnings(!is.na(as.numeric(chunk$c2))), ]
  if (nrow(samp) > 0) {
    samp$Job_ID <- job_id
    sample_rows[[job_id]] <- samp
  }
}

klma <- bind_rows(sample_rows)
klma <- klma %>%
  transmute(
    Original_ID = as.integer(c2),
    d2H = as.numeric(c3),
    d18O = as.numeric(c4),
    WtH = as.numeric(c5),
    WtO = as.numeric(c6),
    OHratio = as.numeric(c7),
    Envelope_ID = str_trim(c8),
    Feather_type_raw = str_trim(c9),
    Job_ID = Job_ID
  )

# ---- 2. Attach per-batch run date and POW reference-standard precision -----
klma$Run_date <- sapply(klma$Job_ID, function(j) batch_info[[j]]$run_date)
klma$POW_d2H_su <- sapply(klma$Job_ID, function(j) batch_info[[j]]$su[["POW"]]["d2H_su"])
klma$POW_d18O_su <- sapply(klma$Job_ID, function(j) batch_info[[j]]$su[["POW"]]["d18O_su"])

# ---- 3. Feather category -----------------------------------------------------
# strip trailing replicate index ("Breast 1" -> "Breast", "Wing covert 1"/"2"
# kept distinct only for within-individual bookkeeping; category
# classification ignores the index) and classify Body vs. Flight from the
# feather-type keyword.
flight_kw <- "rectrix|rectrice|secondary|primary|tertial|covert"
klma <- klma %>%
  mutate(
    Feather_type_clean = str_trim(str_remove(Feather_type_raw, "\\s*(partial|\\d+)$")),
    Feather_cat = ifelse(str_detect(Feather_type_clean, regex(flight_kw, ignore_case = TRUE)),
                          "Flight", "Body"),
    Group = "KLMA"
  )

# ---- 4. Sanity check: every individual should be a KLMA-banded bird -------
# The whole point of switching to this consolidated source file was to stop
# relying on post-hoc filtering to keep non-KLMA individuals out (see file
# header). This is a belt-and-suspenders check, not the primary mechanism.
non_klma <- klma %>% filter(!str_starts(Envelope_ID, "1422-"))
if (nrow(non_klma) > 0) {
  stop("Found ", nrow(non_klma), " sample(s) with a non-KLMA envelope ID (",
       paste(unique(non_klma$Envelope_ID), collapse = ", "), ") in ", path,
       " -- this source file was expected to already be scoped to KLMA",
       " individuals only. Check the file before proceeding.")
}

# ---- 5. QA flags -------------------------------------------------------------
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

# ---- 6. Write cleaned dataset -----------------------------------------------
out <- klma %>%
  select(Original_ID, Envelope_ID, Group, Feather_type_raw, Feather_type_clean,
         Feather_cat, d2H, d18O, WtH, WtO, OHratio, Job_ID, Run_date,
         POW_d2H_su, POW_d18O_su, QC_flag) %>%
  arrange(Original_ID)

write.csv(out, "data/klma_clean.csv", row.names = FALSE)
message("\nWrote data/klma_clean.csv with ", nrow(out), " feather samples from ",
        length(unique(out$Envelope_ID)), " KLMA individuals.")
message("Individuals: ", paste(unique(out$Envelope_ID), collapse = ", "))
message("POW standard uncertainty (d2H) by batch: ",
        paste(sprintf("%s=%.2f permil", names(batch_info),
                       sapply(batch_info, function(b) b$su[["POW"]]["d2H_su"])),
              collapse = "; "))
