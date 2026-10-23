# 04_clean_cn_data.R
# Clean the SIBS carbon/nitrogen isotope results (separate lab, separate run
# from the SIRFER H/O data) and cross-reference against the H/O feather
# records by envelope ID + feather type.
#
# NOTE ON RUN ORDER: this script's first two outputs (cn_clean.csv,
# cn_summary_by_individual.csv, and the raw-d2H/d18O version of
# isotopes_combined_summary.csv) only need data/klma_clean.csv (from
# 00_clean_data.R). But the full C/N-vs-geographic-assignment comparison
# at the end of this script (cn_vs_geography_comparison.csv) needs
# outputs/tables/assignment_summary.csv, which doesn't exist until
# 02_summarize_assignment.R has run. So although this file is numbered
# "04", run it LAST in the full pipeline: 00 -> 05 -> 01 -> 02 -> 03 -> 04.
# If assignment_summary.csv isn't present yet, the last block is skipped
# with a message rather than erroring, so this script still works standalone
# right after 00 for early exploration of the C/N data on its own.

suppressMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})

cn <- read.csv("data/SIBS_CN_230621_Kellerman.csv", stringsAsFactors = FALSE)
names(cn)[1:7] <- c("SIBS_ID", "Sample_ID_raw", "Feather_code", "d15N", "N_pct", "d13C", "C_pct")
cn <- cn %>% filter(!is.na(SIBS_ID) & SIBS_ID != "" & SIBS_ID != "#REF!")
cn <- cn[, 1:7]

# ---- 0. d13C bias correction, from the SIBS lab's own QA standards --------
# data/SIBS_CN_230621_Kellerman_QC.xlsx reports two in-house plant QA
# standards (SIBS-pCo, SIBS-pEc) run alongside the unknown samples, each
# with a known consensus value. After the lab's own isotope-scale
# normalization (via USGS Glutamic Acid 40/41, which check out essentially
# exactly against their known values -- no correction needed there), both
# QA standards still read d13C noticeably enriched (less negative) relative
# to their consensus values -- and strikingly consistently so, despite the
# two standards sitting ~27 permil apart on the d13C scale:
#   SIBS-pCo: measured -10.603 permil vs. consensus -12.7  -> +2.097 permil
#   SIBS-pEc: measured -37.890 permil vs. consensus -39.9  -> +2.010 permil
# That consistency (bias differs by only 0.09 permil across a 27-permil
# range) means a simple additive correction is justified -- a 2-point
# linear calibration through these standards has a slope of 0.997,
# indistinguishable from 1 given the measurement precision, so there's no
# evidence the bias depends on the isotope value itself. d15N shows only a
# small, inconsistent offset between the two standards (+0.31, +0.15
# permil) -- within normal measurement variation, not corrected.
qc_path <- "data/SIBS_CN_230621_Kellerman_QC.xlsx"
qc_raw <- suppressMessages(read_excel(qc_path, sheet = 1, col_names = FALSE))
qc_raw <- as.data.frame(lapply(qc_raw, as.character), stringsAsFactors = FALSE)
qc_val <- function(label_row_text, value_row_text, col) {
  # Each QA standard is a 3-row block: the label row itself holds "Average"
  # (col3) in the same row as the standard's name (col2); "SD (1s)" is the
  # next row; "Consensus value" is two rows below the label. Search the
  # label row and the following 2 rows for the requested col3 value.
  i <- which(qc_raw[[2]] == label_row_text)
  block <- qc_raw[i:(i + 2), ]
  j <- i - 1 + which(block[[3]] == value_row_text)[1]
  as.numeric(qc_raw[[col]][j])
}
pco_bias_d13C <- qc_val("SIBS-pCo", "Average", 5) - qc_val("SIBS-pCo", "Consensus value", 5)
pec_bias_d13C <- qc_val("SIBS-pEc", "Average", 5) - qc_val("SIBS-pEc", "Consensus value", 5)
D13C_CORRECTION <- mean(c(pco_bias_d13C, pec_bias_d13C))
message(sprintf("d13C QA-standard bias: SIBS-pCo = %+.3f permil, SIBS-pEc = %+.3f permil -- applying a %.3f permil correction to all samples (d13C_corrected = d13C_raw - %.3f).",
                pco_bias_d13C, pec_bias_d13C, D13C_CORRECTION, D13C_CORRECTION))

cn$d13C_raw <- cn$d13C
cn$d13C <- cn$d13C_raw - D13C_CORRECTION

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
