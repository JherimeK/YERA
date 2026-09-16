# 04_clean_cn_data.R
# Clean the SIBS carbon/nitrogen isotope results (separate lab, separate run
# from the SIRFER H/O data) and cross-reference against the H/O feather
# records by envelope ID + feather type.

suppressMessages({
  library(dplyr)
  library(stringr)
})

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
    # used elsewhere in this study -- flag rather than assume
    Group = ifelse(str_starts(Envelope_ID, "1422-"), "KLMA", "Non_KLMA_flagged")
  )

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
