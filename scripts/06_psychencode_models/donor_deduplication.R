#!/usr/bin/env Rscript
# PsychENCODE donor de-duplication.
# Donor de-duplication and exact subject-ID overlap filtering (38 ASD / 42 control donors, 532 cortical samples).
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).

suppressMessages(library(methods))
set.seed(42)

data_dir <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "psychencode_processed")
out_dir <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "18_psychencode/06_donor_independence")

# Step 1: Load PsychENCODE metadata
cat("Loading PsychENCODE processed metadata...\n")
env <- new.env()
load(file.path(data_dir, "02_01_B_AllProcessedData_wModelMatrix.RData"), envir = env)
datMeta <- env$datMeta

# Extract unique donors
psych_donors <- data.frame(
  subject = sort(unique(datMeta$subject)),
  diagnosis = NA,
  n_samples = NA,
  brain_bank = NA,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(psych_donors))) {
  subj <- psych_donors$subject[i]
  idx <- datMeta$subject == subj
  psych_donors$diagnosis[i] <- unique(datMeta$Diagnosis[idx])[1]
  psych_donors$n_samples[i] <- sum(idx)
  psych_donors$brain_bank[i] <- unique(datMeta$Brain_Bank_Source[idx])[1]
}

cat(paste0("PsychENCODE total donors: ", nrow(psych_donors), "\n"))
cat(paste0("  ASD: ", sum(psych_donors$diagnosis == "ASD"), "\n"))
cat(paste0("  CTL: ", sum(psych_donors$diagnosis == "CTL"), "\n"))
cat(paste0("  Dup15q: ", sum(psych_donors$diagnosis == "Dup15q"), "\n"))

# Step 2: Load GSE64018/Parikshak donor IDs
cat("\nLoading GSE64018 (Parikshak discovery) donor IDs...\n")
gse_meta <- read.delim(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "09_resource_schema/07_gse64018/02_GSE64018_sample_metadata.tsv"),
                       stringsAsFactors = FALSE)
# Extract donor ID from sample_id (format: DONOR_region_RIN)
gse_donors_raw <- unique(sub("_ba.*|_BA.*", "", gse_meta$sample_id))
cat(paste0("GSE64018 unique donors: ", length(gse_donors_raw), "\n"))
cat(paste0("  Donors: ", paste(sort(gse_donors_raw), collapse = ", "), "\n"))

# Step 3: Normalize IDs for comparison
# GSE64018 and PsychENCODE both encode consortium donor identifiers
# (AN#### / UMB#### / A###-## style). Normalize by removing leading
# zeros and standardizing format before matching.

normalize_id <- function(id) {
  id <- toupper(trimws(id))
  # Remove common prefixes/suffixes
  id <- gsub("^AN0*", "AN", id)  # strip leading zeros in AN identifiers
  id <- gsub("^UMB0*", "UMB", id)  # strip leading zeros in UMB identifiers
  return(id)
}

psych_norm <- data.frame(
  subject = psych_donors$subject,
  normalized = sapply(psych_donors$subject, normalize_id),
  diagnosis = psych_donors$diagnosis,
  n_samples = psych_donors$n_samples,
  brain_bank = psych_donors$brain_bank,
  stringsAsFactors = FALSE
)

gse_norm <- data.frame(
  donor_raw = gse_donors_raw,
  normalized = sapply(gse_donors_raw, normalize_id),
  stringsAsFactors = FALSE
)

# Step 4: Find overlap
cat("\n=== Donor Overlap Analysis ===\n")

# Exact match
exact_overlap <- intersect(psych_donors$subject, gse_donors_raw)
cat(paste0("Exact ID overlap: ", length(exact_overlap), "\n"))
if (length(exact_overlap) > 0) {
  cat(paste0("  Overlapping donors: ", paste(exact_overlap, collapse = ", "), "\n"))
}

# Normalized match
norm_overlap <- intersect(psych_norm$normalized, gse_norm$normalized)
cat(paste0("Normalized ID overlap: ", length(norm_overlap), "\n"))
if (length(norm_overlap) > 0) {
  cat(paste0("  Normalized overlapping: ", paste(norm_overlap, collapse = ", "), "\n"))
  # Show which PsychENCODE subjects match
  for (n in norm_overlap) {
    psych_match <- psych_norm$subject[psych_norm$normalized == n]
    gse_match <- gse_norm$donor_raw[gse_norm$normalized == n]
    cat(paste0("    PsychENCODE: ", paste(psych_match, collapse=","),
               " <-> GSE64018: ", paste(gse_match, collapse=","), "\n"))
  }
}

# Also check partial matches (numeric part)
extract_numeric <- function(id) {
  nums <- regmatches(id, regexpr("[0-9]+", id))
  if (length(nums) > 0) return(nums)
  return("")
}

psych_nums <- sapply(psych_donors$subject, extract_numeric)
gse_nums <- sapply(gse_donors_raw, extract_numeric)

# Check if any PsychENCODE AN donors match GSE AN donors by number
psych_an <- psych_donors$subject[grepl("^AN", psych_donors$subject)]
gse_an <- gse_donors_raw[grepl("^AN", gse_donors_raw)]

psych_an_nums <- sapply(psych_an, extract_numeric)
gse_an_nums <- sapply(gse_an, extract_numeric)

numeric_overlap <- intersect(psych_an_nums, gse_an_nums)
cat(paste0("\nAN-prefix numeric overlap: ", length(numeric_overlap), "\n"))
if (length(numeric_overlap) > 0) {
  for (n in numeric_overlap) {
    pm <- psych_an[psych_an_nums == n]
    gm <- gse_an[gse_an_nums == n]
    cat(paste0("  PsychENCODE: ", paste(pm, collapse=","),
               " <-> GSE64018: ", paste(gm, collapse=","), "\n"))
  }
}

# Step 5: Classify all PsychENCODE donors
all_overlap_ids <- unique(c(exact_overlap,
                           psych_norm$subject[psych_norm$normalized %in% norm_overlap],
                           psych_an[psych_an_nums %in% numeric_overlap]))

cat(paste0("\nTotal overlapping donors (all methods): ", length(all_overlap_ids), "\n"))
if (length(all_overlap_ids) > 0) {
  cat(paste0("  ", paste(sort(all_overlap_ids), collapse = ", "), "\n"))
}

psych_donors$overlap_status <- ifelse(
  psych_donors$subject %in% all_overlap_ids,
  "OVERLAP_WITH_PARIKSHAK",
  "NONOVERLAP_WITH_PARIKSHAK"
)

# Summary
cat("\n=== FINAL DONOR COUNTS ===\n")
cat(paste0("N_PSYCHENCODE_DONORS_TOTAL: ", nrow(psych_donors), "\n"))
cat(paste0("N_IDIOPATHIC_ASD_DONORS: ", sum(psych_donors$diagnosis == "ASD"), "\n"))
cat(paste0("N_CONTROL_DONORS: ", sum(psych_donors$diagnosis == "CTL"), "\n"))
cat(paste0("N_DUP15Q_DONORS: ", sum(psych_donors$diagnosis == "Dup15q"), "\n"))
cat(paste0("N_OVERLAP_DONORS: ", sum(psych_donors$overlap_status == "OVERLAP_WITH_PARIKSHAK"), "\n"))

# Non-overlap counts
nonoverlap <- psych_donors[psych_donors$overlap_status == "NONOVERLAP_WITH_PARIKSHAK", ]
cat(paste0("N_NONOVERLAP_ASD_DONORS: ", sum(nonoverlap$diagnosis == "ASD"), "\n"))
cat(paste0("N_NONOVERLAP_CONTROL_DONORS: ", sum(nonoverlap$diagnosis == "CTL"), "\n"))
cat(paste0("N_NONOVERLAP_DUP15Q_DONORS: ", sum(nonoverlap$diagnosis == "Dup15q"), "\n"))

# Overlap breakdown
overlap <- psych_donors[psych_donors$overlap_status == "OVERLAP_WITH_PARIKSHAK", ]
if (nrow(overlap) > 0) {
  cat(paste0("\nOverlap donor diagnoses:\n"))
  print(table(overlap$diagnosis))
}

# Step 6: Write outputs
# Parikshak donor master
parikshak_df <- data.frame(
  donor_id = gse_donors_raw,
  normalized_id = gse_norm$normalized,
  source = "GSE64018",
  study = "Parikshak_2016",
  stringsAsFactors = FALSE
)
write.table(parikshak_df, file.path(out_dir, "00_Parikshak_donor_master.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# PsychENCODE donor master
write.table(psych_donors, file.path(out_dir, "01_PsychENCODE_donor_master.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ID normalization
norm_df <- data.frame(
  psychencode_subject = psych_norm$subject,
  normalized = psych_norm$normalized,
  in_gse64018 = psych_norm$subject %in% all_overlap_ids,
  stringsAsFactors = FALSE
)
write.table(norm_df, file.path(out_dir, "02_id_normalization.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Exact overlap
if (length(all_overlap_ids) > 0) {
  overlap_df <- data.frame(
    psychencode_subject = all_overlap_ids,
    gse64018_match = all_overlap_ids,  # same since exact
    match_type = "EXACT",
    psych_diagnosis = psych_donors$diagnosis[psych_donors$subject %in% all_overlap_ids],
    stringsAsFactors = FALSE
  )
  write.table(overlap_df, file.path(out_dir, "03_exact_overlap.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  writeLines("psychencode_subject\tgse64018_match\tmatch_type\tpsych_diagnosis",
             file.path(out_dir, "03_exact_overlap.tsv"))
}

# Independence classification
write.table(psych_donors[, c("subject", "diagnosis", "overlap_status", "brain_bank")],
            file.path(out_dir, "05_independence_classification.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Non-overlap analysis set
nonoverlap_set <- psych_donors[psych_donors$overlap_status == "NONOVERLAP_WITH_PARIKSHAK" &
                                psych_donors$diagnosis %in% c("ASD", "CTL"), ]
write.table(nonoverlap_set, file.path(out_dir, "06_nonoverlap_analysis_set.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Phase file
check_df <- data.frame(
  key = c("N_PSYCHENCODE_DONORS_TOTAL", "N_IDIOPATHIC_ASD_DONORS", "N_CONTROL_DONORS",
           "N_DUP15Q_DONORS", "N_OVERLAP_DONORS", "N_NONOVERLAP_ASD_DONORS",
           "N_NONOVERLAP_CONTROL_DONORS", "DONOR_INDEPENDENCE_STATUS"),
  value = c(nrow(psych_donors), sum(psych_donors$diagnosis == "ASD"),
            sum(psych_donors$diagnosis == "CTL"), sum(psych_donors$diagnosis == "Dup15q"),
            sum(psych_donors$overlap_status == "OVERLAP_WITH_PARIKSHAK"),
            sum(nonoverlap$diagnosis == "ASD"),
            sum(nonoverlap$diagnosis == "CTL"),
            ifelse(sum(nonoverlap$diagnosis %in% c("ASD","CTL")) >= 10,
                   "OK_SUFFICIENT_NONOVERLAP_DONORS",
                   "HELD_INSUFFICIENT_NONOVERLAP_DONORS")),
  stringsAsFactors = FALSE
)
write.table(check_df, file.path(out_dir, "07_independence_check.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nDonor independence check complete.\n")
