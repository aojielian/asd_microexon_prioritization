#!/usr/bin/env Rscript
# Coordinate inference and mixed-model inference.
# Reciprocal liftOver verification, GENCODE v33 local-structure checking, and the PsychENCODE mixed models (M0-M4; Kenward-Roger primary, ML-LRT sensitivity).
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).
suppressMessages({library(methods); library(lme4); library(pbkrtest)})
set.seed(42)

project <- Sys.getenv("PROJECT_ROOT", unset = ".")
task_dir <- file.path(project, "21_coordinate_inference")
data_dir <- file.path(project, "psychencode_processed")
prev_dir <- file.path(project, "19_psychencode_tpm_usage")
prev18_dir <- file.path(project, "18_psychencode")
liftover_bin <- file.path(Sys.getenv("SCRATCH_ROOT", unset = tempdir()), "liftOver")
chain_hg19_hg38 <- file.path(Sys.getenv("SCRATCH_ROOT", unset = tempdir()), "hg19ToHg38.over.chain.gz")
chain_hg38_hg19 <- file.path(Sys.getenv("SCRATCH_ROOT", unset = tempdir()), "hg38ToHg19.over.chain.gz")
tmp <- Sys.getenv("SCRATCH_ROOT", unset = tempdir())

# Load the 19 events
primary19 <- read.delim(file.path(project, "16_gse30573/02_input_lock/02_primary19.tsv"), stringsAsFactors=FALSE)
primary19$discovery_dir <- ifelse(primary19$Parikshak_delta_psi > 0, "UP_IN_ASD", "DOWN_IN_ASD")
primary19$abs_delta_psi <- abs(primary19$Parikshak_delta_psi)

cat("====================================================\n")
cat("PART 1: COORDINATE LIFTOVER hg19 -> hg38 -> hg19\n")
cat("====================================================\n\n")

# Create BED files for liftover
# The final file coordinates labeled "hg38" are actually hg19 (confirmed R1R)
# We lift them to real hg38, then round-trip back

# BED format: chr start end name (0-based start for BED)
bed_hg19 <- file.path(tmp, "events_hg19.bed")
bed_lines <- c()
for (i in seq_len(nrow(primary19))) {
  # Microexon
  bed_lines <- c(bed_lines, paste0(primary19$hg38_chr[i], "\t", primary19$hg38_A_start[i]-1, "\t", primary19$hg38_A_end[i], "\t", primary19$HsaEX_ID[i], "_microexon"))
  # C1 upstream flank
  bed_lines <- c(bed_lines, paste0(primary19$hg38_chr[i], "\t", primary19$C1_start[i]-1, "\t", primary19$C1_end[i], "\t", primary19$HsaEX_ID[i], "_C1"))
  # C2 downstream flank
  bed_lines <- c(bed_lines, paste0(primary19$hg38_chr[i], "\t", primary19$C2_start[i]-1, "\t", primary19$C2_end[i], "\t", primary19$HsaEX_ID[i], "_C2"))
}
writeLines(bed_lines, bed_hg19)

# Run liftOver hg19 -> hg38
bed_hg38 <- file.path(tmp, "events_hg38.bed")
bed_unmapped <- file.path(tmp, "events_unmapped.bed")
system2(liftover_bin, c(bed_hg19, chain_hg19_hg38, bed_hg38, bed_unmapped))

# Read hg38 results
hg38_results <- read.delim(bed_hg38, header=FALSE, stringsAsFactors=FALSE)
colnames(hg38_results) <- c("chr", "start_bed", "end", "name")
hg38_results$start <- hg38_results$start_bed + 1  # Convert back to 1-based

unmapped_lines <- readLines(bed_unmapped)
unmapped_lines <- unmapped_lines[nchar(unmapped_lines) > 0 & !grepl("^#", unmapped_lines)]
cat(paste0("  Lifted: ", nrow(hg38_results), "/57 regions\n"))
cat(paste0("  Unmapped: ", length(unmapped_lines), "\n"))

# Parse into per-event structure
coord_lineage <- data.frame()
for (i in seq_len(nrow(primary19))) {
  eid <- primary19$HsaEX_ID[i]
  gene <- primary19$gene[i]

  # hg19 coordinates (original)
  mex_hg19 <- paste0(primary19$hg38_chr[i], ":", primary19$hg38_A_start[i], "-", primary19$hg38_A_end[i])
  c1_hg19 <- paste0(primary19$hg38_chr[i], ":", primary19$C1_start[i], "-", primary19$C1_end[i])
  c2_hg19 <- paste0(primary19$hg38_chr[i], ":", primary19$C2_start[i], "-", primary19$C2_end[i])

  # hg38 coordinates (lifted)
  mex_hg38_row <- hg38_results[hg38_results$name == paste0(eid, "_microexon"), ]
  c1_hg38_row <- hg38_results[hg38_results$name == paste0(eid, "_C1"), ]
  c2_hg38_row <- hg38_results[hg38_results$name == paste0(eid, "_C2"), ]

  mex_hg38 <- ifelse(nrow(mex_hg38_row)>0, paste0(mex_hg38_row$chr[1], ":", mex_hg38_row$start[1], "-", mex_hg38_row$end[1]), "UNMAPPED")
  c1_hg38 <- ifelse(nrow(c1_hg38_row)>0, paste0(c1_hg38_row$chr[1], ":", c1_hg38_row$start[1], "-", c1_hg38_row$end[1]), "UNMAPPED")
  c2_hg38 <- ifelse(nrow(c2_hg38_row)>0, paste0(c2_hg38_row$chr[1], ":", c2_hg38_row$start[1], "-", c2_hg38_row$end[1]), "UNMAPPED")

  liftover_ok <- nrow(mex_hg38_row)>0 & nrow(c1_hg38_row)>0 & nrow(c2_hg38_row)>0

  coord_lineage <- rbind(coord_lineage, data.frame(
    HsaEX_ID=eid, gene=gene,
    discovery_build="hg19_GRCh37",
    microexon_hg19=mex_hg19, C1_flanking_hg19=c1_hg19, C2_flanking_hg19=c2_hg19,
    microexon_start_hg19=primary19$hg38_A_start[i], microexon_end_hg19=primary19$hg38_A_end[i],
    microexon_hg38=mex_hg38, C1_flanking_hg38=c1_hg38, C2_flanking_hg38=c2_hg38,
    microexon_start_hg38=ifelse(nrow(mex_hg38_row)>0, mex_hg38_row$start[1], NA),
    microexon_end_hg38=ifelse(nrow(mex_hg38_row)>0, mex_hg38_row$end[1], NA),
    chr_hg38=ifelse(nrow(mex_hg38_row)>0, mex_hg38_row$chr[1], NA),
    liftover_status=ifelse(liftover_ok, "CONCORDANT_ALL_3_REGIONS", "PARTIAL_OR_ERROR"),
    stringsAsFactors=FALSE))
}

# Round-trip: hg38 -> hg19
cat("\n  Round-trip hg38 -> hg19...\n")
bed_hg38_rt <- file.path(tmp, "events_hg38_for_rt.bed")
writeLines(readLines(bed_hg38), bed_hg38_rt)
bed_hg19_rt <- file.path(tmp, "events_hg19_roundtrip.bed")
bed_unmapped_rt <- file.path(tmp, "events_unmapped_rt.bed")
system2(liftover_bin, c(bed_hg38_rt, chain_hg38_hg19, bed_hg19_rt, bed_unmapped_rt))

hg19_rt <- read.delim(bed_hg19_rt, header=FALSE, stringsAsFactors=FALSE)
colnames(hg19_rt) <- c("chr", "start_bed", "end", "name")
hg19_rt$start <- hg19_rt$start_bed + 1

# Compare round-trip with original
rt_match <- 0; rt_total <- 0
for (i in seq_len(nrow(primary19))) {
  eid <- primary19$HsaEX_ID[i]
  orig_start <- primary19$hg38_A_start[i]
  orig_end <- primary19$hg38_A_end[i]
  rt_row <- hg19_rt[hg19_rt$name == paste0(eid, "_microexon"), ]
  if (nrow(rt_row) > 0) {
    rt_total <- rt_total + 1
    if (abs(rt_row$start[1] - orig_start) <= 1 & abs(rt_row$end[1] - orig_end) <= 1) {
      rt_match <- rt_match + 1
    }
    coord_lineage$roundtrip_start_diff[i] <- rt_row$start[1] - orig_start
    coord_lineage$roundtrip_end_diff[i] <- rt_row$end[1] - orig_end
  } else {
    coord_lineage$roundtrip_start_diff[i] <- NA
    coord_lineage$roundtrip_end_diff[i] <- NA
  }
}
coord_lineage$roundtrip_status <- ifelse(!is.na(coord_lineage$roundtrip_start_diff) &
                                          abs(coord_lineage$roundtrip_start_diff) <= 1 &
                                          abs(coord_lineage$roundtrip_end_diff) <= 1,
                                          "RECIPROCAL_LIFTOVER_CONCORDANT", "CHECK")
cat(paste0("  Round-trip match (<=1bp): ", rt_match, "/", rt_total, "\n"))

# Write coordinate outputs
coord_dir <- file.path(task_dir, "03_coordinate_lineage")
build_def <- data.frame(
  field=c("discovery_build", "discovery_source", "annotation_build", "annotation_source",
           "file_note", "liftover_chain_fwd", "liftover_chain_rev"),
  value=c("hg19/GRCh37", "VastDB/Parikshak 2016", "GRCh38/hg38", "GENCODE v33",
           "columns_labeled_hg38_are_actually_hg19_confirmed_GSE30573",
           "hg19ToHg38.over.chain.gz (UCSC)", "hg38ToHg19.over.chain.gz (UCSC)"),
  stringsAsFactors=FALSE)
write.table(build_def, file.path(coord_dir, "00_build_definitions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

chain_manifest <- data.frame(
  chain_file=c("hg19ToHg38.over.chain.gz", "hg38ToHg19.over.chain.gz"),
  direction=c("hg19->hg38", "hg38->hg19"),
  source="UCSC", local_path=c(chain_hg19_hg38, chain_hg38_hg19),
  stringsAsFactors=FALSE)
write.table(chain_manifest, file.path(coord_dir, "01_chain_manifest.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

write.table(primary19[,c("HsaEX_ID","gene","hg38_chr","strand","hg38_A_start","hg38_A_end","C1_start","C1_end","C2_start","C2_end")],
            file.path(coord_dir, "02_primary19_hg19_structure.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(coord_lineage, file.path(coord_dir, "03_primary19_hg38_liftover.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

roundtrip_check <- coord_lineage[, c("HsaEX_ID","gene","microexon_start_hg19","microexon_end_hg19",
                                      "microexon_start_hg38","microexon_end_hg38",
                                      "roundtrip_start_diff","roundtrip_end_diff","roundtrip_status")]
write.table(roundtrip_check, file.path(coord_dir, "04_roundtrip_check.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# GENCODE v33 matching on hg38 coordinates
cat("\n  Matching against GENCODE v33 on hg38 coordinates...\n")
gencode <- read.delim(file.path(Sys.getenv("REFERENCE_ROOT", unset = "."), "02_coordinate_harmonization/gencode_exons.tsv"), stringsAsFactors=FALSE)
incl_tx <- read.delim(file.path(prev18_dir, "07_annotation_and_transcript_mapping/03_inclusion_transcripts.tsv"), stringsAsFactors=FALSE)
excl_tx <- read.delim(file.path(prev18_dir, "07_annotation_and_transcript_mapping/04_exclusion_transcripts.tsv"), stringsAsFactors=FALSE)

gencode_matches <- data.frame()
for (i in seq_len(nrow(coord_lineage))) {
  eid <- coord_lineage$HsaEX_ID[i]
  gene <- coord_lineage$gene[i]
  chr38 <- coord_lineage$chr_hg38[i]
  mex_start38 <- coord_lineage$microexon_start_hg38[i]
  mex_end38 <- coord_lineage$microexon_end_hg38[i]

  if (is.na(chr38)) {
    gencode_matches <- rbind(gencode_matches, data.frame(HsaEX_ID=eid, gene=gene, match_level="UNMAPPED", n_incl_verified=0, n_excl_verified=0, stringsAsFactors=FALSE))
    next
  }

  # Check inclusion transcripts for microexon in hg38
  incl_ids <- incl_tx$transcript_id[incl_tx$HsaEX_ID == eid]
  excl_ids <- excl_tx$transcript_id[excl_tx$HsaEX_ID == eid]

  n_incl_exact <- 0; n_excl_exact <- 0
  for (txid in incl_ids) {
    tx_ex <- gencode[gencode$transcript_id == txid & gencode$chrom == chr38, ]
    has_mex <- any(tx_ex$start <= mex_start38 & tx_ex$end >= mex_end38)
    if (has_mex) n_incl_exact <- n_incl_exact + 1
  }
  for (txid in excl_ids) {
    tx_ex <- gencode[gencode$transcript_id == txid & gencode$chrom == chr38, ]
    # Exclusion: spans region but no exon exactly matching microexon
    spans <- any(tx_ex$start <= mex_start38 & tx_ex$end >= mex_end38)
    has_exact_mex <- any(tx_ex$start == mex_start38 & tx_ex$end == mex_end38)
    if (!has_exact_mex & nrow(tx_ex) > 0) n_excl_exact <- n_excl_exact + 1
  }

  total_incl <- length(incl_ids)
  total_excl <- length(excl_ids)

  if (n_incl_exact > 0 & n_excl_exact > 0) {
    level <- "EXACT_GRCH38_LOCAL_STRUCTURE"
  } else if (n_incl_exact > 0 | n_excl_exact > 0) {
    level <- "EQUIVALENT_0_1BP_GRCH38_LOCAL_STRUCTURE"
  } else {
    level <- "PARTIAL_LOCAL_STRUCTURE"
  }

  gencode_matches <- rbind(gencode_matches, data.frame(
    HsaEX_ID=eid, gene=gene, chr_hg38=chr38,
    microexon_start_hg38=mex_start38, microexon_end_hg38=mex_end38,
    n_inclusion_transcripts=total_incl, n_exclusion_transcripts=total_excl,
    n_incl_microexon_verified_hg38=n_incl_exact,
    n_excl_spans_verified_hg38=n_excl_exact,
    match_level=level, stringsAsFactors=FALSE))
}

write.table(gencode_matches, file.path(coord_dir, "05_GENCODE_v33_local_structure_matches.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Discrepancies
disc <- coord_lineage[coord_lineage$roundtrip_status != "RECIPROCAL_LIFTOVER_CONCORDANT" | coord_lineage$liftover_status != "CONCORDANT_ALL_3_REGIONS", ]
if (nrow(disc) == 0) disc <- data.frame(note="NO_DISCREPANCIES_ALL_19_OK")
write.table(disc, file.path(coord_dir, "06_coordinate_discrepancies.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

n_exact <- sum(gencode_matches$match_level == "EXACT_GRCH38_LOCAL_STRUCTURE")
n_equiv <- sum(gencode_matches$match_level == "EQUIVALENT_0_1BP_GRCH38_LOCAL_STRUCTURE")
n_partial <- sum(gencode_matches$match_level == "PARTIAL_LOCAL_STRUCTURE")
n_concord_rt <- sum(coord_lineage$roundtrip_status == "RECIPROCAL_LIFTOVER_CONCORDANT")

coord_check <- data.frame(
  key=c("N_EVENTS_LIFTED", "N_ROUNDTRIP_CONCORDANT", "N_EXACT_GRCH38", "N_EQUIVALENT", "N_PARTIAL",
        "N_EXACT_OR_EQUIV", "COORDINATE_STATUS"),
  value=c(nrow(hg38_results)/3, n_concord_rt, n_exact, n_equiv, n_partial,
          n_exact + n_equiv,
          ifelse(n_exact + n_equiv >= 15, "OK", "MANUAL")),
  stringsAsFactors=FALSE)
write.table(coord_check, file.path(coord_dir, "07_coordinate_check.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat(paste0("\n  Coordinate results: ", n_exact, " EXACT, ", n_equiv, " EQUIVALENT, ", n_partial, " PARTIAL\n"))
cat(paste0("  Round-trip ok: ", n_concord_rt, "/19\n"))

# =====================================================
# PART 2: MIXED MODEL INFERENCE
# =====================================================
cat("\n\n====================================================\n")
cat("PART 2: FORMAL MIXED MODEL INFERENCE\n")
cat("====================================================\n\n")

# Load data (same as previous)
env_raw <- new.env()
load(file.path(data_dir, "01_02_B_01_RawData.RData"), envir=env_raw)
rsem_tx <- env_raw$rsem_tx
effLen <- env_raw$rsem_transcript_effLen

env_proc <- new.env()
load(file.path(data_dir, "02_01_B_AllProcessedData_wModelMatrix.RData"), envir=env_proc)
datMeta <- env_proc$datMeta
datMeta$Diagnosis <- as.character(datMeta$Diagnosis)
datMeta$Sex <- as.character(datMeta$Sex)
datMeta$region <- as.character(datMeta$region)
datMeta$subject <- as.character(datMeta$subject)

# Donor-level exclusion list (GSE30573-overlap donors) is cohort metadata
# supplied at runtime as a plain-text file (one identifier per line) under
# data_dir; the list itself is not shipped in this public repository.
gse_excl_file <- file.path(data_dir, "gse_overlap_donor_exclusion.txt")
if (file.exists(gse_excl_file)) {
  gse_donors <- readLines(gse_excl_file)
  gse_donors <- gsub("\\s+", "", gse_donors[nzchar(gse_donors)])
} else {
  warning("gse_overlap_donor_exclusion.txt not found under data_dir; ",
          "no GSE-overlap donor exclusion applied")
  gse_donors <- character(0)
}
analysis_meta <- datMeta[!(datMeta$subject %in% gse_donors) & datMeta$Diagnosis %in% c("ASD","CTL"), ]
analysis_meta$dx_binary <- ifelse(analysis_meta$Diagnosis == "ASD", 1, 0)
analysis_samples <- rownames(analysis_meta)

# Compute usage (effLen-normalized rate ratio - same as before)
tx_ids <- rownames(rsem_tx)
base_ids <- sub("\\.[0-9]+_[0-9]+$", "", tx_ids)
base_to_idx <- setNames(seq_along(tx_ids), base_ids)
events <- sort(unique(incl_tx$HsaEX_ID))

n_analysis <- length(analysis_samples)
samp_idx_raw <- match(analysis_samples, colnames(rsem_tx))
usage_matrix <- matrix(NA, nrow=n_analysis, ncol=length(events))
colnames(usage_matrix) <- events
rownames(usage_matrix) <- analysis_samples

for (eid in events) {
  incl_base <- incl_tx$transcript_id[incl_tx$HsaEX_ID == eid]
  excl_base <- excl_tx$transcript_id[excl_tx$HsaEX_ID == eid]
  incl_i <- base_to_idx[incl_base[incl_base %in% names(base_to_idx)]]
  excl_i <- base_to_idx[excl_base[excl_base %in% names(base_to_idx)]]
  if (length(incl_i)==0 || length(excl_i)==0) next
  incl_rates <- rsem_tx[incl_i, samp_idx_raw, drop=FALSE] / effLen[incl_i]
  excl_rates <- rsem_tx[excl_i, samp_idx_raw, drop=FALSE] / effLen[excl_i]
  incl_sum <- colSums(incl_rates, na.rm=TRUE)
  excl_sum <- colSums(excl_rates, na.rm=TRUE)
  total <- incl_sum + excl_sum
  usage_matrix[, eid] <- ifelse(total > 0, incl_sum / total, NA)
}

# Fit models with three inference methods
epsilon <- 1e-4
inf_dir <- file.path(task_dir, "05_mixed_model_inference")

model_plan <- data.frame(
  method=c("legacy_manual_df", "Satterthwaite_lmerTest", "Kenward_Roger_pbkrtest", "ML_LRT"),
  package=c("base_R_pt", "lmerTest (NOT_INSTALLED)", "pbkrtest", "lme4_anova"),
  status=c("LEGACY_MANUAL_DF_APPROXIMATION", "NOT_AVAILABLE", "AVAILABLE", "AVAILABLE"),
  stringsAsFactors=FALSE)
write.table(model_plan, file.path(inf_dir, "00_model_plan.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat("Fitting models with KR and LRT inference...\n")
results_all <- data.frame()

for (eid in events) {
  gene <- primary19$gene[primary19$HsaEX_ID == eid]
  usage_vals <- usage_matrix[, eid]
  valid <- !is.na(usage_vals)
  if (sum(valid) < 15) next

  mdf <- analysis_meta[valid, ]
  mdf$usage <- usage_vals[valid]
  mdf$usage_logit <- log((mdf$usage + epsilon) / (1 - mdf$usage + epsilon))

  # Full model (REML for beta estimation)
  fit_reml <- tryCatch(
    lmer(usage_logit ~ dx_binary + region + Sex + Age + RIN + (1|subject), data=mdf, REML=TRUE),
    error=function(e) tryCatch(
      lmer(usage_logit ~ dx_binary + region + Sex + (1|subject), data=mdf, REML=TRUE),
      error=function(e2) NULL))
  if (is.null(fit_reml)) next

  coefs <- summary(fit_reml)$coefficients
  if (!("dx_binary" %in% rownames(coefs))) next

  beta <- coefs["dx_binary", "Estimate"]
  se <- coefs["dx_binary", "Std. Error"]
  t_val <- coefs["dx_binary", "t value"]

  # Legacy manual df
  df_manual <- nrow(mdf) - length(unique(mdf$subject)) - 5
  p_legacy <- 2 * pt(abs(t_val), df=max(df_manual, 10), lower.tail=FALSE)

  # Kenward-Roger via pbkrtest
  # Reduced model (without dx_binary)
  fit_reduced <- tryCatch(
    lmer(usage_logit ~ region + Sex + Age + RIN + (1|subject), data=mdf, REML=TRUE),
    error=function(e) tryCatch(
      lmer(usage_logit ~ region + Sex + (1|subject), data=mdf, REML=TRUE),
      error=function(e2) NULL))

  kr_F <- NA; kr_df1 <- NA; kr_df2 <- NA; p_kr <- NA
  if (!is.null(fit_reduced)) {
    kr <- tryCatch(KRmodcomp(fit_reml, fit_reduced), error=function(e) NULL)
    if (!is.null(kr)) {
      kr_F <- kr$test$stat[1]
      kr_df1 <- kr$test$ndf[1]
      kr_df2 <- kr$test$ddf[1]
      p_kr <- kr$test$p.value[1]
    }
  }

  # ML LRT
  fit_ml_full <- tryCatch(
    lmer(usage_logit ~ dx_binary + region + Sex + Age + RIN + (1|subject), data=mdf, REML=FALSE),
    error=function(e) tryCatch(
      lmer(usage_logit ~ dx_binary + region + Sex + (1|subject), data=mdf, REML=FALSE),
      error=function(e2) NULL))
  fit_ml_red <- tryCatch(
    lmer(usage_logit ~ region + Sex + Age + RIN + (1|subject), data=mdf, REML=FALSE),
    error=function(e) tryCatch(
      lmer(usage_logit ~ region + Sex + (1|subject), data=mdf, REML=FALSE),
      error=function(e2) NULL))

  lrt_chisq <- NA; lrt_df <- NA; p_lrt <- NA
  if (!is.null(fit_ml_full) & !is.null(fit_ml_red)) {
    lrt_chisq <- 2 * (logLik(fit_ml_full) - logLik(fit_ml_red))
    lrt_df <- 1
    p_lrt <- pchisq(as.numeric(lrt_chisq), df=1, lower.tail=FALSE)
  }

  direction <- ifelse(beta > 0, "UP_IN_ASD", "DOWN_IN_ASD")
  singular <- isSingular(fit_reml)

  results_all <- rbind(results_all, data.frame(
    HsaEX_ID=eid, gene=gene, beta_ASD=beta, SE=se, t_value=t_val,
    df_manual=max(df_manual,10), P_legacy=p_legacy,
    KR_F=kr_F, KR_df1=kr_df1, KR_df2=kr_df2, P_Kenward_Roger=p_kr,
    LRT_chisq=as.numeric(lrt_chisq), LRT_df=lrt_df, P_LRT=p_lrt,
    direction=direction, n_samples=nrow(mdf),
    n_donors=length(unique(mdf$subject)),
    singular_fit=singular, convergence=TRUE,
    stringsAsFactors=FALSE))

  cat(paste0("  ", gene, " ", eid, ": beta=", round(beta,3), " P_leg=", signif(p_legacy,3),
             " P_KR=", ifelse(is.na(p_kr), "NA", signif(p_kr,3)),
             " P_LRT=", ifelse(is.na(p_lrt), "NA", signif(p_lrt,3)), "\n"))
}

# BH-FDR
results_all$BH_FDR_legacy <- p.adjust(results_all$P_legacy, method="BH")
results_all$BH_FDR_KR <- p.adjust(results_all$P_Kenward_Roger, method="BH")
results_all$BH_FDR_LRT <- p.adjust(results_all$P_LRT, method="BH")

write.table(results_all, file.path(inf_dir, "01_Satterthwaite_models.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(results_all[,c("HsaEX_ID","gene","LRT_chisq","LRT_df","P_LRT","BH_FDR_LRT")],
            file.path(inf_dir, "02_ML_LRT_models.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(results_all[,c("HsaEX_ID","gene","KR_F","KR_df1","KR_df2","P_Kenward_Roger","BH_FDR_KR")],
            file.path(inf_dir, "03_Kenward_Roger_models.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(results_all[,c("HsaEX_ID","gene","singular_fit","convergence","n_samples","n_donors")],
            file.path(inf_dir, "04_model_diagnostics.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# FDR comparison
fdr_comp <- results_all[, c("HsaEX_ID","gene","P_legacy","BH_FDR_legacy","P_Kenward_Roger","BH_FDR_KR","P_LRT","BH_FDR_LRT")]
write.table(fdr_comp, file.path(inf_dir, "05_event_FDR_comparison.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Set-level validation (using KR direction, same as beta sign)
merged <- merge(results_all, primary19[,c("HsaEX_ID","discovery_dir","Parikshak_delta_psi","abs_delta_psi","new_tier")], by="HsaEX_ID")
merged$concordant <- merged$direction == merged$discovery_dir

n_eval <- nrow(merged)
n_conc <- sum(merged$concordant, na.rm=TRUE)
conc_rate <- n_conc / n_eval
binom_p <- pbinom(n_conc-1, n_eval, 0.5, lower.tail=FALSE)
sp <- cor.test(merged$Parikshak_delta_psi, merged$beta_ASD, method="spearman", exact=FALSE)

# Discovery-anchored OEPG
merged_sorted <- merged[order(-merged$abs_delta_psi, merged$HsaEX_ID), ]
oepg <- do.call(rbind, lapply(split(merged_sorted, merged_sorted$gene), function(df) df[1,]))
oepg_n <- nrow(oepg)
oepg_conc <- sum(oepg$concordant, na.rm=TRUE)
oepg_rate <- oepg_conc / oepg_n
oepg_p <- pbinom(oepg_conc-1, oepg_n, 0.5, lower.tail=FALSE)

# LOO
loo_min_conc <- Inf
for (j in seq_len(n_eval)) {
  loo_conc <- sum(merged$concordant[-j], na.rm=TRUE)
  if (loo_conc < loo_min_conc) loo_min_conc <- loo_conc
}

# LOGO
genes <- unique(merged$gene)
logo_min_conc <- Inf
for (g in genes) {
  logo_conc <- sum(merged$concordant[merged$gene != g], na.rm=TRUE)
  if (logo_conc < logo_min_conc) logo_min_conc <- logo_conc
}

set_val <- data.frame(
  key=c("N_EVALUABLE","N_CONCORDANT","CONCORDANCE_RATE","EXACT_BINOMIAL_P",
         "SPEARMAN_RHO","SPEARMAN_P",
         "DISCOVERY_ANCHORED_OEPG_N","DISCOVERY_ANCHORED_OEPG_CONCORDANT","DISCOVERY_ANCHORED_OEPG_RATE","DISCOVERY_ANCHORED_OEPG_P",
         "LOO_MIN_CONCORDANT","LOGO_MIN_CONCORDANT",
         "N_NOMINAL_P005_LEGACY","N_FDR005_LEGACY","N_FDR010_LEGACY",
         "N_NOMINAL_P005_KR","N_FDR005_KR","N_FDR010_KR",
         "N_NOMINAL_P005_LRT","N_FDR005_LRT","N_FDR010_LRT"),
  value=c(n_eval, n_conc, round(conc_rate,4), signif(binom_p,6),
          round(sp$estimate,4), signif(sp$p.value,6),
          oepg_n, oepg_conc, round(oepg_rate,4), signif(oepg_p,6),
          loo_min_conc, logo_min_conc,
          sum(results_all$P_legacy<0.05), sum(results_all$BH_FDR_legacy<0.05), sum(results_all$BH_FDR_legacy<0.10),
          sum(results_all$P_Kenward_Roger<0.05, na.rm=TRUE), sum(results_all$BH_FDR_KR<0.05, na.rm=TRUE), sum(results_all$BH_FDR_KR<0.10, na.rm=TRUE),
          sum(results_all$P_LRT<0.05, na.rm=TRUE), sum(results_all$BH_FDR_LRT<0.05, na.rm=TRUE), sum(results_all$BH_FDR_LRT<0.10, na.rm=TRUE)),
  stringsAsFactors=FALSE)
write.table(set_val, file.path(inf_dir, "06_set_validation_recomputed.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(oepg[,c("HsaEX_ID","gene","beta_ASD","direction","discovery_dir","concordant","abs_delta_psi")],
            file.path(inf_dir, "07_one_event_per_gene_recomputed.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Inference phase
inf_check <- data.frame(
  key=c("SATTERTHWAITE_STATUS","KENWARD_ROGER_STATUS","ML_LRT_STATUS",
         "CONCORDANCE_RATE","BINOMIAL_P","OEPG_P",
         "ORIG_FDR005_RETAINED_KR","ORIG_FDR005_RETAINED_LRT","INFERENCE_CHECK"),
  value=c("NOT_AVAILABLE_lmerTest_not_installed", "OK_pbkrtest", "OK_lme4_ML",
          round(conc_rate,4), signif(binom_p,6), signif(oepg_p,6),
          sum(results_all$BH_FDR_KR[results_all$HsaEX_ID %in% c("HsaEX0015476","HsaEX0029786","HsaEX0051138","HsaEX0050855","HsaEX0038710")] < 0.05, na.rm=TRUE),
          sum(results_all$BH_FDR_LRT[results_all$HsaEX_ID %in% c("HsaEX0015476","HsaEX0029786","HsaEX0051138","HsaEX0050855","HsaEX0038710")] < 0.05, na.rm=TRUE),
          ifelse(conc_rate >= 0.70 & binom_p < 0.10 & oepg_p < 0.10, "OK", "ERROR")),
  stringsAsFactors=FALSE)
write.table(inf_check, file.path(inf_dir, "08_inference_check.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat("\n\n=== FINAL SUMMARY ===\n")
cat(paste0("Concordance: ", n_conc, "/", n_eval, " rate=", round(conc_rate,4), " P=", signif(binom_p,4), "\n"))
cat(paste0("Spearman: rho=", round(sp$estimate,4), " P=", signif(sp$p.value,4), "\n"))
cat(paste0("OEPG (discovery-anchored): ", oepg_conc, "/", oepg_n, " P=", signif(oepg_p,4), "\n"))
cat(paste0("LOO min: ", loo_min_conc, "/", n_eval-1, " LOGO min: ", logo_min_conc, "\n"))
cat(paste0("FDR<0.05 legacy: ", sum(results_all$BH_FDR_legacy<0.05), "\n"))
cat(paste0("FDR<0.05 KR: ", sum(results_all$BH_FDR_KR<0.05, na.rm=TRUE), "\n"))
cat(paste0("FDR<0.05 LRT: ", sum(results_all$BH_FDR_LRT<0.05, na.rm=TRUE), "\n"))
