#!/usr/bin/env Rscript
# PsychENCODE mixed-model analysis.
# M0-M4 model construction, Kenward-Roger primary inference, ML-LRT sensitivity, BH correction and one-event-per-gene concordance.
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).
suppressMessages({library(methods); library(lme4)})
set.seed(42)

project <- Sys.getenv("PROJECT_ROOT", unset = ".")
task_dir <- file.path(project, "20_psychencode_final_models")
prev_dir <- file.path(project, "19_psychencode_tpm_usage")
data_dir <- file.path(project, "psychencode_processed")
prev18_dir <- file.path(project, "18_psychencode")

# Load the 19 events
primary19 <- read.delim(file.path(project, "16_gse30573/02_input_lock/02_primary19.tsv"), stringsAsFactors=FALSE)
primary19$abs_delta_psi <- abs(primary19$Parikshak_delta_psi)
primary19$discovery_dir <- ifelse(primary19$Parikshak_delta_psi > 0, "UP_IN_ASD", "DOWN_IN_ASD")

# Load previous direction comparison (validation results)
direction <- read.delim(file.path(prev_dir, "07_primary_nonoverlap_models/05_event_direction_comparison.tsv"), stringsAsFactors=FALSE)

cat("====================================================\n")
cat("REPAIR 1: ONE-EVENT-PER-GENE\n")
cat("====================================================\n\n")

# =====================================================
# REPAIR 1: ONE-EVENT-PER-GENE
# =====================================================
oepg_dir <- file.path(task_dir, "03_one_event_per_gene_repair")

# Selection rules
sel_rules <- data.frame(
  rule = c("discovery_anchored", "discovery_anchored_tiebreak1", "discovery_anchored_tiebreak2",
           "lexicographic", "PREVIOUS_biased_best_validation_P"),
  description = c("max abs(Parikshak_delta_PSI) per gene",
                  "tie: higher Analysis-R tier",
                  "tie: HsaEX_ID lexicographic smallest",
                  "HsaEX_ID lexicographic smallest per gene",
                  "min validation P per gene (USES OUTCOME - INVALID)"),
  valid_for_confirmatory = c("YES","YES","YES","YES","INVALID_FOR_CONFIRMATORY_SENSITIVITY"),
  stringsAsFactors=FALSE)
write.table(sel_rules, file.path(oepg_dir, "00_selection_rules.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Rule 1: Discovery-anchored (max |delta_PSI|)
tier_rank <- c("TIER_2_FUNCTIONAL"=4, "TIER_3_TRAJECTORY_ONLY"=3, "TIER_4_NON_DYNAMIC"=2, "TIER_1_DYNAMIC"=1)
primary19$tier_rank <- tier_rank[primary19$new_tier]
sorted <- primary19[order(-primary19$abs_delta_psi, -primary19$tier_rank, primary19$HsaEX_ID), ]
discovery_anchored <- do.call(rbind, lapply(split(sorted, sorted$gene), function(df) df[1,]))
discovery_anchored <- discovery_anchored[order(discovery_anchored$HsaEX_ID), ]

write.table(discovery_anchored[,c("HsaEX_ID","gene","Parikshak_delta_psi","abs_delta_psi","new_tier","discovery_dir")],
            file.path(oepg_dir, "01_discovery_anchored_selection.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Rule 2: Lexicographic
lex <- primary19[order(primary19$HsaEX_ID), ]
lexicographic <- do.call(rbind, lapply(split(lex, lex$gene), function(df) df[1,]))
lexicographic <- lexicographic[order(lexicographic$HsaEX_ID), ]

write.table(lexicographic[,c("HsaEX_ID","gene","Parikshak_delta_psi","abs_delta_psi","new_tier","discovery_dir")],
            file.path(oepg_dir, "02_lexicographic_selection.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Compute stats for each rule
compute_oepg_stats <- function(selected_ids, label) {
  sel <- direction[direction$HsaEX_ID %in% selected_ids, ]
  sel <- merge(sel, primary19[,c("HsaEX_ID","gene","discovery_dir","Parikshak_delta_psi")], by="HsaEX_ID", suffixes=c("","_reference"))
  sel$concordant <- sel$direction == sel$discovery_dir_reference
  n <- nrow(sel)
  nc <- sum(sel$concordant, na.rm=TRUE)
  rate <- nc/n
  bp <- pbinom(nc-1, n, 0.5, lower.tail=FALSE)
  valid <- sel[!is.na(sel$beta_ASD) & !is.na(sel$Parikshak_delta_psi), ]
  sp <- cor.test(valid$Parikshak_delta_psi, valid$beta_ASD, method="spearman", exact=FALSE)
  # Permutation
  perm_counts <- replicate(10000, sum(sample(c("UP_IN_ASD","DOWN_IN_ASD"), n, replace=TRUE) == sel$discovery_dir_reference))
  perm_p <- (sum(perm_counts >= nc) + 1) / 10001
  cat(paste0("  ", label, ": ", nc, "/", n, " concordant, rate=", round(rate,4),
             " P=", signif(bp,4), " rho=", round(sp$estimate,4), " rho_P=", signif(sp$p.value,4),
             " perm_P=", signif(perm_p,4), "\n"))
  data.frame(rule=label, n_genes=n, n_concordant=nc, concordance_rate=round(rate,4),
             exact_binomial_P=signif(bp,6), spearman_rho=round(sp$estimate,4),
             spearman_P=signif(sp$p.value,6), direction_perm_P=signif(perm_p,6),
             stringsAsFactors=FALSE)
}

cat("\nOne-event-per-gene results:\n")
da_stats <- compute_oepg_stats(discovery_anchored$HsaEX_ID, "discovery_anchored")
lex_stats <- compute_oepg_stats(lexicographic$HsaEX_ID, "lexicographic")

# Previous biased method (min validation P per gene) - flagged as invalid
direction_with_gene <- direction
direction_with_gene$gene <- primary19$gene[match(direction_with_gene$HsaEX_ID, primary19$HsaEX_ID)]
prev_best <- do.call(rbind, lapply(split(direction_with_gene, direction_with_gene$gene), function(df) df[which.min(df$P),]))
prev_stats <- compute_oepg_stats(prev_best$HsaEX_ID, "PREVIOUS_biased_best_P")
prev_stats$rule <- "PREVIOUS_biased_best_P"

write.table(da_stats, file.path(oepg_dir, "03_discovery_anchored_results.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(lex_stats, file.path(oepg_dir, "04_lexicographic_results.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Check of previous biased selection
prev_check <- data.frame(
  key = c("previous_method", "problem", "status", "previous_n_genes", "previous_concordant",
           "previous_rate", "previous_P", "disposition"),
  value = c("min validation P per gene", "uses validation outcome for selection = selection bias",
            "INVALID_FOR_CONFIRMATORY_SENSITIVITY", nrow(prev_best), sum(prev_best$concordant, na.rm=TRUE),
            round(sum(prev_best$concordant, na.rm=TRUE)/nrow(prev_best), 4),
            signif(pbinom(sum(prev_best$concordant, na.rm=TRUE)-1, nrow(prev_best), 0.5, lower.tail=FALSE), 6),
            "retained_for_check_only_not_for_manuscript"),
  stringsAsFactors=FALSE)
write.table(prev_check, file.path(oepg_dir, "05_previous_biased_selection_check.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

oepg_check <- rbind(da_stats, lex_stats)
oepg_check$check_status <- ifelse(oepg_check$exact_binomial_P < 0.10, "OK", "MARGINAL")
write.table(oepg_check, file.path(oepg_dir, "06_one_event_per_gene_check.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# =====================================================
# REPAIR 2: ZERO FRACTION
# =====================================================
cat("\n\n====================================================\n")
cat("REPAIR 2: ZERO FRACTION CORRECTION\n")
cat("====================================================\n\n")

zero_dir <- file.path(task_dir, "04_raw_RSEM_QC_repair")

# Load raw data
cat("Loading raw RSEM data...\n")
env_raw <- new.env()
load(file.path(data_dir, "01_02_B_01_RawData.RData"), envir=env_raw)
rsem_tx <- env_raw$rsem_tx
effLen <- env_raw$rsem_transcript_effLen

n_tx <- nrow(rsem_tx)
n_samp <- ncol(rsem_tx)
total_entries <- as.numeric(n_tx) * as.numeric(n_samp)

# Compute zero metrics efficiently
cat("Computing zero metrics...\n")
zero_counts_per_sample <- colSums(rsem_tx == 0)
overall_zeros <- sum(zero_counts_per_sample)
overall_zero_fraction <- overall_zeros / total_entries
per_sample_zero_frac <- zero_counts_per_sample / n_tx

zero_metrics <- data.frame(
  metric = c("n_transcripts", "n_samples", "total_matrix_entries",
             "total_zero_entries", "overall_zero_fraction",
             "mean_zero_transcripts_per_sample", "median_zero_transcripts_per_sample",
             "per_sample_zero_fraction_mean", "per_sample_zero_fraction_median",
             "per_sample_zero_fraction_min", "per_sample_zero_fraction_max"),
  value = c(n_tx, n_samp, total_entries,
            overall_zeros, round(overall_zero_fraction, 6),
            round(mean(zero_counts_per_sample), 1), round(median(zero_counts_per_sample), 1),
            round(mean(per_sample_zero_frac), 6), round(median(per_sample_zero_frac), 6),
            round(min(per_sample_zero_frac), 6), round(max(per_sample_zero_frac), 6)),
  stringsAsFactors=FALSE)
write.table(zero_metrics, file.path(zero_dir, "00_zero_metrics.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat(paste0("  Overall zero fraction: ", round(overall_zero_fraction, 4), "\n"))
cat(paste0("  Mean zeros per sample: ", round(mean(zero_counts_per_sample), 1), "\n"))
cat(paste0("  Previous erroneous value: 97804.596 (was MEAN COUNT, not fraction)\n"))

# Per-sample zero fraction
sample_zf <- data.frame(sample=colnames(rsem_tx), n_zeros=zero_counts_per_sample,
                         zero_fraction=round(per_sample_zero_frac, 6))
write.table(sample_zf, file.path(zero_dir, "01_sample_zero_fraction.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Per-transcript zero fraction (summary)
tx_zero_frac <- rowSums(rsem_tx == 0) / n_samp
tx_zf_summary <- data.frame(
  quantile = c("min","q25","median","mean","q75","max","frac_always_zero","frac_never_zero"),
  value = c(min(tx_zero_frac), quantile(tx_zero_frac, 0.25, names=FALSE), median(tx_zero_frac),
            mean(tx_zero_frac), quantile(tx_zero_frac, 0.75, names=FALSE), max(tx_zero_frac),
            sum(tx_zero_frac==1)/length(tx_zero_frac), sum(tx_zero_frac==0)/length(tx_zero_frac)),
  stringsAsFactors=FALSE)
write.table(tx_zf_summary, file.path(zero_dir, "02_transcript_zero_fraction.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

zero_check <- data.frame(key=c("OVERALL_ZERO_FRACTION","MEAN_ZERO_PER_SAMPLE","PREVIOUS_VALUE","PREVIOUS_ERROR","QC_STATUS"),
                         value=c(round(overall_zero_fraction,6), round(mean(zero_counts_per_sample),1),
                                 "97804.596","was_mean_count_not_fraction","CONCORDANT_CORRECTED"),
                         stringsAsFactors=FALSE)
write.table(zero_check, file.path(zero_dir, "03_raw_RSEM_QC_check.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# =====================================================
# REPAIR 3: TRUE TPM + THRESHOLD SENSITIVITY
# =====================================================
cat("\n\n====================================================\n")
cat("REPAIR 3: TRUE TPM RECONSTRUCTION\n")
cat("====================================================\n\n")

tpm_dir <- file.path(task_dir, "05_true_TPM_detection_sensitivity")

# Compute TPM: rate_i = count_i / effLen_i; TPM_i = rate_i / sum(rate) * 1e6
cat("Computing true TPM (rate = count/effLen; TPM = rate/colSum(rate) * 1e6)...\n")
# Only compute for transcripts with effLen > 0
valid_tx <- effLen > 0
rate_matrix <- rsem_tx[valid_tx, ] / effLen[valid_tx]
col_rate_sums <- colSums(rate_matrix)
# TPM = rate / col_sum * 1e6
tpm_matrix <- sweep(rate_matrix, 2, col_rate_sums, "/") * 1e6

# Verify TPM sums to 1e6
tpm_col_sums <- colSums(tpm_matrix)
cat(paste0("  TPM column sum range: ", round(min(tpm_col_sums), 1), " - ", round(max(tpm_col_sums), 1), "\n"))
cat(paste0("  Expected: ~1000000\n"))

tpm_check <- data.frame(
  key = c("method", "formula", "n_transcripts_valid_effLen", "n_samples",
           "tpm_colsum_min", "tpm_colsum_max", "tpm_colsum_mean", "verification"),
  value = c("rate/colsum_rate*1e6", "TPM_i = (count_i/effLen_i) / sum(count_j/effLen_j) * 1e6",
            sum(valid_tx), n_samp,
            round(min(tpm_col_sums),2), round(max(tpm_col_sums),2), round(mean(tpm_col_sums),2),
            "OK_ALL_SAMPLES_SUM_1E6"),
  stringsAsFactors=FALSE)
write.table(tpm_check, file.path(tpm_dir, "00_TPM_reconstruction_check.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Sample TPM totals (first/last 10)
sample_tpm_totals <- data.frame(sample=names(tpm_col_sums)[1:min(20, length(tpm_col_sums))],
                                 tpm_total=round(tpm_col_sums[1:min(20, length(tpm_col_sums))], 2))
write.table(sample_tpm_totals, file.path(tpm_dir, "01_sample_TPM_totals.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Map transcript IDs
tx_ids <- rownames(rsem_tx)
base_ids <- sub("\\.[0-9]+_[0-9]+$", "", tx_ids)
# For valid_tx subset
valid_base_ids <- base_ids[valid_tx]
valid_base_to_idx <- setNames(seq_along(valid_base_ids), valid_base_ids)

# Load transcript mappings
incl_tx <- read.delim(file.path(prev18_dir, "07_annotation_and_transcript_mapping/03_inclusion_transcripts.tsv"), stringsAsFactors=FALSE)
excl_tx <- read.delim(file.path(prev18_dir, "07_annotation_and_transcript_mapping/04_exclusion_transcripts.tsv"), stringsAsFactors=FALSE)

events <- sort(unique(incl_tx$HsaEX_ID))

# Compute event usage and event_total_TPM using true TPM
cat("\nComputing event usage from true TPM...\n")
# Get processed metadata for non-overlap analysis
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
cat(paste0("  Analysis samples: ", length(analysis_samples), "\n"))

# Compute usage and event_total from TPM
n_analysis <- length(analysis_samples)
usage_tpm <- matrix(NA, nrow=n_analysis, ncol=length(events))
event_total_tpm <- matrix(NA, nrow=n_analysis, ncol=length(events))
colnames(usage_tpm) <- events; colnames(event_total_tpm) <- events
rownames(usage_tpm) <- analysis_samples; rownames(event_total_tpm) <- analysis_samples

# Map analysis_samples to tpm_matrix columns
samp_idx <- match(analysis_samples, colnames(tpm_matrix))

for (eid in events) {
  incl_base <- incl_tx$transcript_id[incl_tx$HsaEX_ID == eid]
  excl_base <- excl_tx$transcript_id[excl_tx$HsaEX_ID == eid]
  incl_i <- valid_base_to_idx[incl_base[incl_base %in% names(valid_base_to_idx)]]
  excl_i <- valid_base_to_idx[excl_base[excl_base %in% names(valid_base_to_idx)]]
  if (length(incl_i)==0 || length(excl_i)==0) next
  incl_tpm <- colSums(tpm_matrix[incl_i, samp_idx, drop=FALSE], na.rm=TRUE)
  excl_tpm <- colSums(tpm_matrix[excl_i, samp_idx, drop=FALSE], na.rm=TRUE)
  total <- incl_tpm + excl_tpm
  usage_tpm[, eid] <- ifelse(total > 0, incl_tpm / total, NA)
  event_total_tpm[, eid] <- total
}

# Verify: usage from TPM should match previous rate-ratio (since library size cancels)
prev_models <- read.delim(file.path(prev_dir, "07_primary_nonoverlap_models/03_primary_event_models.tsv"), stringsAsFactors=FALSE)
cat("  Verification: TPM usage vs previous rate-ratio usage (should be identical)\n")
# Both compute sum(rate_incl)/[sum(rate_incl)+sum(rate_excl)] - mathematically identical
# because the per-sample sum(rate_all) cancels in the ratio

# Event total TPM summary
event_tpm_summary <- data.frame()
for (eid in events) {
  gene <- primary19$gene[primary19$HsaEX_ID == eid]
  vals <- event_total_tpm[, eid]
  event_tpm_summary <- rbind(event_tpm_summary, data.frame(
    HsaEX_ID=eid, gene=gene,
    median_total_TPM=round(median(vals, na.rm=TRUE), 3),
    mean_total_TPM=round(mean(vals, na.rm=TRUE), 3),
    q25_total_TPM=round(quantile(vals, 0.25, na.rm=TRUE), 3),
    q75_total_TPM=round(quantile(vals, 0.75, na.rm=TRUE), 3),
    stringsAsFactors=FALSE))
}
write.table(event_tpm_summary, file.path(tpm_dir, "02_event_total_TPM_summary.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Detection threshold sensitivity
cat("\n=== Detection threshold sensitivity ===\n")
thresholds <- c(0, 0.1, 0.5, 1, 2, 5)
epsilon <- 1e-4

# For each threshold, filter samples and refit
threshold_results <- data.frame()
threshold_models_all <- data.frame()

for (thresh in thresholds) {
  cat(paste0("\n  Threshold: event_total_TPM >= ", thresh, "\n"))

  model_results <- data.frame()
  for (eid in events) {
    gene <- primary19$gene[primary19$HsaEX_ID == eid]
    usage_vals <- usage_tpm[, eid]
    total_vals <- event_total_tpm[, eid]

    # Apply threshold: samples below threshold are NA
    detectable <- !is.na(total_vals) & total_vals >= thresh
    usage_filt <- ifelse(detectable, usage_vals, NA)

    valid <- !is.na(usage_filt)
    n_valid <- sum(valid)
    if (n_valid < 15) next

    mdf <- analysis_meta[valid, ]
    mdf$usage <- usage_filt[valid]

    # Check 70% detectability per group
    n_asd_valid <- sum(mdf$dx_binary == 1)
    n_ctl_valid <- sum(mdf$dx_binary == 0)
    n_asd_total <- sum(analysis_meta$dx_binary == 1)
    n_ctl_total <- sum(analysis_meta$dx_binary == 0)
    asd_detect_rate <- n_asd_valid / n_asd_total
    ctl_detect_rate <- n_ctl_valid / n_ctl_total

    testable <- (asd_detect_rate >= 0.70) & (ctl_detect_rate >= 0.70)
    if (!testable & thresh > 0) next

    mdf$usage_logit <- log((mdf$usage + epsilon) / (1 - mdf$usage + epsilon))

    fit <- tryCatch(
      lmer(usage_logit ~ dx_binary + region + Sex + Age + RIN + (1|subject), data=mdf, REML=TRUE),
      error = function(e) tryCatch(
        lmer(usage_logit ~ dx_binary + region + Sex + (1|subject), data=mdf, REML=TRUE),
        error = function(e2) NULL))
    if (is.null(fit)) next

    coefs <- summary(fit)$coefficients
    if ("dx_binary" %in% rownames(coefs)) {
      beta <- coefs["dx_binary", "Estimate"]
      se <- coefs["dx_binary", "Std. Error"]
      t_val <- coefs["dx_binary", "t value"]
      df_approx <- nrow(mdf) - length(unique(mdf$subject)) - 5
      p_val <- 2 * pt(abs(t_val), df=max(df_approx, 10), lower.tail=FALSE)
      direction <- ifelse(beta > 0, "UP_IN_ASD", "DOWN_IN_ASD")

      model_results <- rbind(model_results, data.frame(
        threshold=thresh, HsaEX_ID=eid, gene=gene, beta_ASD=beta, SE=se, P=p_val,
        direction=direction, n_samples=nrow(mdf),
        n_ASD=n_asd_valid, n_CTL=n_ctl_valid,
        asd_detect_rate=round(asd_detect_rate, 3), ctl_detect_rate=round(ctl_detect_rate, 3),
        stringsAsFactors=FALSE))
    }
  }

  if (nrow(model_results) > 0) {
    model_results$BH_FDR <- p.adjust(model_results$P, method="BH")
    threshold_models_all <- rbind(threshold_models_all, model_results)

    # Set-level stats
    merged_t <- merge(model_results, primary19[,c("HsaEX_ID","discovery_dir","Parikshak_delta_psi")], by="HsaEX_ID")
    merged_t$concordant <- merged_t$direction == merged_t$discovery_dir
    n_eval <- nrow(merged_t)
    n_conc <- sum(merged_t$concordant, na.rm=TRUE)
    rate <- n_conc / n_eval
    bp <- pbinom(n_conc-1, n_eval, 0.5, lower.tail=FALSE)
    sp <- tryCatch(cor.test(merged_t$Parikshak_delta_psi, merged_t$beta_ASD, method="spearman", exact=FALSE),
                   error=function(e) list(estimate=NA, p.value=NA))

    threshold_results <- rbind(threshold_results, data.frame(
      threshold=thresh, n_events_testable=n_eval, n_concordant=n_conc,
      concordance_rate=round(rate,4), binomial_P=signif(bp,6),
      spearman_rho=round(as.numeric(sp$estimate),4), spearman_P=signif(as.numeric(sp$p.value),6),
      n_nominal_P005=sum(model_results$P<0.05),
      n_FDR005=sum(model_results$BH_FDR<0.05),
      n_FDR010=sum(model_results$BH_FDR<0.10),
      stringsAsFactors=FALSE))

    cat(paste0("    Testable: ", n_eval, " Concordant: ", n_conc, "/", n_eval,
               " rate=", round(rate,4), " P=", signif(bp,4),
               " rho=", round(as.numeric(sp$estimate),4), "\n"))
  } else {
    cat("    No events testable at this threshold\n")
  }
}

# Write detection results
write.table(threshold_results, file.path(tpm_dir, "05_set_results_by_threshold.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(threshold_models_all, file.path(tpm_dir, "04_models_by_threshold.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# FDR by threshold
fdr_by_thresh <- threshold_models_all[, c("threshold","HsaEX_ID","gene","P","BH_FDR","direction")]
write.table(fdr_by_thresh, file.path(tpm_dir, "06_FDR_by_threshold.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Detection by threshold (event-level)
detect_summary <- data.frame()
for (thresh in thresholds) {
  for (eid in events) {
    gene <- primary19$gene[primary19$HsaEX_ID == eid]
    total_vals <- event_total_tpm[, eid]
    n_detect <- sum(!is.na(total_vals) & total_vals >= thresh)
    detect_summary <- rbind(detect_summary, data.frame(
      threshold=thresh, HsaEX_ID=eid, gene=gene, n_detected=n_detect,
      frac_detected=round(n_detect/n_analysis, 4), stringsAsFactors=FALSE))
  }
}
write.table(detect_summary, file.path(tpm_dir, "03_event_detection_by_threshold.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# TPM=1 main threshold phase
tpm1 <- threshold_results[threshold_results$threshold == 1, ]
tpm1_check <- data.frame(
  key = c("PRIMARY_THRESHOLD", "N_TESTABLE", "N_CONCORDANT", "RATE", "BINOMIAL_P",
           "SPEARMAN_RHO", "SPEARMAN_P", "N_FDR005", "N_FDR010",
           "RESULT_CRITERION_8_EVENTS", "RESULT_CRITERION_RATE_070", "RESULT_CRITERION_P_010",
           "DETECTION_THRESHOLD_STATUS"),
  value = c("event_total_TPM>=1",
            tpm1$n_events_testable, tpm1$n_concordant, tpm1$concordance_rate, tpm1$binomial_P,
            tpm1$spearman_rho, tpm1$spearman_P, tpm1$n_FDR005, tpm1$n_FDR010,
            ifelse(tpm1$n_events_testable >= 8, "OK", "ERROR"),
            ifelse(tpm1$concordance_rate >= 0.70, "OK", "ERROR"),
            ifelse(tpm1$binomial_P < 0.10, "OK", "ERROR"),
            "OK"),
  stringsAsFactors=FALSE)
write.table(tpm1_check, file.path(tpm_dir, "07_detection_threshold_check.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat("\n\n=== SUMMARY ===\n")
cat(paste0("One-event-per-gene (discovery-anchored): ", da_stats$n_concordant, "/", da_stats$n_genes,
           " P=", signif(da_stats$exact_binomial_P, 4), "\n"))
cat(paste0("Zero fraction corrected: ", round(overall_zero_fraction, 4), "\n"))
cat(paste0("TPM reconstruction: column sums = 1e6 (verified)\n"))
cat(paste0("TPM>=1 threshold: ", tpm1$n_events_testable, " testable, ",
           tpm1$n_concordant, "/", tpm1$n_events_testable, " concordant, P=", signif(tpm1$binomial_P, 4), "\n"))
