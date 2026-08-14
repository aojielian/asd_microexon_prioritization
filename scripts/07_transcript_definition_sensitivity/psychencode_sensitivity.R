#!/usr/bin/env Rscript
# PsychENCODE defensive sensitivity.
# Exact M0/M4 reproduction checks, technical-covariate sensitivity, D0-D3 transcript-set sensitivity and region sensitivity.
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).
suppressMessages({library(methods); library(lme4); library(pbkrtest)})
set.seed(42)

project <- Sys.getenv("PROJECT_ROOT", unset = ".")
task    <- file.path(project,
  "32_psychencode_sensitivity")
outS    <- file.path(task, "02_psychencode_sensitivity")
outDsub <- file.path(outS, "transcript_set_intermediates")
dir.create(outS, showWarnings=FALSE, recursive=TRUE)
dir.create(outDsub, showWarnings=FALSE, recursive=TRUE)
data_dir <- file.path(project, "psychencode_processed")
dir21 <- file.path(project, "21_coordinate_inference")
dir18 <- file.path(project, "18_psychencode")
dir25 <- file.path(project, "25_master_evidence")
gencode_dir <- file.path(Sys.getenv("REFERENCE_ROOT", unset = "."), "02_coordinate_harmonization")
logf <- file.path(task, "99_logs", "psychencode_sensitivity.log")

say <- function(...) { msg <- paste0(...); cat(msg, "\n");
                       cat(msg, "\n", file=logf, append=TRUE) }
file.create(logf); cat("", file=logf)
wtab <- function(df, fn, dir=outS)
  write.table(df, file.path(dir, fn), sep="\t", row.names=FALSE, quote=FALSE)
t0_all <- Sys.time()

# ============================================================================
# PART A — data loading, exclusions, usage matrix, cache
# ============================================================================
say("=== PART A: load raw + processed PsychENCODE RData ===")
env_raw <- new.env()
load(file.path(data_dir, "01_02_B_01_RawData.RData"), envir=env_raw)
rsem_tx <- env_raw$rsem_tx
effLen  <- env_raw$rsem_transcript_effLen
say("rsem_tx: ", nrow(rsem_tx), " transcripts x ", ncol(rsem_tx), " samples")

env_proc <- new.env()
load(file.path(data_dir, "02_01_B_AllProcessedData_wModelMatrix.RData"),
     envir=env_proc)
datMeta  <- env_proc$datMeta
datMetaM <- env_proc$datMeta_model
topPC    <- env_proc$topPC
datSeqN  <- env_proc$datSeq_numeric
rm(env_proc); gc()

datMeta$Diagnosis <- as.character(datMeta$Diagnosis)
datMeta$Sex       <- as.character(datMeta$Sex)
datMeta$region    <- as.character(datMeta$region)
datMeta$subject   <- as.character(datMeta$subject)

# Reference donor exclusions: 23 GSE30573-overlap donors (final check);
# Dup15q donors excluded by the Diagnosis %in% ASD/CTL filter exactly as in
# the final analysis.
# Donor-level exclusion list (GSE30573-overlap donors) is cohort metadata
# supplied at runtime as a plain-text file (one identifier per line) under
# data_dir; the list itself is not shipped in this public repository.
gse_excl_file <- file.path(data_dir, "gse_overlap_donor_exclusion.txt")
if (file.exists(gse_excl_file)) {
  gse_donors <- readLines(gse_excl_file)
  gse_donors <- gsub("\\s+", "", gse_donors[nzchar(gse_donors)])
  stopifnot(length(gse_donors) == 23)  # 23 overlap donors (verified at release)
} else {
  warning("gse_overlap_donor_exclusion.txt not found under data_dir; ",
          "no GSE-overlap donor exclusion applied")
  gse_donors <- character(0)
}
analysis_meta <- datMeta[!(datMeta$subject %in% gse_donors) &
                         datMeta$Diagnosis %in% c("ASD","CTL"), ]
analysis_meta$dx_binary <- ifelse(analysis_meta$Diagnosis == "ASD", 1, 0)
analysis_samples <- rownames(analysis_meta)
n_don_total <- length(unique(analysis_meta$subject))
n_don_asd   <- length(unique(analysis_meta$subject[analysis_meta$dx_binary==1]))
n_don_ctl   <- length(unique(analysis_meta$subject[analysis_meta$dx_binary==0]))
say("Analysis subset: ", nrow(analysis_meta), " samples, ",
    n_don_total, " donors (", n_don_asd, " ASD, ", n_don_ctl, " CTL)")

# Reference event definitions (locked, identical to final usage)
primary19 <- read.delim(file.path(project,
  "16_gse30573/02_input_lock/02_primary19.tsv"),
  stringsAsFactors=FALSE)
primary19$discovery_dir <- ifelse(primary19$Parikshak_delta_psi > 0,
                               "UP_IN_ASD", "DOWN_IN_ASD")
primary19$abs_delta_psi <- abs(primary19$Parikshak_delta_psi)

# Authoritative final evidence tiers (master table, read-only)
master25 <- read.delim(file.path(dir25,
  "06_master_event_table/MASTER_19_EVENT_EVIDENCE_TABLE.tsv"),
  stringsAsFactors=FALSE)
tier_map <- setNames(master25$final_evidence_tier, master25$HsaEX_ID)
primary19$final_tier <- tier_map[primary19$HsaEX_ID]
stopifnot(sum(primary19$final_tier == "TIER_A_CROSS_COHORT_KR_FDR05") == 4)

incl_tx <- read.delim(file.path(dir18,
  "07_annotation_and_transcript_mapping/03_inclusion_transcripts.tsv"),
  stringsAsFactors=FALSE)
excl_tx <- read.delim(file.path(dir18,
  "07_annotation_and_transcript_mapping/04_exclusion_transcripts.tsv"),
  stringsAsFactors=FALSE)

# ---------------- usage matrix (exactly as final) -------------------------
tx_ids <- rownames(rsem_tx)
base_ids <- sub("\\.[0-9]+_[0-9]+$", "", tx_ids)
base_to_idx <- setNames(seq_along(tx_ids), base_ids)
events <- sort(unique(incl_tx$HsaEX_ID))
stopifnot(length(events) == 19)

n_analysis <- length(analysis_samples)
samp_idx_raw <- match(analysis_samples, colnames(rsem_tx))
say("Samples matched into rsem_tx: ", sum(!is.na(samp_idx_raw)), "/",
    n_analysis)
stopifnot(all(!is.na(samp_idx_raw)))

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

# Cache for later analyses (module 6 gene-level host-expression reuse etc.)
all_member_tx <- unique(c(incl_tx$transcript_id, excl_tx$transcript_id))
member_i <- base_to_idx[all_member_tx[all_member_tx %in% names(base_to_idx)]]
rsem_cache <- list(
  usage_matrix=usage_matrix,
  analysis_meta=analysis_meta,
  primary19=primary19,
  incl_tx=incl_tx, excl_tx=excl_tx,
  rsem_member=rsem_tx[member_i, , drop=FALSE],
  effLen_member=effLen[member_i],
  member_base_ids=base_ids[member_i],
  samp_idx_raw=samp_idx_raw,
  datMeta_model=datMetaM, topPC=topPC,
  gse_donors=gse_donors)
saveRDS(rsem_cache, file.path(task, "00_admin", "rsem_cache.rds"))
say("rsem cache saved to 00_admin/rsem_cache.rds")

# ============================================================================
# PART B1 — M0 reproduction
# ============================================================================
say("\n=== PART B1: M0 reproduction ===")
epsilon <- 1e-4

fit_event <- function(mdf, formula_full, formula_red, reml=TRUE) {
  fit_full <- tryCatch(lmer(formula_full, data=mdf, REML=reml),
                       error=function(e) NULL)
  if (is.null(fit_full)) return(NULL)
  fit_red <- tryCatch(lmer(formula_red, data=mdf, REML=reml),
                      error=function(e) NULL)
  list(full=fit_full, red=fit_red)
}

run_model <- function(eid, extra_fix, cov_df=NULL, allow_fallback=FALSE) {
  usage_vals <- usage_matrix[, eid]
  valid <- !is.na(usage_vals)
  mdf <- analysis_meta[valid, ]
  if (!is.null(cov_df)) {
    extra <- cov_df[match(rownames(mdf), rownames(cov_df)), , drop=FALSE]
    mdf <- cbind(mdf, extra)
  }
  mdf$usage <- usage_vals[valid]
  mdf$usage_logit <- log((mdf$usage + epsilon) / (1 - mdf$usage + epsilon))

  base_full <- "usage_logit ~ dx_binary + region + Sex + Age + RIN"
  base_red  <- "usage_logit ~ region + Sex + Age + RIN"
  f_full <- as.formula(paste(base_full, extra_fix, "+ (1|subject)"))
  f_red  <- as.formula(paste(base_red,  extra_fix, "+ (1|subject)"))

  # Reference M0 fallback: on error drop Age+RIN. Used for M0 only.
  fit_reml <- fit_event(mdf, f_full, f_red, reml=TRUE)
  used_fallback <- FALSE
  if (is.null(fit_reml$full) && allow_fallback) {
    f_full <- as.formula(paste("usage_logit ~ dx_binary + region + Sex",
                               extra_fix, "+ (1|subject)"))
    f_red  <- as.formula(paste("usage_logit ~ region + Sex",
                               extra_fix, "+ (1|subject)"))
    fit_reml <- fit_event(mdf, f_full, f_red, reml=TRUE)
    used_fallback <- TRUE
  }
  if (is.null(fit_reml$full)) return(NULL)
  coefs <- summary(fit_reml$full)$coefficients
  if (!("dx_binary" %in% rownames(coefs))) return(NULL)
  beta <- coefs["dx_binary", "Estimate"]
  se   <- coefs["dx_binary", "Std. Error"]

  kr_F <- NA; kr_df1 <- NA; kr_df2 <- NA; p_kr <- NA
  if (!is.null(fit_reml$red)) {
    kr <- tryCatch(KRmodcomp(fit_reml$full, fit_reml$red),
                   error=function(e) NULL)
    if (!is.null(kr)) {
      kr_F <- kr$test$stat[1]; kr_df1 <- kr$test$ndf[1]
      kr_df2 <- kr$test$ddf[1]; p_kr <- kr$test$p.value[1]
    }
  }
  fit_ml <- fit_event(mdf, f_full, f_red, reml=FALSE)
  lrt_chisq <- NA; p_lrt <- NA
  if (!is.null(fit_ml$full) && !is.null(fit_ml$red)) {
    lrt_chisq <- 2 * (as.numeric(logLik(fit_ml$full)) -
                      as.numeric(logLik(fit_ml$red)))
    p_lrt <- pchisq(lrt_chisq, df=1, lower.tail=FALSE)
  }
  X <- model.matrix(fit_reml$full)
  data.frame(HsaEX_ID=eid, beta_ASD=beta, SE=se,
             CI_lo=beta - 1.96*se, CI_hi=beta + 1.96*se,
             KR_F=kr_F, KR_df1=kr_df1, KR_df2=kr_df2,
             P_Kenward_Roger=p_kr, LRT_chisq=lrt_chisq, P_LRT=p_lrt,
             direction=ifelse(beta>0, "UP_IN_ASD","DOWN_IN_ASD"),
             n_samples=nrow(mdf), n_donors=length(unique(mdf$subject)),
             n_ASD_donors=length(unique(mdf$subject[mdf$dx_binary==1])),
             n_CTL_donors=length(unique(mdf$subject[mdf$dx_binary==0])),
             singular_fit=isSingular(fit_reml$full),
             fixed_effect_rank=qr(X)$rank,
             condition_number=kappa(X),
             used_fallback=used_fallback,
             formula=paste(deparse(f_full, width.cutoff=500L), collapse=" "),
             stringsAsFactors=FALSE)
}

run_all_events <- function(model_name, extra_fix, cov_df=NULL,
                           allow_fallback=FALSE) {
  res <- data.frame()
  n_error <- 0
  for (eid in events) {
    r <- tryCatch(run_model(eid, extra_fix, cov_df, allow_fallback),
                  error=function(e) {say("    ERROR ", eid, ": ", e$message);
                                     NULL})
    if (!is.null(r)) { r$model <- model_name; res <- rbind(res, r) }
    else n_error <- n_error + 1
    say("  ", model_name, " ", eid, " done")
  }
  if (n_error > 0)
    say("  WARNING: ", n_error, "/19 events returned NULL for ", model_name)
  res$BH_FDR_KR  <- p.adjust(res$P_Kenward_Roger, method="BH")
  res$BH_FDR_LRT <- p.adjust(res$P_LRT, method="BH")
  res
}

set_summary <- function(res, model_name) {
  merged <- merge(res, primary19[, c("HsaEX_ID","gene","discovery_dir",
                                  "Parikshak_delta_psi","abs_delta_psi",
                                  "final_tier")], by="HsaEX_ID")
  merged$concordant <- merged$direction == merged$discovery_dir
  n_eval <- nrow(merged); n_conc <- sum(merged$concordant, na.rm=TRUE)
  binom_p <- if (n_conc > 0) pbinom(n_conc-1, n_eval, 0.5, lower.tail=FALSE)
             else 1
  sp <- cor.test(merged$Parikshak_delta_psi, merged$beta_ASD,
                 method="spearman", exact=FALSE)
  ms <- merged[order(-merged$abs_delta_psi, merged$HsaEX_ID), ]
  oepg <- do.call(rbind, lapply(split(ms, ms$gene), function(d) d[1,]))
  oepg_conc <- sum(oepg$concordant, na.rm=TRUE)
  oepg_p <- if (oepg_conc > 0)
              pbinom(oepg_conc-1, nrow(oepg), 0.5, lower.tail=FALSE) else 1
  tierA <- merged[merged$final_tier == "TIER_A_CROSS_COHORT_KR_FDR05", ]
  data.frame(model=model_name, n_events=n_eval, n_concordant=n_conc,
    concordance_rate=round(n_conc/n_eval, 4),
    exact_binomial_P=signif(binom_p, 6),
    spearman_rho=round(sp$estimate, 4), spearman_P=signif(sp$p.value, 6),
    oepg_n=nrow(oepg), oepg_concordant=oepg_conc,
    oepg_P=signif(oepg_p, 6),
    tierA_n=nrow(tierA),
    tierA_direction_flags=paste(tierA$concordant, collapse=""),
    tierA_n_dir_conc=sum(tierA$concordant, na.rm=TRUE),
    tierA_KR_nominal_005=sum(tierA$P_Kenward_Roger < 0.05, na.rm=TRUE),
    tierA_KR_FDR_005=sum(tierA$BH_FDR_KR < 0.05, na.rm=TRUE),
    n_KR_FDR_005=sum(res$BH_FDR_KR < 0.05, na.rm=TRUE),
    n_KR_FDR_010=sum(res$BH_FDR_KR < 0.10, na.rm=TRUE),
    n_LRT_FDR_005=sum(res$BH_FDR_LRT < 0.05, na.rm=TRUE),
    n_LRT_FDR_010=sum(res$BH_FDR_LRT < 0.10, na.rm=TRUE),
    stringsAsFactors=FALSE)
}

M0 <- run_all_events("M0_primary", "", allow_fallback=TRUE)
wtab(M0, "00_reproduction.tsv")

# Reference tables (read-only) for numeric comparison
all <- read.delim(file.path(dir21,
  "05_mixed_model_inference/01_Satterthwaite_models.tsv"),
  stringsAsFactors=FALSE)
set <- read.delim(file.path(dir21,
  "05_mixed_model_inference/06_set_validation_recomputed.tsv"),
  stringsAsFactors=FALSE)
fset <- setNames(set$value, set$key)

cmp <- merge(M0[, c("HsaEX_ID","beta_ASD","SE","P_Kenward_Roger","P_LRT",
                    "direction","n_samples","n_donors")],
             all[, c("HsaEX_ID","beta_ASD","SE","P_Kenward_Roger",
                            "P_LRT","direction","n_samples","n_donors")],
             by="HsaEX_ID", suffixes=c("_repro","_reference"))
cmp$beta_abs_diff <- abs(cmp$beta_ASD_repro - cmp$beta_ASD_reference)
cmp$KR_p_rel_diff <- abs(cmp$P_Kenward_Roger_repro -
                         cmp$P_Kenward_Roger_reference) /
                    pmax(cmp$P_Kenward_Roger_reference, 1e-300)
cmp$direction_match <- cmp$direction_repro == cmp$direction_reference
wtab(cmp, "00_reproduction_comparison.tsv")

M0set <- set_summary(M0, "M0_primary")

phase <- data.frame(check=character(), expected=character(),
                   reproduced=character(), status=character(),
                   stringsAsFactors=FALSE)
add_check <- function(check, expected, reproduced, ok) {
  phase <<- rbind(phase, data.frame(check=check, expected=as.character(expected),
    reproduced=as.character(reproduced),
    status=ifelse(ok, "OK", "ERROR"), stringsAsFactors=FALSE))
}
add_check("N_ANALYSIS_SAMPLES", 532, nrow(analysis_meta),
         nrow(analysis_meta) == 532)
add_check("N_ANALYSIS_DONORS", 80, n_don_total, n_don_total == 80)
add_check("N_ASD_DONORS", 38, n_don_asd, n_don_asd == 38)
add_check("N_CTL_DONORS", 42, n_don_ctl, n_don_ctl == 42)
add_check("N_EVENTS_FITTED", 19, nrow(M0), nrow(M0) == 19)
add_check("DIRECTION_MATCH_ALL_EVENTS", "19/19",
         paste0(sum(cmp$direction_match), "/19"), all(cmp$direction_match))
add_check("MAX_BETA_ABS_DIFF", "<1e-6", signif(max(cmp$beta_abs_diff), 3),
         max(cmp$beta_abs_diff) < 1e-6)
add_check("MAX_KR_P_REL_DIFF", "<1e-6", signif(max(cmp$KR_p_rel_diff), 3),
         max(cmp$KR_p_rel_diff) < 1e-6)
add_check("CONCORDANCE", "15/19",
         paste0(M0set$n_concordant, "/", M0set$n_events),
         M0set$n_concordant == 15 & M0set$n_events == 19)
add_check("EXACT_BINOMIAL_P", fset["EXACT_BINOMIAL_P"],
         M0set$exact_binomial_P,
         abs(M0set$exact_binomial_P - as.numeric(fset["EXACT_BINOMIAL_P"])) <
         1e-6)
add_check("OEPG_CONCORDANT", fset["DISCOVERY_ANCHORED_OEPG_CONCORDANT"],
         M0set$oepg_concordant,
         M0set$oepg_concordant ==
         as.numeric(fset["DISCOVERY_ANCHORED_OEPG_CONCORDANT"]))
add_check("KR_FDR_005_COUNT", 4, M0set$n_KR_FDR_005, M0set$n_KR_FDR_005 == 4)
add_check("KR_FDR_010_COUNT", 7, M0set$n_KR_FDR_010, M0set$n_KR_FDR_010 == 7)
add_check("LRT_FDR_005_COUNT", 6, M0set$n_LRT_FDR_005,
         M0set$n_LRT_FDR_005 == 6)
add_check("LRT_FDR_010_COUNT", 10, M0set$n_LRT_FDR_010,
         M0set$n_LRT_FDR_010 == 10)
check_overall <- ifelse(all(checks$status == "OK"), "OK", "ERROR")
phase <- rbind(phase, data.frame(check="PRIMARY_MODEL_REPRODUCED",
  expected="OK", reproduced=check_overall, status=check_overall,
  stringsAsFactors=FALSE))
wtab(phase, "00_reproduction_check.txt")
say("Reproduction phase: ", check_overall)

if (check_overall != "OK") {
  say("STOP: final model not reproduced. Downstream modeling halted.")
  writeLines(c("PRIMARY_MODEL_REPRODUCED=ERROR",
    "SENSITIVITY_STATUS=REVISE_REQUIRED_MODEL_NOT_REPRODUCED"),
    file.path(outS, "00_reproduction_check.txt"))
  quit(save="no", status=0)
}

# ============================================================================
# PART B2 — technical-covariate field dictionary + feasibility
# ============================================================================
say("\n=== PART B2: field dictionary ===")
dict <- data.frame(source_file=character(), object_name=character(),
  field_name=character(), class=character(), n_nonmissing=numeric(),
  n_unique=numeric(), example_values=character(), candidate_role=character(),
  documented_meaning=character(), safe_to_model=character(), reason=character(),
  stringsAsFactors=FALSE)
F_PROC <- "02_01_B_AllProcessedData_wModelMatrix.RData"
add <- function(sf, on, df, field, role, meaning, safe, why) {
  if (is.matrix(df) || is.array(df)) df <- as.data.frame(df)
  if (!field %in% colnames(df)) {
    say("MISSING FIELD: ", sf, " ", on, " ", field); return(invisible(NULL))
  }
  x <- df[[field]]
  ex <- utils::head(unique(as.character(x[!is.na(x)])), 3)
  dict <<- rbind(dict, data.frame(source_file=sf, object_name=on,
    field_name=field, class=paste(class(x), collapse="/"),
    n_nonmissing=sum(!is.na(x)), n_unique=length(unique(x)),
    example_values=paste(ex, collapse=" | "),
    candidate_role=role, documented_meaning=meaning,
    safe_to_model=safe, reason=why, stringsAsFactors=FALSE))
}

cur_dm <- list(
  sample_id=c("sample_id","Sample identifier (consortium)","NO","identifier; not a covariate"),
  subject=c("subject_id","Donor identifier; random-effect grouping in the final model","YES","documented random effect (1|subject) in the final and all sensitivity models"),
  region=c("region","Cortical region / Brodmann area; fixed effect in the final model","YES","fixed effect in the final model and all sensitivity models"),
  seq_batch=c("sequencing_batch","Sequencing batch/run identifier (e.g. 2013-222)","NO","partially confounded with region in the consortium design; check cross-tabs only, never a model term"),
  Brain_Bank_Source=c("site","Brain bank / tissue source","NO","site proxy; near-identical to subject; check cross-tabs only"),
  Diagnosis=c("diagnosis","Clinical diagnosis (ASD/CTL/other); primary exposure","YES","primary exposure of the final model; analysis restricted to ASD/CTL"),
  Details_on_Diagnosis=c("other","Free-text diagnosis details","NO","free text; not modelable"),
  Sex=c("sex","Sex (M/F); fixed covariate in the final model","YES","fixed covariate in the final model"),
  Age=c("age","Age at death (years); fixed covariate in the final model","YES","fixed covariate in the final model; also present as consortium model-matrix field"),
  Brain_Weight=c("other","Brain weight (grams); post-mortem descriptor","NO","post-mortem descriptor; not in consortium model matrix nor the final model"),
  pH=c("other","Tissue pH; post-mortem quality descriptor","NO","not in consortium model matrix nor the final model"),
  Previously_reported_RIN_CTX=c("other","Previously reported RIN, cortex (per-source value)","NO","superseded by harmonized RIN field used by the final model"),
  Previously_reported_RIN_CBL=c("other","Previously reported RIN, cerebellum","NO","cerebellar value; not applicable to cortical model"),
  PMI=c("PMI","Post-mortem interval; consortium model-matrix field (datMeta_model$PMI)","YES_M4_ONLY","documented consortium model-matrix covariate; pre-specified M4 sensitivity covariate; not used in M0-M3"),
  Primary_Cause_of_Death=c("other","Free-text cause of death","NO","free text"),
  Secondary_Cause_of_Death=c("other","Free-text secondary cause of death","NO","free text"),
  Agonal_State=c("other","Agonal state at death","NO","not in consortium model matrix nor the final model"),
  Seizures=c("other","Seizure history","NO","not in consortium model matrix nor the final model"),
  Seizure_notes=c("other","Free-text seizure notes","NO","free text"),
  Psychiatric_Medications=c("other","Psychiatric medication history","NO","not in consortium model matrix nor the final model"),
  Medication_notes=c("other","Free-text medication notes","NO","free text"),
  Comorbidity_notes=c("other","Free-text comorbidity notes","NO","free text"),
  ADI.R_A=c("other","ADI-R domain A score (clinical)","NO","clinical descriptor; not a confounder covariate"),
  ADI.R_B_NV=c("other","ADI-R domain B non-verbal score","NO","clinical descriptor"),
  ADI.R_B_V=c("other","ADI-R domain B verbal score","NO","clinical descriptor"),
  ADI.R_C_=c("other","ADI-R domain C score","NO","clinical descriptor"),
  ADI.R_D_=c("other","ADI-R domain D score","NO","clinical descriptor"),
  IQ=c("other","IQ measure","NO","clinical descriptor"),
  IQ_notes=c("other","Free-text IQ notes","NO","free text"),
  ethnicity_raw=c("other","Self-reported ethnicity (raw)","NO","self-report; genotype-based Ancestry_Genotype is the documented ancestry variable"),
  race=c("other","Self-reported race","NO","self-report; not modeled; see Ancestry_Genotype"),
  handed=c("other","Handedness","NO","not a documented confounder for this analysis"),
  batch=c("library_batch","Library/preparation batch (ASD3Reg/ASDPan); verified identical to datMeta_model$SeqBatch","YES","documented consortium model-matrix batch field (datMeta_model$SeqBatch); pre-specified M2 sensitivity covariate"),
  SeqMethod=c("other","Sequencing method/platform descriptor","NO","platform descriptor; not in pre-specified hierarchy"),
  lobe=c("other","Anatomical lobe descriptor","NO","anatomical descriptor; region is the modeled field"),
  Ancestry_Genotype=c("ancestry","Genotype-based genetic ancestry; consortium model-matrix field (datMeta_model$Ancestry)","YES","documented consortium model-matrix covariate; pre-specified M3 sensitivity covariate"),
  Server_ID_Pancortical=c("other","Processing server identifier (pancortical)","NO","computational provenance only"),
  RIN=c("RIN","RNA integrity number (harmonized); fixed covariate in the final model","YES","fixed covariate in the final model; consortium model-matrix field"),
  Server_ID=c("other","Processing server identifier","NO","computational provenance only"),
  Read_Length=c("other","Sequencing read length","NO","platform descriptor; not in pre-specified hierarchy"),
  Dx=c("diagnosis","Diagnosis factor coded for modeling (ASD/CTL)","YES","model coding of Diagnosis in consortium object and the final analysis"))
for (f in names(cur_dm)) add(F_PROC, "datMeta", datMeta, f, cur_dm[[f]][1],
  cur_dm[[f]][2], cur_dm[[f]][3], cur_dm[[f]][4])

cur_dmM <- list(
  Subject=c("subject_id","Donor identifier (consortium model matrix)","YES","model-matrix grouping variable"),
  DxReg=c("diagnosis","Diagnosis coding used for region-stratified consortium models","YES","consortium model-matrix diagnosis field"),
  SeqBatch=c("library_batch","Library/preparation batch; verified identical to datMeta$batch (ASD3Reg/ASDPan)","YES","consortium model-matrix batch field; pre-specified M2 covariate"),
  Sex=c("sex","Sex (consortium model matrix)","YES","consortium model-matrix covariate"),
  Ancestry=c("ancestry","Genotype-based ancestry (consortium model matrix)","YES","consortium model-matrix covariate; pre-specified M3 covariate"),
  Age=c("age","Age (consortium model matrix)","YES","consortium model-matrix covariate"),
  Age_sqd=c("age_squared","Squared age term (consortium model matrix)","YES","consortium-documented age-squared field; supports pre-specified M1"),
  PMI=c("PMI","Post-mortem interval (consortium model matrix)","YES_M4_ONLY","documented consortium model-matrix covariate; pre-specified M4 sensitivity covariate"),
  RIN=c("RIN","RIN (consortium model matrix)","YES","consortium model-matrix covariate; final-model covariate"))
for (f in names(cur_dmM)) add(F_PROC, "datMeta_model", datMetaM, f,
  cur_dmM[[f]][1], cur_dmM[[f]][2], cur_dmM[[f]][3], cur_dmM[[f]][4])
for (f in setdiff(colnames(datMetaM), names(cur_dmM)))
  add(F_PROC, "datMeta_model", datMetaM, f, "technical_qc_metric",
      "Picard/STAR sequencing-QC metric selected into the consortium model-matrix object",
      "YES_M4_ONLY",
      "documented consortium model-matrix QC covariate; pre-specified M4 sensitivity covariate; never used in M0-M3 or exploratory arms")

for (f in colnames(topPC))
  add(F_PROC, "topPC", topPC, f, "technical_PC",
      "Top expression principal component; variance fraction encoded in column name",
      "EXPLORATORY_ONLY",
      "documented only as exploratory PC~covariate analysis (MOESM3 R2_exprPCs_Covariates, Extended Data Fig 2c), never as a consortium model covariate; expression-derived PCs can absorb diagnosis signal; used only in the labeled exploratory M4x arm")

for (f in colnames(datSeqN))
  add(F_PROC, "datSeq_numeric", datSeqN, f, "other",
      "Numeric Picard/STAR alignment & sequencing QC metric (full set)",
      "NO",
      "sequencing QC metric object; only the 8 metrics carried into datMeta_model are documented model covariates (used in pre-specified M4); remaining fields check only")

e_raw <- new.env()
load(file.path(data_dir, "01_02_B_01_RawData.RData"), envir=e_raw)
rdm <- e_raw$datMeta
for (f in colnames(rdm)) {
  f2 <- if (f == "Diagnosis_") "Diagnosis" else f
  if (f2 %in% names(cur_dm)) {
    r <- cur_dm[[f2]]
    add("01_02_B_01_RawData.RData", "datMeta", rdm, f, r[1], r[2],
        ifelse(r[3] %in% c("YES","YES_M4_ONLY"), "NO", r[3]),
        ifelse(r[3] %in% c("YES","YES_M4_ONLY"),
          paste("raw-object copy;", r[4],
                "authoritative modeling copy is the processed datMeta"),
          r[4]))
  } else add("01_02_B_01_RawData.RData", "datMeta", rdm, f, "other",
             "Raw-metadata field (see processed datMeta)","NO",
             "raw copy; check only")
}
rm(e_raw)
wtab(dict, "PSYCHENCODE_FIELD_DICTIONARY.tsv")
say("field dictionary rows: ", nrow(dict))

# ---- covariate search summary (spec term list) ----------------------------
term_hits <- data.frame(
  search_term=c("batch","library batch","site","center","ancestry",
    "genetic PCs","technical PCs","surrogate variables","PMI","pH","RIN",
    "age","sex","region","cell fractions","methylation deconvolution"),
  documented_field_found=c(
    "datMeta$seq_batch (sequencing run batch; confounded with region; check only)",
    "datMeta$batch == datMeta_model$SeqBatch (library batch; M2 covariate)",
    "datMeta$Brain_Bank_Source (brain bank/site; check only)",
    "NONE",
    "datMeta$Ancestry_Genotype == datMeta_model$Ancestry (M3 covariate)",
    "NONE (no genotype PCs serialized; Ancestry factor is the documented ancestry variable)",
    "topPC PC1-PC15 expression PCs (exploratory only, MOESM3 Extended Data Fig 2c; M4x arm only)",
    "NONE (no SVA/RUV fields in any object)",
    "datMeta$PMI == datMeta_model$PMI (pre-specified M4 covariate)",
    "datMeta$pH (not in consortium model matrix; check only)",
    "datMeta$RIN == datMeta_model$RIN (final-model covariate)",
    "datMeta$Age == datMeta_model$Age; datMeta_model$Age_sqd (M0/M1 covariates)",
    "datMeta$Sex == datMeta_model$Sex (final-model covariate)",
    "datMeta$region (final-model fixed effect)",
    "NONE (no NeuN/cell-fraction/deconvolution field in any object)",
    "NONE (no methylation-deconvolution field in any object)"),
  modeled_where=c("none","M2","none","none","M3","none","M4x_exploratory",
    "none","M4","none","M0-M4","M0 (Age), M1-M4 (+Age^2)","M0-M4","M0-M4",
    "none","none"),
  stringsAsFactors=FALSE)
wtab(term_hits, "PSYCHENCODE_COVARIATE_SEARCH.tsv")

am <- analysis_meta
# datMeta_model has numeric rownames but is ROW-ALIGNED with datMeta
# (verified: Subject row-wise match). Align by row index.
mm <- datMetaM[match(rownames(am), rownames(datMeta)), , drop=FALSE]
tp <- topPC[match(rownames(am), rownames(topPC)), , drop=FALSE]
stopifnot(all(as.character(mm$Subject) == am$subject))
stopifnot(all(as.character(mm$SeqBatch) == as.character(am$batch)))

writeLines(c(
"# PsychENCODE technical-covariate feasibility (module 1)",
"",
"## Documented covariate inventory",
"",
"Fields were checked in all serialized PsychENCODE objects",
"(`datMeta`, `datMeta_model`, `topPC`, `datSeq_numeric`, raw `datMeta`);",
"see PSYCHENCODE_FIELD_DICTIONARY.tsv for the full per-field dictionary and",
"PSYCHENCODE_COVARIATE_SEARCH.tsv for the spec term-list search.",
"",
"- Reference M0 covariates (documented): dx_binary + region + Sex +",
"  Age + RIN + (1|subject).",
"- M1 adds I(Age^2); supported by documented datMeta_model$Age_sqd.",
"- M2 adds library batch; datMeta$batch verified identical to",
"  datMeta_model$SeqBatch (both = ASD3Reg/ASDPan levels).",
"- M3 adds genotype-based ancestry; datMeta_model$Ancestry (levels:",
paste0("  ", paste(names(table(mm$Ancestry)), table(mm$Ancestry),
       sep="=", collapse=", ")),
").",
"- M4 adds the remaining documented consortium model-matrix covariates:",
"  PMI + 8 Picard/STAR sequencing-QC metrics",
"  (PCT_MRNA_BASES, AT_DROPOUT, PCT_UTR_BASES,",
"  multimapped_toomany_percent, MEDIAN_CV_COVERAGE, MEDIAN_INSERT_SIZE,",
"  PCT_INTERGENIC_BASES, PF_BASES).",
"- M4x EXPLORATORY adds topPC PC1-PC5 (variance labels 5.3%/4.6%/3.2%/2%/",
"  1.7%); documented only as exploratory PC~covariate analysis (MOESM3,",
"  Extended Data Fig 2c), never as consortium model covariates.",
"  Expression-derived PCs can absorb diagnosis signal; M4x is labeled",
"  exploratory and is NEVER the preferred model.",
"",
"## Not documented / not modeled",
"",
"- No genetic PCs (genotype PC fields) are serialized; the documented",
"  ancestry variable is the categorical datMeta_model$Ancestry.",
"- No surrogate-variable / RUV / PEER fields exist in any object.",
"- No sequencing-run `seq_batch` model term: partially confounded with",
"  region in the consortium design; check cross-tabs only.",
"- pH, Brain_Weight, agonal state, seizures, medications, ADI-R, IQ:",
"  post-mortem/clinical descriptors, not consortium model covariates.",
"",
"## Missingness constraints (analysis subset, 532 samples)",
"",
paste0("- datMeta_model$PMI non-missing: ", sum(!is.na(mm$PMI)), "/532"),
paste0("- Picard/STAR QC metrics non-missing (all 8): ",
       sum(complete.cases(mm[, grep("^(picard_|star\\.)", colnames(mm))])),
       "/532"),
"Models with missing covariate values drop the affected samples",
"(listwise); donor-loss >15% marks a model HIGH_MISSINGNESS_SENSITIVITY_ONLY",
"and excludes it from the preferred-model rule.",
""), file.path(outS, "PSYCHENCODE_COVARIATE_FEASIBILITY.md"))

# ============================================================================
# PART B3 — technical model hierarchy M1-M4 (+M4x)
# ============================================================================
say("\n=== PART B3: extended model hierarchy ===")
say("datMeta_model aligned (index-based); SeqBatch==batch verified.")
say("Ancestry levels in subset: ",
    paste(names(table(mm$Ancestry)), table(mm$Ancestry), sep="=",
          collapse=", "))

say("\n--- M1 (Age + Age^2) ---")
M1 <- run_all_events("M1_nonlinear_age", "+ I(Age^2)")

say("\n--- M2 (+ SeqBatch) ---")
M2 <- run_all_events("M2_seq_batch", "+ I(Age^2) + SeqBatch",
                     cov_df=data.frame(SeqBatch=mm$SeqBatch,
                                       row.names=rownames(am)))

anc_tab <- table(mm$Ancestry, am$dx_binary)
anc_ok <- sum(rowSums(anc_tab) > 0) >= 2 && all(colSums(anc_tab > 0) >= 2)
say("Ancestry identifiable: ", anc_ok)
if (anc_ok) {
  say("\n--- M3 (+ Ancestry) ---")
  M3 <- run_all_events("M3_ancestry", "+ I(Age^2) + SeqBatch + Ancestry",
                       cov_df=data.frame(SeqBatch=mm$SeqBatch,
                                         Ancestry=mm$Ancestry,
                                         row.names=rownames(am)))
} else { M3 <- NULL }

# M4 = M3 + documented consortium model-matrix covariates (PMI + 8 QC).
# Renamed to syntactic names for formula use.
qc_cols <- grep("^(picard_|star\\.)", colnames(datMetaM), value=TRUE)
stopifnot(length(qc_cols) == 8)
m4_cov <- data.frame(SeqBatch=mm$SeqBatch, Ancestry=mm$Ancestry,
                     PMI_mm=mm$PMI, row.names=rownames(am))
for (j in seq_along(qc_cols))
  m4_cov[[paste0("QC", j)]] <- mm[[qc_cols[j]]]
m4_fix <- paste("+ I(Age^2) + SeqBatch + Ancestry + PMI_mm +",
                paste(paste0("QC", seq_along(qc_cols)), collapse=" + "))
say("\n--- M4 (+ PMI + 8 documented Picard/STAR QC metrics) ---")
M4 <- run_all_events("M4_model_matrix_covariates", m4_fix, cov_df=m4_cov)

# M4x EXPLORATORY expression-PC sensitivity (topPC PC1-5; k pre-specified
# from documented variance labels; never preferred).
k_pc <- 5
tp_use <- as.data.frame(tp[, seq_len(k_pc), drop=FALSE])
names(tp_use) <- paste0("PC", seq_len(k_pc))
rownames(tp_use) <- rownames(am)
say("\n--- M4x EXPLORATORY (+ top ", k_pc, " expression PCs) ---")
M4x <- run_all_events("M4x_exploratory_expression_PC",
  paste("+ I(Age^2) + SeqBatch + Ancestry +",
        paste(paste0("PC", seq_len(k_pc)), collapse=" + ")),
  cov_df=cbind(data.frame(SeqBatch=mm$SeqBatch, Ancestry=mm$Ancestry,
                          row.names=rownames(am)), tp_use))

all_models <- M0
for (m in list(M1, M2, M3, M4, M4x))
  if (!is.null(m) && nrow(m) > 0) all_models <- rbind(all_models, m)
wtab(all_models, "TECHNICAL_MODEL_EVENT_RESULTS.tsv")

summ <- M0set
for (m in list(M1, M2, M3, M4, M4x)) {
  if (!is.null(m) && nrow(m) > 0)
    summ <- rbind(summ, set_summary(m, m$model[1]))
  else if (!is.null(m))
    say("WARNING: model object empty, skipped in set summary")
}
wtab(summ, "TECHNICAL_MODEL_SET_VALIDATION.tsv")

# ---- stability vs M0 -------------------------------------------------------
m0_beta <- setNames(M0$beta_ASD, M0$HsaEX_ID)
m0_dir  <- setNames(M0$direction, M0$HsaEX_ID)
m0_krf  <- setNames(M0$BH_FDR_KR, M0$HsaEX_ID)
tierA_ids <- primary19$HsaEX_ID[primary19$final_tier ==
                             "TIER_A_CROSS_COHORT_KR_FDR05"]
stab <- data.frame()
for (m in list(M1, M2, M3, M4, M4x)) {
  if (is.null(m) || nrow(m) == 0) next
  mn <- m$model[1]
  common <- intersect(M0$HsaEX_ID, m$HsaEX_ID)
  b1 <- setNames(m$beta_ASD, m$HsaEX_ID)[common]
  b0 <- m0_beta[common]
  d1 <- setNames(m$direction, m$HsaEX_ID)[common]
  d0 <- m0_dir[common]
  mset <- summ[summ$model == mn, ]
  don0 <- 80
  stab <- rbind(stab, data.frame(model=mn,
    n_events_compared=length(common),
    direction_concordance_vs_M0=paste0(sum(d0==d1), "/", length(common)),
    beta_pearson_r=round(cor(b0, b1), 4),
    beta_spearman_rho=round(cor(b0, b1, method="spearman"), 4),
    median_abs_beta_change=round(median(abs(b1-b0)), 5),
    max_abs_beta_change=round(max(abs(b1-b0)), 5),
    tierA_4of4_direction=paste0(sum(d0[tierA_ids]==d1[tierA_ids]), "/4"),
    tierA_KR_nominal_005=mset$tierA_KR_nominal_005,
    tierA_KR_FDR_005=mset$tierA_KR_FDR_005,
    tierA_FDR_stability_4of4=paste0(
      sum(setNames(m$BH_FDR_KR, m$HsaEX_ID)[tierA_ids] < 0.05,
          na.rm=TRUE), "/4"),
    set_concordance_vs_discovery=paste0(mset$n_concordant, "/",
                                        mset$n_events),
    set_binomial_P=mset$exact_binomial_P,
    oepg_concordance=paste0(mset$oepg_concordant, "/", mset$oepg_n),
    min_donors_per_event=min(m$n_donors),
    max_donors_lost=max(0, don0 - min(m$n_donors)),
    donor_loss_fraction=round((don0 - min(m$n_donors)) / don0, 4),
    high_missingness=ifelse((don0 - min(m$n_donors)) / don0 > 0.15,
      "HIGH_MISSINGNESS_SENSITIVITY_ONLY", "OK"),
    stringsAsFactors=FALSE))
}
wtab(stab, "TECHNICAL_MODEL_STABILITY_VS_M0.tsv")

# ---- pre-specified preferred-model rule (NOT outcome-driven) --------------
# Preferred = most complex documented model (M4 > M3 > M2 > M1 > M0) with:
#   (a) fitted for all 19 events, (b) no rank deficiency,
#   (c) donor loss <= 15%. M4x never eligible.
feasible <- function(m) {
  !is.null(m) && nrow(m) == 19 &&
    all(m$fixed_effect_rank == max(m$fixed_effect_rank)) &&
    (80 - min(m$n_donors)) / 80 <= 0.15
}
preferred <- if (feasible(M4)) "M4_model_matrix_covariates" else
             if (feasible(M3)) "M3_ancestry" else
             if (feasible(M2)) "M2_seq_batch" else
             if (feasible(M1)) "M1_nonlinear_age" else "M0_primary"
pref_res <- switch(preferred,
  M4_model_matrix_covariates=M4, M3_ancestry=M3, M2_seq_batch=M2,
  M1_nonlinear_age=M1, M0_primary=M0)
pref_fix <- switch(preferred,
  M4_model_matrix_covariates=paste("dx_binary + region + Sex + Age +",
    "I(Age^2) + RIN + SeqBatch + Ancestry + PMI_mm +",
    paste(paste0("QC", seq_along(qc_cols)), collapse=" + ")),
  M3_ancestry="dx_binary + region + Sex + Age + I(Age^2) + RIN + SeqBatch + Ancestry",
  M2_seq_batch="dx_binary + region + Sex + Age + I(Age^2) + RIN + SeqBatch",
  M1_nonlinear_age="dx_binary + region + Sex + Age + I(Age^2) + RIN",
  M0_primary="dx_binary + region + Sex + Age + RIN")
pref_cov <- switch(preferred,
  M4_model_matrix_covariates=m4_cov,
  M3_ancestry=data.frame(SeqBatch=mm$SeqBatch, Ancestry=mm$Ancestry,
                         row.names=rownames(am)),
  M2_seq_batch=data.frame(SeqBatch=mm$SeqBatch, row.names=rownames(am)),
  M1_nonlinear_age=NULL, M0_primary=NULL)
writeLines(c(paste0("PREFERRED_TECHNICAL_SENSITIVITY_MODEL=", preferred),
  "RULE=most_complex_documented_model_feasible (all 19 events fitted, no rank deficiency, donor loss <= 15%)",
  "ORDER_TESTED=M4 > M3 > M2 > M1 > M0",
  "M4x_exploratory_expression_PC=EXPLORATORY_NEVER_PREFERRED",
  paste0("M4_feasible=", feasible(M4)),
  paste0("M3_feasible=", feasible(M3)),
  paste0("M2_feasible=", feasible(M2)),
  paste0("M1_feasible=", feasible(M1)),
  paste0("PREFERRED_FIXED_EFFECTS=", pref_fix)),
  file.path(outS, "TECHNICAL_MODEL_PREFERRED_SELECTION.txt"))
say("\nPreferred technical sensitivity model: ", preferred)

# ---- stability report md ----------------------------------------------------
stab_md <- c(
"# Technical model stability report (module 1)",
"",
"## Models",
"",
"| model | fixed effects (beyond dx_binary) |",
"|---|---|",
"| M0_primary | region + Sex + Age + RIN |",
"| M1_nonlinear_age | M0 + I(Age^2) |",
"| M2_seq_batch | M1 + SeqBatch (== datMeta$batch) |",
"| M3_ancestry | M2 + Ancestry |",
"| M4_model_matrix_covariates | M3 + PMI + 8 documented Picard/STAR QC metrics |",
"| M4x_exploratory_expression_PC | M3 + topPC PC1-5 (EXPLORATORY, never preferred) |",
"",
"All models: random intercept (1|subject); logit usage with 1e-4 offset;",
"KR REML primary, ML-LRT sensitivity; BH within each model family.",
"",
"## Set-level validation per model",
"",
"| model | n events | concordance | binomial P | OEPG | KR FDR<.05 | KR FDR<.10 | LRT FDR<.05 | LRT FDR<.10 |",
"|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(summ))) {
  s <- summ[i, ]
  stab_md <- c(stab_md, sprintf(
    "| %s | %d | %d/%d | %s | %d/%d | %d | %d | %d | %d |",
    s$model, s$n_events, s$n_concordant, s$n_events,
    format(s$exact_binomial_P, digits=4), s$oepg_concordant, s$oepg_n,
    s$n_KR_FDR_005, s$n_KR_FDR_010, s$n_LRT_FDR_005, s$n_LRT_FDR_010))
}
stab_md <- c(stab_md, "",
"## Stability vs M0 (per-event betas/directions)",
"",
"| model | direction vs M0 | beta Pearson r | max |delta beta| | Tier A 4/4 direction | Tier A KR FDR<.05 | donor loss |",
"|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(stab))) {
  s <- stab[i, ]
  stab_md <- c(stab_md, sprintf(
    "| %s | %s | %.4f | %.5f | %s | %s | %s (%.0f%%) |",
    s$model, s$direction_concordance_vs_M0, s$beta_pearson_r,
    s$max_abs_beta_change, s$tierA_4of4_direction, s$tierA_KR_FDR_005,
    s$high_missingness, 100*s$donor_loss_fraction))
}
stab_md <- c(stab_md, "",
paste0("## Preferred model: ", preferred),
"",
"Pre-specified rule: most complex documented model feasible for all 19",
"events with no rank deficiency and <=15% donor loss. M4x is exploratory",
"and never preferred. Feasibility flags:",
paste0("M4=", feasible(M4), " M3=", feasible(M3), " M2=", feasible(M2),
       " M1=", feasible(M1)),
"",
"Interpretation rules: direction stability and 15/19 set concordance are",
"the primary robustness readouts; FDR movement under additional documented",
"covariates is reported, not used for model shopping. Any event whose",
"direction flips under a documented-covariate model is flagged",
"TECHNICAL_COVARIATE_SENSITIVE in module 8.",
"")
writeLines(stab_md, file.path(outS, "TECHNICAL_MODEL_STABILITY_REPORT.md"))
# NOTE for downstream: also expose the per-event stability table
wtab(stab, "TECHNICAL_MODEL_SET_SUMMARY.tsv")

# ============================================================================
# PART C — cell-composition feasibility + region sensitivity
# ============================================================================
say("\n=== PART C: cell-composition check + region sensitivity ===")

# C.1 composition field check (term scan of all metadata columns)
comp_terms <- c("NeuN","neuron","neuronal","excitatory","inhibitory",
  "astrocyte","oligodendrocyte","OPC","microglia","endothelial","cell_prop",
  "cell_fraction","methylation","deconvolution","composition","glia",
  "proportion","CIBERSORT","BRETIGEA","MuSiC","xCell","CIBERSORTx")
scan_cols <- c(colnames(datMeta), colnames(datMetaM), colnames(topPC),
               colnames(datSeqN), colnames(rsem_tx)[0])
hits <- grep(paste(comp_terms, collapse="|"), scan_cols,
             ignore.case=TRUE, value=TRUE)
comp_check <- data.frame(
  source=c("RData objects (datMeta, datMeta_model, topPC, datSeq_numeric): column-name scan",
           "rsem_tx / expression objects: no per-cell metadata (bulk sample x transcript matrix)"),
  method=c(paste0("grep of all colnames for: ",
                  paste(comp_terms, collapse="/")),
           "object structure inspection"),
  hits=c(length(hits), 0),
  hit_fields=c(if (length(hits)) paste(hits, collapse="; ") else "NONE",
               "NONE"),
  conclusion=c(
    if (length(hits)) "candidate fields require manual documentation check"
    else "NO cell-composition variable in any serialized object",
    "bulk-tissue matrices carry no cell-level information"),
  stringsAsFactors=FALSE)
wtab(comp_check, "CELL_COMPOSITION_FIELD_CHECK.tsv")
comp_class <- if (length(hits) == 0)
  "NOT_AVAILABLE_IN_EXISTING_PROCESSED_DATA" else "PARTIAL_PROXY_ONLY"
say("composition classification: ", comp_class)

# C.2 host-gene expression matrix (gene-level mean of datExpr log2 values)
e1 <- new.env()
load(file.path(data_dir, "02_01_B_AllProcessedData_wModelMatrix.RData"),
     envir=e1)
datExpr <- e1$datExpr
rm(e1); gc()
gt <- read.delim(file.path(gencode_dir, "gencode_transcripts.tsv"),
                 stringsAsFactors=FALSE)
base_expr <- sub("\\.[0-9]+_[0-9]+$", "", rownames(datExpr))
tx2gene <- setNames(gt$gene_name, gt$transcript_id)
gene_of_expr <- tx2gene[base_expr]
say("datExpr transcripts mapped to GENCODE v33 genes: ",
    sum(!is.na(gene_of_expr)), "/", length(base_expr))
host_genes <- unique(primary19$gene)
hge <- matrix(NA, nrow=nrow(am), ncol=length(host_genes),
              dimnames=list(rownames(am), host_genes))
expr_samp_idx <- match(rownames(am), colnames(datExpr))
for (g in host_genes) {
  i <- which(gene_of_expr == g)
  if (length(i) == 0) next
  hge[, g] <- colMeans(datExpr[i, expr_samp_idx, drop=FALSE], na.rm=TRUE)
}
rm(datExpr); gc()

fit_beta <- function(eid, mdf, fix) {
  f <- as.formula(paste("usage_logit ~", fix, "+ (1|subject)"))
  fit <- tryCatch(lmer(f, data=mdf, REML=TRUE), error=function(e) NULL)
  model_type <- "lmer_REML"
  if (is.null(fit)) {
    # strata with one sample per donor: random intercept not identifiable;
    # documented fixed-effects-only fallback
    f2 <- as.formula(paste("usage_logit ~", fix))
    fit <- tryCatch(lm(f2, data=mdf), error=function(e) NULL)
    model_type <- "lm_no_random_effect_single_sample_donors"
  }
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients
  if (!("dx_binary" %in% rownames(cf))) return(NULL)
  data.frame(HsaEX_ID=eid, beta_ASD=cf["dx_binary","Estimate"],
             SE=cf["dx_binary","Std. Error"],
             n_samples=nrow(mdf), n_donors=length(unique(mdf$subject)),
             singular=if (model_type == "lmer_REML") isSingular(fit) else NA,
             model_type=model_type, stringsAsFactors=FALSE)
}
prep <- function(sub_meta, eid, cov_df=NULL) {
  u <- usage_matrix[rownames(sub_meta), eid]
  ok <- !is.na(u)
  mdf <- sub_meta[ok, ]
  if (!is.null(cov_df)) {
    extra <- cov_df[match(rownames(mdf), rownames(cov_df)), , drop=FALSE]
    mdf <- cbind(mdf, extra)
  }
  mdf$usage <- u[ok]
  mdf$usage_logit <- log((mdf$usage + epsilon) / (1 - mdf$usage + epsilon))
  mdf
}

# C.3 Alternative A: region-stratified direction (>=5 ASD & >=5 CTL donors)
say("--- Alt A: region-stratified direction analysis ---")
don_ct <- do.call(rbind, lapply(split(am, am$region), function(d)
  data.frame(region=d$region[1],
             asd_donors=length(unique(d$subject[d$dx_binary==1])),
             ctl_donors=length(unique(d$subject[d$dx_binary==0])),
             n_samples=nrow(d), stringsAsFactors=FALSE)))
don_ct$eligible <- don_ct$asd_donors >= 5 & don_ct$ctl_donors >= 5
wtab(don_ct, "REGION_DONOR_ELIGIBILITY.tsv")
say("eligible regions: ", sum(don_ct$eligible), "/", nrow(don_ct))

reg_res <- data.frame()
for (rg in don_ct$region[don_ct$eligible]) {
  sub_meta <- am[am$region == rg, ]
  for (eid in events) {
    mdf <- prep(sub_meta, eid)
    r <- fit_beta(eid, mdf, "dx_binary + Sex + Age + RIN")
    if (!is.null(r)) { r$region <- rg; reg_res <- rbind(reg_res, r) }
  }
  say("  region ", rg, " done (", nrow(sub_meta), " samples)")
}
het <- do.call(rbind, lapply(events, function(eid) {
  d <- reg_res[reg_res$HsaEX_ID == eid, ]
  if (nrow(d) < 2) return(NULL)
  w <- 1 / d$SE^2
  mu <- sum(w * d$beta_ASD) / sum(w)
  Q <- sum(w * (d$beta_ASD - mu)^2)
  data.frame(HsaEX_ID=eid, n_regions=nrow(d), Cochran_Q=round(Q, 3),
             df=nrow(d)-1, P_het=round(pchisq(Q, nrow(d)-1,
             lower.tail=FALSE), 4),
             I2=round(max(0, (Q - (nrow(d)-1))) / max(Q, 1e-9) * 100, 1),
             stringsAsFactors=FALSE)
}))
reg_out <- merge(reg_res, primary19[, c("HsaEX_ID","gene","discovery_dir")],
                 by="HsaEX_ID")
reg_out$direction <- ifelse(reg_out$beta_ASD > 0, "UP_IN_ASD", "DOWN_IN_ASD")
reg_out$concordant <- reg_out$direction == reg_out$discovery_dir
wtab(reg_out, "REGION_SENSITIVITY_RESULTS.tsv")
wtab(het, "REGION_HETEROGENEITY.tsv")
say("Alt A done.")

# C.4 Alternative B: leave-one-region-out refits of the preferred model
say("--- Alt B: leave-one-region-out (preferred model: ", preferred, ") ---")
loro <- data.frame()
for (rg in unique(am$region)) {
  sub_meta <- am[am$region != rg, ]
  for (eid in events) {
    mdf <- prep(sub_meta, eid, cov_df=pref_cov)
    r <- fit_beta(eid, mdf, pref_fix)
    if (!is.null(r)) { r$region_excluded <- rg; loro <- rbind(loro, r) }
  }
  say("  LORO ", rg, " done")
}
loro <- merge(loro, pref_res[, c("HsaEX_ID","beta_ASD","direction")],
              by="HsaEX_ID", suffixes=c("_loro","_full"))
loro$beta_change <- loro$beta_ASD_loro - loro$beta_ASD_full
loro <- merge(loro, primary19[, c("HsaEX_ID","discovery_dir","final_tier")],
              by="HsaEX_ID")
loro$concordant_discovery <- ifelse(loro$beta_ASD_loro > 0, "UP_IN_ASD",
                                    "DOWN_IN_ASD") == loro$discovery_dir
loro$tierA <- loro$final_tier == "TIER_A_CROSS_COHORT_KR_FDR05"
loro_sum <- do.call(rbind, lapply(split(loro, loro$region_excluded),
  function(d) data.frame(region_excluded=d$region_excluded[1],
    n_events=nrow(d),
    discovery_concordance=paste0(sum(d$concordant_discovery), "/", nrow(d)),
    tierA_direction_stable=paste0(sum(sign(d$beta_ASD_loro[d$tierA]) ==
                                      sign(d$beta_ASD_full[d$tierA])), "/4"),
    max_abs_beta_change=round(max(abs(d$beta_change)), 4),
    min_donors=min(d$n_donors), stringsAsFactors=FALSE)))
wtab(loro, "LEAVE_ONE_REGION_OUT_RESULTS.tsv")
wtab(loro_sum, "LEAVE_ONE_REGION_OUT_SUMMARY.tsv")
say("Alt B done.")

# C.5 Alternative C: host-gene expression sensitivity (preferred covariate set)
say("--- Alt C: host-gene expression sensitivity ---")
altC <- data.frame()
for (eid in events) {
  g <- primary19$gene[primary19$HsaEX_ID == eid]
  mdf <- prep(am, eid, cov_df=pref_cov)
  mdf$host_expr <- hge[rownames(mdf), g]
  have <- !is.na(mdf$host_expr)
  if (sum(have) < nrow(mdf) * 0.9) {
    say("  ", eid, ": host expr missing; skip"); next
  }
  mdf <- mdf[have, ]
  r <- fit_beta(eid, mdf, paste("dx_binary + host_expr +",
    gsub("^dx_binary \\+ ", "", pref_fix)))
  if (!is.null(r)) { r$host_gene <- g
                     r$label <- "HOST_GENE_EXPRESSION_SENSITIVITY"
                     altC <- rbind(altC, r) }
  say("  Alt C ", eid, " done")
}
altC <- merge(altC, pref_res[, c("HsaEX_ID","beta_ASD","direction")],
              by="HsaEX_ID", suffixes=c("_C","_pref"))
altC <- merge(altC, primary19[, c("HsaEX_ID","discovery_dir")], by="HsaEX_ID")
altC$beta_change <- altC$beta_ASD_C - altC$beta_ASD_pref
altC$concordant_discovery <- ifelse(altC$beta_ASD_C > 0, "UP_IN_ASD",
                                    "DOWN_IN_ASD") == altC$discovery_dir
wtab(altC, "HOST_GENE_EXPRESSION_SENSITIVITY_RESULTS.tsv")
say("Alt C done: ", nrow(altC), " events fitted.")

# C.6 Region interaction sensitivity (dx_binary:region; ML-LRT vs preferred)
say("--- region interaction sensitivity ---")
inter <- data.frame()
for (eid in events) {
  mdf <- prep(am, eid, cov_df=pref_cov)
  f_full <- as.formula(paste("usage_logit ~", pref_fix,
                             "+ dx_binary:region + (1|subject)"))
  f_red  <- as.formula(paste("usage_logit ~", pref_fix, "+ (1|subject)"))
  fml <- tryCatch(lmer(f_full, data=mdf, REML=FALSE), error=function(e) NULL)
  rml <- tryCatch(lmer(f_red,  data=mdf, REML=FALSE), error=function(e) NULL)
  if (is.null(fml) || is.null(rml)) {
    inter <- rbind(inter, data.frame(HsaEX_ID=eid, LRT_chisq=NA, df=NA,
      P_interaction=NA, convergence="ERROR", stringsAsFactors=FALSE))
    next
  }
  df_dif <- length(fixef(fml)) - length(fixef(rml))
  lrt <- 2 * (as.numeric(logLik(fml)) - as.numeric(logLik(rml)))
  inter <- rbind(inter, data.frame(HsaEX_ID=eid, LRT_chisq=round(lrt, 4),
    df=df_dif, P_interaction=pchisq(lrt, df_dif, lower.tail=FALSE),
    convergence="OK", stringsAsFactors=FALSE))
  say("  interaction ", eid, " done")
}
inter$BH_FDR_interaction <- p.adjust(inter$P_interaction, method="BH")
inter <- merge(inter, primary19[, c("HsaEX_ID","gene")], by="HsaEX_ID")
wtab(inter, "REGION_INTERACTION_RESULTS.tsv")
say("region interaction done.")

# C.7 feasibility md + results table
cc_rows <- rbind(
  data.frame(analysis="CELL_COMPOSITION_CLASSIFICATION",
    result=comp_class, detail="no documented cell-fraction variable in any PsychENCODE processed object; see CELL_COMPOSITION_FIELD_CHECK.tsv",
    stringsAsFactors=FALSE),
  data.frame(analysis="CELL_COMPOSITION_ADJUSTMENT",
    result="NOT_FEASIBLE_FROM_AVAILABLE_PROCESSED_DATA",
    detail="no composition-adjusted model fitted or claimed",
    stringsAsFactors=FALSE),
  data.frame(analysis="REGION_STRATIFIED_DIRECTION",
    result=paste0(sum(don_ct$eligible), "/", nrow(don_ct), " regions eligible; see REGION_SENSITIVITY_RESULTS.tsv"),
    detail=paste0("eligible-region discovery concordance: ",
      paste(tapply(reg_out$concordant, reg_out$HsaEX_ID, sum), collapse=",")
      ), stringsAsFactors=FALSE),
  data.frame(analysis="LEAVE_ONE_REGION_OUT",
    result=paste0("preferred model ", preferred, "; see LEAVE_ONE_REGION_OUT_SUMMARY.tsv"),
    detail=paste0("Tier A direction stability across LORO folds: ",
      paste(loro_sum$tierA_direction_stable, collapse=", ")),
    stringsAsFactors=FALSE),
  data.frame(analysis="REGION_INTERACTION",
    result=paste0(sum(inter$P_interaction < 0.05, na.rm=TRUE),
                  " nominal dx:region interactions"),
    detail="ML-LRT of dx_binary:region added to preferred model; BH reported",
    stringsAsFactors=FALSE),
  data.frame(analysis="HOST_GENE_EXPRESSION_SENSITIVITY",
    result=paste0(nrow(altC), "/19 events fitted; see HOST_GENE_EXPRESSION_SENSITIVITY_RESULTS.tsv"),
    detail="gene-level mean log2 expression of host gene added to preferred covariate set; NOT a deconvolution substitute",
    stringsAsFactors=FALSE))
wtab(cc_rows, "CELL_COMPOSITION_RESULTS.tsv")

writeLines(c(
"# Cell-type composition feasibility (module 1)",
"",
"## Classification",
"",
"```text",
paste0("CLASSIFICATION=", comp_class),
"CELL_COMPOSITION_ADJUSTMENT=NOT_FEASIBLE_FROM_AVAILABLE_PROCESSED_DATA",
"```",
"",
"No sample-level or donor-level cell-type composition variable exists in",
"any PsychENCODE RData object (datMeta, datMeta_model, topPC,",
"datSeq_numeric; column-name scan for neuron/glial/astrocyte/OPC/microglia/",
"endothelial/fraction/deconvolution/methylation terms: 0 hits; see",
"CELL_COMPOSITION_FIELD_CHECK.tsv). No NeuN fraction, CIBERSORT/BRETIGEA/",
"MuSiC or methylation-deconvolution estimate is available locally, and",
"controlled-data re-derivation is out of scope. Therefore no",
"composition-adjusted model was fit and none is claimed.",
"",
"## Defensible alternatives performed",
"",
"- Region-stratified direction analysis (>=5 ASD and >=5 CTL donors per",
"  region): REGION_SENSITIVITY_RESULTS.tsv + REGION_DONOR_ELIGIBILITY.tsv +",
"  Cochran-Q/I2 heterogeneity in REGION_HETEROGENEITY.tsv. Low-powered",
"  strata are reported, not overinterpreted.",
"- Leave-one-region-out refits of the preferred technical model",
paste0("  (", preferred, "): LEAVE_ONE_REGION_OUT_RESULTS.tsv /"),
"  LEAVE_ONE_REGION_OUT_SUMMARY.tsv.",
"- Region interaction sensitivity: dx_binary:region ML-LRT against the",
"  preferred model (REGION_INTERACTION_RESULTS.tsv, BH within family).",
"- Host-gene expression sensitivity (HOST_GENE_EXPRESSION_SENSITIVITY):",
"  gene-level mean log2 expression of the host gene added to the preferred",
"  covariate set. This adjusts for overall host-gene expression level; it",
"  is NOT a cell-type deconvolution substitute and does not use the usage",
"  ratio or its denominator.",
"",
"Cells are not replicates: all models keep (1|subject) donor grouping.",
""), file.path(outS, "CELL_COMPOSITION_FEASIBILITY.md"))
say("PART C done.")

# ============================================================================
# PART D — transcript-set sensitivity D0-D4
# ============================================================================
say("\n=== PART D: transcript membership (GENCODE v33) ===")
# D1 strict-local membership rules (reference coordinates):
#   exon match = |start diff| <= 1 bp AND |end diff| <= 1 bp
#   inclusion clean-local = microexon present, immediate transcription-order
#     neighbours match reference up+down flanks, no other exon strictly
#     inside the inter-flank gap
#   exclusion clean-local = reference flanks adjacent (clean skip), no other
#     exon inside the gap
# D2 representative pair = highest median effective-length-normalized
#   expression over the 532 analysis samples; MANE Select / GENCODE basic /
#   APPRIS principal metadata NOT available locally (documented);
#   ties -> lexicographically smallest transcript ID.
m21 <- read.delim(file.path(dir21,
  "04_transcript_membership_repair",
  "00_event_transcript_membership_GRCh38.tsv"), stringsAsFactors=FALSE)
gex <- read.delim(file.path(gencode_dir, "gencode_exons.tsv"),
                  stringsAsFactors=FALSE)
gtr <- read.delim(file.path(gencode_dir, "gencode_transcripts.tsv"),
                  stringsAsFactors=FALSE)

parse_coord <- function(x) {
  p1 <- strsplit(x, ":")[[1]]; p2 <- strsplit(p1[2], "-")[[1]]
  list(s=as.numeric(p2[1]), e=as.numeric(p2[2]))
}
match_ex <- function(st, en, ref, tol=1) {
  abs(st - ref$s) <= tol && abs(en - ref$e) <= tol
}

rsem_member <- rsem_cache$rsem_member
member_base <- rsem_cache$member_base_ids
effLen_member <- rsem_cache$effLen_member
rates <- rsem_member[, samp_idx_raw] / effLen_member
med_rate <- apply(rates, 1, function(z) median(z, na.rm=TRUE))
names(med_rate) <- member_base

master <- data.frame()
set_def <- data.frame()
set_lists <- list()
for (eid in events) {
  gene <- primary19$gene[primary19$HsaEX_ID == eid]
  ref <- m21[m21$HsaEX_ID == eid, ][1, ]
  me <- parse_coord(ref$microexon_hg38)
  up <- parse_coord(ref$upstream_flanking_exon_hg38)
  dn <- parse_coord(ref$downstream_flanking_exon_hg38)
  gap_lo <- min(up$e, dn$e); gap_hi <- max(up$s, dn$s)
  d0_incl <- incl_tx$transcript_id[incl_tx$HsaEX_ID == eid]
  d0_excl <- excl_tx$transcript_id[excl_tx$HsaEX_ID == eid]

  tex <- gex[gex$gene_name == gene, ]
  tex <- tex[order(tex$transcript_id, tex$exon_number), ]
  tids <- unique(tex$transcript_id)
  ttype <- setNames(gtr$transcript_type[match(tids, gtr$transcript_id)], tids)

  d1_incl <- character(); d1_excl <- character()
  for (t in tids) {
    ex <- tex[tex$transcript_id == t, ]
    n <- nrow(ex)
    is_me <- vapply(seq_len(n), function(i)
      match_ex(ex$start[i], ex$end[i], me), TRUE)
    is_up <- vapply(seq_len(n), function(i)
      match_ex(ex$start[i], ex$end[i], up), TRUE)
    is_dn <- vapply(seq_len(n), function(i)
      match_ex(ex$start[i], ex$end[i], dn), TRUE)
    has_me <- any(is_me)
    in_gap <- ex$start > gap_lo & ex$end < gap_hi & !(is_me | is_up | is_dn)
    extra <- sum(in_gap)
    if (has_me) {
      i <- which(is_me)[1]
      a <- if (i > 1) ex[i-1, ] else ex[NA_integer_, ]
      b <- if (i < n) ex[i+1, ] else ex[NA_integer_, ]
      adj <- ((i-1) %in% which(is_up) && (i+1) %in% which(is_dn)) ||
             ((i-1) %in% which(is_dn) && (i+1) %in% which(is_up))
    } else {
      sp <- which((is_up & c(FALSE, is_dn[-n])) |
                  (is_dn & c(FALSE, is_up[-n])))
      if (length(sp) >= 1) { i0 <- sp[1]; a <- ex[i0, ]; b <- ex[i0+1, ]
                             adj <- TRUE }
      else { a <- ex[NA_integer_, ]; b <- ex[NA_integer_, ]; adj <- FALSE }
    }
    if (is.na(a$start[1])) { up_own <- a; dn_own <- b }
    else if (is.na(b$start[1])) { up_own <- a; dn_own <- b }
    else if (a$start <= b$start) { up_own <- a; dn_own <- b }
    else { up_own <- b; dn_own <- a }
    fmt <- function(e) if (is.na(e$start[1])) NA_character_
                       else paste0(e$chrom[1], ":", e$start[1], "-",
                                   e$end[1])
    same_flanks <- adj && extra == 0
    if (has_me && t %in% d0_incl && same_flanks) d1_incl <- c(d1_incl, t)
    if (!has_me && t %in% d0_excl && same_flanks) d1_excl <- c(d1_excl, t)

    pc <- ttype[t] == "protein_coding"
    det <- t %in% member_base && !is.na(med_rate[t]) && med_rate[t] > 0
    cls <- character()
    if (t %in% d0_incl) cls <- c(cls, "D0_INCLUSION")
    if (t %in% d0_excl) cls <- c(cls, "D0_EXCLUSION")
    if (t %in% d1_incl) cls <- c(cls, "D1_STRICT_LOCAL_INCLUSION")
    if (t %in% d1_excl) cls <- c(cls, "D1_STRICT_LOCAL_EXCLUSION")
    if (pc && t %in% d0_incl) cls <- c(cls, "D3_PROTEIN_CODING_INCLUSION")
    if (pc && t %in% d0_excl) cls <- c(cls, "D3_PROTEIN_CODING_EXCLUSION")
    if (length(cls) == 0) cls <- "NOT_IN_PRIMARY_SETS"
    reason <- if (has_me && same_flanks)
        "same immediate flanks; local-only microexon presence; no extra exon in window"
      else if (!has_me && same_flanks)
        "clean skip of microexon between reference flanks; no extra exon in window"
      else if (has_me && extra > 0)
        "additional alternative exon within local window"
      else if (has_me)
        "microexon present but immediate flanks differ from reference"
      else if (extra > 0)
        "lacks microexon; extra alternative exon within local window"
      else if (t %in% d0_excl)
        "lacks microexon; flank adjacency not clean under 1bp rule"
      else "not a member; structural context recorded for check"
    struct_diff <- if (!(t %in% c(d0_incl, d0_excl))) "YES_non_member"
                   else if (extra > 0 || !same_flanks) "LOCAL_DIFFERENCES"
                   else "NO"
    master <- rbind(master, data.frame(event_id=eid, gene=gene,
      transcript_id=t,
      transcript_version=sub("^ENST[0-9]+(\\.[0-9]+)$", "\\1", t),
      transcript_type=ttype[t],
      protein_coding_status=ifelse(pc, "protein_coding",
                                   paste0("non_coding:", ttype[t])),
      expression_detected=ifelse(det, "YES", "NO"),
      includes_microexon=ifelse(has_me, "YES", "NO"),
      excludes_microexon=ifelse(!has_me, "YES", "NO"),
      upstream_flanking_exon=fmt(up_own),
      downstream_flanking_exon=fmt(dn_own),
      same_local_flanks=ifelse(same_flanks, "YES", "NO"),
      other_alternative_exons_within_local_window=ifelse(extra > 0,
        "YES", "NO"),
      other_transcript_structure_differences=struct_diff,
      membership_class=paste(cls, collapse=";"),
      ambiguity_reason=reason, stringsAsFactors=FALSE))
  }
  pick <- function(cand) {
    if (length(cand) == 0) return(NA_character_)
    m <- med_rate[cand]; m[is.na(m)] <- -1
    cand[order(-m, cand)][1]
  }
  d1_ok <- length(d1_incl) >= 1 && length(d1_excl) >= 1
  set_def <- rbind(set_def, data.frame(event_id=eid, gene=gene,
    n_D0_incl=length(d0_incl), n_D0_excl=length(d0_excl),
    n_D1_incl=length(d1_incl), n_D1_excl=length(d1_excl),
    D1_analyzable=ifelse(d1_ok, "YES", "NOT_ANALYZABLE_STRICT_LOCAL"),
    D2_incl=pick(d0_incl), D2_excl=pick(d0_excl),
    D3_incl=pick(d0_incl[ttype[d0_incl] == "protein_coding"]),
    D3_excl=pick(d0_excl[ttype[d0_excl] == "protein_coding"]),
    n_gene_transcripts=length(tids),
    D4_eligible=ifelse(length(d0_incl) >= 2 || length(d0_excl) >= 2,
                       "YES", "NO"),
    stringsAsFactors=FALSE))
  set_lists[[eid]] <- list(d1_incl=d1_incl, d1_excl=d1_excl,
    d2_incl=pick(d0_incl), d2_excl=pick(d0_excl),
    d3_incl=d0_incl[ttype[d0_incl] == "protein_coding"],
    d3_excl=d0_excl[ttype[d0_excl] == "protein_coding"])
  say("  membership ", eid, " done (", length(tids), " transcripts)")
}
wtab(master, "TRANSCRIPT_MEMBERSHIP_MASTER.tsv")
wtab(set_def, "TRANSCRIPT_SET_DEFINITIONS.tsv", dir=outDsub)
for (eid in events)
  saveRDS(set_lists[[eid]], file.path(outDsub, paste0("sets_", eid, ".rds")))

say("\n=== PART D: transcript-set models (preferred model: ",
    preferred, ") ===")
idx_by_base <- setNames(seq_along(member_base), member_base)
usage_for <- function(eid, incl, excl) {
  ii <- idx_by_base[incl[incl %in% names(idx_by_base)]]
  ee <- idx_by_base[excl[excl %in% names(idx_by_base)]]
  if (length(ii) == 0 || length(ee) == 0) return(rep(NA, nrow(am)))
  ir <- colSums(rsem_member[ii, samp_idx_raw, drop=FALSE] /
                effLen_member[ii], na.rm=TRUE)
  er <- colSums(rsem_member[ee, samp_idx_raw, drop=FALSE] /
                effLen_member[ee], na.rm=TRUE)
  tot <- ir + er
  ifelse(tot > 0, ir / tot, NA)
}
fit_set_event <- function(u, cov_df=NULL) {
  ok <- !is.na(u)
  mdf <- am[ok, ]
  if (!is.null(cov_df)) {
    extra <- cov_df[match(rownames(mdf), rownames(cov_df)), , drop=FALSE]
    mdf <- cbind(mdf, extra)
  }
  if (nrow(mdf) < 20) return(NULL)
  mdf$usage_logit <- log((u[ok] + epsilon) / (1 - u[ok] + epsilon))
  f_full <- as.formula(paste("usage_logit ~", pref_fix, "+ (1|subject)"))
  f_red  <- as.formula(paste("usage_logit ~",
    sub("dx_binary \\+ ", "", pref_fix), "+ (1|subject)"))
  fr <- tryCatch(lmer(f_full, data=mdf, REML=TRUE), error=function(e) NULL)
  if (is.null(fr)) return(NULL)
  cf <- summary(fr)$coefficients
  if (!("dx_binary" %in% rownames(cf))) return(NULL)
  beta <- cf["dx_binary","Estimate"]; se <- cf["dx_binary","Std. Error"]
  kr <- tryCatch(KRmodcomp(fr, tryCatch(lmer(f_red, data=mdf, REML=TRUE),
             error=function(e) NULL)), error=function(e) NULL)
  pkr <- if (!is.null(kr)) kr$test$p.value[1] else NA
  fml <- tryCatch(lmer(f_full, data=mdf, REML=FALSE), error=function(e) NULL)
  rml <- tryCatch(lmer(f_red,  data=mdf, REML=FALSE), error=function(e) NULL)
  plrt <- NA
  if (!is.null(fml) && !is.null(rml))
    plrt <- pchisq(2*(as.numeric(logLik(fml)) - as.numeric(logLik(rml))),
                   1, lower.tail=FALSE)
  data.frame(beta_ASD=beta, SE=se, CI_lo=beta-1.96*se, CI_hi=beta+1.96*se,
             P_KR=pkr, P_LRT=plrt,
             direction=ifelse(beta>0,"UP_IN_ASD","DOWN_IN_ASD"),
             n_samples=nrow(mdf), n_donors=length(unique(mdf$subject)),
             singular=isSingular(fr), stringsAsFactors=FALSE)
}

d0_incl <- setNames(lapply(events, function(e)
  incl_tx$transcript_id[incl_tx$HsaEX_ID == e]), events)
d0_excl <- setNames(lapply(events, function(e)
  excl_tx$transcript_id[excl_tx$HsaEX_ID == e]), events)

run_def <- function(def_name, get_lists) {
  res <- data.frame()
  for (eid in events) {
    L <- get_lists(eid)
    if (is.null(L)) next
    u <- usage_for(eid, L$incl, L$excl)
    r <- fit_set_event(u, cov_df=pref_cov)
    if (!is.null(r)) { r$HsaEX_ID <- eid; r$definition <- def_name
                       res <- rbind(res, r) }
    say("  ", def_name, " ", eid, " done")
  }
  if (nrow(res) > 1) {
    res$BH_FDR_KR  <- p.adjust(res$P_KR,  method="BH")
    res$BH_FDR_LRT <- p.adjust(res$P_LRT, method="BH")
  }
  say("  ", def_name, ": ", nrow(res), " events fitted")
  res
}

say("--- D0 final sets ---")
D0s <- run_def("D0", function(e) list(incl=d0_incl[[e]], excl=d0_excl[[e]]))
say("--- D1 strict local ---")
D1s <- run_def("D1", function(e) {
  s <- set_lists[[e]]
  if (length(s$d1_incl) == 0 || length(s$d1_excl) == 0) return(NULL)
  list(incl=s$d1_incl, excl=s$d1_excl)})
say("--- D2 representative pair ---")
D2s <- run_def("D2", function(e) {
  s <- set_lists[[e]]
  if (is.na(s$d2_incl) || is.na(s$d2_excl)) return(NULL)
  list(incl=s$d2_incl, excl=s$d2_excl)})
say("--- D3 protein-coding only ---")
D3s <- run_def("D3", function(e) {
  s <- set_lists[[e]]
  if (length(s$d3_incl) == 0 || length(s$d3_excl) == 0) return(NULL)
  list(incl=s$d3_incl, excl=s$d3_excl)})

allres <- do.call(rbind, list(D0s, D1s, D2s, D3s))
allres <- merge(allres, primary19[, c("HsaEX_ID","gene","discovery_dir",
                                   "final_tier")], by="HsaEX_ID")
wtab(allres, "TRANSCRIPT_SET_EVENT_RESULTS.tsv")

# ---- D4 leave-one-transcript-out (REML beta/SE; no KR for runtime) -------
say("--- D4 leave-one-transcript-out ---")
d4 <- data.frame()
for (eid in events) {
  ci <- d0_incl[[eid]]; ce <- d0_excl[[eid]]
  if (length(ci) < 2 && length(ce) < 2) next
  for (drop_tx in c(ci, ce)) {
    u <- usage_for(eid, setdiff(ci, drop_tx), setdiff(ce, drop_tx))
    ok <- !is.na(u)
    mdf <- am[ok, ]
    if (!is.null(pref_cov)) {
      extra <- pref_cov[match(rownames(mdf), rownames(pref_cov)), ,
                        drop=FALSE]
      mdf <- cbind(mdf, extra)
    }
    mdf$usage_logit <- log((u[ok]+epsilon)/(1-u[ok]+epsilon))
    fr <- tryCatch(lmer(as.formula(paste("usage_logit ~", pref_fix,
      "+ (1|subject)")), data=mdf, REML=TRUE), error=function(e) NULL)
    if (is.null(fr)) next
    cf <- summary(fr)$coefficients
    if (!("dx_binary" %in% rownames(cf))) next
    d4 <- rbind(d4, data.frame(HsaEX_ID=eid, dropped=drop_tx,
      side=ifelse(drop_tx %in% ci, "inclusion", "exclusion"),
      beta_ASD=cf["dx_binary","Estimate"], SE=cf["dx_binary","Std. Error"],
      stringsAsFactors=FALSE))
  }
  say("  D4 ", eid, " done")
}
d4 <- merge(d4, D0s[, c("HsaEX_ID","beta_ASD","direction")], by="HsaEX_ID",
            suffixes=c("_loo","_D0"))
d4$dir_flip <- ifelse(d4$beta_ASD_loo > 0, "UP_IN_ASD", "DOWN_IN_ASD") !=
                     d4$direction
wtab(d4, "LEAVE_ONE_TRANSCRIPT_OUT_RESULTS.tsv")

# ---- summaries ------------------------------------------------------------
sum_def <- function(res, name) {
  if (is.null(res) || nrow(res) == 0)
    return(data.frame(definition=name, n_analyzable=0,
      n_analyzable_Tier_A=0, note="none", stringsAsFactors=FALSE))
  m <- merge(res, primary19[, c("HsaEX_ID","discovery_dir","final_tier",
                             "abs_delta_psi","gene")], by="HsaEX_ID")
  m$conc <- m$direction == m$discovery_dir
  bp <- pbinom(sum(m$conc)-1, nrow(m), 0.5, lower.tail=FALSE)
  ms <- m[order(-m$abs_delta_psi, m$HsaEX_ID), ]
  oepg <- do.call(rbind, lapply(split(ms, ms$gene), function(d) d[1,]))
  tA <- m[m$final_tier == "TIER_A_CROSS_COHORT_KR_FDR05", ]
  data.frame(definition=name, n_analyzable=nrow(m),
    n_analyzable_Tier_A=nrow(tA),
    direction_concordance=paste0(sum(m$conc), "/", nrow(m)),
    exact_binomial_P=signif(bp, 6),
    oepg_concordance=paste0(sum(oepg$conc), "/", nrow(oepg)),
    tierA_direction_flags=paste(tA$conc, collapse=""),
    tierA_KR_FDR05=sum(tA$BH_FDR_KR < 0.05, na.rm=TRUE),
    n_KR_FDR05=sum(m$BH_FDR_KR < 0.05, na.rm=TRUE),
    n_LRT_FDR05=sum(m$BH_FDR_LRT < 0.05, na.rm=TRUE),
    stringsAsFactors=FALSE)
}
summD <- do.call(rbind, list(sum_def(D0s,"D0"), sum_def(D1s,"D1"),
                             sum_def(D2s,"D2"), sum_def(D3s,"D3")))
bc <- function(a, b) {
  common <- intersect(a$HsaEX_ID, b$HsaEX_ID)
  if (length(common) < 2) return(NA)
  round(cor(a$beta_ASD[match(common, a$HsaEX_ID)],
            b$beta_ASD[match(common, b$HsaEX_ID)]), 4)
}
flips <- function(a, b) {
  common <- intersect(a$HsaEX_ID, b$HsaEX_ID)
  sum(a$direction[match(common, a$HsaEX_ID)] !=
      b$direction[match(common, b$HsaEX_ID)])
}
sens <- data.frame(metric=c("D0_vs_D1_beta_correlation",
  "D0_vs_D2_beta_correlation", "D0_vs_D3_beta_correlation",
  "D0_vs_D1_direction_flips", "D0_vs_D2_direction_flips",
  "D0_vs_D3_direction_flips", "D1_non_analyzable_events",
  "D3_non_analyzable_events", "D4_eligible_events",
  "D4_max_abs_beta_change_vs_D0", "D4_direction_flips_total"),
  value=c(bc(D0s,D1s), bc(D0s,D2s), bc(D0s,D3s),
          flips(D0s,D1s), flips(D0s,D2s), flips(D0s,D3s),
          sum(!events %in% D1s$HsaEX_ID), sum(!events %in% D3s$HsaEX_ID),
          length(unique(d4$HsaEX_ID)),
          if (nrow(d4)) round(max(abs(d4$beta_ASD_loo - d4$beta_ASD_D0)), 4)
          else NA,
          sum(d4$dir_flip)), stringsAsFactors=FALSE)
wtab(summD, "TRANSCRIPT_SET_SUMMARY.tsv")
wtab(sens, "TRANSCRIPT_SET_SENSITIVITY_METRICS.tsv")

# ---- sensitivity report md -------------------------------------------------
d1_error <- set_def$event_id[set_def$D1_analyzable ==
                            "NOT_ANALYZABLE_STRICT_LOCAL"]
rep_md <- c(
"# Transcript-set sensitivity report (module 1)",
"",
"## Definitions tested",
"",
"- D0: final all-transcript inclusion/exclusion sets (locked).",
"- D1: strict local structure - same immediate upstream and downstream",
"  exons (1 bp tolerance), local difference only the target microexon,",
"  no additional local alternative exon within the inter-flank window;",
"  membership restricted to final D0 members satisfying the local rule.",
"- D2: deterministic representative pair per side: highest median",
"  effective-length-normalized expression over the 532 analysis samples.",
"  Tie-break chain per spec: expression > MANE Select > GENCODE basic >",
"  APPRIS principal > lexical transcript ID. MANE/basic/APPRIS metadata",
"  are NOT available locally (no MANE/APPRIS/basic annotation files in the",
"  project); ties after expression therefore resolve lexically.",
"  DOCUMENTED_LIMITATION.",
"- D3: protein-coding-only members (GENCODE v33 transcript_type).",
"- D4: leave-one-transcript-out over D0 sets (REML betas; KR omitted for",
"  runtime).",
"",
paste0("All definition models use the preferred technical model ("),
paste0(preferred, "); no selection by ASD effect at any step."),
"",
"## Analyzability",
"",
paste0("- D1 strict-local analyzable: ", sum(set_def$D1_analyzable=="YES"),
       "/19; errors: ",
       if (length(d1_error)) paste(d1_error, collapse=", ") else "none"),
paste0("- D3 protein-coding analyzable: ", nrow(D3s), "/19"),
paste0("- D4 eligible events: ", length(unique(d4$HsaEX_ID))),
"",
"## Set summaries",
"",
"| definition | n | concordance | binomial P | OEPG | Tier A KR FDR<.05 |",
"|---|---|---|---|---|---|")
for (i in seq_len(nrow(summD))) {
  s <- summD[i, ]
  rep_md <- c(rep_md, sprintf("| %s | %d | %s | %s | %s | %s |",
    s$definition, s$n_analyzable, s$direction_concordance,
    format(s$exact_binomial_P, digits=4), s$oepg_concordance,
    s$tierA_KR_FDR05))
}
rep_md <- c(rep_md, "",
"## Sensitivity metrics",
"",
paste0("- beta correlations vs D0: D1=", sens$value[1],
       " D2=", sens$value[2], " D3=", sens$value[3]),
paste0("- direction flips vs D0: D1=", sens$value[4],
       " D2=", sens$value[5], " D3=", sens$value[6]),
paste0("- D4 max |beta change| vs D0: ", sens$value[10],
       "; total direction flips: ", sens$value[11]),
"",
"## Interpretation rules",
"",
"- Events whose direction flips under D1/D2/D3 are flagged",
"  ANNOTATION_SENSITIVE in module 8 (no tier changes).",
"- Strict-definition errors (D1 not analyzable) are reported as",
"  ANNOTATION_SENSITIVE with NA results, never as discordant.",
"- MANE/basic/APPRIS absence is a documented limitation of the",
"  representative-pair tie-break, not a result.",
"")
writeLines(rep_md, file.path(outS, "TRANSCRIPT_SET_SENSITIVITY_REPORT.md"))

# ---- input/output manifest + QC check for the first analysis phase --------------------------
in_manifest <- data.frame(
  path=c(file.path(data_dir, "01_02_B_01_RawData.RData"),
         file.path(data_dir, "02_01_B_AllProcessedData_wModelMatrix.RData"),
         file.path(project, "16_gse30573/02_input_lock/02_primary19.tsv"),
         file.path(dir18, "07_annotation_and_transcript_mapping/03_inclusion_transcripts.tsv"),
         file.path(dir18, "07_annotation_and_transcript_mapping/04_exclusion_transcripts.tsv"),
         file.path(dir21, "04_transcript_membership_repair/00_event_transcript_membership_GRCh38.tsv"),
         file.path(dir21, "05_mixed_model_inference/01_Satterthwaite_models.tsv"),
         file.path(dir21, "05_mixed_model_inference/06_set_validation_recomputed.tsv"),
         file.path(dir25, "06_master_event_table/MASTER_19_EVENT_EVIDENCE_TABLE.tsv"),
         file.path(gencode_dir, "gencode_exons.tsv"),
         file.path(gencode_dir, "gencode_transcripts.tsv")),
  role=c("psychencode_raw_rsem_counts+effLen",
         "psychencode_processed_metadata+model_matrix+topPC",
         "primary19_events", "inclusion_transcript_sets",
         "exclusion_transcript_sets",
         "dir21_GRCh38_reference_local_coordinates",
         "M0_reference_comparison", "set_validation_keys",
         "evidence_tiers", "GENCODE_v33_exons",
         "GENCODE_v33_transcripts"),
  exists=NA, stringsAsFactors=FALSE)
in_manifest$exists <- file.exists(in_manifest$path)
wtab(in_manifest, "INPUT_MANIFEST.tsv")
out_files <- sort(list.files(outS, full.names=FALSE, recursive=TRUE))
out_manifest <- data.frame(
  path=file.path("02_psychencode_sensitivity", out_files),
  size_bytes=file.info(file.path(outS, out_files))$size,
  stringsAsFactors=FALSE)
wtab(out_manifest, "OUTPUT_MANIFEST.tsv")

qc1 <- data.frame(phase=c(
  "PRIMARY_MODEL_REPRODUCED",
  "PSYCHENCODE_COVARIATE_CHECK_COMPLETE",
  "PSYCHENCODE_TECHNICAL_HIERARCHY_COMPLETE",
  "CELL_COMPOSITION_CLASSIFIED",
  "REGION_SENSITIVITY_COMPLETE",
  "TRANSCRIPT_MEMBERSHIP_REBUILT",
  "TRANSCRIPT_SET_MODELS_COMPLETE",
  "NO_OUTCOME_DRIVEN_SELECTION",
  "SEED"),
  status=c(check_overall,
    ifelse(nrow(dict) > 0, "OK", "ERROR"),
    ifelse(all(c(!is.null(M1), !is.null(M2), !is.null(M4))), "OK", "ERROR"),
    "OK",
    ifelse(nrow(reg_out) > 0 && nrow(loro) > 0, "OK", "ERROR"),
    ifelse(nrow(master) > 0, "OK", "ERROR"),
    ifelse(nrow(allres) > 0, "OK", "ERROR"),
    "OK_BY_CONSTRUCTION (pre-specified hierarchy, rule-based preferred model, no post-hoc covariate/transcript choice)",
    "42"),
  stringsAsFactors=FALSE)
wtab(qc1, "QC_CHECKS.txt")
say("\nStep 1 checks written.")
say("ANALYSIS DONE in ", round(as.numeric(difftime(Sys.time(), t0_all,
    units="mins")), 1), " min.")
