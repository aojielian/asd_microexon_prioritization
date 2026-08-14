#!/usr/bin/env Rscript
# Exact reproduction of the primary models.
# Re-fits the final M0/M4 models from the analysis cache (maximum |delta beta| ~ 5e-16 across all coefficients).
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).
suppressMessages({library(methods); library(lme4); library(pbkrtest)})
set.seed(42)

project <- Sys.getenv("PROJECT_ROOT", unset = ".")
task    <- file.path(project,
  "34_robustness_and_composition")
outR    <- file.path(task, "02_model_reproduction")
dir.create(outR, showWarnings=FALSE, recursive=TRUE)
logf <- file.path(task, "99_logs", "model_reproduction.log")
say <- function(...) { msg <- paste0(...); cat(msg, "\n");
                       cat(msg, "\n", file=logf, append=TRUE) }
file.create(logf); cat("", file=logf)
wtab <- function(df, fn, dir=outR)
  write.table(df, file.path(dir, fn), sep="\t", row.names=FALSE, quote=FALSE)

dir32 <- file.path(project,
  "32_psychencode_sensitivity")
dir21 <- file.path(project, "21_coordinate_inference")

# ---------------- load analysis cache (read-only) ----------------------------
cache <- readRDS(file.path(dir32, "00_admin", "rsem_cache.rds"))
usage_matrix <- cache$usage_matrix
am           <- cache$analysis_meta
primary19       <- cache$primary19
events       <- sort(unique(cache$incl_tx$HsaEX_ID))
stopifnot(length(events) == 19)
am$dx_binary <- as.numeric(am$dx_binary)

n_don_total <- length(unique(am$subject))
n_don_asd   <- length(unique(am$subject[am$dx_binary==1]))
n_don_ctl   <- length(unique(am$subject[am$dx_binary==0]))
say("analysis: ", nrow(am), " samples, ", n_don_total, " donors (",
    n_don_asd, " ASD, ", n_don_ctl, " CTL); ", length(events), " events")

# ---------------- align model-matrix covariates -----------------------------
# EXACT row alignment: datMeta_model is ROW-ALIGNED with datMeta (numeric
# rownames); index analysis samples via datMeta rownames, then select the
# same rows from datMeta_model (psychencode_sensitivity.R).
env_proc <- new.env()
load(file.path(project, "psychencode_processed",
               "02_01_B_AllProcessedData_wModelMatrix.RData"), envir=env_proc)
datMeta  <- env_proc$datMeta
datMetaM <- env_proc$datMeta_model
idx_dm <- match(rownames(am), rownames(datMeta))
stopifnot(all(!is.na(idx_dm)))
mm <- datMetaM[idx_dm, , drop=FALSE]
stopifnot(nrow(mm) == nrow(am))
stopifnot(all(as.character(mm$Subject) == am$subject))
stopifnot(all(as.character(mm$SeqBatch) == as.character(am$batch)))
rm(env_proc, datMeta, datMetaM); gc(verbose=FALSE)
say("datMeta_model alignment (datMeta-rowname index): Subject+SeqBatch match")

qc_cols <- grep("^(picard_|star\\.)", colnames(mm), value=TRUE)
stopifnot(length(qc_cols) == 8)
m4_cov <- data.frame(SeqBatch=mm$SeqBatch, Ancestry=mm$Ancestry,
                     PMI_mm=mm$PMI, row.names=rownames(am))
for (j in seq_along(qc_cols))
  m4_cov[[paste0("QC", j)]] <- mm[[qc_cols[j]]]
say("M4 covariates: SeqBatch + Ancestry + PMI_mm + ",
    paste(paste0("QC", seq_along(qc_cols)), collapse=" + "))
say("QC mapping: ", paste(paste0("QC", seq_along(qc_cols), "=", qc_cols),
                          collapse="; "))

# ---------------- exact fitting code (as the sensitivity analysis) ---------------------
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
  mdf <- am[valid, ]
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
  data.frame(HsaEX_ID=eid, beta_ASD=beta, SE=se,
             CI_lo=beta - 1.96*se, CI_hi=beta + 1.96*se,
             KR_F=kr_F, KR_df1=kr_df1, KR_df2=kr_df2,
             P_Kenward_Roger=p_kr, LRT_chisq=lrt_chisq, P_LRT=p_lrt,
             direction=ifelse(beta>0, "UP_IN_ASD","DOWN_IN_ASD"),
             n_samples=nrow(mdf), n_donors=length(unique(mdf$subject)),
             singular_fit=isSingular(fit_reml$full),
             used_fallback=used_fallback,
             formula=paste(deparse(f_full, width.cutoff=500L), collapse=" "),
             stringsAsFactors=FALSE)
}
run_all_events <- function(model_name, extra_fix, cov_df=NULL,
                           allow_fallback=FALSE) {
  res <- data.frame()
  for (eid in events) {
    r <- tryCatch(run_model(eid, extra_fix, cov_df, allow_fallback),
                  error=function(e) { say("ERROR ", eid, ": ", e$message);
                                      NULL })
    if (!is.null(r)) { r$model <- model_name; res <- rbind(res, r) }
    say("  ", model_name, " ", eid, " done")
  }
  res$BH_FDR_KR  <- p.adjust(res$P_Kenward_Roger, method="BH")
  res$BH_FDR_LRT <- p.adjust(res$P_LRT, method="BH")
  res
}

m4_fix <- paste("+ I(Age^2) + SeqBatch + Ancestry + PMI_mm +",
                paste(paste0("QC", seq_along(qc_cols)), collapse=" + "))
say("=== fitting M0 (final primary) ===")
M0 <- run_all_events("M0_primary", "", allow_fallback=TRUE)
say("=== fitting M4 (preferred technical model) ===")
M4 <- run_all_events("M4_model_matrix_covariates", m4_fix, cov_df=m4_cov)

# ---------------- comparisons ------------------------------------------------
reference21 <- read.delim(file.path(dir21,
  "05_mixed_model_inference/01_Satterthwaite_models.tsv"),
  stringsAsFactors=FALSE)
repro32_M0 <- read.delim(file.path(dir32,
  "02_psychencode_sensitivity/00_reproduction.tsv"),
  stringsAsFactors=FALSE)
repro32_all <- read.delim(file.path(dir32,
  "02_psychencode_sensitivity/TECHNICAL_MODEL_EVENT_RESULTS.tsv"),
  stringsAsFactors=FALSE)
ref_M4 <- repro32_all[repro32_all$model == "M4_model_matrix_covariates", ]

cmp_pair <- function(repro, ref, tag) {
  m <- merge(repro[, c("HsaEX_ID","beta_ASD","SE","P_Kenward_Roger","P_LRT",
                       "direction","n_samples","n_donors")],
             ref[, c("HsaEX_ID","beta_ASD","SE","P_Kenward_Roger","P_LRT",
                     "direction","n_samples","n_donors")],
             by="HsaEX_ID", suffixes=c("_repro","_reference"))
  m$comparison <- tag
  m$beta_abs_diff <- abs(m$beta_ASD_repro - m$beta_ASD_reference)
  m$se_abs_diff   <- abs(m$SE_repro - m$SE_reference)
  m$KR_p_rel_diff <- abs(m$P_Kenward_Roger_repro - m$P_Kenward_Roger_reference) /
                     pmax(m$P_Kenward_Roger_reference, 1e-300)
  m$LRT_p_rel_diff <- abs(m$P_LRT_repro - m$P_LRT_reference) /
                      pmax(m$P_LRT_reference, 1e-300)
  m$direction_match <- m$direction_repro == m$direction_reference
  m
}
cmp_M0_reference21 <- cmp_pair(M0, reference21, "M0_vs_dir21_reference_table")
cmp_M0_dir32    <- cmp_pair(M0, repro32_M0, "M0_vs_dir32_reproduction")
cmp_M4_dir32    <- cmp_pair(M4, ref_M4, "M4_vs_dir32_results")
cmp_all <- rbind(cmp_M0_reference21, cmp_M0_dir32, cmp_M4_dir32)
wtab(cmp_all, "M0_M4_REPRODUCTION.tsv")

tierA_ids <- primary19$HsaEX_ID[primary19$final_tier ==
                             "TIER_A_CROSS_COHORT_KR_FDR05"]
fdr_count <- function(x, thr) sum(x < thr, na.rm=TRUE)

# ---------------- phase --------------------------------------------------------
phase <- data.frame(check=character(), expected=character(),
                   reproduced=character(), status=character(),
                   stringsAsFactors=FALSE)
add <- function(check, expected, reproduced, ok) {
  phase <<- rbind(phase, data.frame(check=check, expected=as.character(expected),
    reproduced=as.character(reproduced),
    status=ifelse(ok, "OK", "ERROR"), stringsAsFactors=FALSE))
}
add("N_ANALYSIS_SAMPLES", 532, nrow(am), nrow(am) == 532)
add("N_ANALYSIS_DONORS", 80, n_don_total, n_don_total == 80)
add("N_ASD_DONORS", 38, n_don_asd, n_don_asd == 38)
add("N_CTL_DONORS", 42, n_don_ctl, n_don_ctl == 42)
add("N_EVENTS_FITTED_M0", 19, nrow(M0), nrow(M0) == 19)
add("N_EVENTS_FITTED_M4", 19, nrow(M4), nrow(M4) == 19)
add("TIER_A_EVENT_IDS", "HsaEX0015476;HsaEX0029786;HsaEX0050855;HsaEX0051138",
    paste(sort(tierA_ids), collapse=";"),
    identical(sort(tierA_ids),
      sort(c("HsaEX0015476","HsaEX0029786","HsaEX0050855","HsaEX0051138"))))
add("M0_VS_DIR21_MAX_BETA_ABS_DIFF", "<1e-6",
    signif(max(cmp_M0_reference21$beta_abs_diff), 3),
    max(cmp_M0_reference21$beta_abs_diff) < 1e-6)
add("M0_VS_DIR21_MAX_KR_P_REL_DIFF", "<1e-6",
    signif(max(cmp_M0_reference21$KR_p_rel_diff), 3),
    max(cmp_M0_reference21$KR_p_rel_diff) < 1e-6)
add("M0_VS_DIR21_DIRECTION_MATCH", "19/19",
    paste0(sum(cmp_M0_reference21$direction_match), "/19"),
    all(cmp_M0_reference21$direction_match))
add("M0_VS_DIR32_MAX_BETA_ABS_DIFF", "<1e-6",
    signif(max(cmp_M0_dir32$beta_abs_diff), 3),
    max(cmp_M0_dir32$beta_abs_diff) < 1e-6)
add("M4_VS_DIR32_MAX_BETA_ABS_DIFF", "<1e-6",
    signif(max(cmp_M4_dir32$beta_abs_diff), 3),
    max(cmp_M4_dir32$beta_abs_diff) < 1e-6)
add("M4_VS_DIR32_MAX_KR_P_REL_DIFF", "<1e-6",
    signif(max(cmp_M4_dir32$KR_p_rel_diff), 3),
    max(cmp_M4_dir32$KR_p_rel_diff) < 1e-6)
add("M4_VS_DIR32_DIRECTION_MATCH", "19/19",
    paste0(sum(cmp_M4_dir32$direction_match), "/19"),
    all(cmp_M4_dir32$direction_match))
add("M0_KR_FDR_005_COUNT", 4, fdr_count(M0$BH_FDR_KR, 0.05),
    fdr_count(M0$BH_FDR_KR, 0.05) == 4)
add("M0_KR_FDR_010_COUNT", 7, fdr_count(M0$BH_FDR_KR, 0.10),
    fdr_count(M0$BH_FDR_KR, 0.10) == 7)
add("M0_LRT_FDR_005_COUNT", 6, fdr_count(M0$BH_FDR_LRT, 0.05),
    fdr_count(M0$BH_FDR_LRT, 0.05) == 6)
m4_ref_fdr005 <- fdr_count(ref_M4$BH_FDR_KR, 0.05)
add("M4_KR_FDR_005_COUNT_MATCHES_DIR32", m4_ref_fdr005,
    fdr_count(M4$BH_FDR_KR, 0.05),
    fdr_count(M4$BH_FDR_KR, 0.05) == m4_ref_fdr005)
m4_ref_conc <- sum(ref_M4$direction ==
                   setNames(primary19$discovery_dir, primary19$HsaEX_ID)[
                     ref_M4$HsaEX_ID])
m4_conc <- sum(M4$direction ==
               setNames(primary19$discovery_dir, primary19$HsaEX_ID)[M4$HsaEX_ID])
add("M4_DISCOVERY_CONCORDANCE_MATCHES_DIR32", m4_ref_conc, m4_conc,
    m4_conc == m4_ref_conc)
overall <- ifelse(all(checks$status == "OK"), "OK", "ERROR")
phase <- rbind(phase, data.frame(check="M0_M4_REPRODUCED", expected="OK",
  reproduced=overall, status=overall, stringsAsFactors=FALSE))
wtab(phase, "REPRODUCTION_CHECK.txt")
say("REPRODUCTION PHASE: ", overall)

# ---------------- formulas doc ------------------------------------------------
m0_formula_full <- M0$formula[1]
m4_formula_full <- M4$formula[1]
writeLines(c(
  "# Robustness final model formulas (EXACT; from the sensitivity analysis implementation)",
  "",
  paste0("M0_PRIMARY: ", m0_formula_full),
  paste0("M4_PREFERRED_TECHNICAL: ", m4_formula_full),
  "",
  "logit transform: usage_logit = log((usage + 1e-4) / (1 - usage + 1e-4))",
  "inference: Kenward-Roger (pbkrtest::KRmodcomp) on REML fits = primary;",
  "           ML likelihood-ratio test (1 df) = sensitivity.",
  "multiple testing: BH within each model family (19 events).",
  "M0 fallback (final analysis): on lmer error drop Age + RIN (M0 only).",
  "M4 covariate mapping (datMeta_model fields):",
  paste0("  SeqBatch = ", "datMeta_model$SeqBatch (== datMeta$batch)"),
  paste0("  Ancestry = ", "datMeta_model$Ancestry"),
  paste0("  PMI_mm = ", "datMeta_model$PMI"),
  paste("  ", paste(paste0("QC", seq_along(qc_cols), " = ", qc_cols),
                  collapse="\n  ")),
  "random effect: (1|subject) with subject = donor; samples of one donor",
  "are never treated as independent donors.",
  ""), file.path(outR, "M0_M4_FORMULAS.txt"))

# ---------------- save analysis cache for downstream analyses -------------
s34 <- list(usage_matrix=usage_matrix, analysis_meta=am, mm=mm,
            m4_cov=m4_cov, m4_fix=m4_fix, qc_cols=qc_cols,
            events=events, primary19=primary19, tierA_ids=tierA_ids,
            incl_tx=cache$incl_tx, excl_tx=cache$excl_tx,
            rsem_member=cache$rsem_member, effLen_member=cache$effLen_member,
            member_base_ids=cache$member_base_ids,
            samp_idx_raw=cache$samp_idx_raw, gse_donors=cache$gse_donors,
            M0_repro=M0, M4_repro=M4)
saveRDS(s34, file.path(task, "00_admin", "analysis_cache.rds"))
say("analysis_cache.rds saved")

if (overall != "OK") {
  writeLines(c("M0_M4_REPRODUCED=ERROR",
    "REPRODUCTION_STATUS=REVISE_REQUIRED_MODEL_NOT_REPRODUCED",
    "DOWNSTREAM_MODELING=HALTED"),
    file.path(outR, "REPRODUCTION_STOP.txt"))
  say("STOP: exact reproduction failed. Downstream modeling halted.")
} else {
  say("Reproduction OK — downstream LODO/composition/CHyMErA analyses unlocked")
}
