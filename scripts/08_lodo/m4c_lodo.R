#!/usr/bin/env Rscript
# Composition-adjusted model M4C and its LODO.
# Fits the cell-composition-adjusted model M4C and LODO for the four Tier A events.
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).
suppressMessages({
  library(methods); library(lme4); library(pbkrtest); library(parallel)
})
set.seed(42)

task <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "34_robustness_and_composition")
dir32 <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "32_psychencode_sensitivity")
compD <- file.path(task, "04_cell_composition")
outD  <- file.path(task, "05_composition_adjusted_models")
dir.create(outD, showWarnings = FALSE, recursive = TRUE)
setD  <- file.path(dir32, "02_psychencode_sensitivity", "transcript_set_intermediates")
logf  <- file.path(task, "99_logs", "m4c_lodo.log")
say <- function(...) { msg <- paste0(...); cat(msg, "\n");
                       cat(msg, "\n", file = logf, append = TRUE) }
file.create(logf); cat("", file = logf)
wtab <- function(df, fn, dir = outD)
  write.table(df, file.path(dir, fn), sep = "\t", row.names = FALSE,
              quote = FALSE)

cache <- readRDS(file.path(task, "00_admin", "analysis_cache.rds"))
usage_matrix <- cache$usage_matrix
am           <- cache$analysis_meta
m4_cov       <- cache$m4_cov
m4_fix       <- cache$m4_fix
events       <- cache$events
primary19       <- cache$primary19
tierA_ids    <- cache$tierA_ids
gene_of      <- setNames(primary19$gene, primary19$HsaEX_ID)
disc_dir     <- setNames(primary19$discovery_dir, primary19$HsaEX_ID)
abs_dpsi     <- setNames(primary19$abs_delta_psi, primary19$HsaEX_ID)
final_tier   <- setNames(primary19$final_tier, primary19$HsaEX_ID)
full_M4      <- cache$M4_repro
epsilon      <- 1e-4
stopifnot(length(events) == 19, length(tierA_ids) == 4)

# ---------------- composition PCs (composition module) -------------------------
pc <- read.delim(file.path(compD, "COMPOSITION_PC_SCORES.tsv"),
                 stringsAsFactors = FALSE)
pc_cols <- grep("^CompPC", colnames(pc), value = TRUE)
k <- length(pc_cols)
say("composition PCs: k=", k, " (", paste(pc_cols, collapse = ", "),
    ")")
stopifnot(all(pc$sample_id %in% rownames(am)),
          all(rownames(am) %in% pc$sample_id))
pc_mat <- as.matrix(pc[match(rownames(am), pc$sample_id), pc_cols,
                       drop = FALSE])
rownames(pc_mat) <- rownames(am)
colnames(pc_mat) <- pc_cols

# M4C covariate block = M4 model-matrix covariates + composition PCs
m4c_cov <- cbind(m4_cov, pc_mat)
m4c_extra_fix <- paste(m4_fix, "+", paste(pc_cols, collapse = " + "))
m4c_formula_str <- paste0("usage_logit ~ dx_binary + region + Sex + Age + ",
                          "RIN", " ", m4c_extra_fix, " + (1 | subject)")
say("M4C formula: ", m4c_formula_str)

donors <- sort(unique(am$subject))
stopifnot(length(donors) == 80)

# ---------------- fitting machinery (verbatim from Robustness-1 LODO) --------
fit_event <- function(mdf, formula_full, formula_red, reml = TRUE) {
  fit_full <- tryCatch(lmer(formula_full, data = mdf, REML = reml),
                       error = function(e) e)
  if (inherits(fit_full, "error")) return(list(full = NULL, red = NULL,
                                               err_full = fit_full$message))
  fit_red <- tryCatch(lmer(formula_red, data = mdf, REML = reml),
                      error = function(e) NULL)
  list(full = fit_full, red = fit_red, err_full = NULL)
}
conv_of <- function(fit) {
  msgs <- tryCatch(fit@optinfo$conv$lme4$messages, error = function(e) NULL)
  if (is.null(msgs)) msgs <- character(0)
  data.frame(converged = length(msgs) == 0,
             singular  = isSingular(fit),
             conv_messages = paste(msgs, collapse = "; "),
             stringsAsFactors = FALSE)
}
run_model_sub <- function(eid, am_sub, extra_fix, cov_df_sub,
                          allow_fallback, removed_donor, removed_dx) {
  base <- data.frame(removed_donor = removed_donor,
                     removed_diagnosis = removed_dx,
                     model_id = NA_character_, event_id = eid,
                     gene = gene_of[[eid]],
                     n_samples = NA_integer_, n_donors = NA_integer_,
                     beta = NA_real_, SE = NA_real_,
                     CI_low = NA_real_, CI_high = NA_real_,
                     KR_P = NA_real_, KR_BH_FDR = NA_real_,
                     direction = NA_character_,
                     converged = NA, singular = NA, rank_deficient = NA,
                     used_fallback = NA,
                     error_reason = NA_character_,
                     stringsAsFactors = FALSE)
  usage_vals <- usage_matrix[, eid]
  row_keep <- rownames(am_sub)
  valid <- !is.na(usage_vals[row_keep])
  mdf <- am_sub[valid, ]
  if (!is.null(cov_df_sub)) {
    extra <- cov_df_sub[match(rownames(mdf), rownames(cov_df_sub)), ,
                        drop = FALSE]
    mdf <- cbind(mdf, extra)
  }
  if (nrow(mdf) < 10 || length(unique(mdf$subject)) < 3) {
    base$error_reason <- "INSUFFICIENT_SAMPLES_AFTER_DELETION"
    return(base)
  }
  mdf$usage <- usage_vals[row_keep][valid]
  mdf$usage_logit <- log((mdf$usage + epsilon) /
                         (1 - mdf$usage + epsilon))
  base_full <- "usage_logit ~ dx_binary + region + Sex + Age + RIN"
  base_red  <- "usage_logit ~ region + Sex + Age + RIN"
  f_full <- as.formula(paste(base_full, extra_fix, "+ (1|subject)"))
  f_red  <- as.formula(paste(base_red,  extra_fix, "+ (1|subject)"))
  fit_reml <- fit_event(mdf, f_full, f_red, reml = TRUE)
  used_fallback <- FALSE
  if (!is.null(fit_reml$err_full) || is.null(fit_reml$full)) {
    if (allow_fallback) {
      f_full <- as.formula(paste("usage_logit ~ dx_binary + region + Sex",
                                 extra_fix, "+ (1|subject)"))
      f_red  <- as.formula(paste("usage_logit ~ region + Sex",
                                 extra_fix, "+ (1|subject)"))
      fit_reml2 <- fit_event(mdf, f_full, f_red, reml = TRUE)
      if (is.null(fit_reml2$full) && !is.null(fit_reml$err_full)) {
        base$error_reason <- paste0("LMER_ERROR: ", fit_reml$err_full)
        base$rank_deficient <- grepl("rank", fit_reml$err_full,
                                     ignore.case = TRUE)
        return(base)
      }
      if (!is.null(fit_reml2$full)) {
        fit_reml <- fit_reml2; used_fallback <- TRUE
      } else if (!is.null(fit_reml$full)) {
        used_fallback <- FALSE
      } else {
        base$error_reason <- paste0("LMER_ERROR: ", fit_reml$err_full)
        base$rank_deficient <- grepl("rank", fit_reml$err_full,
                                     ignore.case = TRUE)
        return(base)
      }
    } else if (!is.null(fit_reml$err_full)) {
      base$error_reason <- paste0("LMER_ERROR: ", fit_reml$err_full)
      base$rank_deficient <- grepl("rank", fit_reml$err_full,
                                   ignore.case = TRUE)
      return(base)
    }
  }
  if (is.null(fit_reml$full)) {
    base$error_reason <- "FULL_MODEL_NULL"
    return(base)
  }
  ci <- conv_of(fit_reml$full)
  coefs <- summary(fit_reml$full)$coefficients
  if (!("dx_binary" %in% rownames(coefs))) {
    base$error_reason <- "DX_COEF_ALIASED"
    base$converged <- ci$converged; base$singular <- ci$singular
    return(base)
  }
  beta <- coefs["dx_binary", "Estimate"]
  se   <- coefs["dx_binary", "Std. Error"]
  p_kr <- NA_real_
  if (!is.null(fit_reml$red)) {
    kr <- tryCatch(KRmodcomp(fit_reml$full, fit_reml$red),
                   error = function(e) NULL)
    if (!is.null(kr)) p_kr <- kr$test$p.value[1]
  }
  base$n_samples <- nrow(mdf)
  base$n_donors <- length(unique(mdf$subject))
  base$beta <- beta; base$SE <- se
  base$CI_low <- beta - 1.96 * se; base$CI_high <- beta + 1.96 * se
  base$KR_P <- p_kr
  base$direction <- ifelse(beta > 0, "UP_IN_ASD", "DOWN_IN_ASD")
  base$converged <- ci$converged; base$singular <- ci$singular
  base$rank_deficient <- FALSE
  base$used_fallback <- used_fallback
  base
}

# ============================================================================
# SECTION 15: primary M4C (all 19 D0 events), full sample
# ============================================================================
say("\n=== SECTION 15: primary M4C (all 19 D0 events) ===")
m4c_rows <- list()
for (eid in events) {
  r <- tryCatch(run_model_sub(eid, am, m4c_extra_fix, m4c_cov, FALSE,
                              "NONE_full_sample", "NONE"),
                error = function(e) {
                  d <- data.frame(removed_donor = "NONE_full_sample",
                                  removed_diagnosis = "NONE",
                                  model_id = NA_character_,
                                  event_id = eid, gene = gene_of[[eid]],
                                  n_samples = NA_integer_,
                                  n_donors = NA_integer_,
                                  beta = NA_real_, SE = NA_real_,
                                  CI_low = NA_real_, CI_high = NA_real_,
                                  KR_P = NA_real_, KR_BH_FDR = NA_real_,
                                  direction = NA_character_,
                                  converged = NA, singular = NA,
                                  rank_deficient = NA, used_fallback = NA,
                                  error_reason = paste0("UNEXPECTED: ",
                                                          e$message),
                                  stringsAsFactors = FALSE)
                  d })
  r$model_id <- "M4C_primary"
  m4c_rows[[length(m4c_rows) + 1]] <- r
  say("  M4C ", eid, " (", gene_of[[eid]], "): beta=",
      ifelse(is.na(r$beta), "NA", signif(r$beta, 4)),
      ifelse(is.na(r$error_reason), "", paste0(" [", r$error_reason, "]")))
}
m4c <- do.call(rbind, m4c_rows)
rownames(m4c) <- NULL
m4c$KR_BH_FDR <- p.adjust(m4c$KR_P, method = "BH")
m4c$formula <- m4c_formula_str
m4c$complete_case_loss <- nrow(am) - m4c$n_samples
wtab(m4c, "M4C_EVENT_RESULTS.tsv")
n_m4c_ok <- sum(is.na(m4c$error_reason))
say("primary M4C fitted: ", n_m4c_ok, "/19")
m4c_results_complete <- ifelse(n_m4c_ok == 19, "OK", "PARTIAL_OR_NOT_FEASIBLE")

# ---------------- Section 15: M4 vs M4C comparison -------------------------
say("\n=== SECTION 15: M4 vs M4C comparison ===")
m4_beta <- setNames(full_M4$beta_ASD, full_M4$HsaEX_ID)
m4_dir  <- setNames(full_M4$direction, full_M4$HsaEX_ID)
m4_fdr  <- setNames(full_M4$BH_FDR_KR, full_M4$HsaEX_ID)
m4_se   <- setNames(full_M4$SE, full_M4$HsaEX_ID)
m4c_beta <- setNames(m4c$beta, m4c$event_id)
m4c_dir  <- setNames(m4c$direction, m4c$event_id)
m4c_fdr  <- setNames(m4c$KR_BH_FDR, m4c$event_id)

common <- intersect(events, m4c$event_id[!is.na(m4c$beta)])
dir_retention_all <- mean(m4c_dir[common] == m4_dir[common])
b4 <- m4_beta[common]; b4c <- m4c_beta[common]
pearson_r  <- if (length(common) >= 3) cor(b4, b4c) else NA_real_
spearman_r <- if (length(common) >= 3)
  cor(b4, b4c, method = "spearman") else NA_real_
med_retention <- median(abs(b4c) / abs(b4))
tA_common <- intersect(tierA_ids, common)
tierA_dir_4of4 <- sum(m4c_dir[tA_common] == m4_dir[tA_common])
tierA_beta_retention <- abs(m4c_beta[tA_common]) / abs(m4_beta[tA_common])
tierA_fdr05_m4  <- sum(m4_fdr[tierA_ids] < 0.05, na.rm = TRUE)
tierA_fdr05_m4c <- sum(m4c_fdr[tA_common] < 0.05, na.rm = TRUE)
# discovery-direction concordance (15/19 reference is M0-based;
# here we report M4C vs the discovery direction for all 19)
disc_conc_m4c <- sum(m4c_dir[common] == disc_dir[common])
# one-event-per-gene (OEPG): per gene pick event with max abs_delta_psi
ord <- order(-abs_dpsi[common], common)
oepg_events <- common[ord][!duplicated(gene_of[common[ord]])]
oepg_conc <- sum(m4c_dir[oepg_events] == disc_dir[oepg_events])

say("all-event direction retention: ", signif(dir_retention_all, 4),
    " (", sum(m4c_dir[common] == m4_dir[common]), "/", length(common), ")")
say("M4-M4C beta Pearson r=", signif(pearson_r, 4),
    " Spearman rho=", signif(spearman_r, 4))
say("median effect retention: ", signif(med_retention, 4))
say("Tier A direction: ", tierA_dir_4of4, "/4")
say("Tier A beta retention: ", paste(signif(tierA_beta_retention, 3),
                                     collapse = ", "))
say("Tier A FDR<0.05: M4=", tierA_fdr05_m4, " M4C=", tierA_fdr05_m4c)
say("discovery-direction concordance (M4C): ", disc_conc_m4c, "/19")
say("OEPG concordance (M4C): ", oepg_conc, "/", length(oepg_events))

# ---------------- per-event interpretation (pre-specified rules) -----------
# Pre-specified magnitude threshold for "attenuated" classification.
RETAIN_MAGNITUDE <- 0.8
interp <- character(length(common))
names(interp) <- common
for (eid in common) {
  b0 <- m4_beta[eid]; b1 <- m4c_beta[eid]
  if (is.na(b0) || is.na(b1)) { interp[eid] <- "NOT_ASSESSABLE"; next }
  d0 <- m4_dir[eid]; d1 <- m4c_dir[eid]
  sig0 <- (!is.na(m4_fdr[eid])) && m4_fdr[eid] < 0.05
  sig1 <- (!is.na(m4c_fdr[eid])) && m4c_fdr[eid] < 0.05
  reten <- abs(b1) / abs(b0)
  if (d1 != d0) {
    interp[eid] <- "COMPOSITION_SENSITIVE_DIRECTION"
  } else if (sig0 && !sig1) {
    interp[eid] <- "COMPOSITION_SENSITIVE_SIGNIFICANCE_ONLY"
  } else if (reten < RETAIN_MAGNITUDE) {
    interp[eid] <- "COMPOSITION_ROBUST_DIRECTION_ATTENUATED_MAGNITUDE"
  } else {
    interp[eid] <- "COMPOSITION_ROBUST_DIRECTION_AND_MAGNITUDE"
  }
}
# events that failed M4C fit
for (eid in setdiff(events, common)) interp[eid] <- "NOT_ASSESSABLE"

cmp <- data.frame(
  event_id = events, gene = gene_of[events],
  M4_beta = m4_beta[events], M4_SE = m4_se[events],
  M4_direction = m4_dir[events], M4_KR_FDR = m4_fdr[events],
  M4C_beta = m4c_beta[events], M4C_SE = m4c$SE[match(events, m4c$event_id)],
  M4C_direction = m4c_dir[events],
  M4C_KR_FDR = m4c_fdr[events],
  direction_retained = m4c_dir[events] == m4_dir[events],
  abs_effect_retention = abs(m4c_beta[events]) / abs(m4_beta[events]),
  discovery_dir = disc_dir[events],
  M4C_matches_discovery = m4c_dir[events] == disc_dir[events],
  final_tier = final_tier[events],
  interpretation = interp[events],
  stringsAsFactors = FALSE)
wtab(cmp, "M4_VS_M4C_SUMMARY.tsv")

# ============================================================================
# SECTION 16: explicit-fraction sensitivity (secondary) M4F
#   M4F = M4 + 6 of 7 broad fractions (one class omitted as reference to
#   avoid the sum-to-1 collinearity). OMITTED class = OPC (documented; the
#   model column space is invariant to which single class is omitted).
# ============================================================================
say("\n=== SECTION 16: explicit-fraction sensitivity (M4F, secondary) ===")
frac <- read.delim(file.path(compD, "COMPOSITION_FRACTIONS_HARMONIZED.tsv"),
                   stringsAsFactors = FALSE)
CLASSES <- c("Excitatory_neuron", "Inhibitory_neuron", "Astrocyte",
             "Oligodendrocyte", "OPC", "Microglia_immune",
             "Endothelial_mural")
OMIT <- "OPC"
frac_mat <- as.matrix(frac[match(rownames(am), frac$sample_id),
                           setdiff(CLASSES, OMIT), drop = FALSE])
rownames(frac_mat) <- rownames(am)
colnames(frac_mat) <- paste0("frac_", setdiff(CLASSES, OMIT))
m4f_cov <- cbind(m4_cov, frac_mat)
m4f_extra_fix <- paste(m4_fix, "+", paste(colnames(frac_mat),
                                          collapse = " + "))
say("M4F omitted reference class: ", OMIT, " ; included fractions: ",
    paste(colnames(frac_mat), collapse = ", "))
m4f_rows <- list()
for (eid in events) {
  r <- tryCatch(run_model_sub(eid, am, m4f_extra_fix, m4f_cov, FALSE,
                              "NONE_full_sample", "NONE"),
                error = function(e) NULL)
  if (is.null(r)) {
    r <- data.frame(removed_donor = "NONE_full_sample",
                    removed_diagnosis = "NONE", model_id = NA_character_,
                    event_id = eid, gene = gene_of[[eid]],
                    n_samples = NA_integer_, n_donors = NA_integer_,
                    beta = NA_real_, SE = NA_real_, CI_low = NA_real_,
                    CI_high = NA_real_, KR_P = NA_real_,
                    KR_BH_FDR = NA_real_, direction = NA_character_,
                    converged = NA, singular = NA, rank_deficient = NA,
                    used_fallback = NA,
                    error_reason = "UNEXPECTED_ERROR",
                    stringsAsFactors = FALSE)
  }
  r$model_id <- "M4F_explicit_fractions"
  m4f_rows[[length(m4f_rows) + 1]] <- r
}
m4f <- do.call(rbind, m4f_rows); rownames(m4f) <- NULL
m4f$KR_BH_FDR <- p.adjust(m4f$KR_P, method = "BH")
m4f$formula <- paste0("usage_logit ~ dx_binary + region + Sex + Age + RIN ",
                      m4f_extra_fix, " + (1 | subject)")
m4f$omitted_reference_class <- OMIT
m4f$note <- "SECONDARY exploratory explicit-fraction model; column space invariant to omitted class"
wtab(m4f, "M4F_EVENT_RESULTS.tsv")
say("M4F fitted: ", sum(is.na(m4f$error_reason)), "/19")

# ============================================================================
# Amendment Sec 5/6: primary- vs sensitivity-reference M4C comparison
# ============================================================================
say("\n=== Amendment Sec 5/6: primary vs sensitivity-reference M4C ===")
sens_frac <- read.delim(file.path(compD,
  "COMPOSITION_FRACTIONS_SENSITIVITY_OVERLAP_EXCLUDED.tsv"),
  stringsAsFactors = FALSE)
Fs <- as.matrix(sens_frac[match(rownames(am), sens_frac$sample_id), CLASSES])
rownames(Fs) <- rownames(am)
pos_s <- Fs[Fs > 0]
delta_s <- 0.5 * min(pos_s)
Frs <- Fs; Frs[Frs == 0] <- delta_s
gm_s <- exp(rowMeans(log(Frs)))
clr_s <- log(Frs) - log(gm_s)
pca_s <- prcomp(clr_s, center = TRUE, scale. = FALSE)
cve_s <- pca_s$sdev^2 / sum(pca_s$sdev^2)
# retain the SAME k as the final primary decision for direct comparability
pc_sens <- pca_s$x[, seq_len(k), drop = FALSE]
colnames(pc_sens) <- paste0("CompPC", seq_len(k), "_sens")
say("sensitivity PC variance explained: ",
    paste(round(cve_s[1:min(6, length(cve_s))], 3), collapse = ", "))
write.table(data.frame(sample_id = rownames(pc_sens), pc_sens,
                       check.names = FALSE),
            file.path(compD, "COMPOSITION_PC_SCORES_SENSITIVITY.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
m4cs_cov <- cbind(m4_cov, pc_sens)
m4cs_extra_fix <- paste(m4_fix, "+", paste(colnames(pc_sens),
                                           collapse = " + "))
m4cs_rows <- list()
for (eid in events) {
  r <- tryCatch(run_model_sub(eid, am, m4cs_extra_fix, m4cs_cov, FALSE,
                              "NONE_full_sample", "NONE"),
                error = function(e) NULL)
  if (is.null(r)) {
    r <- data.frame(removed_donor = "NONE_full_sample",
                    removed_diagnosis = "NONE", model_id = NA_character_,
                    event_id = eid, gene = gene_of[[eid]],
                    n_samples = NA_integer_, n_donors = NA_integer_,
                    beta = NA_real_, SE = NA_real_, CI_low = NA_real_,
                    CI_high = NA_real_, KR_P = NA_real_,
                    KR_BH_FDR = NA_real_, direction = NA_character_,
                    converged = NA, singular = NA, rank_deficient = NA,
                    used_fallback = NA, error_reason = "UNEXPECTED_ERROR",
                    stringsAsFactors = FALSE)
  }
  r$model_id <- "M4C_sensitivity_reference"
  m4cs_rows[[length(m4cs_rows) + 1]] <- r
}
m4cs <- do.call(rbind, m4cs_rows); rownames(m4cs) <- NULL
m4cs$KR_BH_FDR <- p.adjust(m4cs$KR_P, method = "BH")
pv <- data.frame(
  event_id = events, gene = gene_of[events],
  M4C_primary_beta = m4c_beta[events],
  M4C_sens_beta = setNames(m4cs$beta, m4cs$event_id)[events],
  M4C_primary_direction = m4c_dir[events],
  M4C_sens_direction = setNames(m4cs$direction, m4cs$event_id)[events],
  M4C_primary_KR_FDR = m4c_fdr[events],
  M4C_sens_KR_FDR = setNames(m4cs$KR_BH_FDR, m4cs$event_id)[events],
  stringsAsFactors = FALSE)
pv$direction_agree <- pv$M4C_primary_direction == pv$M4C_sens_direction
both_ok <- !is.na(pv$M4C_primary_beta) & !is.na(pv$M4C_sens_beta)
pv_beta_cor <- if (sum(both_ok) >= 3)
  cor(pv$M4C_primary_beta[both_ok], pv$M4C_sens_beta[both_ok]) else NA
pv$note <- ""
wtab(pv, "M4C_PRIMARY_VS_SENSITIVITY.tsv")
say("primary vs sensitivity M4C: beta Pearson r=",
    ifelse(is.na(pv_beta_cor), "NA", signif(pv_beta_cor, 5)),
    " ; direction agreement ", sum(pv$direction_agree, na.rm = TRUE),
    "/", nrow(pv))

# ============================================================================
# SECTION 17: D0-D3 under composition adjustment (reuse final sets)
# ============================================================================
say("\n=== SECTION 17: D0-D3 under M4C covariates ===")
rsem_member    <- cache$rsem_member
member_base    <- cache$member_base_ids
effLen_member  <- cache$effLen_member
samp_idx_raw   <- cache$samp_idx_raw
incl_tx        <- cache$incl_tx
excl_tx        <- cache$excl_tx
idx_by_base <- setNames(seq_along(member_base), member_base)
usage_for <- function(eid, incl, excl) {
  ii <- idx_by_base[incl[incl %in% names(idx_by_base)]]
  ee <- idx_by_base[excl[excl %in% names(idx_by_base)]]
  if (length(ii) == 0 || length(ee) == 0) return(rep(NA, nrow(am)))
  ir <- colSums(as.matrix(rsem_member[ii, samp_idx_raw, drop = FALSE]) /
                effLen_member[ii], na.rm = TRUE)
  er <- colSums(as.matrix(rsem_member[ee, samp_idx_raw, drop = FALSE]) /
                effLen_member[ee], na.rm = TRUE)
  tot <- ir + er
  ifelse(tot > 0, ir / tot, NA)
}
# M4C fixed-effects string for set models (dx_binary included)
m4c_set_fix_full <- paste("dx_binary + region + Sex + Age + RIN",
                          m4c_extra_fix)
fit_set_event_m4c <- function(u, cov_df) {
  ok <- !is.na(u)
  mdf <- am[ok, ]
  extra <- cov_df[match(rownames(mdf), rownames(cov_df)), , drop = FALSE]
  mdf <- cbind(mdf, extra)
  if (nrow(mdf) < 20) return(NULL)
  mdf$usage_logit <- log((u[ok] + epsilon) / (1 - u[ok] + epsilon))
  f_full <- as.formula(paste("usage_logit ~", m4c_set_fix_full,
                             "+ (1|subject)"))
  f_red  <- as.formula(paste("usage_logit ~",
    sub("dx_binary \\+ ", "", m4c_set_fix_full), "+ (1|subject)"))
  fr <- tryCatch(lmer(f_full, data = mdf, REML = TRUE),
                 error = function(e) NULL)
  if (is.null(fr)) return(NULL)
  cf <- summary(fr)$coefficients
  if (!("dx_binary" %in% rownames(cf))) return(NULL)
  beta <- cf["dx_binary", "Estimate"]; se <- cf["dx_binary", "Std. Error"]
  kr <- tryCatch(KRmodcomp(fr, tryCatch(lmer(f_red, data = mdf, REML = TRUE),
             error = function(e) NULL)), error = function(e) NULL)
  pkr <- if (!is.null(kr)) kr$test$p.value[1] else NA_real_
  data.frame(beta_ASD = beta, SE = se,
             CI_lo = beta - 1.96 * se, CI_hi = beta + 1.96 * se,
             P_KR = pkr, direction = ifelse(beta > 0, "UP_IN_ASD",
                                            "DOWN_IN_ASD"),
             n_samples = nrow(mdf), n_donors = length(unique(mdf$subject)),
             singular = isSingular(fr), stringsAsFactors = FALSE)
}
d0_incl <- setNames(lapply(events, function(e)
  incl_tx$transcript_id[incl_tx$HsaEX_ID == e]), events)
d0_excl <- setNames(lapply(events, function(e)
  excl_tx$transcript_id[excl_tx$HsaEX_ID == e]), events)
set_lists <- list()
for (eid in events) {
  sf <- file.path(setD, paste0("sets_", eid, ".rds"))
  set_lists[[eid]] <- if (file.exists(sf)) readRDS(sf) else NULL
}
run_def_m4c <- function(def_name, get_lists) {
  res <- data.frame()
  for (eid in events) {
    L <- get_lists(eid)
    if (is.null(L)) next
    u <- usage_for(eid, L$incl, L$excl)
    r <- fit_set_event_m4c(u, cov_df = m4c_cov)
    if (!is.null(r)) { r$HsaEX_ID <- eid; r$definition <- def_name
                       res <- rbind(res, r) }
  }
  if (nrow(res) > 1) res$BH_FDR_KR <- p.adjust(res$P_KR, method = "BH")
  else res$BH_FDR_KR <- NA_real_
  say("  ", def_name, " (M4C): ", nrow(res), " events fitted")
  res
}
D0s <- run_def_m4c("D0", function(e)
  list(incl = d0_incl[[e]], excl = d0_excl[[e]]))
D1s <- run_def_m4c("D1", function(e) {
  s <- set_lists[[e]]; if (is.null(s)) return(NULL)
  if (length(s$d1_incl) == 0 || length(s$d1_excl) == 0) return(NULL)
  list(incl = s$d1_incl, excl = s$d1_excl) })
D2s <- run_def_m4c("D2", function(e) {
  s <- set_lists[[e]]; if (is.null(s)) return(NULL)
  if (is.na(s$d2_incl) || is.na(s$d2_excl)) return(NULL)
  list(incl = s$d2_incl, excl = s$d2_excl) })
D3s <- run_def_m4c("D3", function(e) {
  s <- set_lists[[e]]; if (is.null(s)) return(NULL)
  if (length(s$d3_incl) == 0 || length(s$d3_excl) == 0) return(NULL)
  list(incl = s$d3_incl, excl = s$d3_excl) })
allres <- do.call(rbind, list(D0s, D1s, D2s, D3s))
allres <- merge(allres, primary19[, c("HsaEX_ID", "gene", "discovery_dir",
                                   "final_tier")], by = "HsaEX_ID")
wtab(allres, "M4C_D0_D3_EVENT_RESULTS.tsv")

# ---- correlation of M4C D0-D3 betas with non-composition results ----
d32 <- read.delim(file.path(dir32, "02_psychencode_sensitivity",
                            "TRANSCRIPT_SET_EVENT_RESULTS.tsv"),
                  stringsAsFactors = FALSE)
corr_rows <- list()
for (dn in c("D0", "D1", "D2", "D3")) {
  a <- allres[allres$definition == dn, ]
  b <- d32[d32$definition == dn, ]
  ce <- intersect(a$HsaEX_ID, b$HsaEX_ID)
  if (length(ce) >= 3) {
    r <- cor(a$beta_ASD[match(ce, a$HsaEX_ID)],
             b$beta_ASD[match(ce, b$HsaEX_ID)])
    fl <- sum(a$direction[match(ce, a$HsaEX_ID)] !=
              b$direction[match(ce, b$HsaEX_ID)])
  } else { r <- NA_real_; fl <- NA_integer_ }
  corr_rows[[length(corr_rows) + 1]] <- data.frame(
    definition = dn, n_events_m4c = nrow(a), n_events_dir32 = nrow(b),
    n_compared = length(ce), beta_correlation_vs_non_composition = r,
    n_direction_flips_vs_non_composition = fl, stringsAsFactors = FALSE)
}
corr_tab <- do.call(rbind, corr_rows)

# ---- per-definition set summary (M4C) -------------------------------------
sum_def_m4c <- function(res, name) {
  if (is.null(res) || nrow(res) == 0)
    return(data.frame(definition = name, n_analyzable = 0,
      n_analyzable_Tier_A = 0, direction_concordance = "0/0",
      exact_binomial_P = NA_real_, oepg_concordance = "0/0",
      tierA_direction_flags = "", tierA_KR_FDR05 = 0, n_KR_FDR05 = 0,
      note = "none_fitted", stringsAsFactors = FALSE))
  m <- merge(res, primary19[, c("HsaEX_ID", "discovery_dir", "final_tier",
                             "abs_delta_psi", "gene")], by = "HsaEX_ID")
  m$conc <- m$direction == m$discovery_dir
  bp <- pbinom(sum(m$conc) - 1, nrow(m), 0.5, lower.tail = FALSE)
  ms <- m[order(-m$abs_delta_psi, m$HsaEX_ID), ]
  oepg <- do.call(rbind, lapply(split(ms, ms$gene), function(d) d[1, ]))
  tA <- m[m$final_tier == "TIER_A_CROSS_COHORT_KR_FDR05", ]
  data.frame(definition = name, n_analyzable = nrow(m),
    n_analyzable_Tier_A = nrow(tA),
    direction_concordance = paste0(sum(m$conc), "/", nrow(m)),
    exact_binomial_P = signif(bp, 6),
    oepg_concordance = paste0(sum(oepg$conc), "/", nrow(oepg)),
    tierA_direction_flags = paste(tA$conc, collapse = ""),
    tierA_KR_FDR05 = sum(tA$BH_FDR_KR < 0.05, na.rm = TRUE),
    n_KR_FDR05 = sum(m$BH_FDR_KR < 0.05, na.rm = TRUE),
    note = "", stringsAsFactors = FALSE)
}
summ_m4c <- do.call(rbind, list(sum_def_m4c(D0s, "D0"),
                                sum_def_m4c(D1s, "D1"),
                                sum_def_m4c(D2s, "D2"),
                                sum_def_m4c(D3s, "D3")))
set_summary_out <- merge(summ_m4c, corr_tab, by = "definition",
                         all = TRUE, sort = FALSE)
wtab(set_summary_out, "M4C_D0_D3_SET_SUMMARY.tsv")

# ============================================================================
# SECTION 18: Tier A M4C LODO (80 donors x 4 Tier A events = 320 fits)
# ============================================================================
say("\n=== SECTION 18: Tier A M4C LODO ===")
dx_of_donor <- setNames(
  tapply(am$dx_binary, am$subject, function(x) unique(x)[1]), donors)
lodo_donor <- function(did) {
  am_sub  <- am[am$subject != did, ]
  cov_sub <- m4c_cov[rownames(m4c_cov) %in% rownames(am_sub), , drop = FALSE]
  rdx <- ifelse(dx_of_donor[[did]] == 1, "ASD", "CTL")
  rows <- list()
  for (eid in tierA_ids) {
    r <- tryCatch(run_model_sub(eid, am_sub, m4c_extra_fix, cov_sub, FALSE,
                                did, rdx),
                  error = function(e) {
                    data.frame(removed_donor = did, removed_diagnosis = rdx,
                               model_id = NA_character_, event_id = eid,
                               gene = gene_of[[eid]],
                               n_samples = NA_integer_, n_donors = NA_integer_,
                               beta = NA_real_, SE = NA_real_,
                               CI_low = NA_real_, CI_high = NA_real_,
                               KR_P = NA_real_, KR_BH_FDR = NA_real_,
                               direction = NA_character_,
                               converged = NA, singular = NA,
                               rank_deficient = NA, used_fallback = NA,
                               error_reason = paste0("UNEXPECTED: ",
                                                       e$message),
                               stringsAsFactors = FALSE) })
    r$model_id <- "M4C_primary_LODO"
    rows[[length(rows) + 1]] <- r
  }
  do.call(rbind, rows)
}
t0 <- Sys.time()
N_CORES <- 10
lodo_list <- mclapply(donors, lodo_donor, mc.cores = N_CORES,
                      mc.preschedule = TRUE)
lodo_long <- do.call(rbind, lodo_list)
rownames(lodo_long) <- NULL
say("Tier A M4C LODO fits done in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    " min ; rows=", nrow(lodo_long), " (expect 320)")
# BH within each donor-deletion family
lodo_long$KR_BH_FDR <- NA_real_
for (key in unique(lodo_long$removed_donor)) {
  w <- which(lodo_long$removed_donor == key)
  lodo_long$KR_BH_FDR[w] <- p.adjust(lodo_long$KR_P[w], method = "BH")
}
wtab(lodo_long, "TIER_A_M4C_LODO_RESULTS.tsv")
n_lodo_error <- sum(!is.na(lodo_long$error_reason))
say("Tier A M4C LODO failed fits: ", n_lodo_error)

# ---- Tier A M4C LODO summary (influence metrics, Section 9 criteria) -----
m4c_full_beta <- setNames(m4c$beta, m4c$event_id)[tierA_ids]
m4c_full_se   <- setNames(m4c$SE, m4c$event_id)[tierA_ids]
m4c_full_dir  <- setNames(m4c$direction, m4c$event_id)[tierA_ids]
lodo_summ <- list()
for (eid in tierA_ids) {
  bf <- m4c_full_beta[eid]; sef <- m4c_full_se[eid]
  dfull <- m4c_full_dir[eid]
  s <- lodo_long[lodo_long$event_id == eid &
                 is.na(lodo_long$error_reason), ]
  if (nrow(s) == 0 || is.na(bf) || abs(bf) < 1e-6) {
    lodo_summ[[length(lodo_summ) + 1]] <- data.frame(
      model_id = "M4C_primary_LODO", event_id = eid, gene = gene_of[[eid]],
      full_beta = bf, full_SE = sef, full_direction = dfull,
      n_deletions_ok = nrow(s), NA_REASON = "FULL_BETA_NEAR_ZERO_OR_ERROR",
      stringsAsFactors = FALSE)
    next
  }
  retention <- abs(s$beta) / abs(bf)
  dfb <- (bf - s$beta) / sef
  dirp <- s$direction == dfull
  lodo_summ[[length(lodo_summ) + 1]] <- data.frame(
    model_id = "M4C_primary_LODO", event_id = eid, gene = gene_of[[eid]],
    full_beta = bf, full_SE = sef, full_direction = dfull,
    n_deletions_ok = nrow(s), NA_REASON = "",
    LODO_beta_min = min(s$beta), LODO_beta_max = max(s$beta),
    direction_preservation_fraction = mean(dirp),
    direction_reversal_count = sum(!dirp),
    min_abs_effect_retention = min(retention),
    median_abs_effect_retention = median(retention),
    max_abs_beta_change = max(abs(s$beta - bf)),
    max_rel_beta_change = max(abs(s$beta - bf) / abs(bf)),
    max_abs_DFBETA = max(abs(dfb)),
    most_influential_donor = s$removed_donor[which.max(abs(dfb))],
    most_influential_donor_dx =
      s$removed_diagnosis[which.max(abs(dfb))],
    nominal_p_lt_005_fraction = mean(s$KR_P < 0.05, na.rm = TRUE),
    bh_fdr_lt_005_fraction = mean(s$KR_BH_FDR < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE)
}
lodo_summ <- do.call(rbind, lodo_summ)
lodo_summ$direction_preserved <-
  lodo_summ$direction_preservation_fraction >= 0.95
lodo_summ$no_reversal <- lodo_summ$direction_reversal_count == 0
lodo_summ$min_retention_ok <- lodo_summ$min_abs_effect_retention >= 0.50
lodo_summ$no_high_influence <- lodo_summ$max_abs_DFBETA < 1.0
lodo_summ$TIER_A_M4C_LODO_CONFIRMED <- lodo_summ$direction_preserved &
  lodo_summ$no_reversal & lodo_summ$min_retention_ok &
  lodo_summ$no_high_influence
wtab(lodo_summ, "TIER_A_M4C_LODO_SUMMARY.tsv")
tierA_m4c_direction_stable_n <-
  sum(lodo_summ$direction_preservation_fraction >= 0.95, na.rm = TRUE)
tierA_m4c_lodo_overall <-
  ifelse(nrow(lodo_summ) == 4 && all(lodo_summ$TIER_A_M4C_LODO_CONFIRMED),
         "OK", "ERROR")
say("Tier A M4C LODO overall: ", tierA_m4c_lodo_overall)

# ---------------- phase values for Robustness-7 -------------------------------
tierA_m4c_fdr05_n <- sum(m4c_fdr[tierA_ids] < 0.05, na.rm = TRUE)
phase <- c(
  "# Robustness-3 composition-adjusted model phase values",
  paste0("M4C_RESULTS_COMPLETE=", m4c_results_complete),
  paste0("M4C_EVENTS_FITTED=", n_m4c_ok, "/19"),
  paste0("M4C_FORMULA=", m4c_formula_str),
  paste0("COMPOSITION_PCS_USED_k=", k),
  paste0("M4_VS_M4C_ALL_EVENT_DIRECTION_RETENTION=",
         signif(dir_retention_all, 4)),
  paste0("M4_VS_M4C_BETA_PEARSON_R=",
         ifelse(is.na(pearson_r), "NA", signif(pearson_r, 4))),
  paste0("M4_VS_M4C_BETA_SPEARMAN_RHO=",
         ifelse(is.na(spearman_r), "NA", signif(spearman_r, 4))),
  paste0("M4_VS_M4C_MEDIAN_EFFECT_RETENTION=",
         signif(med_retention, 4)),
  paste0("TIER_A_M4C_DIRECTION_4OF4=", tierA_dir_4of4, "/4"),
  paste0("TIER_A_M4C_FDR05_N=", tierA_m4c_fdr05_n),
  paste0("TIER_A_M4C_DIRECTION_STABLE_N=", tierA_m4c_direction_stable_n,
         "/4"),
  paste0("TIER_A_M4C_LODO_OVERALL=", tierA_m4c_lodo_overall),
  paste0("M4C_DISCOVERY_DIRECTION_CONCORDANCE=", disc_conc_m4c, "/19"),
  paste0("M4C_OEPG_CONCORDANCE=", oepg_conc, "/", length(oepg_events)),
  paste0("M4C_PRIMARY_VS_SENSITIVITY_BETA_R=",
         ifelse(is.na(pv_beta_cor), "NA", signif(pv_beta_cor, 5))))
writeLines(phase, file.path(outD, "m4c_check_values.txt"))
say("\nM4C analysis complete. M4C_RESULTS_COMPLETE=", m4c_results_complete)
