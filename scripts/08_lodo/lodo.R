#!/usr/bin/env Rscript
# Leave-one-donor-out robustness.
# M0/M4 LODO across all 80 donors (3,040 fits), Tier A influence metrics and DFBETA diagnostics.
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).
suppressMessages({
  library(methods); library(lme4); library(pbkrtest); library(parallel)
})
set.seed(42)

task <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "34_robustness_and_composition")
outD <- file.path(task, "03_lodo")
ckpt <- file.path(outD, "checkpoints")
dir.create(ckpt, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(task, "99_logs", "lodo.log")
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
full_M0      <- cache$M0_repro
full_M4      <- cache$M4_repro
epsilon      <- 1e-4

donors <- sort(unique(am$subject))
dx_of_donor <- setNames(
  tapply(am$dx_binary, am$subject, function(x) unique(x)[1]), donors)
stopifnot(length(donors) == 80, length(events) == 19)
say("LODO: ", length(donors), " donors x ", length(events),
    " events x 2 models = ", length(donors) * length(events) * 2, " fits")

# ---------------- exact fitting code (as the sensitivity analysis_0b) --------
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
  base$model_id <- NA_character_  # filled by caller
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

donor_job <- function(did) {
  ck <- file.path(ckpt, paste0("donor_", did, ".rds"))
  if (file.exists(ck)) return(invisible("skip"))
  am_sub  <- am[am$subject != did, ]
  cov_sub <- m4_cov[rownames(m4_cov) %in% rownames(am_sub), , drop = FALSE]
  rdx <- ifelse(dx_of_donor[[did]] == 1, "ASD", "CTL")
  rows <- list()
  for (eid in events) {
    r0 <- tryCatch(run_model_sub(eid, am_sub, "", NULL, TRUE, did, rdx),
                   error = function(e) {
                     d <- data.frame(removed_donor = did,
                                     removed_diagnosis = rdx,
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
    r0$model_id <- "M0_primary"
    r4 <- tryCatch(run_model_sub(eid, am_sub, m4_fix, cov_sub, FALSE,
                                 did, rdx),
                   error = function(e) {
                     d <- data.frame(removed_donor = did,
                                     removed_diagnosis = rdx,
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
    r4$model_id <- "M4_model_matrix_covariates"
    rows[[length(rows) + 1]] <- r0
    rows[[length(rows) + 1]] <- r4
  }
  saveRDS(do.call(rbind, rows), ck)
  invisible("done")
}

# ---------------- run with resumable parallel checkpointing -----------------
done0 <- length(list.files(ckpt, pattern = "^donor_.*\\.rds$"))
say("checkpoints already present: ", done0, " / ", length(donors))
t0 <- Sys.time()
N_CORES <- 10
if (done0 < length(donors)) {
  res <- mclapply(donors, donor_job, mc.cores = N_CORES,
                  mc.preschedule = TRUE)
  n_err <- sum(vapply(res, function(x) identical(x, "skip"), logical(1))
               == FALSE & vapply(res, is.character, logical(1)) == FALSE)
  say("donor jobs finished in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      " min")
}
ck_files <- sort(list.files(ckpt, pattern = "^donor_.*\\.rds$",
                            full.names = TRUE))
stopifnot(length(ck_files) == length(donors))
long <- do.call(rbind, lapply(ck_files, readRDS))
rownames(long) <- NULL
say("long table: ", nrow(long), " fit rows (expect 3040)")

# ---------------- BH within donor-deletion x model family -------------------
long$KR_BH_FDR <- NA_real_
for (key in unique(paste(long$removed_donor, long$model_id))) {
  w <- which(paste(long$removed_donor, long$model_id) == key)
  long$KR_BH_FDR[w] <- p.adjust(long$KR_P[w], method = "BH")
}

# ---------------- job manifest ----------------------------------------------
long$status <- ifelse(is.na(long$error_reason), "DONE", "ERROR")
manifest <- long[, c("removed_donor", "removed_diagnosis", "model_id",
                     "event_id", "gene", "n_samples", "n_donors", "status",
                     "error_reason", "converged", "singular",
                     "rank_deficient", "used_fallback")]
wtab(manifest, "LODO_JOB_MANIFEST.tsv")

errors <- long[!is.na(long$error_reason),
              c("removed_donor", "removed_diagnosis", "model_id",
                "event_id", "gene", "error_reason", "n_samples")]
wtab(errors, "LODO_MODEL_ERRORS.tsv")
say("failed fits: ", nrow(errors))

conv <- data.frame(model_id = character(), n_fits = integer(),
                   n_done = integer(), n_converged_ok = integer(),
                   n_convergence_warning = integer(),
                   n_singular = integer(), n_rank_deficient = integer(),
                   n_other_error = integer(), n_used_fallback = integer(),
                   stringsAsFactors = FALSE)
for (mid in c("M0_primary", "M4_model_matrix_covariates")) {
  s <- long[long$model_id == mid, ]
  conv <- rbind(conv, data.frame(
    model_id = mid, n_fits = nrow(s),
    n_done = sum(is.na(s$error_reason)),
    n_converged_ok = sum(s$converged == TRUE, na.rm = TRUE),
    n_convergence_warning = sum(s$converged == FALSE, na.rm = TRUE),
    n_singular = sum(s$singular == TRUE, na.rm = TRUE),
    n_rank_deficient = sum(s$rank_deficient == TRUE, na.rm = TRUE),
    n_other_error = sum(!is.na(s$error_reason) &
                          s$rank_deficient != TRUE, na.rm = TRUE),
    n_used_fallback = sum(s$used_fallback == TRUE, na.rm = TRUE),
    stringsAsFactors = FALSE))
}
wtab(conv, "LODO_CONVERGENCE_SUMMARY.tsv")

# ---------------- influence metrics (Section 8) -----------------------------
full_beta <- rbind(
  data.frame(model_id = "M0_primary",
             event_id = full_M0$HsaEX_ID, beta_full = full_M0$beta_ASD,
             SE_full = full_M0$SE, direction_full = full_M0$direction),
  data.frame(model_id = "M4_model_matrix_covariates",
             event_id = full_M4$HsaEX_ID, beta_full = full_M4$beta_ASD,
             SE_full = full_M4$SE, direction_full = full_M4$direction))
summ <- list()
for (i in seq_len(nrow(full_beta))) {
  mid <- full_beta$model_id[i]; eid <- full_beta$event_id[i]
  bf <- full_beta$beta_full[i]; sef <- full_beta$SE_full[i]
  dfull <- full_beta$direction_full[i]
  s <- long[long$model_id == mid & long$event_id == eid &
            is.na(long$error_reason), ]
  if (nrow(s) == 0) {
    summ[[length(summ) + 1]] <- data.frame(
      model_id = mid, event_id = eid, gene = gene_of[[eid]],
      full_beta = bf, full_SE = sef, full_direction = dfull,
      n_deletions_ok = 0, NA_REASON = "ALL_DELETIONS_ERROR",
      stringsAsFactors = FALSE)
    next
  }
  near_zero <- abs(bf) < 1e-6
  if (near_zero) {
    summ[[length(summ) + 1]] <- data.frame(
      model_id = mid, event_id = eid, gene = gene_of[[eid]],
      full_beta = bf, full_SE = sef, full_direction = dfull,
      n_deletions_ok = nrow(s), NA_REASON = "FULL_BETA_NEAR_ZERO",
      stringsAsFactors = FALSE)
    next
  }
  retention <- abs(s$beta) / abs(bf)
  dfb <- (bf - s$beta) / sef
  dirp <- s$direction == dfull
  asd_del <- s$removed_diagnosis == "ASD"
  summ[[length(summ) + 1]] <- data.frame(
    model_id = mid, event_id = eid, gene = gene_of[[eid]],
    full_beta = bf, full_SE = sef, full_direction = dfull,
    n_deletions_ok = nrow(s),
    NA_REASON = "",
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
    bh_fdr_lt_010_fraction = mean(s$KR_BH_FDR < 0.10, na.rm = TRUE),
    asd_deletion_direction_fraction =
      if (any(asd_del)) mean(dirp[asd_del]) else NA_real_,
    ctl_deletion_direction_fraction =
      if (any(!asd_del)) mean(dirp[!asd_del]) else NA_real_,
    asd_deletion_beta_min = if (any(asd_del)) min(s$beta[asd_del]) else NA,
    asd_deletion_beta_max = if (any(asd_del)) max(s$beta[asd_del]) else NA,
    ctl_deletion_beta_min = if (any(!asd_del)) min(s$beta[!asd_del]) else NA,
    ctl_deletion_beta_max = if (any(!asd_del)) max(s$beta[!asd_del]) else NA,
    stringsAsFactors = FALSE)
}
summ <- do.call(rbind, summ)
for (cn in setdiff(colnames(summ), c("model_id", "event_id", "gene",
                                     "NA_REASON"))) {
  if (is.numeric(summ[[cn]]))
    summ[[cn]][is.na(summ[[cn]])] <- NA_real_
}
wtab(summ, "LODO_EVENT_SUMMARY.tsv")
wtab(long, "LODO_EVENT_RESULTS_LONG.tsv")

# ---------------- Tier A criteria (Section 9) -------------------------------
tierA <- summ[summ$event_id %in% tierA_ids & summ$NA_REASON == "", ]
tierA$direction_preserved <-
  tierA$direction_preservation_fraction >= 0.95
tierA$no_reversal <- tierA$direction_reversal_count == 0
tierA$min_retention_ok <- tierA$min_abs_effect_retention >= 0.50
tierA$no_high_influence <- tierA$max_abs_DFBETA < 1.0
tierA$TIER_A_LODO_CONFIRMED <- tierA$direction_preserved &
  tierA$no_reversal & tierA$min_retention_ok & tierA$no_high_influence
wtab(tierA, "TIER_A_LODO_SUMMARY.tsv")

check_lines <- c(
  "# Tier A LODO confirmed phase (Robustness module 1, Section 9)",
  paste0("TIER_A_EVENTS: ", paste(tierA_ids, collapse = ";")),
  "criteria: direction_preservation>=0.95; no direction reversal;",
  "          min_abs_effect_retention>=0.50; no |DFBETA|>=1.0",
  "FDR retention is auxiliary only and does not phase.")
for (i in seq_len(nrow(tierA))) {
  check_lines <- c(check_lines, paste0(
    tierA$model_id[i], " ", tierA$event_id[i], ": ",
    ifelse(tierA$TIER_A_LODO_CONFIRMED[i], "OK", "ERROR"),
    " (dir_pres=", signif(tierA$direction_preservation_fraction[i], 3),
    ", reversals=", tierA$direction_reversal_count[i],
    ", min_ret=", signif(tierA$min_abs_effect_retention[i], 3),
    ", maxDFBETA=", signif(tierA$max_abs_DFBETA[i], 3), ")"))
}
expected_tierA_rows <- length(tierA_ids) * 2
tierA_complete <- nrow(tierA) == expected_tierA_rows
overall <- ifelse(tierA_complete && nrow(tierA) > 0 &&
                    all(tierA$TIER_A_LODO_CONFIRMED), "OK", "ERROR")
check_lines <- c(check_lines,
  paste0("TIER_A_ROWS_FOUND=", nrow(tierA), " expected=",
         expected_tierA_rows, " complete=", tierA_complete))
check_lines <- c(check_lines, paste0("TIER_A_LODO_OVERALL=", overall))
writeLines(check_lines, file.path(outD, "tier_a_lodo_summary.txt"))

# ---------------- report -----------------------------------------------------
rep <- c(
  "# LODO influence report (Robustness module 1)", "",
  paste0("- 80 donors x 19 events x 2 models (M0 primary, M4 preferred) = ",
         nrow(long), " deletion refits."),
  paste0("- Failed fits: ", nrow(errors),
         " (see LODO_MODEL_ERRORS.tsv)."),
  paste0("- Convergence: see LODO_CONVERGENCE_SUMMARY.tsv."),
  paste0("- BH adjustment: within each donor-deletion x model family ",
         "(19 events)."), "",
  "## Tier A robustness", "")
for (i in seq_len(nrow(tierA))) {
  rep <- c(rep, paste0(
    "- ", tierA$event_id[i], " (", tierA$gene[i], ") ", tierA$model_id[i],
    ": ", ifelse(tierA$TIER_A_LODO_CONFIRMED[i], "OK", "ERROR"),
    " | direction preserved in ",
    round(100 * tierA$direction_preservation_fraction[i], 1), "% of ",
    tierA$n_deletions_ok[i], " deletions",
    " | min effect retention ",
    signif(tierA$min_abs_effect_retention[i], 3),
    " | max |DFBETA| ", signif(tierA$max_abs_DFBETA[i], 3),
    " (most influential donor: ", tierA$most_influential_donor[i], ", ",
    tierA$most_influential_donor_dx[i], ")"))
}
rep <- c(rep, "", paste0("Exact deletion refits are authoritative; ",
                         "approximate influence metrics are supplemental."),
         paste0("TIER_A_LODO_OVERALL=", overall))
writeLines(rep, file.path(outD, "LODO_INFLUENCE_REPORT.md"))

say("LODO complete. TIER_A_LODO_OVERALL=", overall)
