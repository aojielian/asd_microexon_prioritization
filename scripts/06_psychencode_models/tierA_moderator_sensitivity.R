#!/usr/bin/env Rscript
# Tier A diagnosis x sex and diagnosis x age sensitivity.
# Moderator sensitivity models for the four Tier A events (final Supplementary Table S12E).
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).

suppressMessages({library(methods); library(lme4); library(pbkrtest)})
set.seed(42)

project <- Sys.getenv("PROJECT_ROOT", unset = ".")
root    <- file.path(project,
                     "40_moderator_and_direction")
out_dir <- file.path(root, "02_tierA_moderator")
dir.create(out_dir, showWarnings = FALSE)
loglines <- character(0)
say <- function(...) {
  m <- paste0(...)
  cat(m, "\n")
  loglines <<- c(loglines, m)
}

sc <- readRDS(file.path(
  project,
  "34_robustness_and_composition",
  "00_admin", "analysis_cache.rds"))
usage_matrix <- sc$usage_matrix
am           <- sc$analysis_meta
m4_cov       <- sc$m4_cov
m4_fix       <- sc$m4_fix          # "+ I(Age^2) + SeqBatch + Ancestry +
                                   #  PMI_mm + QC1 + ... + QC8"
tierA_ids    <- sc$tierA_ids
gene_of      <- setNames(as.character(sc$primary19$gene),
                         as.character(sc$primary19$HsaEX_ID))
am$dx_binary <- as.numeric(am$dx_binary)
stopifnot(all(rownames(m4_cov) == rownames(am)))
stopifnot(length(tierA_ids) == 4)
say("Loaded analysis cache: ", nrow(am), " samples, ",
    length(unique(am$subject)), " donors")
say("Tier A events: ", paste(tierA_ids, collapse = ", "))
say("Reference M4 RHS extension (m4_fix): ", m4_fix)

epsilon <- 1e-4
wtab <- function(df, fn)
  write.table(df, file.path(out_dir, fn), sep = "\t",
              row.names = FALSE, quote = FALSE)
fmt <- function(x) formatC(x, digits = 17, format = "g")

rhs_fix <- paste("dx_binary + region + Sex + Age + RIN", m4_fix,
                 collapse = " ")
say("Base (final M4) RHS: ", rhs_fix)

mods <- data.frame(moderator = c("sex", "age"),
                   term = c("dx_binary:Sex", "dx_binary:Age"),
                   stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# build per-event model frame (complete cases of the final metadata)
# ---------------------------------------------------------------------------
build_mdf <- function(eid) {
  usage_vals <- usage_matrix[, eid]
  valid <- !is.na(usage_vals)
  mdf <- am[valid, ]
  extra <- m4_cov[match(rownames(mdf), rownames(m4_cov)), , drop = FALSE]
  mdf <- cbind(mdf, extra)
  mdf$usage <- usage_vals[valid]
  mdf$usage_logit <- log((mdf$usage + epsilon) /
                         (1 - mdf$usage + epsilon))
  # same-sample filters as final M4; complete-case requirement for the
  # interaction models (final metadata have no missingness, so zero rows
  # are dropped here)
  cc_cols <- c("dx_binary", "region", "Sex", "Age", "RIN", "subject",
               names(m4_cov))
  keep <- complete.cases(mdf[, cc_cols])
  attr(mdf, "n_dropped_cc") <- sum(!keep)
  mdf[keep, ]
}

# ---------------------------------------------------------------------------
# feasibility table (spec section 6)
# ---------------------------------------------------------------------------
donor_tab <- am[!duplicated(am$subject), c("subject", "Sex", "dx_binary",
                                           "Age")]
feas_rows <- list()
fit_store <- list()
res_rows <- list()
note_rows <- list()
lrt_rows <- list()

for (eid in tierA_ids) {
  mdf <- build_mdf(eid)
  say("COMPLETE-CASE filter ", eid, ": dropped ",
      attr(mdf, "n_dropped_cc"), " of 532 samples")
  g <- gene_of[eid]
  dsub <- mdf[!duplicated(mdf$subject), ]
  frow <- data.frame(
    event_id = eid, gene = g,
    n_samples_base_M4 = nrow(mdf),
    n_donors_base_M4 = length(unique(mdf$subject)),
    n_ASD_donors = sum(dsub$dx_binary == 1),
    n_CTL_donors = sum(dsub$dx_binary == 0),
    n_male_donors = sum(dsub$Sex == "M"),
    n_female_donors = sum(dsub$Sex == "F"),
    n_ASD_male_donors = sum(dsub$dx_binary == 1 & dsub$Sex == "M"),
    n_ASD_female_donors = sum(dsub$dx_binary == 1 & dsub$Sex == "F"),
    n_CTL_male_donors = sum(dsub$dx_binary == 0 & dsub$Sex == "M"),
    n_CTL_female_donors = sum(dsub$dx_binary == 0 & dsub$Sex == "F"),
    age_n_nonmissing = sum(!is.na(mdf$Age)),
    age_min = min(mdf$Age), age_q1 = as.numeric(quantile(mdf$Age, 0.25)),
    age_median = as.numeric(quantile(mdf$Age, 0.5)),
    age_q3 = as.numeric(quantile(mdf$Age, 0.75)),
    age_max = max(mdf$Age),
    sex_interaction_full_rank = NA, age_interaction_full_rank = NA,
    sex_model_converged = NA, age_model_converged = NA,
    sex_model_singular = NA, age_model_singular = NA,
    stringsAsFactors = FALSE)

  for (k in seq_len(nrow(mods))) {
    mod <- mods$moderator[k]; term <- mods$term[k]
    rhs_full <- paste(rhs_fix, "+", term)
    X_full <- model.matrix(as.formula(paste0("~ ", rhs_full)),
                           data = mdf)
    full_rank <- qr(X_full)$rank == ncol(X_full)
    f_full <- as.formula(paste("usage_logit ~", rhs_full,
                               "+ (1|subject)"))
    f_red  <- as.formula(paste("usage_logit ~", rhs_fix,
                               "+ (1|subject)"))
    warns <- character(0)
    fit_full <- tryCatch(
      withCallingHandlers(
        lmer(f_full, data = mdf, REML = TRUE),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
      error = function(e) {
        warns <<- c(warns, paste("ERROR:", conditionMessage(e)))
        NULL
      })
    fit_red <- tryCatch(lmer(f_red, data = mdf, REML = TRUE),
                        error = function(e) NULL)
    converged <- (!is.null(fit_full)) && (!is.null(fit_red))
    singular <- if (converged) isSingular(fit_full) else NA
    kr_F <- NA; ndf <- NA; ddf <- NA; p_kr <- NA
    if (converged) {
      kr <- tryCatch(KRmodcomp(fit_full, fit_red),
                     error = function(e) NULL)
      if (!is.null(kr)) {
        kr_F <- kr$test$stat[1]; ndf <- kr$test$ndf[1]
        ddf <- kr$test$ddf[1]; p_kr <- kr$test$p.value[1]
      }
    }
    estimable <- full_rank && converged && !is.na(p_kr)
    reason <- ""
    if (!full_rank) reason <- paste0("interaction design rank deficient: ",
                                     "rank ", qr(X_full)$rank, "/",
                                     ncol(X_full))
    if (!converged && reason == "") reason <- paste(
      "model did not converge:", paste(unique(warns), collapse = " | "))
    if (converged && is.na(p_kr) && reason == "")
      reason <- "KRmodcomp failed"

    # coefficient-level extraction from the full model
    est <- NA; se <- NA; ci_lo <- NA; ci_hi <- NA
    coef_name <- NA
    if (converged) {
      cf <- coef(summary(fit_full))
      nm <- rownames(cf)
      hit <- if (mod == "sex") grep("^dx_binary:Sex", nm) else
             grep("^dx_binary:Age$", nm)
      if (length(hit) == 1) {
        coef_name <- nm[hit]
        est <- cf[hit, "Estimate"]; se <- cf[hit, "Std. Error"]
        if (!is.na(ddf)) {
          ci_lo <- est - qt(0.975, ddf) * se
          ci_hi <- est + qt(0.975, ddf) * se
        }
      } else {
        reason <- paste0(reason, if (nzchar(reason)) "; " else "",
                         "interaction coefficient not found in ",
                         paste(nm, collapse = ","))
        estimable <- FALSE
      }
    }

    res_rows[[length(res_rows) + 1]] <- data.frame(
      gene = g, event_id = eid, moderator = mod,
      interaction_term = ifelse(is.na(coef_name), term, coef_name),
      estimate = fmt(est), SE = fmt(se), df = fmt(ddf),
      CI_low = fmt(ci_lo), CI_high = fmt(ci_hi),
      KR_P = fmt(p_kr), model_n = nrow(mdf),
      donor_n = length(unique(mdf$subject)),
      converged = toupper(as.character(converged)),
      singular = toupper(as.character(singular)),
      estimable = ifelse(estimable, "YES", "NOT_ESTIMABLE"),
      stringsAsFactors = FALSE)
    note_rows[[length(note_rows) + 1]] <- data.frame(
      event_id = eid, gene = g, moderator = mod,
      status = ifelse(estimable, "ESTIMABLE", "NOT_ESTIMABLE"),
      full_rank = full_rank, converged = converged,
      singular_fit = singular, KR_F = fmt(kr_F),
      reason_not_estimable = reason, stringsAsFactors = FALSE)
    if (mod == "sex") {
      frow$sex_interaction_full_rank <- full_rank
      frow$sex_model_converged <- converged
      frow$sex_model_singular <- singular
    } else {
      frow$age_interaction_full_rank <- full_rank
      frow$age_model_converged <- converged
      frow$age_model_singular <- singular
    }
    fit_store[[paste(eid, mod, sep = "|")]] <-
      list(full = fit_full, red = fit_red)
    say("FIT ", g, " ", eid, " moderator=", mod,
        " full_rank=", full_rank, " converged=", converged,
        " singular=", singular, " estimable=", estimable,
        " KR_P=", fmt(p_kr))
  }
  feas_rows[[eid]] <- frow
}

feas <- do.call(rbind, feas_rows); rownames(feas) <- NULL
wtab(feas, "TIERA_MODERATOR_FEASIBILITY.tsv")
say("WROTE TIERA_MODERATOR_FEASIBILITY.tsv (", nrow(feas), " rows)")

res <- do.call(rbind, res_rows); rownames(res) <- NULL
wtab(res, "TIERA_MODERATOR_INTERACTION_RESULTS.tsv")
say("WROTE TIERA_MODERATOR_INTERACTION_RESULTS.tsv (", nrow(res),
    " rows = 4 events x 2 moderators)")

notes <- do.call(rbind, note_rows); rownames(notes) <- NULL
wtab(notes, "TIERA_MODERATOR_MODEL_STATUS_NOTES.tsv")
say("WROTE TIERA_MODERATOR_MODEL_STATUS_NOTES.tsv (", nrow(notes),
    " rows)")

# ---------------------------------------------------------------------------
# secondary sensitivity: nested ML likelihood-ratio test (part of the
# authoritative M4 framework as sensitivity only)
# ---------------------------------------------------------------------------
for (eid in tierA_ids) {
  mdf <- build_mdf(eid)
  for (k in seq_len(nrow(mods))) {
    mod <- mods$moderator[k]; term <- mods$term[k]
    rhs_full <- paste(rhs_fix, "+", term)
    lf <- tryCatch(lmer(as.formula(paste("usage_logit ~", rhs_full,
                                         "+ (1|subject)")),
                        data = mdf, REML = FALSE),
                   error = function(e) NULL)
    lr <- tryCatch(lmer(as.formula(paste("usage_logit ~", rhs_fix,
                                         "+ (1|subject)")),
                        data = mdf, REML = FALSE),
                   error = function(e) NULL)
    chisq <- NA; p_lrt <- NA
    if (!is.null(lf) && !is.null(lr)) {
      an <- anova(lr, lf)
      chisq <- an$Chisq[2]; p_lrt <- an$`Pr(>Chisq)`[2]
    }
    lrt_rows[[length(lrt_rows) + 1]] <- data.frame(
      gene = gene_of[eid], event_id = eid, moderator = mod,
      role = "SECONDARY_SENSITIVITY_ONLY",
      LRT_chisq = fmt(chisq), LRT_P = fmt(p_lrt),
      stringsAsFactors = FALSE)
  }
}
lrt <- do.call(rbind, lrt_rows); rownames(lrt) <- NULL
wtab(lrt, "TIERA_MODERATOR_LRT_SECONDARY.tsv")
say("WROTE TIERA_MODERATOR_LRT_SECONDARY.tsv (secondary sensitivity)")

# ---------------------------------------------------------------------------
# BH across the prespecified family of exactly 8 tests
# ---------------------------------------------------------------------------
est_rows <- res[res$estimable == "YES", ]
pvals <- suppressWarnings(as.numeric(res$KR_P))
pvals[is.na(pvals)] <- 1.0           # non-estimable cannot be significant
stopifnot(length(pvals) == 8)
bh8 <- p.adjust(pvals, method = "BH")
bh_tab <- data.frame(gene = res$gene, event_id = res$event_id,
                     moderator = res$moderator,
                     raw_P = res$KR_P,
                     BH_FDR_8 = fmt(bh8),
                     significant_FDR05 =
                       ifelse(bh8 < 0.05, "YES", "NO"),
                     stringsAsFactors = FALSE)
wtab(bh_tab, "TIERA_MODERATOR_BH8.tsv")
say("WROTE TIERA_MODERATOR_BH8.tsv (family fixed at 8 tests)")
say("MODERATOR_TESTS_FDR05=", sum(bh8 < 0.05), "/8")
say("NOMINAL_RAW_P_LT_0.05=", sum(pvals < 0.05), "/8")
for (i in seq_len(nrow(res))) {
  say("  ", res$gene[i], " ", res$moderator[i],
      " raw=", res$KR_P[i], " BH8=", fmt(bh8[i]))
}
say("PHASE A_models_run=4 events x 2 moderators = 8 fits")

writeLines(loglines, file.path(root, "06_logs",
                               "tierA_moderator_sensitivity.log"))
