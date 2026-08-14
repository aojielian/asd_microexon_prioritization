#!/usr/bin/env Rscript
# Probability-scale Tier A effects.
# Marginal standardized fixed-effect predictions of adjusted usage for the four Tier A events with parametric simulation confidence intervals (seed 42).
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).

suppressPackageStartupMessages({
  library(lme4)
  library(pbkrtest)
  library(MASS)
})
set.seed(42)
options(stringsAsFactors = FALSE)

cache_path <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "34_robustness_and_composition/00_admin/analysis_cache.rds")
pc_path    <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "34_robustness_and_composition/04_cell_composition/COMPOSITION_PC_SCORES.tsv")
out_dir    <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "36_probability_scale_and_protein/02_adjusted_usage_effects")

eps <- 1e-4
B   <- 1000  # parametric simulation draws

cache <- readRDS(cache_path)
usage  <- cache$usage_matrix          # 532 x 19 untransformed inclusion fraction
ameta  <- cache$analysis_meta         # 532 x 42
m4cov  <- cache$m4_cov                # SeqBatch Ancestry PMI_mm QC1..QC8
events <- cache$events
tierA  <- cache$tierA_ids
pc     <- read.delim(pc_path, stringsAsFactors = FALSE)

stopifnot(nrow(usage) == 532, length(events) == 19, length(tierA) == 4)
stopifnot(all(rownames(usage) == ameta$sample_id),
          all(rownames(usage) == rownames(m4cov)),
          all(rownames(usage) == pc$sample_id))

# ---- build shared covariate frame ----
md0 <- data.frame(
  sample_id  = ameta$sample_id,
  subject    = factor(ameta$subject),
  dx_binary  = ameta$dx_binary,
  region     = factor(ameta$region),
  Sex        = factor(ameta$Sex),
  Age        = ameta$Age,
  RIN        = ameta$RIN,
  SeqBatch   = factor(m4cov$SeqBatch, levels = levels(m4cov$SeqBatch)),
  Ancestry   = factor(m4cov$Ancestry, levels = levels(m4cov$Ancestry)),
  PMI_mm     = m4cov$PMI_mm,
  QC1 = m4cov$QC1, QC2 = m4cov$QC2, QC3 = m4cov$QC3, QC4 = m4cov$QC4,
  QC5 = m4cov$QC5, QC6 = m4cov$QC6, QC7 = m4cov$QC7, QC8 = m4cov$QC8,
  CompPC1 = pc$CompPC1, CompPC2 = pc$CompPC2
)

f_m0  <- usage_logit ~ dx_binary + region + Sex + Age + RIN + (1 | subject)
f_m4  <- usage_logit ~ dx_binary + region + Sex + Age + RIN + I(Age^2) + SeqBatch +
                       Ancestry + PMI_mm + QC1 + QC2 + QC3 + QC4 + QC5 + QC6 + QC7 + QC8 +
                       (1 | subject)
f_m4c <- update(f_m4, . ~ . + CompPC1 + CompPC2)
flist <- list(M0 = f_m0, M4 = f_m4, M4C = f_m4c)

# ---- fit all 19 events x 3 models; KR P ----
res_all <- list()
fits    <- list()
for (mid in names(flist)) {
  for (ev in events) {
    u <- usage[, ev]
    mdf <- md0
    mdf$usage_logit <- log((u + eps) / (1 - u + eps))
    ff_nb <- nobars(update(flist[[mid]], usage_logit ~ .))
    cc <- complete.cases(model.frame(ff_nb, mdf, na.action = na.ok))
    mdf <- mdf[cc, ]
    key <- paste(mid, ev, sep = "|")
    # plain lmer call, exactly as model_reproduction.R fit_event()
    fit_full <- tryCatch(lmer(flist[[mid]], data = mdf, REML = TRUE),
                         error = function(e) {
                           cat(sprintf("[probability-scale-4] FIT ERROR %s %s: %s\n", mid, ev, conditionMessage(e)))
                           NULL
                         })
    if (is.null(fit_full)) {
      res_all[[key]] <- data.frame(model = mid, event_id = ev, gene = NA, n_samples = nrow(mdf),
                                   beta = NA, SE = NA, KR_P = NA, converged = FALSE)
      next
    }
    fit_red <- lmer(update(flist[[mid]], . ~ . - dx_binary), data = mdf, REML = TRUE)
    kr <- KRmodcomp(fit_full, fit_red)
    cf <- summary(fit_full)$coefficients
    res_all[[key]] <- data.frame(
      model = mid, event_id = ev,
      gene = cache$primary19$gene[match(ev, cache$primary19$HsaEX_ID)],
      n_samples = nrow(mdf), n_donors = length(unique(mdf$subject)),
      beta = cf["dx_binary", "Estimate"], SE = cf["dx_binary", "Std. Error"],
      KR_P = kr$test$p.value[1],
      converged = TRUE, stringsAsFactors = FALSE)
    fits[[key]] <- list(fit = fit_full, data = mdf)
  }
  cat(sprintf("[probability-scale-4] %s fits done (%d events)\n", mid, sum(sapply(res_all, function(r) r$model[1] == mid))))
}
res_df <- do.call(rbind, res_all)
res_df$KR_BH_FDR <- ave(res_df$KR_P, res_df$model, FUN = function(p) p.adjust(p, method = "BH"))
write.table(res_df, file = file.path(out_dir, "REFIT_ALL19_M0_M4_M4C.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[probability-scale-4] wrote REFIT_ALL19_M0_M4_M4C.tsv\n")

# ---- verify final Tier A values ----
m0 <- cache$M0_repro
m4 <- cache$M4_repro
m4c_file <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "34_robustness_and_composition/05_composition_adjusted_models/M4C_EVENT_RESULTS.tsv")
m4c <- read.delim(m4c_file, stringsAsFactors = FALSE)
m4c <- m4c[m4c$removed_donor == "NONE_full_sample" & m4c$model_id == "M4C_primary", ]

check_rows <- list()
for (ev in tierA) {
  r0 <- res_df[res_df$model == "M0" & res_df$event_id == ev, ]
  f0 <- m0[m0$HsaEX_ID == ev, ]
  r4 <- res_df[res_df$model == "M4" & res_df$event_id == ev, ]
  f4 <- m4[m4$HsaEX_ID == ev, ]
  rc <- res_df[res_df$model == "M4C" & res_df$event_id == ev, ]
  fc <- m4c[m4c$event_id == ev, ]
  check_rows[[length(check_rows) + 1]] <- data.frame(
    event_id = ev,
    M0_beta_refit = r0$beta, M0_beta_reference = f0$beta_ASD, M0_beta_abs_diff = abs(r0$beta - f0$beta_ASD),
    M0_KR_P_refit = r0$KR_P, M0_KR_P_reference = f0$P_Kenward_Roger,
    M4_beta_refit = r4$beta, M4_beta_reference = f4$beta_ASD, M4_beta_abs_diff = abs(r4$beta - f4$beta_ASD),
    M4_KR_P_refit = r4$KR_P, M4_KR_P_reference = f4$P_Kenward_Roger,
    M4C_beta_refit = rc$beta, M4C_beta_reference = fc$beta, M4C_beta_abs_diff = abs(rc$beta - fc$beta),
    M4C_KR_P_refit = rc$KR_P, M4C_KR_P_reference = fc$KR_P, stringsAsFactors = FALSE)
}
check_df <- do.call(rbind, check_rows)
write.table(check_df, file = file.path(out_dir, "REFIT_VS_PRIMARY_TIERA.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
max_beta_diff <- max(c(check_df$M0_beta_abs_diff, check_df$M4_beta_abs_diff, check_df$M4C_beta_abs_diff))
cat(sprintf("[probability-scale-4] max |refit - final| beta over Tier A x 3 models = %.3e\n", max_beta_diff))
if (max_beta_diff > 1e-6) {
  cat("ADJUSTED_USAGE_DIRECTION_MATCHES_BETA=ERROR (refit mismatch, stop)\n")
  quit(status = 2)
}

# ---- marginal adjusted predictions + parametric-simulation CI (Tier A only) ----
plogis_safe <- plogis
out_rows <- list()
for (mid in names(flist)) {
  for (ev in tierA) {
    key <- paste(mid, ev, sep = "|")
    fit <- fits[[key]]$fit
    mdf <- fits[[key]]$data
    n_samples <- nrow(mdf)
    n_donors  <- length(unique(mdf$subject))

    # prediction frames: keep covariates, force dx
    nd0 <- mdf; nd0$dx_binary <- 0
    nd1 <- mdf; nd1$dx_binary <- 1

    # fixed-effect design matrices (drop the random-effect term from the formula)
    ff <- nobars(formula(fit))
    X0 <- model.matrix(ff, nd0)[, names(fixef(fit)), drop = FALSE]
    X1 <- model.matrix(ff, nd1)[, names(fixef(fit)), drop = FALSE]
    beta_hat <- fixef(fit)
    b <- beta_hat["dx_binary"]
    se <- sqrt(vcov(fit)["dx_binary", "dx_binary"])
    kr_p <- res_df$KR_P[res_df$model == mid & res_df$event_id == ev]
    fdr  <- res_df$KR_BH_FDR[res_df$model == mid & res_df$event_id == ev]

    m_ctl <- mean(plogis_safe(X0 %*% beta_hat))
    m_asd <- mean(plogis_safe(X1 %*% beta_hat))
    diff  <- m_asd - m_ctl

    # parametric simulation from fixed-effect vcov (eigen-based, PD-safe)
    V  <- as.matrix(vcov(fit))
    mu <- beta_hat
    eig <- eigen(V, symmetric = TRUE)
    half <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(mu))
    Z <- matrix(rnorm(B * length(mu)), nrow = B)
    draws <- sweep(Z %*% t(half), 2, mu, "+")
    h0 <- plogis_safe(X0 %*% t(draws))   # n_samples x B
    h1 <- plogis_safe(X1 %*% t(draws))
    d_draws <- colMeans(h1) - colMeans(h0)
    ci <- quantile(d_draws, probs = c(0.025, 0.975), names = FALSE)

    direction_check <- ifelse(sign(diff) == sign(b), "OK", "ERROR")

    out_rows[[length(out_rows) + 1]] <- data.frame(
      gene = cache$primary19$gene[match(ev, cache$primary19$HsaEX_ID)],
      event_id = ev, model = mid,
      n_samples = n_samples, n_donors = n_donors,
      beta_logit = b,
      beta_95CI_low = b - 1.96 * se, beta_95CI_high = b + 1.96 * se,
      KR_P = kr_p, KR_BH_FDR = fdr,
      adjusted_control_usage = m_ctl, adjusted_ASD_usage = m_asd,
      adjusted_ASD_minus_control = diff,
      adjusted_difference_percentage_points = diff * 100,
      adjusted_difference_95CI_low = ci[1], adjusted_difference_95CI_high = ci[2],
      prediction_method = "marginal_standardized_fixed_effect_prediction_over_reference_532_sample_covariate_distribution",
      ci_method = sprintf("parametric_simulation_fixed_effect_vcov_B%d_seed42", B),
      direction_check = direction_check,
      source_model_file = "sensitivity refit from analysis_cache.rds (formula per M0_M4_FORMULAS.txt)",
      source_data_file = cache_path,
      stringsAsFactors = FALSE)
    cat(sprintf("[probability-scale-4] %s %s: beta=%.4f adjCtl=%.4f adjASD=%.4f diff_pp=%.3f [%s]\n",
                mid, ev, b, m_ctl, m_asd, diff * 100, direction_check))
  }
}
out_df <- do.call(rbind, out_rows)
write.table(out_df, file = file.path(out_dir, "TIER_A_ADJUSTED_TRANSCRIPT_USAGE_M0_M4_M4C.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[probability-scale-4] wrote TIER_A_ADJUSTED_TRANSCRIPT_USAGE_M0_M4_M4C.tsv\n")
cat(sprintf("ADJUSTED_USAGE_DIRECTION_MATCHES_BETA=%s\n",
            ifelse(all(out_df$direction_check == "OK"), "OK", "ERROR")))
