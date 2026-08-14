#!/usr/bin/env Rscript
# Neuron-merged composition sensitivity.
# 6-class to neuron-merged sensitivity with k = 3 PCs and M4C refits.
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).

suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(pbkrtest))
options(stringsAsFactors = FALSE)

dir34   <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "34_robustness_and_composition")
out_dir <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "36_probability_scale_and_protein/04_neuron_merged_composition")

cache  <- readRDS(file.path(dir34, "00_admin", "analysis_cache.rds"))
frac   <- read.delim(file.path(dir34, "04_cell_composition",
                               "COMPOSITION_FRACTIONS_HARMONIZED.tsv"),
                     stringsAsFactors = FALSE)
usage  <- cache$usage_matrix
ameta  <- cache$analysis_meta
m4cov  <- cache$m4_cov
events <- cache$events
tierA  <- cache$tierA_ids
m4_frz <- cache$M4_repro
m4c <- read.delim(file.path(dir34, "05_composition_adjusted_models",
                                   "M4C_EVENT_RESULTS.tsv"),
                         stringsAsFactors = FALSE)
m4c <- m4c[m4c$removed_donor == "NONE_full_sample" &
                         m4c$model_id == "M4C_primary", ]

stopifnot(nrow(frac) == 532, nrow(usage) == 532)
stopifnot(all(frac$sample_id == rownames(usage)))

# ---- build 6-class neuron-merged fractions ----
CLS6 <- c("Neuron", "Astrocyte", "Oligodendrocyte", "OPC",
          "Microglia_immune", "Endothelial_mural")
F6 <- data.frame(sample_id    = frac$sample_id,
                 subject_id   = frac$subject_id,
                 region       = frac$region,
                 diagnosis    = frac$diagnosis,
                 Neuron       = frac$Excitatory_neuron + frac$Inhibitory_neuron,
                 Astrocyte    = frac$Astrocyte,
                 Oligodendrocyte = frac$Oligodendrocyte,
                 OPC          = frac$OPC,
                 Microglia_immune = frac$Microglia_immune,
                 Endothelial_mural = frac$Endothelial_mural)
stopifnot(all(is.finite(as.matrix(F6[, CLS6]))),
          all(as.matrix(F6[, CLS6]) >= 0))
rs <- rowSums(as.matrix(F6[, CLS6]))
# the 7-class TSV is written at limited decimal precision; row sums deviate
# from 1 by <= 3e-6 in the source file, so tolerance 1e-4 is used (same logic
# as the analysis QC: all fractions finite and non-negative, sums ~ 1)
stopifnot(max(abs(rs - 1)) < 1e-4)
cat(sprintf("[probability-scale-6] neuron-merged fractions: %d samples x %d classes; row sums max deviation = %.3e\n",
            nrow(F6), length(CLS6), max(abs(rs - 1))))
write.table(F6, file.path(out_dir, "NEURON_MERGED_FRACTIONS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- final PC pipeline on 6-class matrix ----
F <- as.matrix(F6[, CLS6])
rownames(F) <- F6$sample_id
pos <- F[F > 0]
delta <- 0.5 * min(pos)
Fr <- F; Fr[Fr == 0] <- delta
n_zeros <- sum(F == 0)
gm  <- exp(rowMeans(log(Fr)))
clr <- log(Fr) - log(gm)
pc  <- prcomp(clr, center = TRUE, scale. = FALSE)
ev  <- pc$sdev^2
cve <- ev / sum(ev)
cumv <- cumsum(cve)
k80 <- which(cumv >= 0.80)[1]
k <- min(3, k80); if (k < 2) k <- 2
cat(sprintf("[probability-scale-6] zero replacement delta=%.5g (n zeros=%d); PC variance: %s\n",
            delta, n_zeros, paste(round(cve, 3), collapse = ", ")))
cat(sprintf("[probability-scale-6] k80=%d -> retained k=%d; cumulative variance=%.4f\n",
            k80, k, cumv[k]))

pc_sum <- data.frame(pc = paste0("PC", seq_along(cve)),
                     sdev = pc$sdev,
                     variance_explained = cve,
                     cumulative_variance = cumv,
                     retained = seq_along(cve) <= k)
pc_sum <- cbind(pc_sum, as.data.frame(pc$rotation))
write.table(pc_sum, file.path(out_dir, "NEURON_MERGED_CLR_PC_SUMMARY.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

scores <- pc$x[, seq_len(k), drop = FALSE]
colnames(scores) <- paste0("NMCompPC", seq_len(k))
stopifnot(all(rownames(scores) == rownames(usage)))

# ---- model frame (final M4 covariates + neuron-merged PCs) ----
md0 <- data.frame(
  sample_id = ameta$sample_id,
  subject   = factor(ameta$subject),
  dx_binary = ameta$dx_binary,
  region    = factor(ameta$region),
  Sex       = factor(ameta$Sex),
  Age       = ameta$Age,
  RIN       = ameta$RIN,
  SeqBatch  = factor(m4cov$SeqBatch, levels = levels(m4cov$SeqBatch)),
  Ancestry  = factor(m4cov$Ancestry, levels = levels(m4cov$Ancestry)),
  PMI_mm    = m4cov$PMI_mm,
  QC1 = m4cov$QC1, QC2 = m4cov$QC2, QC3 = m4cov$QC3, QC4 = m4cov$QC4,
  QC5 = m4cov$QC5, QC6 = m4cov$QC6, QC7 = m4cov$QC7, QC8 = m4cov$QC8)
md0 <- cbind(md0, as.data.frame(scores))

pc_terms <- paste(paste0("NMCompPC", seq_len(k)), collapse = " + ")
f_nm <- as.formula(paste("usage_logit ~ dx_binary + region + Sex + Age + RIN +",
                         "I(Age^2) + SeqBatch + Ancestry + PMI_mm +",
                         "QC1 + QC2 + QC3 + QC4 + QC5 + QC6 + QC7 + QC8 +",
                         pc_terms, "+ (1 | subject)"))
eps <- 1e-4
rows <- list()
for (ev in events) {
  u <- usage[, ev]
  mdf <- md0
  mdf$usage_logit <- log((u + eps) / (1 - u + eps))
  fit_full <- tryCatch(lmer(f_nm, data = mdf, REML = TRUE),
                       error = function(e) {
                         cat(sprintf("[probability-scale-6] FIT ERROR %s: %s\n", ev, conditionMessage(e)))
                         NULL })
  if (is.null(fit_full)) {
    rows[[length(rows) + 1]] <- data.frame(
      event_id = ev, gene = NA, n_samples = nrow(mdf), beta = NA, SE = NA,
      KR_P = NA, converged = FALSE)
    next
  }
  fit_red <- lmer(update(f_nm, . ~ . - dx_binary), data = mdf, REML = TRUE)
  kr <- KRmodcomp(fit_full, fit_red)
  cf <- summary(fit_full)$coefficients
  rows[[length(rows) + 1]] <- data.frame(
    event_id = ev,
    gene = cache$primary19$gene[match(ev, cache$primary19$HsaEX_ID)],
    n_samples = nrow(mdf), n_donors = length(unique(mdf$subject)),
    beta = cf["dx_binary", "Estimate"], SE = cf["dx_binary", "Std. Error"],
    KR_P = kr$test$p.value[1], converged = TRUE)
}
res <- do.call(rbind, rows)
res$KR_BH_FDR <- p.adjust(res$KR_P, method = "BH")
res$direction_negative <- res$beta < 0
res$model_id <- "M4C_NeuronMerged"
res$formula <- paste(deparse(f_nm, width.cutoff = 500), collapse = " ")
write.table(res, file.path(out_dir, "M4C_NEURON_MERGED_ALL19.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[probability-scale-6] M4C_NeuronMerged fitted for %d/%d events\n",
            sum(res$converged), length(events)))

# ---- Tier A summary ----
ta <- res[match(tierA, res$event_id), ]
ta <- ta[, c("gene", "event_id", "n_samples", "n_donors", "beta", "SE",
             "KR_P", "KR_BH_FDR", "direction_negative")]
write.table(ta, file.path(out_dir, "M4C_NEURON_MERGED_TIERA_SUMMARY.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- comparison vs M4 and vs primary M4C ----
m4_beta  <- m4_frz$beta_ASD[match(res$event_id, m4_frz$HsaEX_ID)]
m4_dir   <- m4_frz$direction[match(res$event_id, m4_frz$HsaEX_ID)]
m4c_beta <- m4c$beta[match(res$event_id, m4c$event_id)]
m4c_fdr  <- m4c$KR_BH_FDR[match(res$event_id, m4c$event_id)]

pear  <- cor(m4_beta, res$beta, method = "pearson")
spear <- cor(m4_beta, res$beta, method = "spearman")
dir_retained <- sum(sign(res$beta) == sign(m4_beta))
ta_dir_retained <- sum(sign(ta$beta) == sign(m4_frz$beta_ASD[match(tierA, m4_frz$HsaEX_ID)]))
retention <- abs(res$beta) / abs(m4_beta)
median_retention <- median(retention, na.rm = TRUE)

vs <- data.frame(
  event_id = res$event_id, gene = res$gene,
  M4_beta_reference = m4_beta, M4_direction_reference = m4_dir,
  M4C_primary_beta_reference = m4c_beta, M4C_primary_KR_BH_FDR_reference = m4c_fdr,
  M4C_NeuronMerged_beta = res$beta,
  M4C_NeuronMerged_KR_P = res$KR_P,
  M4C_NeuronMerged_KR_BH_FDR = res$KR_BH_FDR,
  direction_matches_M4 = sign(res$beta) == sign(m4_beta),
  direction_matches_M4C_primary = sign(res$beta) == sign(m4c_beta),
  abs_effect_retention_vs_M4 = retention,
  is_TierA = res$event_id %in% tierA)
write.table(vs, file.path(out_dir, "M4C_NEURON_MERGED_VS_PRIMARY_M4C.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ta_sig <- ta$KR_BH_FDR[ta$event_id %in% c("HsaEX0015476", "HsaEX0050855", "HsaEX0051138")] < 0.05
herc4 <- ta[ta$event_id == "HsaEX0029786", ]

cat(sprintf("[probability-scale-6] Pearson M4 vs M4C_NeuronMerged beta = %.4f\n", pear))
cat(sprintf("[probability-scale-6] Spearman M4 vs M4C_NeuronMerged beta = %.4f\n", spear))
cat(sprintf("[probability-scale-6] direction retention vs M4 = %d/19\n", dir_retained))
cat(sprintf("[probability-scale-6] Tier A direction retention = %d/4\n", ta_dir_retained))
cat(sprintf("[probability-scale-6] median absolute effect retention vs M4 = %.4f\n", median_retention))
cat(sprintf("[probability-scale-6] CLASP1/PTK2/PTPRF FDR<0.05 = %d/3\n", sum(ta_sig)))
cat(sprintf("[probability-scale-6] HERC4 beta=%.4f KR_P=%.4f FDR=%.4f direction_negative=%s\n",
            herc4$beta, herc4$KR_P, herc4$KR_BH_FDR, herc4$direction_negative))
cat("NEURON_MERGED_COMPOSITION_COMPLETED=OK\n")
cat(sprintf("NEURON_MERGED_ALL19_MODEL_COMPLETED=%s\n",
            ifelse(sum(res$converged) == 19, "OK", "ERROR")))
cat(sprintf("NEURON_MERGED_TIERA_DIRECTION_STATUS_REPORTED=%s\n",
            ifelse(ta_dir_retained == 4, "OK_4_OF_4", "REPORTED")))
saveRDS(list(res = res, vs = vs, ta = ta, k = k, cve = cve, pear = pear,
             spear = spear, dir_retained = dir_retained,
             ta_dir_retained = ta_dir_retained,
             median_retention = median_retention),
        file.path(out_dir, "NEURON_MERGED_RDATA.rds"))
cat("[probability-scale-6] DONE\n")
