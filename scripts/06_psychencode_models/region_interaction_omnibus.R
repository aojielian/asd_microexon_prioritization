#!/usr/bin/env Rscript
# Diagnosis x region interaction omnibus.
# Omnibus tests of diagnosis x region interaction for the four Tier A events (final Supplementary Table S12C).
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).

suppressMessages({library(methods); library(lme4); library(pbkrtest)})
set.seed(42)

project <- Sys.getenv("PROJECT_ROOT", unset = ".")
root    <- file.path(project, "39_rule_and_numeric_verification")
out_dir <- file.path(root, "05_D_region_interaction")
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
m4_fix       <- sc$m4_fix
events       <- sc$events
tierA_ids    <- sc$tierA_ids
mm17         <- sc$mm
am$dx_binary <- as.numeric(am$dx_binary)
stopifnot(all(rownames(m4_cov) == rownames(am)))
stopifnot(all(as.character(mm17$Subject) == am$subject))

wtab <- function(df, fn)
  write.table(df, file.path(out_dir, fn), sep = "\t",
              row.names = FALSE, quote = FALSE)
fmt <- function(x) formatC(x, digits = 15, format = "g")

# ---------------------------------------------------------------------------
# 1. Donor-by-region counts (exact 11 region groupings)
# ---------------------------------------------------------------------------
regions <- levels(factor(am$region))
cnt <- data.frame(region = regions, stringsAsFactors = FALSE)
cnt$n_ASD_donors <- vapply(regions, function(r)
  length(unique(am$subject[am$region == r & am$dx_binary == 1])), integer(1))
cnt$n_CTL_donors <- vapply(regions, function(r)
  length(unique(am$subject[am$region == r & am$dx_binary == 0])), integer(1))
cnt$n_total_donors <- vapply(regions, function(r)
  length(unique(am$subject[am$region == r])), integer(1))
cnt$n_ASD_samples <- vapply(regions, function(r)
  sum(am$region == r & am$dx_binary == 1), integer(1))
cnt$n_CTL_samples <- vapply(regions, function(r)
  sum(am$region == r & am$dx_binary == 0), integer(1))
cnt$n_total_samples <- vapply(regions, function(r)
  sum(am$region == r), integer(1))
rownames(cnt) <- NULL
wtab(cnt, "region_by_diagnosis_donor_counts.tsv")
say("WROTE region_by_diagnosis_donor_counts.tsv (",
    nrow(cnt), " regions)")

mat <- rbind(ASD_donors = cnt$n_ASD_donors, CTL_donors = cnt$n_CTL_donors)
colnames(mat) <- cnt$region
mat <- cbind(diagnosis = rownames(mat), as.data.frame(mat,
              check.names = FALSE))
wtab(mat, "donor_count_matrix_2x11.tsv")
say("WROTE donor_count_matrix_2x11.tsv (2 x ", length(regions), ")")

# ---------------------------------------------------------------------------
# 2. Feasibility criteria C1-C4
# ---------------------------------------------------------------------------
feas <- data.frame(criterion = character(0), status = character(0),
                   evidence = character(0), stringsAsFactors = FALSE)
add_feas <- function(cr, st, ev) {
  feas <<- rbind(feas, data.frame(criterion = cr, status = st,
                                  evidence = ev, stringsAsFactors = FALSE))
  say("FEASIBILITY ", cr, ": ", st, " — ", ev)
}

# C1: >= 5 ASD and >= 5 CTL donors per region
c1_ok <- all(cnt$n_ASD_donors >= 5 & cnt$n_CTL_donors >= 5)
add_feas("C1_min_donors_per_region",
         ifelse(c1_ok, "OK", "ERROR"),
         paste0("min ASD donors/region = ", min(cnt$n_ASD_donors),
                "; min CTL donors/region = ", min(cnt$n_CTL_donors),
                " (threshold >= 5 each)"))

# C3: metadata issues (NA covariates, donor-diagnosis consistency)
cov_names <- c("SeqBatch", "Ancestry", "PMI_mm",
               paste0("QC", 1:8))
na_base <- sum(is.na(am$dx_binary) + is.na(am$region) + is.na(am$Sex) +
               is.na(am$Age) + is.na(am$RIN) + is.na(am$subject)) > 0
na_cov  <- any(vapply(m4_cov[, cov_names],
                      function(x) any(is.na(x)), logical(1)))
diag_per_donor <- tapply(am$dx_binary, am$subject,
                         function(x) length(unique(x)))
c3_ok <- (!na_base) && (!na_cov) && all(diag_per_donor == 1)
add_feas("C3_no_metadata_issues",
         ifelse(c3_ok, "OK", "ERROR"),
         paste0("NA in base covariates = ", na_base,
                "; NA in M4 model-matrix covariates = ", na_cov,
                "; donors with >1 diagnosis value = ",
                sum(diag_per_donor > 1)))

# C4: no region merging — exact 11 regions, names equal the prespecified list
elig <- read.delim(file.path(
  project,
  "32_psychencode_sensitivity",
  "02_psychencode_sensitivity", "REGION_DONOR_ELIGIBILITY.tsv"),
  stringsAsFactors = FALSE)
c4_ok <- length(regions) == 11 &&
  setequal(regions, as.character(elig$region))
add_feas("C4_no_region_merging",
         ifelse(c4_ok, "OK", "ERROR"),
         paste0("n_regions = ", length(regions),
                "; identical to the prespecified region list = ",
                setequal(regions, as.character(elig$region))))

# C2: full-rank M4 + dx_binary:region design (checked per Tier A event below;
# placeholder set after fits). Position recorded here, evidence filled later.

feasibility_ok <- c1_ok && c3_ok && c4_ok  # C2 added after design check

# ---------------------------------------------------------------------------
# 3. Tier A omnibus interaction fits (only if C1/C3/C4 already ok)
# ---------------------------------------------------------------------------
epsilon <- 1e-4
gene_of <- setNames(as.character(sc$primary19$gene),
                    as.character(sc$primary19$HsaEX_ID))

fit_rows <- list()
c2_ok <- TRUE
c2_evidence <- character(0)
if (feasibility_ok) {
  for (eid in tierA_ids) {
    usage_vals <- usage_matrix[, eid]
    valid <- !is.na(usage_vals)
    mdf <- am[valid, ]
    extra <- m4_cov[match(rownames(mdf), rownames(m4_cov)), , drop = FALSE]
    mdf <- cbind(mdf, extra)
    mdf$usage <- usage_vals[valid]
    mdf$usage_logit <- log((mdf$usage + epsilon) /
                           (1 - mdf$usage + epsilon))

    rhs_fix <- paste("dx_binary + region + Sex + Age + RIN", m4_fix,
                     collapse = " ")
    rhs_full <- paste(rhs_fix, "+ dx_binary:region")
    X_full <- model.matrix(as.formula(paste0("~ ", rhs_full)), data = mdf)
    fr <- qr(X_full)$rank == ncol(X_full)
    c2_ok <- c2_ok && fr
    c2_evidence <- c(c2_evidence,
                     paste0(eid, ":rank=", qr(X_full)$rank, "/",
                            ncol(X_full)))

    f_full <- as.formula(paste("usage_logit ~", rhs_full, "+ (1|subject)"))
    f_red  <- as.formula(paste("usage_logit ~", rhs_fix,  "+ (1|subject)"))
    fit_full <- tryCatch(lmer(f_full, data = mdf, REML = TRUE),
                         error = function(e) NULL)
    fit_red  <- tryCatch(lmer(f_red,  data = mdf, REML = TRUE),
                         error = function(e) NULL)
    converged <- (!is.null(fit_full)) && (!is.null(fit_red))
    kr_F <- NA; ndf <- NA; ddf <- NA; p_kr <- NA
    if (converged) {
      kr <- tryCatch(KRmodcomp(fit_full, fit_red),
                     error = function(e) NULL)
      if (!is.null(kr)) {
        kr_F <- kr$test$stat[1]; ndf <- kr$test$ndf[1]
        ddf <- kr$test$ddf[1]; p_kr <- kr$test$p.value[1]
      }
    }
    fit_rows[[eid]] <- data.frame(
      gene = gene_of[eid], event_id = eid,
      n_samples = nrow(mdf), n_donors = length(unique(mdf$subject)),
      interaction_test_method =
        "KRmodcomp_REML(M4+dx_binary:region vs M4), pbkrtest",
      interaction_df = ndf, KR_F = kr_F, KR_denom_df = ddf,
      interaction_p_raw = p_kr,
      model_converged = converged, full_rank = fr,
      singular_fit = if (converged) isSingular(fit_full) else NA,
      stringsAsFactors = FALSE)
    say("FIT ", eid, " converged=", converged, " full_rank=", fr,
        " df=", ndf, " P_KR=", fmt(p_kr))
  }
} else {
  say("FEASIBILITY C1/C3/C4 ERROR — Tier A interaction fits NOT RUN")
}

feasibility_ok <- feasibility_ok && c2_ok
add_feas("C2_full_rank_interaction_design",
         ifelse(c2_ok, "OK", "ERROR"),
         paste(c2_evidence, collapse = "; "))
feas <- rbind(feas, data.frame(
  criterion = "INTERACTION_FEASIBILITY",
  status = ifelse(feasibility_ok, "OK", "NOT_RUN"),
  evidence = paste0("C1=", c1_ok, " C2=", c2_ok, " C3=", c3_ok,
                    " C4=", c4_ok),
  stringsAsFactors = FALSE))
wtab(feas, "region_interaction_feasibility.tsv")
say("WROTE region_interaction_feasibility.tsv (",
    nrow(feas), " rows)")
say("PHASE D_interaction_feasibility=",
    ifelse(feasibility_ok, "OK", "NOT_RUN"))

# ---------------------------------------------------------------------------
# 4. Omnibus table + BH across the 4 Tier A events + prior evidence
# ---------------------------------------------------------------------------
if (feasibility_ok) {
  omn <- do.call(rbind, fit_rows)
  rownames(omn) <- NULL
  omn$interaction_p_bh4 <- p.adjust(omn$interaction_p_raw, method = "BH")
  omn$interpretation <- ifelse(
    omn$interaction_p_bh4 < 0.05,
    "significant diagnosis x region interaction after BH across 4 Tier A events",
    "no significant diagnosis x region interaction after BH across 4 Tier A events")

  prior <- read.delim(file.path(
    project,
    "32_psychencode_sensitivity",
    "02_psychencode_sensitivity", "REGION_INTERACTION_RESULTS.tsv"),
    stringsAsFactors = FALSE)
  pc <- colnames(prior)
  id_col <- intersect(c("HsaEX_ID", "event_id"), pc)[1]
  p_col  <- intersect(c("p_lrt", "P_LRT", "p", "P_interaction"), pc)[1]
  ch_col <- intersect(c("chisq", "LRT_chisq"), pc)[1]
  bh_col <- intersect(c("BH_FDR_interaction"), pc)[1]
  pm <- match(omn$event_id, as.character(prior[[id_col]]))
  omn$prior_dir32_ML_LRT_chisq <- if (!is.na(ch_col)) prior[[ch_col]][pm] else NA
  omn$prior_dir32_ML_LRT_p     <- if (!is.na(p_col))  prior[[p_col]][pm]  else NA
  omn$prior_dir32_ML_LRT_bh    <- if (!is.na(bh_col)) prior[[bh_col]][pm] else NA
  omn$prior_evidence_source <-
    "32_.../02_psychencode_sensitivity/REGION_INTERACTION_RESULTS.tsv"
  wtab(omn, "tier_a_diagnosis_region_omnibus.tsv")
  say("WROTE tier_a_diagnosis_region_omnibus.tsv (",
      nrow(omn), " rows)")
  say("PHASE D_region_interaction_run=YES")
  say("OMNIBUS BH4 significant: ",
      sum(omn$interaction_p_bh4 < 0.05), "/4")
  for (i in seq_len(nrow(omn))) {
    say("  ", omn$gene[i], " ", omn$event_id[i],
        " P_KR=", fmt(omn$interaction_p_raw[i]),
        " BH4=", fmt(omn$interaction_p_bh4[i]))
  }
} else {
  say("PHASE D_region_interaction_run=NO (feasibility ",
      ifelse(feasibility_ok, "OK", "NOT_RUN"), ")")
}

writeLines(loglines, file.path(root, "08_logs",
                               "region_interaction_omnibus.log"))
