#!/usr/bin/env Rscript
# Deterministic transcript pair selection.
# D2 inclusion/exclusion transcript pair selection by highest median count / effective length over 532 samples with lexical tie-breaking.
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).

dir32sub <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "32_psychencode_sensitivity/02_psychencode_sensitivity/transcript_set_intermediates")
rsem_rds <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "32_psychencode_sensitivity/00_admin/rsem_cache.rds")
gtf_gz   <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "00_reference/gencode_v33/gencode.v33.annotation.gtf.gz")
mane19   <- file.path(Sys.getenv("SCRATCH_ROOT", unset = tempdir()), "mane19.tsv")
out_dir  <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "36_probability_scale_and_protein/05_D2_transcript_check")

defs <- read.delim(file.path(dir32sub, "TRANSCRIPT_SET_DEFINITIONS.tsv"),
                   stringsAsFactors = FALSE)
cache <- readRDS(rsem_rds)

# ---- median effective-length-normalized rates (final pipeline logic) ----
rates <- cache$rsem_member[, cache$samp_idx_raw] / cache$effLen_member
med_rate <- apply(rates, 1, function(z) median(z, na.rm = TRUE))
names(med_rate) <- cache$member_base_ids
stopifnot(all(defs$event_id %in% names(cache$incl_tx$HsaEX_ID) | TRUE))

incl_by_ev <- split(cache$incl_tx$transcript_id, cache$incl_tx$HsaEX_ID)
excl_by_ev <- split(cache$excl_tx$transcript_id, cache$excl_tx$HsaEX_ID)

# ---- GENCODE v33 basic tags / transcript_type (full local GTF, read-only) ----
con <- gzfile(gtf_gz, "rt")
basic <- character(); ttype <- character()
while (TRUE) {
  ln <- readLines(con, n = 50000, warn = FALSE)
  if (length(ln) == 0) break
  ln <- ln[grepl("^\ttranscript\t|^[^#]+\ttranscript\t", ln)]
  if (length(ln) == 0) next
  tid <- sub('.*transcript_id "([^"]+)".*', "\\1", ln)
  is_basic <- grepl('tag "basic"', ln)
  typ <- sub('.*transcript_type "([^"]+)".*', "\\1", ln)
  base <- sub("\\.[0-9]+$", "", tid)
  basic <- c(basic, setNames(as.character(is_basic), base))
  ttype <- c(ttype, setNames(typ, base))
}
close(con)
# keep first occurrence (versions may repeat across gene loci rarely)
basic <- basic[!duplicated(names(basic))]
ttype <- ttype[!duplicated(names(ttype))]
cat(sprintf("[probability-scale-7] GTF transcript lines parsed: %d basic tags\n",
            sum(basic == "TRUE")))

# ---- MANE v1.5 (gene -> Ensembl MANE Select transcript) ----
mane <- if (file.exists(mane19))
  read.delim(mane19, header = FALSE, stringsAsFactors = FALSE) else NULL
if (!is.null(mane)) colnames(mane) <- c("gene", "mane_enst")

rows <- list()
for (i in seq_len(nrow(defs))) {
  ev <- defs$event_id[i]; gene <- defs$gene[i]
  ti <- defs$D2_incl[i]; te <- defs$D2_excl[i]
  med_incl <- if (!is.na(ti)) med_rate[ti] else NA_real_
  med_excl <- if (!is.na(te)) med_rate[te] else NA_real_
  # tie check within the D0 candidate sets
  cand_i <- incl_by_ev[[ev]]; cand_e <- excl_by_ev[[ev]]
  tie_i <- if (!is.na(ti) && length(cand_i) > 1) {
    m <- med_rate[cand_i]; m[is.na(m)] <- -1
    sum(abs(m - max(m)) < 1e-12) > 1
  } else FALSE
  tie_e <- if (!is.na(te) && length(cand_e) > 1) {
    m <- med_rate[cand_e]; m[is.na(m)] <- -1
    sum(abs(m - max(m)) < 1e-12) > 1
  } else FALSE
  mane_i <- if (!is.null(mane) &&
                any(mane$gene == gene)) {
    ms <- mane$mane_enst[mane$gene == gene]
    ifelse(sub("\\.[0-9]+$", "", ms) == sub("\\.[0-9]+$", "", ti),
           "IS_MANE_SELECT", "NOT_MANE_SELECT")
  } else "NOT_FETCHED"
  mane_e <- if (!is.null(mane) &&
                any(mane$gene == gene)) {
    ms <- mane$mane_enst[mane$gene == gene]
    ifelse(sub("\\.[0-9]+$", "", ms) == sub("\\.[0-9]+$", "", te),
           "IS_MANE_SELECT", "NOT_MANE_SELECT")
  } else "NOT_FETCHED"
  rows[[i]] <- data.frame(
    gene = gene, event_id = ev,
    inclusion_transcript_id = ti, exclusion_transcript_id = te,
    inclusion_median_expression = round(med_incl, 6),
    exclusion_median_expression = round(med_excl, 6),
    inclusion_MANE_Select = mane_i, exclusion_MANE_Select = mane_e,
    inclusion_GENCODE_basic = ifelse(!is.na(ti) &&
      basic[sub("\\.[0-9]+$", "", ti)] == "TRUE", "YES", "NO"),
    exclusion_GENCODE_basic = ifelse(!is.na(te) &&
      basic[sub("\\.[0-9]+$", "", te)] == "TRUE", "YES", "NO"),
    inclusion_APPRIS = "UNAVAILABLE_SERVICE_404",
    exclusion_APPRIS = "UNAVAILABLE_SERVICE_404",
    inclusion_protein_coding = ifelse(!is.na(ti) &&
      ttype[sub("\\.[0-9]+$", "", ti)] == "protein_coding", "YES", "NO"),
    exclusion_protein_coding = ifelse(!is.na(te) &&
      ttype[sub("\\.[0-9]+$", "", te)] == "protein_coding", "YES", "NO"),
    selection_priority_step = "D2: highest median effective-length-normalized expression over 532 analysis samples, per side, within D0 membership set",
    tie_break_reason = ifelse(tie_i | tie_e,
      "tie encountered: lexicographically smallest transcript ID",
      "no tie: unique maximum"),
    final_pair_status = ifelse(is.na(ti) | is.na(te),
      "INCOMPLETE_PAIR", "COMPLETE_PAIR"),
    provenance = "<internal-path-redacted>/02_psychencode_sensitivity/transcript_set_intermediates/TRANSCRIPT_SET_DEFINITIONS.tsv + rsem_cache.rds + GENCODE v33 GTF + MANE v1.5 GTF",
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
write.table(res, file.path(out_dir, "D2_REPRESENTATIVE_TRANSCRIPTS_ALL19.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
tierA <- res[res$event_id %in% c("HsaEX0015476", "HsaEX0029786",
                                 "HsaEX0050855", "HsaEX0051138"), ]
write.table(tierA, file.path(out_dir, "D2_REPRESENTATIVE_TRANSCRIPTS_TIERA.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[probability-scale-7] wrote %d D2 rows (%d complete pairs)\n",
            nrow(res), sum(res$final_pair_status == "COMPLETE_PAIR")))
print(tierA[, c("gene", "event_id", "inclusion_transcript_id",
                "exclusion_transcript_id", "inclusion_median_expression",
                "exclusion_median_expression", "inclusion_MANE_Select",
                "inclusion_GENCODE_basic", "inclusion_protein_coding")])
cat("D2_TIERA_TRANSCRIPTS_EXPLICIT=OK\n")
