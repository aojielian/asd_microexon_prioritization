#!/usr/bin/env Rscript
# Reference gene TPM preprocessing.
# Prepares reference cell-type gene TPM data for deconvolution.
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).
task <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "34_robustness_and_composition")
proj <- Sys.getenv("PROJECT_ROOT", unset = ".")
outD <- file.path(task, "04_cell_composition")
dir.create(outD, showWarnings = FALSE)
logf <- file.path(task, "99_logs", "gene_tpm.log")
say <- function(...) { msg <- paste0(...); cat(msg, "\n");
                       cat(msg, "\n", file = logf, append = TRUE) }
file.create(logf); cat("", file = logf)

cache <- readRDS(file.path(task, "00_admin", "analysis_cache.rds"))
am <- cache$analysis_meta
samp_idx_raw <- cache$samp_idx_raw
say("analysis samples: ", nrow(am))

say("loading raw RData ...")
env <- new.env()
load(file.path(proj, "psychencode_processed", "01_02_B_01_RawData.RData"),
     envir = env)
rsem_tx <- env$rsem_tx
effLen  <- env$rsem_transcript_effLen
say("rsem_tx: ", nrow(rsem_tx), " x ", ncol(rsem_tx))
stopifnot(length(effLen) == nrow(rsem_tx))

# ---- restrict to analysis samples immediately (memory) --------------------
rsem_tx <- as.matrix(rsem_tx[, samp_idx_raw, drop = FALSE])
colnames(rsem_tx) <- rownames(am)
say("subset to analysis samples: ", nrow(rsem_tx), " x ", ncol(rsem_tx))

# ---- transcript -> gene map -----------------------------------------------
t2g <- read.delim(file.path(task, "01_source_inventory", "derived",
                            "GENCODE_V33_TX2GENE_FROM_GTF.tsv"),
                  stringsAsFactors = FALSE)
base_ids <- sub("\\.[0-9]+_[0-9]+$", "", rownames(rsem_tx))
gene_of_tx <- t2g$gene_id[match(base_ids, t2g$transcript_id)]
gene_name_of_tx <- t2g$gene_name[match(base_ids, t2g$transcript_id)]
n_mapped <- sum(is.na(gene_of_tx) == FALSE)
say("transcripts mapped to a GENCODE v33 gene: ", n_mapped, " / ",
    nrow(rsem_tx))

# ---- TPM -------------------------------------------------------------------
say("computing TPM ...")
n_bad_el <- sum(effLen <= 0)
say("transcripts with non-positive effective length (masked to rate 0): ",
    n_bad_el)
safe_el <- effLen
safe_el[safe_el <= 0] <- 1                   # avoid div-by-zero
rate <- rsem_tx / safe_el                    # recycle effLen per row
rate[effLen <= 0, ] <- 0                     # zero-effLen contribute nothing
norm <- colSums(rate)                        # per-sample rate sum
stopifnot(all(norm > 0), all(is.finite(norm)))
tpm <- t(t(rate) / norm) * 1e6
rm(rate, rsem_tx); gc()
say("TPM done; col-sum check (should be ~1e6): ",
    signif(mean(colSums(tpm)), 7),
    " ; any NaN: ", any(is.nan(tpm)))

# ---- aggregate to gene -----------------------------------------------------
say("aggregating transcripts -> gene ...")
keep <- is.na(gene_of_tx) == FALSE
gene_ids <- gene_of_tx[keep]
gene_tpm <- rowsum(tpm[keep, , drop = FALSE], group = gene_ids,
                   reorder = TRUE)
# gene symbol lookup
gid2sym <- t2g$gene_name[match(rownames(gene_tpm), t2g$gene_id)]
say("gene-level matrix: ", nrow(gene_tpm), " genes x ", ncol(gene_tpm),
    " samples")

saveRDS(list(gene_tpm = gene_tpm, gene_symbol = gid2sym,
             sample_id = colnames(gene_tpm),
             n_transcripts_mapped = n_mapped,
             n_transcripts_zero_efflen_masked = n_bad_el,
             method = "RSEM expected counts / effLen -> TPM, summed per gene; zero-effLen transcripts masked to rate 0",
             annotation = "GENCODE_V33_TX2GENE_FROM_GTF.tsv"),
        file.path(outD, "PSYCHENCODE_GENE_TPM.rds"))
say("saved PSYCHENCODE_GENE_TPM.rds")
say("GENE_TPM_MATRIX Genes=", nrow(gene_tpm), " Samples=", ncol(gene_tpm))
