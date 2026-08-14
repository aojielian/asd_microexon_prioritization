#!/usr/bin/env Rscript
# Marker validation and composition PCs.
# Marker-gene validation and CLR + PCA of composition fractions (k = 2, 85.4% cumulative variance).
# Paths: configured via environment variables PROJECT_ROOT, DATA_ROOT, REFERENCE_ROOT, LIFTOVER_PATH (see config/paths_template.yaml).
task <- file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "34_robustness_and_composition")
outD <- file.path(task, "04_cell_composition")
logf <- file.path(task, "99_logs", "marker_pc.log")
say <- function(...) { msg <- paste0(...); cat(msg, "\n");
                       cat(msg, "\n", file = logf, append = TRUE) }
file.create(logf); cat("", file = logf)

CLASSES <- c("Excitatory_neuron","Inhibitory_neuron","Astrocyte",
             "Oligodendrocyte","OPC","Microglia_immune","Endothelial_mural")

# ---------------- pre-specified canonical marker sets -----------------------
MARKERS <- list(
  Excitatory_neuron = c("SLC17A7","CAMK2A","NEUROD6","SATB2","TBR1","RORB"),
  Inhibitory_neuron = c("GAD1","GAD2","SLC32A1","DLX1","DLX2"),
  Astrocyte         = c("GFAP","AQP4","SLC1A3","ALDH1L1","GJA1","ALDOC"),
  Oligodendrocyte   = c("MBP","PLP1","MOG","MAG","CNP"),
  OPC               = c("PDGFRA","CSPG4","OLIG1","VCAN"),
  Microglia_immune  = c("AIF1","C1QA","C1QB","C1QC","TMEM119","P2RY12","CX3CR1"),
  Endothelial_mural = c("CLDN5","FLT1","PECAM1","ACTA2","PDGFRB","NOTCH3"))

# ---------------- load gene TPM + fractions --------------------------------
g <- readRDS(file.path(outD, "PSYCHENCODE_GENE_TPM.rds"))
gene_tpm <- g$gene_tpm; gene_symbol <- g$gene_symbol
say("gene TPM: ", nrow(gene_tpm), " x ", ncol(gene_tpm))

frac <- read.delim(file.path(outD, "COMPOSITION_FRACTIONS_HARMONIZED.tsv"),
                   stringsAsFactors = FALSE)
samples <- frac$sample_id
keep <- match(samples, colnames(gene_tpm))
stopifnot(all(is.na(keep) == FALSE))
expr <- gene_tpm[, keep, drop = FALSE]          # genes x samples (aligned)
colnames(expr) <- samples
say("aligned expression to ", ncol(expr), " samples")

# ---------------- marker scores --------------------------------------------
logtpm <- log2(expr + 1)
# z-score each gene across samples
gz <- t(scale(t(as.matrix(logtpm))))
score <- matrix(NA_real_, nrow = ncol(expr), ncol = length(CLASSES),
                dimnames = list(samples, CLASSES))
for (bc in CLASSES) {
  mg <- MARKERS[[bc]]
  idx <- which(gene_symbol %in% mg)
  used <- gene_symbol[idx]
  if (length(idx) == 0) {
    say("  ", bc, ": no marker genes present in bulk; score NA")
    next
  }
  score[, bc] <- colMeans(gz[idx, , drop = FALSE], na.rm = TRUE)
  say("  ", bc, ": markers used ", paste(used, collapse = ","))
}

# ---------------- Spearman validation --------------------------------------
val <- data.frame(marker_class = character(), fraction_class = character(),
                  rho = numeric(), P = numeric(), n = integer(),
                  match_type = character(), stringsAsFactors = FALSE)
for (mc in CLASSES) {
  for (fc in CLASSES) {
    ok <- is.finite(score[, mc]) & is.finite(frac[[fc]])
    x <- score[ok, mc]; y <- frac[ok, fc]
    if (length(unique(x)) > 1 && length(unique(y)) > 1) {
      st <- cor.test(x, y, method = "spearman", exact = FALSE)
      rho <- unname(st$estimate); p <- st$p.value
    } else { rho <- NA_real_; p <- NA_real_ }
    val <- rbind(val, data.frame(marker_class = mc, fraction_class = fc,
                                 rho = rho, P = p, n = sum(ok),
                                 match_type = ifelse(mc == fc, "MATCHING",
                                                     "OFF_TARGET"),
                                 stringsAsFactors = FALSE))
  }
}
write.table(val, file.path(outD, "COMPOSITION_MARKER_VALIDATION.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ok/error per class
class_ok <- c()
for (mc in CLASSES) {
  m <- val[val$marker_class == mc & val$match_type == "MATCHING", "rho"]
  ot <- val[val$marker_class == mc & val$match_type == "OFF_TARGET", "rho"]
  class_ok[mc] <- (!is.na(m)) && m >= 0.3 && (!is.na(m)) &&
                    m > max(ot, na.rm = TRUE)
}
n_ok <- sum(class_ok)
marker_ok <- n_ok >= 6
say("marker validation: ", n_ok, "/7 classes ok; overall=",
    ifelse(marker_ok, "OK", "ERROR"))

# ---------------- coverage / region / diagnosis (descriptive) --------------
# fractions file already carries sample_id / subject_id / region / diagnosis
fracm <- frac
region_tab <- aggregate(fracm[, CLASSES], by = list(region = fracm$region),
                        FUN = function(z) round(mean(z, na.rm = TRUE), 4))
dx_tab <- aggregate(fracm[, CLASSES], by = list(dx = fracm$diagnosis),
                    FUN = function(z) round(mean(z, na.rm = TRUE), 4))
write.table(region_tab, file.path(outD, "COMPOSITION_REGION_MEANS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dx_tab, file.path(outD, "COMPOSITION_DX_MEANS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------- PC lock -------------------------------------------------
F <- as.matrix(frac[, CLASSES])
rownames(F) <- frac$sample_id
# documented zero replacement: zeros -> 0.5 * min positive fraction
pos <- F[F > 0]
delta <- 0.5 * min(pos)
Fr <- F; Fr[Fr == 0] <- delta
say("zero replacement delta=", signif(delta, 5),
    " ; n zeros replaced=", sum(F == 0))
# CLR transform (per sample across classes)
gm <- exp(rowMeans(log(Fr)))
clr <- log(Fr) - log(gm)
# PCA
pc <- prcomp(clr, center = TRUE, scale. = FALSE)
ev <- pc$sdev^2
cve <- ev / sum(ev)
cumv <- cumsum(cve)
k80 <- which(cumv >= 0.80)[1]
k <- min(3, k80)
if (k < 2) k <- 2
say("PC variance explained: ", paste(round(cve[1:6], 3), collapse = ", "))
say("k for >=80%=", k80, " -> retained k=", k, " (cap 3, min 2)")

scores <- pc$x[, seq_len(k), drop = FALSE]
colnames(scores) <- paste0("CompPC", seq_len(k))
loadings <- as.data.frame(pc$rotation[, seq_len(k), drop = FALSE])
loadings$class <- CLASSES
write.table(data.frame(sample_id = rownames(scores), scores,
                       check.names = FALSE),
            file.path(outD, "COMPOSITION_PC_SCORES.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(loadings, file.path(outD, "COMPOSITION_PC_LOADINGS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# decision doc
cell_status <- ifelse(marker_ok, "FEASIBLE_LOCAL_REFERENCE",
                      "NOT_FEASIBLE_FROM_LOCAL_REFERENCE")
dl <- c(
  "# Composition PC decision (Robustness module)", "",
  paste0("- Marker validation: ", n_ok, "/7 classes ok (criterion: ",
         "matching-class Spearman rho >= 0.3 AND greater than all ",
         "off-target rhos). Overall: ", ifelse(marker_ok, "OK", "ERROR")),
  paste0("- CELL_COMPOSITION_STATUS=", cell_status), "",
  "- Zero replacement: zeros -> 0.5 x min positive fraction (documented).",
  "- CLR transform per sample across the 7 broad classes.",
  "- PCA on CLR fractions.",
  paste0("- Variance explained per PC: ",
         paste(round(cve[1:min(6, length(cve))], 3), collapse = ", ")),
  paste0("- PCs retained for M4C: ", k,
         " (minimum PCs explaining >=80% variance, capped at 3, minimum 2)."),
  paste0("- Cumulative variance of retained PCs: ",
         round(cumv[k], 3)), "",
  "This decision is PRIMARY before any event-level M4C modeling.")
writeLines(dl, file.path(outD, "COMPOSITION_PC_DECISION.md"))

# QC report
qr <- c(
  "# Composition QC report (Robustness module)", "",
  paste0("- Samples deconvolved: ", nrow(frac)),
  paste0("- Donors: ", length(unique(frac$subject_id))), "",
  "## Marker validation (Spearman, matching vs off-target)", "")
for (mc in CLASSES) {
  m <- val[val$marker_class == mc & val$match_type == "MATCHING", ]
  qr <- c(qr, paste0("- ", mc, ": matching rho=", signif(m$rho, 3),
                     " (n=", m$n, ") ",
                     ifelse(class_ok[mc], "OK", "ERROR")))
}
qr <- c(qr, "", paste0("OVERALL_MARKER_VALIDATION=",
                       ifelse(marker_ok, "OK", "ERROR")),
        paste0("CELL_COMPOSITION_STATUS=", cell_status))
writeLines(qr, file.path(outD, "COMPOSITION_QC_REPORT.md"))
say("DONE. CELL_COMPOSITION_STATUS=", cell_status)
