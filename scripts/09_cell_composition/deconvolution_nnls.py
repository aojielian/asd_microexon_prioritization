import os
#!/usr/bin/env python3
"""Robustness module: reference-based cell-composition deconvolution of the
532 PsychENCODE samples.

Method (documented, amendment Section 6 hierarchy):
  The local brainSCOPE reference is log2-scale pseudo-bulk (NOT raw counts),
  so Bisque's raw-counts assumption (method A) is not met. Primary method is
  non-negative least squares (NNLS, method C) on:
    - reference: brainSCOPE pseudo-bulk linearized (2^x - 1), averaged across
      donors per subtype, then equal-weight averaged across subtypes per broad
      class, then each class profile scaled to sum 1e6 (TPM-like);
    - signature genes: top class-specific genes selected from the reference
      (documented, pre-specified rule: per-class mean / (per-class mean +
      mean of all other classes), top N_TOP_PER_CLASS per class, union);
    - bulk: PsychENCODE gene-level TPM (RSEM expected counts / effLen -> TPM,
      summed per GENCODE v33 gene; exported binary from PSYCHENCODE_GENE_TPM.rds).
  Fractions = NNLS coefficients renormalized to sum 1.
  A SENSITIVITY reference excluding the 6 overlapping retained donors is also
  built and compared (amendment Section 5).

Outputs (04_cell_composition/):
  COMPOSITION_FRACTIONS_RAW.tsv
  COMPOSITION_FRACTIONS_HARMONIZED.tsv
  COMPOSITION_DECONVOLUTION_QC.tsv
  COMPOSITION_FRACTIONS_SENSITIVITY_OVERLAP_EXCLUDED.tsv
  COMPOSITION_PRIMARY_VS_SENSITIVITY.tsv
  COMPOSITION_SIGNATURE_GENES.tsv
"""
import gzip, os
import numpy as np
from scipy.optimize import nnls

BS = os.path.join(os.environ.get("REFERENCE_ROOT", "."), "step02_standardize/brainSCOPE")
TASK = os.path.join(os.environ.get("PROJECT_ROOT", "."), "34_robustness_and_composition")
OUT = os.path.join(TASK, "04_cell_composition")
EXP = os.path.join(os.environ.get("SCRATCH_ROOT", "/tmp"), "gene_tpm_export")
N_TOP_PER_CLASS = 150

MAPPING = {
    "Excitatory_neuron": ["L2.3.IT", "L4.IT", "L5.6.NP", "L5.ET", "L5.IT",
                          "L6.CT", "L6.IT", "L6.IT.Car3", "L6b"],
    "Inhibitory_neuron": ["Chandelier__Pvalb", "Lamp5", "Lamp5.Lhx6", "Pax6",
                          "Sncg", "Sst__Sst.Chodl", "Vip"],
    "Astrocyte": ["Astro"],
    "Oligodendrocyte": ["Oligo"],
    "OPC": ["OPC"],
    "Microglia_immune": ["Micro.PVM", "Immune"],
    "Endothelial_mural": ["Endo__VLMC", "PC", "SMC"],
}
CLASSES = list(MAPPING.keys())
# The overlap-retained donor list (sensitivity reference) is cohort
# metadata supplied at runtime as a plain-text file (one identifier
# per line) under TASK; the list itself is not shipped in this
# public repository.
_overlap_file = os.path.join(TASK, "gse_overlap_donor_exclusion.txt")
if os.path.exists(_overlap_file):
    OVERLAP_RETAINED = set(filter(None, map(str.strip,
                                        open(_overlap_file).read().splitlines())))
else:
    print("WARNING: gse_overlap_donor_exclusion.txt not found under TASK; "
          "no overlap-retained donor exclusion applied", flush=True)
    OVERLAP_RETAINED = set()

# ------------------------------------------------ 1. gene_name -> ENSG map
print("building gene_name -> ENSG map ...", flush=True)
name2ensg = {}
ambig = set()
t2g_path = os.path.join(TASK, "01_source_inventory", "derived",
                        "GENCODE_V33_TX2GENE_FROM_GTF.tsv")
with open(t2g_path) as fh:
    next(fh)
    for line in fh:
        tx, chrom, gid, gname = line.rstrip("\n").split("\t")
        if gname in ambig:
            continue
        if gname in name2ensg and name2ensg[gname] != gid:
            ambig.add(gname)
            del name2ensg[gname]
        else:
            name2ensg[gname] = gid
print(f"  unique gene_name->ENSG: {len(name2ensg)}, ambiguous dropped: "
      f"{len(ambig)}", flush=True)

# ------------------------------------------------ 2. subtype profiles (both refs in one ok)
import pandas as pd

def read_ct_profiles(ct):
    base = os.path.join(BS, ct)
    df = pd.read_csv(os.path.join(base, f"{ct}.expression.tsv.gz"),
                     sep="\t", index_col=0)
    cols = df.columns.tolist()
    sens_idx = [i for i, c in enumerate(cols) if c not in OVERLAP_RETAINED]
    vals = df.values.astype(np.float64)
    lin = np.power(2.0, vals) - 1.0                 # linearize log2 scale
    mean_all = lin.mean(axis=1)
    mean_sens = lin[:, sens_idx].mean(axis=1) if sens_idx else None
    gene_names = df.index.tolist()
    # aggregate by ENSG (mean over duplicated ENSG rows)
    acc_all, acc_sens = {}, {}
    for gi, gname in enumerate(gene_names):
        ensg = name2ensg.get(gname)
        if ensg is None:
            continue
        a = acc_all.get(ensg)
        acc_all[ensg] = (a[0] + mean_all[gi], a[1] + 1) if a else (mean_all[gi], 1)
        if mean_sens is not None:
            s = acc_sens.get(ensg)
            acc_sens[ensg] = ((s[0] + mean_sens[gi], s[1] + 1) if s
                              else (mean_sens[gi], 1))
    prof_all = {k: v[0] / v[1] for k, v in acc_all.items()}
    prof_sens = {k: v[0] / v[1] for k, v in acc_sens.items()}
    return prof_all, prof_sens, len(cols), len(sens_idx)

print("building subtype profiles (primary + sensitivity) ...", flush=True)
sub_primary, sub_sens = {}, {}
for bc in CLASSES:
    for ct in MAPPING[bc]:
        pa, ps, n_all, n_sens = read_ct_profiles(ct)
        sub_primary[ct] = pa
        sub_sens[ct] = ps
        print(f"  {ct}: {len(pa)} genes; donors primary={n_all} "
              f"sensitivity={n_sens}", flush=True)

# ------------------------------------------------ 3. broad class profiles (equal subtype weight)
def build_class_profile(subtype_dict):
    allg = set()
    for ct in subtype_dict:
        allg |= set(subtype_dict[ct])
    prof = {}
    for bc in CLASSES:
        vals = {}
        for g in allg:
            xs = [subtype_dict[ct][g] for ct in MAPPING[bc]
                  if g in subtype_dict[ct]]
            if xs:
                vals[g] = sum(xs) / len(xs)
        tot = sum(vals.values())
        if tot > 0:
            vals = {g: v / tot * 1e6 for g, v in vals.items()}
        prof[bc] = vals
    return prof

print("building broad-class profiles ...", flush=True)
cls_primary = build_class_profile(sub_primary)
cls_sens = build_class_profile(sub_sens)

# ------------------------------------------------ 4. signature gene selection
def select_signature(class_prof):
    genes = set()
    for bc in CLASSES:
        genes |= set(class_prof[bc])
    genes = sorted(genes)
    sig = set()
    for bc in CLASSES:
        in_c = class_prof[bc]
        scores = []
        for g in genes:
            x = in_c.get(g, 0.0)
            others = [class_prof[o].get(g, 0.0) for o in CLASSES if o != bc]
            denom = x + (sum(others) / len(others))
            scores.append((x / denom if denom > 0 else 0.0, g))
        scores.sort(reverse=True)
        for sc, g in scores[:N_TOP_PER_CLASS]:
            if class_prof[bc].get(g, 0.0) > 0:
                sig.add(g)
    return sorted(sig)

print("selecting signature genes ...", flush=True)
sig_genes = select_signature(cls_primary)
print(f"  signature genes (union): {len(sig_genes)}", flush=True)
with open(os.path.join(OUT, "COMPOSITION_SIGNATURE_GENES.tsv"), "w") as fh:
    fh.write("gene_id\tin_classes\n")
    for g in sig_genes:
        inc = [bc for bc in CLASSES if cls_primary[bc].get(g, 0.0) > 0]
        fh.write(f"{g}\t{';'.join(inc)}\n")

# ------------------------------------------------ 5. load bulk gene TPM
print("loading bulk gene TPM ...", flush=True)
genes = [l.strip() for l in open(os.path.join(EXP, "genes.txt"))]
samples = [l.strip() for l in open(os.path.join(EXP, "samples.txt"))]
mat = np.fromfile(os.path.join(EXP, "matrix.bin"), dtype=np.float64)
bulk = mat.reshape((len(genes), len(samples)), order="F")
gene_idx = {g: i for i, g in enumerate(genes)}

# restrict signature genes to those present in bulk
sig_present = [g for g in sig_genes if g in gene_idx]
print(f"  signature genes present in bulk: {len(sig_present)}/{len(sig_genes)}",
      flush=True)
rows = [gene_idx[g] for g in sig_present]

# ------------------------------------------------ 6. reference matrix (sig genes x classes)
def make_ref_matrix(class_prof, sig_present):
    R = np.zeros((len(sig_present), len(CLASSES)))
    for j, bc in enumerate(CLASSES):
        for i, g in enumerate(sig_present):
            R[i, j] = class_prof[bc].get(g, 0.0)
    return R

R_primary = make_ref_matrix(cls_primary, sig_present)
R_sens = make_ref_matrix(cls_sens, sig_present)

# ------------------------------------------------ 7. NNLS per sample
def deconvolve(bulk_sub, R):
    ns = bulk_sub.shape[1]
    coefs = np.zeros((ns, R.shape[1]))
    r2 = np.zeros(ns)
    resnorm = np.zeros(ns)
    for s in range(ns):
        y = bulk_sub[:, s]
        x, r = nnls(R, y)
        coefs[s] = x
        resnorm[s] = r
        ss_tot = np.sum((y - y.mean()) ** 2)
        rss = r ** 2
        r2[s] = 1 - rss / ss_tot if ss_tot > 0 else np.nan
    return coefs, r2, resnorm

bulk_sub = bulk[rows, :]
print("running NNLS (primary) ...", flush=True)
C_p, r2_p, res_p = deconvolve(bulk_sub, R_primary)
print("running NNLS (sensitivity) ...", flush=True)
C_s, r2_s, res_s = deconvolve(bulk_sub, R_sens)

# normalize to fractions
def normalize(C):
    s = C.sum(axis=1, keepdims=True)
    s[s == 0] = np.nan
    return C / s

F_p = normalize(C_p)
F_s = normalize(C_s)

# ------------------------------------------------ 8. sample metadata
# pec_532_sample_meta.txt: sample_id, subject, region, dx_binary
sample_meta = {}
donor_of_sample = {}
with open(os.path.join(os.environ.get("SCRATCH_ROOT", "/tmp"), "pec_532_sample_meta.txt")) as fh:
    for line in fh:
        sid, subj, region, dx = line.rstrip("\n").split("\t")
        sample_meta[sid] = (region, dx)
        donor_of_sample[sid] = subj

# ------------------------------------------------ 9. write outputs
def write_fractions(path, F, label):
    with open(path, "w") as fh:
        fh.write("sample_id\tsubject_id\tregion\tdiagnosis\t" +
                 "\t".join(CLASSES) + "\n")
        for si, sid in enumerate(samples):
            region, dx = sample_meta.get(sid, ("", ""))
            subj = donor_of_sample.get(sid, "")
            vals = "\t".join(f"{F[si, j]:.6f}" if not np.isnan(F[si, j])
                             else "NA" for j in range(len(CLASSES)))
            fh.write(f"{sid}\t{subj}\t{region}\t{dx}\t{vals}\n")

write_fractions(os.path.join(OUT, "COMPOSITION_FRACTIONS_HARMONIZED.tsv"),
                F_p, "primary")
# RAW = unnormalized NNLS coefficients
with open(os.path.join(OUT, "COMPOSITION_FRACTIONS_RAW.tsv"), "w") as fh:
    fh.write("sample_id\t" + "\t".join(f"{c}_coef" for c in CLASSES) +
             "\tcoef_sum\n")
    for si, sid in enumerate(samples):
        vals = "\t".join(f"{C_p[si, j]:.6g}" for j in range(len(CLASSES)))
        fh.write(f"{sid}\t{vals}\t{C_p[si].sum():.6g}\n")

write_fractions(os.path.join(OUT,
                "COMPOSITION_FRACTIONS_SENSITIVITY_OVERLAP_EXCLUDED.tsv"),
                F_s, "sensitivity")

# QC
with open(os.path.join(OUT, "COMPOSITION_DECONVOLUTION_QC.tsv"), "w") as fh:
    fh.write("sample_id\tsubject_id\tregion\tdiagnosis\tfraction_sum_primary\t"
             "r_squared_primary\tnnls_residual_primary\tnegative_coef_count\t"
             "missing_fraction_count\tmapping_confidence\tfit_quality\n")
    for si, sid in enumerate(samples):
        region, dx = sample_meta.get(sid, ("", ""))
        subj = donor_of_sample.get(sid, "")
        n_neg = int((C_p[si] < 0).sum())
        n_na = int(np.isnan(F_p[si]).sum())
        mc = len(sig_present) / len(sig_genes)
        fq = "GOOD" if (not np.isnan(r2_p[si]) and r2_p[si] > 0.5) else "MODERATE" \
             if (not np.isnan(r2_p[si]) and r2_p[si] > 0.2) else "POOR"
        fh.write(f"{sid}\t{subj}\t{region}\t{dx}\t{C_p[si].sum():.6g}\t"
                 f"{r2_p[si]:.4f}\t{res_p[si]:.6g}\t{n_neg}\t{n_na}\t"
                 f"{mc:.4f}\t{fq}\n")

# primary vs sensitivity comparison
with open(os.path.join(OUT, "COMPOSITION_PRIMARY_VS_SENSITIVITY.tsv"), "w") as fh:
    fh.write("class\tspearman_primary_vs_sensitivity\tmean_abs_difference\n")
    from scipy.stats import spearmanr
    for j, bc in enumerate(CLASSES):
        a, b = F_p[:, j], F_s[:, j]
        m = ~(np.isnan(a) | np.isnan(b))
        rho = spearmanr(a[m], b[m]).correlation if m.sum() > 2 else np.nan
        mad = np.nanmean(np.abs(a - b))
        fh.write(f"{bc}\t{rho:.4f}\t{mad:.6f}\n")

print("DECONVOLUTION DONE", flush=True)
print(f"  median R^2 primary: {np.nanmedian(r2_p):.4f}", flush=True)
print(f"  mean coef sum primary: {np.nanmean(C_p.sum(axis=1)):.4g}", flush=True)
