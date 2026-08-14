#!/usr/bin/env python3
"""Final publication Main Tables 1 and 2 (4-row Tier A source data).

Generates the publication-exact source-data files for manuscript
Main Table 1 "Primary Tier A microexon events and model stability" and
Main Table 2 "Tier A robustness, probability-scale effects and protein
context".  All numeric values are read from the prespecified derived analysis
files below (no statistic is recomputed); the protein-context strings are
the curated manuscript annotations recorded in the Tier A protein-context
analysis (UniProt features and AlphaFold pLDDT context).

Prespecified derived inputs (read-only):
  - master table                  GRCh38 coordinates + discovery delta PSI
  - adjusted usage M0/M4/M4C      betas, KR FDRs, probability-scale
                                  adjusted differences and 95% CIs
  - Tier A LODO summary (M4)      minimum effect retention, max |DFBETA|
  - neuron-merged Tier A summary  KR BH FDR
  - Tier A protein context        curated context strings

Outputs (into 41_submission_figures_and_tables/03_final_tables/):
  Table_1_TierA_model_stability.tsv
  Table_2_TierA_robustness_probability_protein.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supplementary_common import (ROOT, DIR34, D36_ADJ, D36_NM, MASTER, rd)

OUT = os.path.join(ROOT, "41_submission_figures_and_tables",
                   "03_final_tables")
os.makedirs(OUT, exist_ok=True)

LODO_SUMMARY = os.path.join(DIR34, "03_lodo", "TIER_A_LODO_SUMMARY.tsv")
ADJ_USAGE = os.path.join(D36_ADJ,
                         "TIER_A_ADJUSTED_TRANSCRIPT_USAGE_M0_M4_M4C.tsv")
NEURON_TIER_A = os.path.join(D36_NM, "M4C_NEURON_MERGED_TIERA_SUMMARY.tsv")

# Publication order of the four Tier A events (manuscript order).
TIER_A = [
    ("CLASP1", "HsaEX0015476"),
    ("HERC4", "HsaEX0029786"),
    ("PTK2", "HsaEX0050855"),
    ("PTPRF", "HsaEX0051138"),
]

# Curated protein-context annotations (manuscript Table 2 wording),
# derived from the Tier A protein-features and AlphaFold context analyses.
PROTEIN_CONTEXT = {
    "HsaEX0015476":
        "9-aa insertion, residues 673-682; disordered MAPRE1/3 "
        "interaction region",
    "HsaEX0029786":
        "8-aa segment, residues 643-650; annotated alternative sequence",
    "HsaEX0050855":
        "6-aa insertion before E393; 4 aa from Y397",
    "HsaEX0051138":
        "9-aa insertion, residues 772-780; fibronectin type-III repeat 5",
}

F4 = "{:.4f}"      # beta / FDR / retention display precision
F2 = "{:.2f}"      # percentage-point display precision


def load_master():
    rows = rd(MASTER)
    out = {}
    for r in rows:
        eid = r.get("HsaEX_ID")
        if eid in {e for _, e in TIER_A}:
            out[eid] = r
    return out


def fmt_coord(r):
    start = int(r["start_hg38"])
    end = int(r["end_hg38"])
    return "{chr}:{s:,}-{e:,}".format(chr=r["chr_hg38"], s=start, e=end)


def load_adjusted():
    """(event_id, model) -> dict of probability-scale values."""
    out = {}
    for r in rd(ADJ_USAGE):
        out[(r["event_id"], r["model"])] = {
            "beta": float(r["beta_logit"]),
            "fdr": float(r["KR_BH_FDR"]),
            "diff_pp": float(r["adjusted_difference_percentage_points"]),
            "ci_low": float(r["adjusted_difference_95CI_low"]),
            "ci_high": float(r["adjusted_difference_95CI_high"]),
        }
    return out


def load_lodo():
    """event_id -> (min retention, max |DFBETA|) under the M4 model."""
    out = {}
    for r in rd(LODO_SUMMARY):
        if r["model_id"].startswith("M4"):
            out[r["event_id"]] = (
                float(r["min_abs_effect_retention"]),
                float(r["max_abs_DFBETA"]))
    return out


def load_neuron():
    return {r["event_id"]: float(r["KR_BH_FDR"]) for r in rd(NEURON_TIER_A)}


def write(path, header, rows):
    with open(path, "w", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    print("wrote", path, "(%d rows)" % len(rows))


def main():
    master = load_master()
    adj = load_adjusted()
    lodo = load_lodo()
    neuron = load_neuron()

    t1 = []
    t2 = []
    for gene, eid in TIER_A:
        assert eid in master and eid in lodo and eid in neuron, eid
        m = master[eid]
        beta = {}
        fdr = {}
        diff = {}
        for mod in ("M0", "M4", "M4C"):
            v = adj[(eid, mod)]
            beta[mod] = F4.format(v["beta"])
            fdr[mod] = F4.format(v["fdr"])
            diff[mod] = v
        t1.append([
            gene, eid, fmt_coord(m),
            F4.format(float(m["Parikshak_delta_PSI"])),
            beta["M0"], fdr["M0"],
            beta["M4"], fdr["M4"],
            beta["M4C"], fdr["M4C"],
        ])
        ret, dfb = lodo[eid]
        t2.append([
            gene, eid,
            F4.format(ret), F4.format(dfb),
            F2.format(diff["M4C"]["diff_pp"]),
            F2.format(diff["M4C"]["ci_low"] * 100.0),
            F2.format(diff["M4C"]["ci_high"] * 100.0),
            F4.format(neuron[eid]),
            PROTEIN_CONTEXT[eid],
        ])

    write(os.path.join(OUT, "Table_1_TierA_model_stability.tsv"),
          ["gene", "event_id", "grch38_coordinates", "discovery_delta_PSI",
           "M0_beta", "M0_KR_FDR", "M4_beta", "M4_KR_FDR",
           "M4C_beta", "M4C_KR_FDR"], t1)
    write(os.path.join(OUT, "Table_2_TierA_robustness_probability_protein.tsv"),
          ["gene", "event_id",
           "M4_LODO_min_retention", "M4_LODO_max_abs_DFBETA",
           "M4C_adjusted_difference_percentage_points",
           "M4C_adjusted_difference_95CI_low_pp",
           "M4C_adjusted_difference_95CI_high_pp",
           "neuron_merged_KR_BH_FDR",
           "protein_context"], t2)

    # publication-exact value assertions (round-trip verification)
    assert t1[0][3:] == ["-0.1169", "-0.3008", "0.0112", "-0.2404",
                         "0.0145", "-0.1317", "0.0290"], t1[0]
    assert t1[1][3:] == ["-0.0039", "-0.2869", "0.0112", "-0.2297",
                         "0.0145", "-0.1215", "0.1573"], t1[1]
    assert t1[2][3:] == ["-0.0795", "-0.1707", "0.0112", "-0.1396",
                         "0.0088", "-0.0952", "0.0222"], t1[2]
    assert t1[3][3:] == ["-0.0165", "-0.2904", "0.0112", "-0.2521",
                         "0.0088", "-0.1633", "0.0290"], t1[3]
    assert t2[0][2:8] == ["0.8659", "0.4149", "-1.62", "-2.73", "-0.63",
                          "0.0335"], t2[0]
    assert t2[1][2:8] == ["0.8483", "0.4703", "-1.94", "-3.73", "-0.01",
                          "0.0586"], t2[1]
    assert t2[2][2:8] == ["0.8056", "0.6769", "-2.19", "-3.40", "-1.00",
                          "0.0148"], t2[2]
    assert t2[3][2:8] == ["0.7549", "0.8587", "-1.06", "-1.74", "-0.39",
                          "0.0026"], t2[3]
    print("FINAL_MAIN_TABLES_1_2=OK rows=4/4")
    print("OUT=" + OUT)


if __name__ == "__main__":
    main()
