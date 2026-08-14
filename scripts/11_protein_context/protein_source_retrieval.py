#!/usr/bin/env python3
"""Probability-scale Phase 5 — external annotation integration (per-protein, lightweight).

Sources (all fetched 2026-08-07, small per-protein queries only):
- UniProt reviewed entries (rest.uniprot.org): Q7Z460 CLASP1, Q5GLZ8 HERC4,
  Q05397 PTK2, P10586 PTPRF — features, InterPro/Pfam cross-references,
  Ensembl isoform mapping.
- MANE release 1.5 (ftp.ncbi.nlm.nih.gov, stream-filtered to the 4 genes).
- AlphaFold DB (alphafold.ebi.ac.uk): per-residue pLDDT from the latest
  monomer model CIF.
- APPRIS: service returned HTTP 404 (unavailable) -> reported as missing,
  not substituted.

No whole proteome, no whole GENCODE download, no BAM/FASTQ.
"""
import json
import os
import re

TMP = os.environ.get("SCRATCH_ROOT", "/tmp")
OUT = os.path.join(os.environ.get("PROJECT_ROOT", "."), "36_probability_scale_and_protein/03_tierA_protein_mapping")

GENES = {
    "CLASP1": {"event": "HsaEX0015476", "acc": "Q7Z460",
               "window": "RSRSANPAGA", "n_insert": 9,
               "d3_incl_tx": "ENST00000397587"},
    "HERC4":  {"event": "HsaEX0029786", "acc": "Q5GLZ8",
               "window": "DVNHGLTE", "n_insert": 8,
               "d3_incl_tx": "ENST00000395198"},
    "PTK2":   {"event": "HsaEX0050855", "acc": "Q05397",
               "window": "DEISGDE", "n_insert": 6,
               "d3_incl_tx": "ENST00000521986"},
    "PTPRF":  {"event": "HsaEX0051138", "acc": "P10586",
               "window": "WRPEESEDY", "n_insert": 9,
               "d3_incl_tx": "ENST00000372407"},
}
# extended exclusion junction contexts (translated from VastDB Seq_C1/Seq_C2 in
# the final CDS frame; computed below from raw sequences)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(seq):
    seq = seq.upper().replace('U', 'T')
    return ''.join(CODON_TABLE.get(seq[i:i + 3], 'X')
                   for i in range(0, len(seq) - len(seq) % 3, 3))


# VastDB flanking sequences (EVENT_INFO, transcript orientation)
FLANKS = {
    "CLASP1": ("CACAGTCCCAGC", "CTGGCAGCCGGT", 1),
    "HERC4":  ("GCCTATGGAATG", "TTGGCAGATATC", 0),
    "PTK2":   ("TCTCTGTGTCAG", "AAACAGATGATT", 1),
    "PTPRF":  ("GCCGAGGCCCAG", "GAAACCACTATC", 0),
}

up = {}
for gene, g in GENES.items():
    up[gene] = json.load(open(f"{TMP}/up_{g['acc']}.json"))


def feat_loc(f):
    s = f["location"]["start"]["value"]
    e = f["location"]["end"]["value"]
    return s, e


# ---- locate insertion site on the UniProt canonical sequence ----
sites = {}
for gene, g in GENES.items():
    seq = up[gene]["sequence"]["value"]
    acc = g["acc"]
    w = g["window"]
    pos = seq.find(w)
    if pos >= 0:
        lo, hi = pos + 1, pos + len(w)
        canonical_state = "canonical_isoform_INCLUDES_microexon"
        site_basis = f"junction window {w} found in {acc} canonical sequence"
    else:
        # canonical = exclusion isoform: build the exclusion-junction peptide
        # in the correct reading frame and locate it on the canonical sequence.
        # frame = number of bases of the junction-spanning codon contributed
        # by the upstream exon (GENCODE phase at the junction).
        c1, c2, frame = FLANKS[gene]
        if frame:
            up_nt = c1[:-frame]
            up_nt = up_nt[len(up_nt) % 3:]                 # align to boundary
            jun_nt = c1[len(c1) - frame:] + c2[:3 - frame]  # spanning codon
            dn_nt = c2[3 - frame:]
            dn_nt = dn_nt[:len(dn_nt) - len(dn_nt) % 3]
        else:
            up_nt = c1[len(c1) % 3:]
            jun_nt = ""
            dn_nt = c2[:len(c2) - len(c2) % 3]
        ctx = translate(up_nt) + translate(jun_nt) + translate(dn_nt)
        n_up = len(up_nt) // 3
        pos = seq.find(ctx)
        if pos >= 0:
            # the junction lies after the upstream in-frame codons, plus the
            # junction-spanning codon if one exists (frame > 0)
            split_res = pos + n_up + (1 if jun_nt else 0)
            lo, hi = split_res, split_res
            canonical_state = "canonical_isoform_EXCLUDES_microexon"
            if jun_nt:
                split_word = (f"insertion splits canonical residue {split_res} "
                              f"(frame-{frame} codon split; exclusion residue "
                              f"{translate(jun_nt)} at this position)")
            else:
                split_word = f"insertion after canonical residue {split_res}"
            site_basis = (f"exclusion junction context {ctx} found in {acc}; "
                          + split_word)
        else:
            lo = hi = None
            canonical_state = "UNRESOLVED"
            site_basis = "neither inclusion window nor exclusion context found"
    sites[gene] = {"lo": lo, "hi": hi, "state": canonical_state, "basis": site_basis}
    print(f"[protein-mapping] {gene}: {site_basis} -> site {lo}-{hi} ({canonical_state})")

# ---- parse AlphaFold CIF per-residue pLDDT ----
plddt = {}
for gene, g in GENES.items():
    acc = g["acc"]
    seq_id_to_res = {}
    res_b = {}
    in_site = False
    cols = []
    with open(f"{TMP}/af_{acc}.cif") as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split()[0].split(".", 1)[1])
                in_site = True
                continue
            if in_site:
                if line.startswith("#") or line.startswith("_"):
                    in_site = False
                    continue
                parts = line.split()
                if len(parts) < len(cols):
                    continue
                d = dict(zip(cols, parts))
                if d.get("label_atom_id") != "CA":
                    continue
                try:
                    resid = int(d["label_seq_id"])
                    b = float(d["B_iso_or_equiv"])
                except (ValueError, KeyError):
                    continue
                res_b[resid] = b
    plddt[gene] = res_b
    print(f"[protein-mapping] {gene} AlphaFold residues with pLDDT: {len(res_b)}")


def plddt_class(v):
    if v >= 90:
        return "very_high"
    if v >= 70:
        return "confident"
    if v >= 50:
        return "low"
    return "very_low"


# ---- build feature proximity table ----
FEAT_TYPES = ["Domain", "Motif", "Region", "Binding site", "Active site", "Site",
              "Modified residue", "Cross-link", "Coiled coil", "Compositional bias",
              "Transmembrane", "Intramembrane", "Topological domain", "DNA binding",
              "Zinc finger", "Alternative sequence"]
feat_rows = []
for gene, g in GENES.items():
    acc = g["acc"]
    ev = g["event"]
    lo, hi = sites[gene]["lo"], sites[gene]["hi"]
    seqlen = up[gene]["sequence"]["length"]
    for f in up[gene].get("features", []):
        if f["type"] not in FEAT_TYPES:
            continue
        s, e = feat_loc(f)
        if lo is None:
            dist, prox = None, "UNRESOLVED"
        elif not (e < lo or s > hi):
            dist, prox = 0, "overlap"
        else:
            dist = max(lo - e, s - hi)
            prox = "proximal_le_10aa" if dist <= 10 else \
                   ("nearby_11_30aa" if dist <= 30 else "distant_gt_30aa")
        evd = ";".join(sorted({(x.get("source") or "curated") for x in f.get("evidences", [])})) or "curated"
        feat_rows.append({
            "gene": gene, "event_id": ev, "uniprot_acc": acc,
            "uniprot_seqlen": seqlen, "feature_type": f["type"],
            "feature_description": f.get("description", ""),
            "feature_start": s, "feature_end": e,
            "insertion_site_start": lo, "insertion_site_end": hi,
            "site_state_on_canonical": sites[gene]["state"],
            "distance_residues": dist, "proximity": prox,
            "evidence_source": evd, "evidence_class": "DIRECT_DATABASE_ANNOTATION"})
with open(os.path.join(OUT, "TIER_A_PROTEIN_FEATURES.tsv"), "w") as fh:
    cols = list(feat_rows[0].keys())
    fh.write("\t".join(cols) + "\n")
    for r in feat_rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"[protein-mapping] wrote {len(feat_rows)} feature rows")

# ---- UniProt / InterPro / Pfam mapping ----
map_rows = []
for gene, g in GENES.items():
    acc = g["acc"]
    lo, hi = sites[gene]["lo"], sites[gene]["hi"]
    xrefs = up[gene].get("uniProtKBCrossReferences", [])
    ip = [x for x in xrefs if x["database"] == "InterPro"]
    pf = [x for x in xrefs if x["database"] == "Pfam"]
    for x in ip:
        name = next((p["value"] for p in x["properties"] if p["key"] == "EntryName"), "")
        map_rows.append({"gene": gene, "event_id": g["event"], "uniprot_acc": acc,
                         "db": "InterPro", "db_id": x["id"], "entry_name": name,
                         "note": "entry-level cross-reference (UniProt)"})
    for x in pf:
        name = next((p["value"] for p in x["properties"] if p["key"] == "EntryName"), "")
        map_rows.append({"gene": gene, "event_id": g["event"], "uniprot_acc": acc,
                         "db": "Pfam", "db_id": x["id"], "entry_name": name,
                         "note": "entry-level cross-reference (UniProt)"})
    # domain features with positions (for proximity context)
    for f in up[gene].get("features", []):
        if f["type"] != "Domain":
            continue
        s, e = feat_loc(f)
        evs = f.get("evidences", [])
        ipid = next((x["id"] for x in evs if x.get("source") == "InterPro"), "")
        prox = "NA"
        if lo is not None:
            if not (e < lo or s > hi):
                prox = "overlap"
            else:
                dist = max(lo - e, s - hi)
                prox = "proximal_le_10aa" if dist <= 10 else \
                       ("nearby_11_30aa" if dist <= 30 else "distant_gt_30aa")
        map_rows.append({"gene": gene, "event_id": g["event"], "uniprot_acc": acc,
                         "db": "Domain_feature", "db_id": ipid or "NA",
                         "entry_name": f"{f.get('description','')} [{s}-{e}] proximity={prox}",
                         "note": "positioned domain feature on canonical sequence"})
with open(os.path.join(OUT, "TIER_A_UNIPROT_INTERPRO_MAPPING.tsv"), "w") as fh:
    cols = ["gene", "event_id", "uniprot_acc", "db", "db_id", "entry_name", "note"]
    fh.write("\t".join(cols) + "\n")
    for r in map_rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"[protein-mapping] wrote {len(map_rows)} InterPro/Pfam mapping rows")

# ---- AlphaFold context table ----
af_rows = []
for gene, g in GENES.items():
    acc = g["acc"]
    lo, hi = sites[gene]["lo"], sites[gene]["hi"]
    pl = plddt[gene]
    api = json.load(open(f"{TMP}/af_api_{acc}.json")) if os.path.exists(f"{TMP}/af_api_{acc}.json") else None
    if lo is not None and pl:
        win = [pl[r] for r in range(lo, hi + 1) if r in pl]
        ext = [pl[r] for r in range(max(1, lo - 5), hi + 6) if r in pl]
        mean_win = sum(win) / len(win) if win else None
        mean_ext = sum(ext) / len(ext) if ext else None
    else:
        win, ext, mean_win, mean_ext = [], [], None, None
    allv = list(pl.values())
    af_rows.append({
        "gene": gene, "event_id": g["event"], "uniprot_acc": acc,
        "model_id": f"AF-{acc}-F1", "model_version": "v6 (latest; fetched 2026-08-07)",
        "method": "AlphaFold Monomer v2.0 pipeline",
        "site_start": lo, "site_end": hi,
        "site_state_on_canonical": sites[gene]["state"],
        "n_residues_modeled": len(allv),
        "mean_plddt_site": round(mean_win, 2) if mean_win is not None else "NA",
        "mean_plddt_site_pm5": round(mean_ext, 2) if mean_ext is not None else "NA",
        "site_plddt_class": plddt_class(mean_win) if mean_win is not None else "NA",
        "interpretation": ("site lies in a high-confidence (likely structured) region"
                           if (mean_win or 0) >= 70 else
                           "site lies in a low-confidence (likely disordered/flexible) region")
                          if mean_win is not None else "UNRESOLVED"})
with open(os.path.join(OUT, "TIER_A_ALPHAFOLD_CONTEXT.tsv"), "w") as fh:
    cols = list(af_rows[0].keys())
    fh.write("\t".join(cols) + "\n")
    for r in af_rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"[protein-mapping] wrote {len(af_rows)} AlphaFold context rows")
for r in af_rows:
    print("  ", r["gene"], "site", r["site_start"], "-", r["site_end"],
          "pLDDT", r["mean_plddt_site"], r["site_plddt_class"])
