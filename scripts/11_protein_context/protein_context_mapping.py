#!/usr/bin/env python3
"""Probability-scale Phase 5 — Tier A coding/protein consequence mapping (local part).

Maps microexon -> transcript -> CDS consequence -> peptide consequence ->
protein residue coordinates for the 4 Tier A events, using ONLY local
authoritative sources:

- final 19-event master table for hg38 event coordinates;
- final transcript-membership table for D0 membership;
- TRANSCRIPT_SET_DEFINITIONS.tsv for D2/D3 representative pairs;
- TRANSCRIPT_MEMBERSHIP_MASTER.tsv for protein-coding status;
- GENCODE v33 GTF (00_reference) for exon/CDS structure, reading frame,
  transcript_type and the 'basic' tag;
- local VastDB hg38 PROT_IMPACT (protein-impact ontology) and EVENT_INFO
  (microexon + flanking exon sequences).

Evidence classes used: DIRECT_DATABASE_ANNOTATION (VastDB), GENCODE_DERIVED
(exon/CDS/frame), SEQUENCE_PREDICTION (peptide translation), and
UNAVAILABLE_* rows for external resources that cannot be reached.
"""
import gzip
import os
import re

ROOT = os.environ.get("PROJECT_ROOT", ".")
OUT = os.path.join(ROOT, "36_probability_scale_and_protein",
                   "03_tierA_protein_mapping")
GTF = os.path.join(ROOT, "00_reference/gencode_v33/gencode.v33.annotation.gtf.gz")
MEMBERSHIP = os.path.join(ROOT, "20_psychencode_final_models",
                          "06_transcript_structure_final_check",
                          "00_event_transcript_membership.tsv")
SET_DEF = os.path.join(ROOT, "32_psychencode_sensitivity",
                       "02_psychencode_sensitivity", "transcript_set_intermediates",
                       "TRANSCRIPT_SET_DEFINITIONS.tsv")
MEM_MASTER = os.path.join(ROOT, "32_psychencode_sensitivity",
                          "02_psychencode_sensitivity", "TRANSCRIPT_MEMBERSHIP_MASTER.tsv")
MASTER = os.path.join(ROOT, "25_master_evidence",
                      "06_master_event_table", "MASTER_19_EVENT_EVIDENCE_TABLE.tsv")
PROT_IMPACT = os.path.join(ROOT, "05_vastdb/hg38/PROT_IMPACT-hg38-v3.tab.gz")
EVENT_INFO_LOCAL = os.path.join(os.environ.get("SCRATCH_ROOT", "/tmp"), "tiera_event_info.tsv")

TIERA = {
    "HsaEX0015476": "CLASP1",
    "HsaEX0029786": "HERC4",
    "HsaEX0050855": "PTK2",
    "HsaEX0051138": "PTPRF",
}
GENES = sorted(set(TIERA.values()))

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


def attr(s, key):
    m = re.search(key + r' "([^"]+)"', s)
    return m.group(1) if m else None


def attrs_tags(s):
    return re.findall(r'tag "([^"]+)"', s)


# ---- load final master coordinates ----
master = {}
with open(MASTER) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        r = dict(zip(hdr, line.rstrip('\n').split('\t')))
        if r["HsaEX_ID"] in TIERA:
            master[r["HsaEX_ID"]] = r
print("[probability-scale-5] master rows:", {k: (v["chr_hg38"], v["start_hg38"], v["end_hg38"])
                                   for k, v in master.items()})

# ---- load D0 membership ----
d0 = {ev: {"incl": [], "excl": []} for ev in TIERA}
with open(MEMBERSHIP) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        r = dict(zip(hdr, line.rstrip('\n').split('\t')))
        if r["HsaEX_ID"] in TIERA:
            role = "incl" if r["inclusion_or_exclusion"] == "inclusion" else "excl"
            d0[r["HsaEX_ID"]][role].append(r["transcript_id"])
for ev in TIERA:
    d0[ev]["incl"] = sorted(set(d0[ev]["incl"]))
    d0[ev]["excl"] = sorted(set(d0[ev]["excl"]))

# ---- load D2/D3 definitions ----
dsets = {}
with open(SET_DEF) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        r = dict(zip(hdr, line.rstrip('\n').split('\t')))
        if r["event_id"] in TIERA:
            dsets[r["event_id"]] = r

# ---- load transcript membership master (protein-coding status) ----
mem_master = {}
with open(MEM_MASTER) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        r = dict(zip(hdr, line.rstrip('\n').split('\t')))
        if r.get("event_id") in TIERA:
            mem_master[(r["event_id"], r["transcript_id"].split('.')[0])] = r

# ---- load VastDB PROT_IMPACT ----
onto = {}
with gzip.open(PROT_IMPACT, "rt") as fh:
    fh.readline()
    for line in fh:
        ev, o = line.rstrip('\n').split('\t')
        if ev in TIERA:
            onto[ev] = o
print("[probability-scale-5] VastDB ONTO:", onto)

# ---- load VastDB EVENT_INFO sequences ----
vast_seq = {}
with open(EVENT_INFO_LOCAL) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip('\n').split('\t')
        vast_seq[f[1]] = {"coord_A": f[9], "Seq_C1": f[11], "Seq_A": f[12],
                          "Seq_C2": f[13], "LE_n": int(f[7])}

# ---- parse GENCODE v33 GTF for the 4 genes ----
tx = {}   # tid -> dict
gene_pattern = re.compile('gene_name "(' + '|'.join(GENES) + ')"')
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if not gene_pattern.search(line):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        chrom, feat, start, end, strand, a = f[0], f[2], int(f[3]), int(f[4]), f[6], f[8]
        tid = attr(a, 'transcript_id')
        if tid is None:
            continue
        tid_nv = tid.split('.')[0]
        if tid_nv not in tx:
            tx[tid_nv] = {
                "tid": tid_nv, "gene": attr(a, 'gene_name'), "chrom": chrom,
                "strand": strand, "ttype": attr(a, 'transcript_type'),
                "tags": attrs_tags(a), "exons": [], "cds": [],
                "protein_id": attr(a, 'protein_id')}
        if feat == 'exon':
            tx[tid_nv]["exons"].append((start, end))
        elif feat == 'CDS':
            tx[tid_nv]["cds"].append((start, end, int(a.split('frame ')[1].split(';')[0]) if 'frame ' in a else None))
print(f"[probability-scale-5] parsed {len(tx)} transcripts for {GENES}")

for t in tx.values():
    rev = (t["strand"] == '-')
    t["exons"].sort(reverse=rev)
    t["cds"].sort(reverse=rev)
    t["cds_len"] = sum(e - s + 1 for s, e, _ in t["cds"])


def spliced_cds_segments(t):
    """CDS segments in transcript order with cumulative offsets."""
    segs, off = [], 0
    for s, e, ph in t["cds"]:
        segs.append((s, e, off))
        off += e - s + 1
    return segs


rows = []
for ev, gene in TIERA.items():
    m = master[ev]
    fs, fe = int(m["start_hg38"]), int(m["end_hg38"])
    chrom = m["chr_hg38"]
    len = fe - fs + 1
    sd = dsets[ev]
    seq = vast_seq[ev]
    # GENCODE exon containing the microexon: search inclusion transcripts
    gencode_exon = None
    for tid in d0[ev]["incl"]:
        t = tx.get(tid)
        if not t:
            continue
        for s, e in t["exons"]:
            if s <= fs + 1 and e >= fe and (e - s + 1) <= len:
                cand = (s, e)
                if gencode_exon is None or cand == gencode_exon:
                    gencode_exon = cand
    gs, ge = gencode_exon
    glen = ge - gs + 1
    boundary_shift = (fs - gs, fe - ge)
    print(f"[probability-scale-5] {gene} {ev}: final {chrom}:{fs}-{fe} ({len} nt); "
          f"GENCODE v33 exon {chrom}:{gs}-{ge} ({glen} nt); boundary shift {boundary_shift}")

    for role, tids in (("inclusion", d0[ev]["incl"]), ("exclusion", d0[ev]["excl"])):
        for tid in tids:
            t = tx.get(tid)
            row = {
                "gene": gene, "event_id": ev, "transcript_id": tid,
                "transcript_set": role,
                "transcript_type": t["ttype"] if t else "NOT_IN_GENCODE_V33",
                "gencode_basic": "basic" in (t["tags"] if t else []),
                "strand": t["strand"] if t else "",
                "microexon_reference_hg38": f"{chrom}:{fs}-{fe}",
                "microexon_reference_length_nt": len,
                "microexon_gencode_v33_exon": f"{chrom}:{gs}-{ge}",
                "microexon_gencode_v33_length_nt": glen,
                "vs_gencode_boundary_shift_bp": f"start{boundary_shift[0]:+d}/end{boundary_shift[1]:+d}",
            }
            set_role = []
            if sd.get("D2_incl", "").split('.')[0] == tid: set_role.append("D2_incl")
            if sd.get("D2_excl", "").split('.')[0] == tid: set_role.append("D2_excl")
            if sd.get("D3_incl", "").split('.')[0] == tid: set_role.append("D3_incl")
            if sd.get("D3_excl", "").split('.')[0] == tid: set_role.append("D3_excl")
            row["representative_set"] = ";".join(set_role) if set_role else "D0_only"
            mmr = mem_master.get((ev, tid))
            row["dir32_protein_coding_status"] = mmr["protein_coding_status"] if mmr else "NA"

            if not t or t["cds_len"] == 0:
                row.update({
                    "cds_overlap_nt": 0, "coding_status": "noncoding_transcript_no_CDS",
                    "inclusion_frame": "NA", "in_frame": "NA", "inserted_aa_count": 0,
                    "inserted_peptide": "NA", "junction_context": "NA",
                    "protein_residue_positions": "NA", "vastdb_onto": onto[ev],
                    "consequence": "no coding consequence in this transcript",
                    "ambiguity_note": "transcript has no annotated CDS",
                    "evidence_class": "GENCODE_DERIVED"})
                rows.append(row)
                continue

            segs = spliced_cds_segments(t)
            if role == "inclusion":
                # find CDS segment overlapping the GENCODE microexon exon
                hit = None
                for s, e, off in segs:
                    ov = min(e, ge) - max(s, gs) + 1
                    if ov > 0:
                        hit = (s, e, off, ov)
                        break
                if hit is None:
                    # microexon exon not in CDS -> locate position relative to CDS
                    # count CDS bases upstream of the exon in transcript order
                    if t["strand"] == '-':
                        up = sum(e2 - s2 + 1 for s2, e2, _ in segs if s2 > ge)
                        downstream_cds = any(s2 <= gs for s2, e2, _ in segs)
                        upstream_cds = any(e2 >= ge for s2, e2, _ in segs)
                    else:
                        up = sum(e2 - s2 + 1 for s2, e2, _ in segs if e2 < gs)
                        downstream_cds = any(s2 >= ge for s2, e2, _ in segs)
                        upstream_cds = any(e2 <= gs for s2, e2, _ in segs)
                    if not upstream_cds:
                        loc = "upstream_of_CDS_5UTR"
                    elif not downstream_cds:
                        loc = "downstream_of_CDS_3UTR"
                    else:
                        loc = "intronic_relative_to_CDS?"
                    row.update({
                        "cds_overlap_nt": 0, "coding_status": loc,
                        "inclusion_frame": "NA", "in_frame": "NA", "inserted_aa_count": 0,
                        "inserted_peptide": "NA", "junction_context": "NA",
                        "protein_residue_positions": "NA", "vastdb_onto": onto[ev],
                        "consequence": "microexon outside annotated CDS of this transcript",
                        "ambiguity_note": "event is UTR-level in this transcript",
                        "evidence_class": "GENCODE_DERIVED"})
                    rows.append(row)
                    continue
                s, e, off, ov = hit
                fully = (ov == glen) and (s <= gs) and (e >= ge)
                frame = off % 3
                # CDS bases preceding the microexon = off + (genomic CDS bases in this
                # segment before the microexon start, transcript orientation)
                if t["strand"] == '-':
                    before_in_seg = e - ge
                else:
                    before_in_seg = gs - s
                cds_before = off + before_in_seg
                f0 = cds_before % 3
                if fully and glen % 3 == 0:
                    in_frame = "YES"
                    n_aa = glen // 3
                    seqA = seq["Seq_A"]
                    assert len(seqA) == glen, (gene, len(seqA), glen)
                    c1, c2 = seq["Seq_C1"], seq["Seq_C2"]
                    if f0 == 0:
                        ins_pep = translate(seqA)
                        jctx = (f"insertion between codons; no codon split; "
                                f"inserted peptide = {ins_pep}")
                        ins_window = ins_pep
                        exc_window = ""
                        pos_txt = (f"inserted after residue {cds_before // 3} of the "
                                   f"exclusion isoform; occupies residues "
                                   f"{cds_before // 3 + 1}-{cds_before // 3 + n_aa} "
                                   f"of the inclusion isoform")
                        res_lo = cds_before // 3 + 1
                        res_hi = cds_before // 3 + n_aa
                    else:
                        tail = c1[-f0:] if f0 else ""
                        head = c2[:3 - f0] if f0 else ""
                        window = tail + seqA + head
                        assert len(window) % 3 == 0
                        ins_window = translate(window)
                        exc_window = translate(tail + head)
                        n_new = len(ins_window) - 1  # net residues added
                        jctx = (f"insertion splits a codon (frame {f0}); junction window "
                                f"inclusion={ins_window} exclusion={exc_window}; net +{n_new} residues")
                        ins_pep = ins_window
                        res_lo = cds_before // 3 + 1
                        res_hi = cds_before // 3 + n_aa
                        pos_txt = (f"splits residue {cds_before // 3 + 1} of the exclusion "
                                   f"isoform; insertion window occupies residues "
                                   f"{res_lo}-{res_hi} of the inclusion isoform")
                    row.update({
                        "cds_overlap_nt": ov, "coding_status": "coding_fully_within_CDS",
                        "inclusion_frame": f0, "in_frame": in_frame,
                        "inserted_aa_count": n_aa,
                        "inserted_peptide": seqA if f0 == 0 else f"{seqA} (frame {f0})",
                        "translated_window_inclusion": ins_window,
                        "translated_window_exclusion": exc_window if exc_window else "NA",
                        "junction_context": jctx,
                        "protein_residue_positions": pos_txt,
                        "residue_start_inclusion": res_lo,
                        "residue_end_inclusion": res_hi,
                        "vastdb_onto": onto[ev],
                        "consequence": f"in-frame insertion of {n_aa} aa",
                        "ambiguity_note": "none",
                        "evidence_class": "GENCODE_DERIVED;SEQUENCE_PREDICTION"})
                elif fully and glen % 3 != 0:
                    row.update({
                        "cds_overlap_nt": ov, "coding_status": "coding_fully_within_CDS",
                        "inclusion_frame": f0, "in_frame": "NO",
                        "inserted_aa_count": "NA",
                        "inserted_peptide": "NA",
                        "junction_context": f"frameshift ({glen} nt mod 3 = {glen % 3})",
                        "protein_residue_positions": "NA",
                        "vastdb_onto": onto[ev],
                        "consequence": "FRAMESHIFT upon inclusion",
                        "ambiguity_note": "none",
                        "evidence_class": "GENCODE_DERIVED"})
                else:
                    row.update({
                        "cds_overlap_nt": ov, "coding_status": "partial_CDS_overlap",
                        "inclusion_frame": f0, "in_frame": "PARTIAL",
                        "inserted_aa_count": "NA", "inserted_peptide": "NA",
                        "junction_context": f"microexon straddles a CDS boundary (overlap {ov}/{glen})",
                        "protein_residue_positions": "NA",
                        "vastdb_onto": onto[ev],
                        "consequence": "microexon straddles CDS boundary",
                        "ambiguity_note": "requires boundary-specific interpretation",
                        "evidence_class": "GENCODE_DERIVED"})
                rows.append(row)
            else:
                # exclusion transcript: the microexon is skipped, so no CDS segment
                # overlaps the skipped interval; count CDS bases strictly 5' of it
                if t["strand"] == '-':
                    cds_before = sum(e2 - s2 + 1 for s2, e2, _ in segs if s2 > fe)
                    # add portion of any CDS segment overlapping (ge..fe) region 5' side: none (skipped)
                else:
                    cds_before = sum(e2 - s2 + 1 for s2, e2, _ in segs if e2 < fs)
                f_excl = cds_before % 3
                row.update({
                    "cds_overlap_nt": 0,
                    "coding_status": "coding_skipped_exon_reference",
                    "inclusion_frame": f_excl,
                    "in_frame": "reference",
                    "inserted_aa_count": 0,
                    "inserted_peptide": "none (exclusion isoform)",
                    "junction_context": f"CDS bases 5' of skipped microexon = {cds_before}",
                    "protein_residue_positions": "NA",
                    "vastdb_onto": onto[ev],
                    "consequence": "exclusion isoform lacks the microexon peptide",
                    "ambiguity_note": "none",
                    "evidence_class": "GENCODE_DERIVED"})
                rows.append(row)

# ---- write ----
cols = ["gene", "event_id", "transcript_id", "transcript_set", "representative_set",
        "transcript_type", "dir32_protein_coding_status", "gencode_basic", "strand",
        "microexon_reference_hg38", "microexon_reference_length_nt",
        "microexon_gencode_v33_exon", "microexon_gencode_v33_length_nt",
        "vs_gencode_boundary_shift_bp", "cds_overlap_nt", "coding_status",
        "inclusion_frame", "in_frame", "inserted_aa_count", "inserted_peptide",
        "translated_window_inclusion", "translated_window_exclusion",
        "junction_context", "protein_residue_positions",
        "residue_start_inclusion", "residue_end_inclusion",
        "vastdb_onto", "consequence", "ambiguity_note", "evidence_class"]
for r in rows:
    for c in cols:
        r.setdefault(c, "NA")
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "TIER_A_CODING_CONSEQUENCE.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"[probability-scale-5] wrote {len(rows)} rows to TIER_A_CODING_CONSEQUENCE.tsv")
for r in rows:
    if r["representative_set"] != "D0_only" or r["transcript_set"] == "inclusion":
        print("  ", r["gene"], r["transcript_id"], r["transcript_set"],
              r["representative_set"], "->", r["consequence"],
              "| pep:", r.get("inserted_peptide"))
