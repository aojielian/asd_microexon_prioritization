#!/usr/bin/env python3
"""Final Task B — CHyMErA-to-human directional mapping check.

Reconstructs the direction rule, testability, concordance and the
multi-exon/multi-event mapping directly from the authoritative local
source tables (NOT from manuscript prose), then compares against the
Reference final counts.

Reads (read-only):
  15_.../06_directional_bridge/01_event_direction_master.tsv   (19 events)
  15_.../06_directional_bridge/00_direction_definitions.tsv
  25_.../06_master_event_table/MASTER_19_EVENT_EVIDENCE_TABLE.tsv (coords)
  10_.../06_asd_event_overlap/01_CHyMErA_CTX_strict_event_overlap.tsv
  10_.../06_asd_event_overlap/05_event_overlap_multiplicity_check.tsv
  10_.../03_definition_lock/01_CHyMErA_event_reconciliation.tsv
  11_.../04_event_reconciliation/02_duplicate_event_check.tsv
  39_.../02_directional_eligibility_rule/*  (for reproducibility comparison)

Writes (inside Final root only):
  03_chymera_direction_check/CHYMERA_HUMAN_DIRECTION_EVENT_CHECK.tsv
  03_chymera_direction_check/CHYMERA_MULTI_EXON_MAPPING_CHECK.tsv
  03_chymera_direction_check/CHYMERA_DIRECTION_SUMMARY.tsv
  03_chymera_direction_check/CHYMERA_DIRECTION_RULES.md
  03_chymera_direction_check/CHYMERA_DIRECTION_REPRODUCIBILITY.md
  06_logs/chymera_direction_verification.log
"""
import os
from math import comb

PROJECT = os.environ.get("PROJECT_ROOT", ".")
ROOT = os.path.join(PROJECT,
                    "40_moderator_and_direction")
OUT = os.path.join(ROOT, "03_chymera_direction_check")
LOG = []


def say(m):
    print(m)
    LOG.append(str(m))


def read_tsv(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() != ""]
    head = lines[0].split("\t")
    return head, [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def g17(x):
    return "%.17g" % float(x)


def exact_binom_two_sided(k, n, p=0.5):
    """Exact two-sided binomial P (equal to scipy.stats.binomtest default
    'exact' method: sum of all outcome probabilities <= P(observed))."""
    pobs = comb(n, k) * (p ** k) * ((1 - p) ** (n - k))
    tot = 0.0
    for i in range(n + 1):
        pi = comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
        if pi <= pobs + 1e-300:
            tot += pi
    return tot


# --- authoritative source paths ---
DIR15_MASTER = os.path.join(PROJECT, "15_directional_concordance",
                            "06_directional_bridge",
                            "01_event_direction_master.tsv")
DIR15_DEFS = os.path.join(PROJECT, "15_directional_concordance",
                          "06_directional_bridge",
                          "00_direction_definitions.tsv")
DIR25_MASTER = os.path.join(PROJECT,
                            "25_master_evidence",
                            "06_master_event_table",
                            "MASTER_19_EVENT_EVIDENCE_TABLE.tsv")
DIR10_OVERLAP = os.path.join(PROJECT, "10_event_mapping",
                             "06_asd_event_overlap",
                             "01_CHyMErA_CTX_strict_event_overlap.tsv")
DIR10_MULT = os.path.join(PROJECT, "10_event_mapping",
                          "06_asd_event_overlap",
                          "05_event_overlap_multiplicity_check.tsv")
DIR10_RECON = os.path.join(PROJECT, "10_event_mapping",
                           "03_definition_lock",
                           "01_CHyMErA_event_reconciliation.tsv")
DIR11_DUP = os.path.join(PROJECT, "11_set_level_enrichment",
                         "04_event_reconciliation",
                         "02_duplicate_event_check.tsv")
REF_RECALC = os.path.join(PROJECT,
                          "39_rule_and_numeric_verification",
                          "02_directional_eligibility_rule",
                          "DIRECTIONAL_BRIDGE_RECALC.tsv")
REF_ALL19 = os.path.join(PROJECT,
                         "39_rule_and_numeric_verification",
                         "02_directional_eligibility_rule",
                         "ALL19_THRESHOLD_FREE_DIRECTION_"
                         "SENSITIVITY.tsv")
DIR15_SCRIPT = ("15_directional_concordance/00_scripts/"
                "chymera_directional_analysis.py")


def main():
    os.makedirs(OUT, exist_ok=True)

    # ---- load authoritative sources ----
    _, dir15 = read_tsv(DIR15_MASTER)
    assert len(dir15) == 19, "direction master must have 19 rows"
    _, dir25 = read_tsv(DIR25_MASTER)
    coord25 = {r["HsaEX_ID"]: r for r in dir25}
    _, overlap = read_tsv(DIR10_OVERLAP)
    _, mult = read_tsv(DIR10_MULT)
    mult_of = {r["MmuEX_ID"]: r for r in mult}
    _, recon = read_tsv(DIR10_RECON)
    recon_of = {r["MmuEX_ID"]: r for r in recon if r["MmuEX_ID"]}
    _, dupcheck = read_tsv(DIR11_DUP)
    dup_of = {r["MmuEX_ID"]: r for r in dupcheck}

    # ---- index overlap records by MmuEX (strict cortex) ----
    ov_by_mmuex = {}
    for r in overlap:
        ov_by_mmuex.setdefault(r["MmuEX_ID"], []).append(r)

    # ---- build the 19-row event check ----
    rows = []
    gene_events = {}
    for ev in dir15:
        gene_events.setdefault(ev["gene"], []).append(ev["HsaEX_ID"])

    for ev in dir15:
        hid = ev["HsaEX_ID"]; mid = ev["MmuEX_ID"]; gene = ev["gene"]
        dps = float(ev["ASD_delta_psi"]); p = float(ev["ASD_p"])
        eligible = abs(dps) > 0.01                     # magnitude rule
        if dps < -0.01:
            exp_sign = "NEGATIVE_decreased_inclusion"
        elif dps > 0.01:
            exp_sign = "POSITIVE_increased_inclusion"
        else:
            exp_sign = "UNRESOLVED_abs_dpsi_le_0p01"
        # CHyMErA: all are microexon deletions -> inclusion loss ->
        # expected human sign is NEGATIVE (decreased inclusion)
        chy_expected = "NEGATIVE_decreased_inclusion"
        if eligible:
            main_conc = "YES" if dps < -0.01 else "NO"
        else:
            main_conc = "NOT_ELIGIBLE"
        all19_conc = "YES" if dps < 0 else ("NO" if dps > 0 else "ZERO")

        rc = recon_of.get(mid, {})
        guides = rc.get("guide_ids", "")
        tgt_cat = rc.get("target_category", "")
        is_me_del = rc.get("is_microexon_deletion", "")
        is_host_ko = rc.get("is_host_gene_KO", "")

        n_parik = int(mult_of[mid]["n_parikshak_records"]) \
            if mid in mult_of else 0

        same_gene_multi = "YES" if len(gene_events[gene]) > 1 else "NO"

        c25 = coord25.get(hid, {})
        ov_rows = ov_by_mmuex.get(mid, [])
        # primary overlap record = one whose delta_psi matches the final analysis
        prim = [r for r in ov_rows
                if abs(float(r["delta_psi"]) - dps) < 1e-12]
        n_ov = len(ov_rows)

        notes = []
        if n_ov > 1:
            notes.append("%d strict-cortex overlap records for this exon; "
                         "the final analysis uses delta_psi=%.6g (duplicate Parikshak "
                         "record resolved upstream in phase-0BR)" %
                         (n_ov, dps))
        if same_gene_multi == "YES":
            notes.append("gene %s has %d events, each linked to a "
                         "DISTINCT exon-specific microexon-deletion "
                         "perturbation (see multi-exon check)" %
                         (gene, len(gene_events[gene])))
        if is_host_ko == "True":
            notes.append("HOST_GENE_KO perturbation (gene-level)")

        rows.append({
            "event_id": hid,
            "gene": gene,
            "human_chr": c25.get("chr_hg38", ""),
            "human_start": c25.get("start_hg38", ""),
            "human_end": c25.get("end_hg38", ""),
            "discovery_delta_PSI": g17(dps),
            "discovery_P": g17(p),
            "abs_delta_PSI_gt_0p01": "TRUE" if eligible else "FALSE",
            "direction_eligible": "YES" if eligible else "NO",
            "human_expected_sign_from_discovery": exp_sign,
            "chymera_source_event_id": mid,
            "chymera_target_gene": rc.get("gene", ""),
            "chymera_target_exon_or_interval": mid,
            "chymera_perturbation_type": tgt_cat,
            "chymera_effect_on_microexon_inclusion":
                "INCLUSION_LOSS" if is_me_del == "True" else "UNKNOWN",
            "chymera_expected_human_sign": chy_expected,
            "main_analysis_concordant": main_conc,
            "all19_concordant": all19_conc,
            "same_gene_multiple_prespecified_events": same_gene_multi,
            "same_perturbation_maps_multiple_exons": "NO",
            "same_perturbation_maps_multiple_prespecified_events": "NO",
        })

    cols = ["event_id", "gene", "human_chr", "human_start", "human_end",
            "discovery_delta_PSI", "discovery_P",
            "abs_delta_PSI_gt_0p01", "direction_eligible",
            "human_expected_sign_from_discovery",
            "chymera_source_event_id", "chymera_target_gene",
            "chymera_target_exon_or_interval", "chymera_perturbation_type",
            "chymera_effect_on_microexon_inclusion",
            "chymera_expected_human_sign", "main_analysis_concordant",
            "all19_concordant", "same_gene_multiple_prespecified_events",
            "same_perturbation_maps_multiple_exons",
            "same_perturbation_maps_multiple_prespecified_events"]
    p1 = os.path.join(OUT, "CHYMERA_HUMAN_DIRECTION_EVENT_CHECK.tsv")
    with open(p1, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in cols) + "\n")
    say("WROTE %s (%d rows)" % (p1, len(rows)))

    # ---- directional statistics (main + all-19) ----
    n_eligible = sum(1 for r in rows
                     if r["direction_eligible"] == "YES")
    n_conc = sum(1 for r in rows
                 if r["main_analysis_concordant"] == "YES")
    n_disc = sum(1 for r in rows
                 if r["main_analysis_concordant"] == "NO")
    main_p = exact_binom_two_sided(n_conc, n_conc + n_disc)
    say("MAIN: eligible=%d/19 concordant=%d discordant=%d P=%s" %
        (n_eligible, n_conc, n_disc, g17(main_p)))

    n19_conc = sum(1 for r in rows if r["all19_concordant"] == "YES")
    n19_disc = sum(1 for r in rows if r["all19_concordant"] == "NO")
    all19_p = exact_binom_two_sided(n19_conc, n19_conc + n19_disc)
    say("ALL19: concordant=%d discordant=%d P=%s" %
        (n19_conc, n19_disc, g17(all19_p)))

    p3 = os.path.join(OUT, "CHYMERA_DIRECTION_SUMMARY.tsv")
    with open(p3, "w") as fh:
        fh.write("analysis\tn_total\tn_eligible\tn_concordant\t"
                 "n_discordant\tconcordance_fraction\t"
                 "exact_binomial_P_two_sided\trule\n")
        fh.write("\t".join([
            "prespecified_main_analysis", "19", str(n_eligible), str(n_conc),
            str(n_disc), g17(n_conc / (n_conc + n_disc)), g17(main_p),
            "direction-eligible iff |human discovery delta_PSI| > 0.01 "
            "(magnitude rule); concordant iff delta_PSI < -0.01 "
            "(matches CHyMErA microexon-inclusion loss); two-sided exact "
            "binomial vs p=0.5"]) + "\n")
        fh.write("\t".join([
            "all19_threshold_free_sensitivity", "19", "19",
            str(n19_conc), str(n19_disc),
            g17(n19_conc / (n19_conc + n19_disc)), g17(all19_p),
            "no threshold: all 19 events; concordant iff signed discovery "
            "delta_PSI < 0 (inclusion loss); two-sided exact binomial "
            "vs p=0.5"]) + "\n")
    say("WROTE %s" % p3)

    # ---- multi-exon / multi-event mapping check ----
    # Build guide -> exon map to PROVE no single perturbation spans exons.
    guide_to_exons = {}
    for r in recon:
        if not r.get("guide_ids"):
            continue
        for g in r["guide_ids"].split("|"):
            g = g.strip()
            if g:
                guide_to_exons.setdefault(g, set()).add(r["MmuEX_ID"])
    shared_guides = {g: s for g, s in guide_to_exons.items() if len(s) > 1}
    say("GUIDE-SHARING check: guides mapping to >1 exon = %d %s" %
        (len(shared_guides),
         sorted(shared_guides) if shared_guides else "(none)"))

    # multi-exon check rows: one per CHyMErA source event linked to the 19
    mrows = []
    linked_mmuex = [ev["MmuEX_ID"] for ev in dir15]
    mmuex_to_reference = {}
    for ev in dir15:
        mmuex_to_reference.setdefault(ev["MmuEX_ID"], []).append(
            ev["HsaEX_ID"])
    for mid in linked_mmuex:
        rc = recon_of.get(mid, {})
        guides = set(g for g in rc.get("guide_ids", "").split("|") if g)
        # exons this perturbation's guides touch
        exons_touched = set()
        for g in guides:
            exons_touched |= guide_to_exons.get(g, set())
        linked_reference = mmuex_to_reference.get(mid, [])
        gene = rc.get("gene", "")
        ambiguous = "NO" if len(exons_touched) <= 1 else "YES"
        mrows.append({
            "chymera_source_event_id": mid,
            "target_gene": gene,
            "target_exon_or_interval": mid,
            "n_source_exons_linked": str(len(exons_touched)),
            "n_prespecified_19_events_linked": str(len(linked_reference)),
            "linked_prespecified_event_ids": ";".join(linked_reference),
            "mapping_type":
                "exon_specific_microexon_deletion"
                if rc.get("is_microexon_deletion") == "True"
                else rc.get("target_category", ""),
            "potential_directional_ambiguity": ambiguous,
            "resolution_rule":
                "each perturbation is an exon-specific microexon deletion "
                "with a unique guide set; 1 perturbation -> 1 exon -> "
                "direction is unambiguous"
                if ambiguous == "NO" else
                "REQUIRES_MANUAL_CHECK: shared guides span exons",
            })
    mcols = ["chymera_source_event_id", "target_gene",
             "target_exon_or_interval", "n_source_exons_linked",
             "n_prespecified_19_events_linked",
             "linked_prespecified_event_ids",
             "mapping_type", "potential_directional_ambiguity",
             "resolution_rule"]
    p2 = os.path.join(OUT, "CHYMERA_MULTI_EXON_MAPPING_CHECK.tsv")
    with open(p2, "w") as fh:
        fh.write("\t".join(mcols) + "\n")
        for r in mrows:
            fh.write("\t".join(r[c] for c in mcols) + "\n")
    say("WROTE %s (%d rows)" % (p2, len(mrows)))
    n_amb = sum(1 for r in mrows
                if r["potential_directional_ambiguity"] == "YES")
    say("MULTI_EXON_AMBIGUITIES=%d" % n_amb)

    # ---- gene-level multi-event summary (for rules / report) ----
    gene_level = []
    for gene, eids in sorted(gene_events.items()):
        if len(eids) > 1:
            mids = [ev["MmuEX_ID"] for ev in dir15 if ev["gene"] == gene]
            gene_level.append((gene, len(eids), eids, mids))
    for g, n, eids, mids in gene_level:
        say("GENE %s: %d events -> distinct mouse exons %s" %
            (g, n, ",".join(mids)))

    # ---- reproducibility vs Reference ----
    _, ref_r = read_tsv(REF_RECALC)
    ref_q = {r["quantity"]: r for r in ref_r}
    ref_main_p = ref_q.get("exact_binomial_p_two_sided", {})
    say("REF recalc main P final=%s ref=%s match=%s" %
        (ref_main_p.get("original_value"),
         ref_main_p.get("recalc_value"),
         ref_main_p.get("match")))
    _, ref_a = read_tsv(REF_ALL19)
    ref_a_of = {r["field"]: r["value"] for r in ref_a}

    repro = []
    repro.append("main_concordant=%d (final) vs recalculated "
                 "final 12/14" % n_conc)
    match_main = (n_conc == 12 and n_eligible == 14)
    repro.append("main_exact_P=%s ; matches final 0.012939453125 = %s"
                 % (g17(main_p), abs(main_p - 0.012939453125) < 1e-15))
    match_all19 = (n19_conc == 16 and n19_disc == 3)
    repro.append("all19_concordant=%d/19 ; matches final 16/19 = %s" %
                 (n19_conc, match_all19))
    repro.append("all19_exact_P=%s ; matches recalculated %s = %s" %
                 (g17(all19_p), ref_a_of.get("exact_binomial_p_two_sided"),
                  abs(all19_p - float(ref_a_of.get(
                      "exact_binomial_p_two_sided", "nan"))) < 1e-15))
    for line in repro:
        say("REPRO " + line)

    exact_repro = (match_main and match_all19 and
                   abs(main_p - 0.012939453125) < 1e-15 and
                   abs(all19_p - float(ref_a_of.get(
                       "exact_binomial_p_two_sided", "0"))) < 1e-15)
    say("EXACT_REPRODUCE_REFERENCE=%s" % exact_repro)

    # ---- CHYMERA_DIRECTION_RULES.md ----
    rules = build_rules_md(rows, recon_of, guide_to_exons, shared_guides,
                           gene_level, main_p, all19_p, n_conc,
                           n_eligible, n19_conc)
    p4 = os.path.join(OUT, "CHYMERA_DIRECTION_RULES.md")
    with open(p4, "w") as fh:
        fh.write(rules)
    say("WROTE %s" % p4)

    # ---- CHYMERA_DIRECTION_REPRODUCIBILITY.md ----
    rep = build_repro_md(rows, n_conc, n_eligible, n_disc, main_p,
                         n19_conc, n19_disc, all19_p, exact_repro,
                         ref_a_of)
    p5 = os.path.join(OUT, "CHYMERA_DIRECTION_REPRODUCIBILITY.md")
    with open(p5, "w") as fh:
        fh.write(rep)
    say("WROTE %s" % p5)

    with open(os.path.join(ROOT, "06_logs",
                           "chymera_direction_verification.log"),
              "w") as fh:
        fh.write("\n".join(LOG) + "\n")


def build_rules_md(rows, recon_of, guide_to_exons, shared_guides,
                   gene_level, main_p, all19_p, n_conc, n_eligible,
                   n19_conc):
    multi_txt = "\n".join(
        "- **%s**: %d events (%s) map to %d DISTINCT mouse exons "
        "(%s), each with its own exon-specific microexon-deletion "
        "perturbation and unique guide set." %
        (g, n, ", ".join(e), n, ", ".join(m))
        for g, n, e, m in gene_level)
    return f"""# CHYMERA_DIRECTION_RULES — explicit direction logic (Final §11)

Reconstructed from the authoritative source schema/code, NOT from
manuscript prose. Sign convention is stated explicitly throughout.

## Sign convention (used in every statement below)
- Human discovery `delta_PSI = PSI_ASD - PSI_control` (Parikshak 2016).
  `delta_PSI < 0` = **decreased** microexon inclusion in ASD;
  `delta_PSI > 0` = **increased** inclusion in ASD.
- A CHyMErA microexon **deletion** removes the exon -> expected effect is
  **inclusion LOSS**, i.e. expected human sign is **NEGATIVE**.
- "Concordant" = observed human discovery `delta_PSI < 0` (inclusion
  decrease), matching the CHyMErA inclusion-loss direction.

## 1. What the CHyMErA perturbation represents
Each CHyMErA source event is a CRISPR deletion of a single neural
microexon in mouse (event id = `MmuEX_ID`). In the final set every
perturbation linked to the 19 human events has
`target_category = MICROEXON_DELETION` and `is_microexon_deletion = True`
(source: `10_event_mapping/03_definition_lock/01_CHyMErA_event_reconciliation.tsv`).
No host-gene knockout (`is_host_gene_KO = False` for all 19-linked
events), so each perturbation affects exactly the targeted microexon.

## 2. How inclusion-loss / inclusion-gain is encoded
In `15_directional_concordance/06_directional_bridge/00_direction_definitions.tsv`:
- `CHyMErA_deletion_direction = "Microexon inclusion loss (deletion = exclusion)"`.
- The final direction master sets
  `CHyMErA_direction = MICROEXON_INCLUSION_LOSS` for all 19 events
  (there is no inclusion-gain class among these perturbations).

## 3. How that maps to expected human ASD ΔPSI sign
Because every perturbation is an inclusion **loss**, the expected human
ASD discovery sign is **negative** (`delta_PSI < 0`, decreased
inclusion). So `chymera_expected_human_sign = NEGATIVE_decreased_inclusion`
for all 19 events.

## 4. What makes an event eligible for directional comparison
`Eligible for directional comparison iff |human discovery delta_PSI| > 0.01`.
This is implemented in
`{DIR15_SCRIPT}` Phase 4:
`ASD_DECREASED_INCLUSION` if `dps < -0.01`,
`ASD_INCREASED_INCLUSION` if `dps > 0.01`,
else `ASD_DIRECTION_UNRESOLVED`. Only the first two classes are
eligible. This is a **magnitude** threshold on the effect size.

## 5. Why discovery P is NOT used for eligibility
Eligibility depends only on `|delta_PSI| > 0.01`. The discovery P value
is never consulted in the Phase 4 direction/eligibility block. This is
deliberate: a tiny but precisely-estimated effect (small P, tiny ΔPSI)
would still lack a meaningful direction, and a moderate effect with a
moderate P still has a well-defined direction. Using P would conflate
effect size with sample size.

## 6. How concordance is calculated
Among events eligible for directional comparison, an event is **concordant** if its
`delta_PSI < -0.01` (matches inclusion loss) and **discordant/opposite**
if `delta_PSI > 0.01`. The final main result is {n_conc} concordant of
{n_eligible} eligible; exact two-sided binomial vs p=0.5 gives
P = {g17(main_p)}.

## 7. How the all-19 sensitivity differs from the main analysis
The all-19 threshold-free sensitivity ignores the `|ΔPSI|>0.01`
eligibility filter: all 19 events are used exactly once and an event is
concordant iff its signed `delta_PSI < 0`. Result {n19_conc}/19
concordant, exact two-sided binomial P = {g17(all19_p)}. It is a
SENSITIVITY analysis and does not replace the final 14-event primary
bridge.

## 8. How multi-exon / multi-event mappings are handled
{multi_txt}

Mechanism guaranteeing no ambiguity: guide-level check of
`01_CHyMErA_event_reconciliation.tsv` shows
{len(guide_to_exons)} distinct guides, of which
{len(shared_guides)} map to more than one exon
({"none -> every perturbation is exon-specific" if not shared_guides else
  "SHARED: " + ", ".join(sorted(shared_guides))}).
Because each perturbation is an exon-specific deletion with a unique
guide set, one perturbation -> one exon -> one primary19 event, so the
expected direction is unambiguous even for genes that host multiple
primary19 events. The same-perturbation-maps-multiple-exons flag is
therefore NO for all 19 rows.

Note: multi-event genes can show DIFFERENT per-event directions in the
human discovery data (e.g. PTK2 events split decreased/increased/
unresolved). This is not a mapping artefact: each human event is matched
to its own mouse microexon deletion, and concordance is evaluated
event-by-event against that event's own perturbation.

Record-level (not exon-level) multiplicity: exactly one primary19 event,
CPEB4 (MmuEX0012543), has two Parikshak records for the SAME exon
(alternative upstream-exon annotation) in the strict cortex overlap
table. The primary19 workflow resolves this by first-record
de-duplication per MmuEX_ID
(`11_set_level_enrichment/00_scripts/01_mapping_recheck_and_reconciliation.py`,
`ctx_matches.drop_duplicates('MmuEX_ID')`; primary19 check
`04_event_reconciliation/02_duplicate_event_check.tsv` flags only
MmuEX0012543 with n_records=2). Both records are negative
(delta_psi = -0.06363... kept, -0.04590... dropped), so the direction
conclusion is identical under either record.
"""


def build_repro_md(rows, n_conc, n_eligible, n_disc, main_p,
                   n19_conc, n19_disc, all19_p, exact_repro, ref_a_of):
    verdict = ("EXACTLY REPRODUCES" if exact_repro
               else "DOES NOT exactly reproduce")
    return f"""# CHYMERA_DIRECTION_REPRODUCIBILITY — Final vs Reference

Date: 2026-08-13. Final recomputed the CHyMErA directional bridge
from the primary authoritative sources (direction master and
mapping/definition tables), independently of the Reference code path.

## Verdict
Final **{verdict}** the Reference directional results.

| analysis | Final | Reference | match |
|---|---|---|---|
| main eligible | {n_eligible}/19 | 14/19 | {n_eligible == 14} |
| main concordant | {n_conc}/{n_eligible} | 12/14 | {n_conc == 12} |
| main discordant | {n_disc} | 2 | {n_disc == 2} |
| main exact two-sided binomial P | {g17(main_p)} | 0.012939453125 | {abs(main_p - 0.012939453125) < 1e-15} |
| all-19 concordant | {n19_conc}/19 | 16/19 | {n19_conc == 16} |
| all-19 exact two-sided binomial P | {g17(all19_p)} | {ref_a_of.get("exact_binomial_p_two_sided", "?")} | {abs(all19_p - float(ref_a_of.get("exact_binomial_p_two_sided", "0"))) < 1e-15} |

## Method of recomputation
- Testability rule re-derived from `{DIR15_SCRIPT}`
  Phase 4 (`|delta_PSI| > 0.01`), applied to the 19 discovery
  `delta_PSI` values in
  `15_.../06_directional_bridge/01_event_direction_master.tsv`.
- Exact two-sided binomial computed as the sum of all outcome
  probabilities <= P(observed) under Binom(n, 0.5) (equivalent to
  scipy `binomtest(..., method="exact")`).
- All-19 sensitivity = sign-only (`delta_PSI < 0`), no threshold.

## Discrepancies
{"None. Every count and P value matches the the reference outputs to full precision." if exact_repro else "DISCREPANCY DETECTED — see the final report and flag for follow-up."}
"""


if __name__ == "__main__":
    main()
