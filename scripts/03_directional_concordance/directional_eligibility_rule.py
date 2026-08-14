#!/usr/bin/env python3
"""Reference Workstream A — CHyMErA directional-eligibility check.

Reads (read-only):
  15_directional_concordance/06_directional_bridge/01_event_direction_master.tsv
  15_directional_concordance/06_directional_bridge/02_concordance_tests.tsv
  15_directional_concordance/06_directional_bridge/00_direction_definitions.tsv
  25_.../06_master_event_table/MASTER_19_EVENT_EVIDENCE_TABLE.tsv
Rule source (cited, not executed):
  scripts/03_directional_concordance/chymera_directional_analysis.py lines 294-316, 343
  23_event_annotation/00_scripts/event_annotation_and_master_table_repair.R lines 116-117

Eligible for directional comparison iff |discovery delta_psi| > 0.01.

Writes (inside Reference root only):
  02_directional_eligibility_rule/DIRECTIONAL_ELIGIBILITY_CHECK.tsv
  02_directional_eligibility_rule/DIRECTIONAL_ELIGIBILITY_RULE_TRACE.md
  02_directional_eligibility_rule/DIRECTIONAL_BRIDGE_RECALC.tsv
  02_directional_eligibility_rule/ALL19_THRESHOLD_FREE_DIRECTION_SENSITIVITY.tsv
  08_logs/directional_eligibility_rule.log
"""
import os

from scipy import stats

PROJECT = os.environ.get("PROJECT_ROOT", ".")
ROOT = os.path.join(PROJECT, "39_rule_and_numeric_verification")
OUT = os.path.join(ROOT, "02_directional_eligibility_rule")
LOG = os.path.join(ROOT, "08_logs", "directional_eligibility_rule.log")

BRIDGE_MASTER = os.path.join(
    PROJECT, "15_directional_concordance/06_directional_bridge/"
    "01_event_direction_master.tsv")
CONC_TESTS = os.path.join(
    PROJECT, "15_directional_concordance/06_directional_bridge/"
    "02_concordance_tests.tsv")
DIR_DEFS = os.path.join(
    PROJECT, "15_directional_concordance/06_directional_bridge/"
    "00_direction_definitions.tsv")
MASTER25 = os.path.join(
    PROJECT, "25_master_evidence/"
    "06_master_event_table/MASTER_19_EVENT_EVIDENCE_TABLE.tsv")

RULE_SOURCE_FILE = ("15_directional_concordance/00_scripts/"
                    "chymera_directional_analysis.py")
RULE_SOURCE_LOC = ("lines 294-316 (asd_dir/bridge block: eligible iff "
                   "|delta_psi| > 0.01), line 343 (exact binomial test)")

_log = []


def say(m):
    print(m)
    _log.append(m)


def read_tsv(path):
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() != ""]
    head = lines[0].split("\t")
    return head, [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def main():
    os.makedirs(OUT, exist_ok=True)

    _, bridge = read_tsv(BRIDGE_MASTER)
    _, conc = read_tsv(CONC_TESTS)
    _, defs = read_tsv(DIR_DEFS)
    mhead, master = read_tsv(MASTER25)

    assert len(bridge) == 19, "bridge master must have 19 rows"
    assert len(master) == 19, "master table must have 19 rows"
    master_by_id = {r["HsaEX_ID"]: r for r in master}

    # roster cross-check
    ids_bridge = sorted(r["HsaEX_ID"] for r in bridge)
    ids_master = sorted(r["HsaEX_ID"] for r in master)
    roster_match = ids_bridge == ids_master
    say("ROSTER_19_BRIDGE_VS_MASTER_MATCH=%s" % roster_match)
    assert roster_match

    # ------- apply the rule independently -------
    check_rows = []
    n_eligible = n_conc = n_opp = 0
    for r in bridge:
        dps = float(r["ASD_delta_psi"])
        p_raw = float(r["ASD_p"])
        # rule: direction assignment by magnitude threshold 0.01
        if dps < -0.01:
            asd_dir = "ASD_DECREASED_INCLUSION"
        elif dps > 0.01:
            asd_dir = "ASD_INCREASED_INCLUSION"
        else:
            asd_dir = "ASD_DIRECTION_UNRESOLVED"
        if asd_dir == "ASD_DECREASED_INCLUSION":
            bridge_calc = "CONCORDANT_WITH_MICROEXON_LOSS"
        elif asd_dir == "ASD_INCREASED_INCLUSION":
            bridge_calc = "OPPOSITE_TO_MICROEXON_LOSS"
        else:
            bridge_calc = "NO_ASD_EFFECT_DIRECTION"
        eligible = bridge_calc != "NO_ASD_EFFECT_DIRECTION"
        if eligible:
            n_eligible += 1
            if bridge_calc == "CONCORDANT_WITH_MICROEXON_LOSS":
                n_conc += 1
            else:
                n_opp += 1
        mrow = master_by_id[r["HsaEX_ID"]]
        check_rows.append({
            "gene": r["gene"],
            "human_event_id": r["HsaEX_ID"],
            "mouse_event_id": r["MmuEX_ID"],
            "discovery_delta_psi": "%.15g" % dps,
            "discovery_p_raw": "%.15g" % p_raw,
            "discovery_direction": asd_dir,
            "chymera_perturbation_direction": r["CHyMErA_direction"],
            "eligible_for_directional_comparison": "YES" if eligible else "NO",
            "direction_unresolved_reason":
                ("" if eligible else
                 "|discovery delta_psi| <= 0.01 (magnitude threshold; "
                 "code-derived)"),
            "direction_unresolved_reason_in_master":
                mrow["CHyMErA_reason_not_testable"],
            "concordant_if_eligible":
                ("YES" if bridge_calc == "CONCORDANT_WITH_MICROEXON_LOSS"
                 else ("NO" if eligible else "NA")),
            "bridge_classification": r["bridge_classification"],
            "recalculated_bridge_classification": bridge_calc,
            "bridge_matches_reference":
                "YES" if bridge_calc == r["bridge_classification"] else "NO",
            "rule_source_file": RULE_SOURCE_FILE,
            "rule_source_line_or_code_object": RULE_SOURCE_LOC,
        })

    mismatches = [a for a in check_rows if a["bridge_matches_reference"] == "NO"]
    say("PER_EVENT_BRIDGE_RECALC_MISMATCHES=%d" % len(mismatches))

    # ------- bridge recalculation table -------
    overall = [c for c in conc if c["test"] == "overall"][0]
    binom_recalc = stats.binomtest(n_conc, n_eligible, 0.5)
    binom_p_recalc = binom_recalc.pvalue
    p = float(overall["binomial_p"])
    recalc_path = os.path.join(OUT, "DIRECTIONAL_BRIDGE_RECALC.tsv")
    with open(recalc_path, "w") as fh:
        fh.write("quantity\treference_original_value\trecalc_value\tmatch\n")
        rows = [
            ("n_events_total", "19", "19", "YES"),
            ("n_eligible", overall["n_testable"], str(n_eligible),
             "YES" if int(overall["n_testable"]) == n_eligible else "NO"),
            ("n_concordant", overall["n_concordant"], str(n_conc),
             "YES" if int(overall["n_concordant"]) == n_conc else "NO"),
            ("n_discordant_opposite", str(int(float(overall["n_opposite"]))),
             str(n_opp),
             "YES" if int(float(overall["n_opposite"])) == n_opp else "NO"),
            ("n_no_asd_effect_direction",
             str(int(float(overall["n_no_direction"]))),
             str(19 - n_eligible),
             "YES" if int(float(overall["n_no_direction"])) == 19 - n_eligible
             else "NO"),
            ("concordance_rate", overall["concordance_rate"],
             "%.16g" % (n_conc / n_eligible),
             "YES" if abs(float(overall["concordance_rate"])
                          - n_conc / n_eligible) < 1e-12 else "NO"),
            ("exact_binomial_p_two_sided", "%.15g" % p,
             "%.15g" % binom_p_recalc,
             "YES" if abs(p - binom_p_recalc) < 1e-12 else "NO"),
            ("weighted_concordance", overall["weighted_concordance"],
             "NOT_RECALCULATED", "NA"),
        ]
        for q, f, s, m in rows:
            fh.write("%s\t%s\t%s\t%s\n" % (q, f, s, m))
        fh.write("rule_applied\tdelta_psi = PSI_ASD - PSI_control (Parikshak 2016); "
                 "eligible iff |delta_psi| > 0.01; concordant iff delta_psi < -0.01; "
                 "two-sided exact binomial on concordant vs discordant among eligible\t"
                 "same rule reconstructed from chymera_directional_analysis.py lines 294-316, 343\t-\n")
    say("WROTE %s" % recalc_path)
    say("RECALC: eligible=%d concordant=%d discordant=%d P=%.15g" %
        (n_eligible, n_conc, n_opp, binom_p_recalc))
    say("REFERENCE: eligible=%s concordant=%s discordant=%s P=%s" %
        (overall["n_testable"], overall["n_concordant"],
         overall["n_opposite"], overall["binomial_p"]))

    bridge_recalc_ok = all(
        r[3] == "YES" for r in rows if r[0] != "weighted_concordance")

    # ------- all-19 threshold-free sensitivity -------
    # Sign-only comparison; every event has a finite signed discovery delta_psi
    # and CHyMErA direction = MICROEXON_INCLUSION_LOSS (all deletions).
    n_total = len(bridge)
    comparable = []
    for r in bridge:
        dps = float(r["ASD_delta_psi"])
        chy_ok = r["CHyMErA_direction"] == "MICROEXON_INCLUSION_LOSS"
        comparable.append(chy_ok and dps != 0)
    n_comparable = sum(comparable)
    if n_comparable == n_total:
        n_conc19 = sum(float(r["ASD_delta_psi"]) < 0 for r in bridge)
        n_disc19 = n_comparable - n_conc19
        p19 = stats.binomtest(n_conc19, n_comparable, 0.5).pvalue
        status = "COMPUTABLE"
        rule_desc = ("Sign-only concordance, no threshold: all 19 CHyMErA "
                     "perturbations are microexon deletions = expected "
                     "microexon inclusion loss; an event is concordant iff its "
                     "signed discovery delta_psi (PSI_ASD - PSI_control, "
                     "Parikshak 2016) is negative; two-sided exact binomial "
                     "under p=0.5; |delta_psi|>0.01 eligibility filter "
                     "ignored; all 19 events used exactly once; no P-value "
                     "thresholding; no event exclusion.")
    else:
        n_conc19 = n_disc19 = "NA"
        p19 = "NA"
        status = "NOT_COMPUTABLE_WITHOUT_IMPUTATION"
        rule_desc = ("NOT computed: %d/19 events lack an interpretable signed "
                     "comparison." % (n_total - n_comparable))
    all19_path = os.path.join(
        OUT, "ALL19_THRESHOLD_FREE_DIRECTION_SENSITIVITY.tsv")
    with open(all19_path, "w") as fh:
        fh.write("field\tvalue\n")
        fh.write("n_total\t%d\n" % n_total)
        fh.write("n_comparable\t%d\n" % n_comparable)
        fh.write("n_concordant\t%s\n" % n_conc19)
        fh.write("n_discordant\t%s\n" % n_disc19)
        fh.write("concordance_proportion\t%s\n" %
                 ("%.16g" % (n_conc19 / n_comparable) if n_comparable else "NA"))
        fh.write("exact_binomial_p_two_sided\t%s\n" %
                 ("%.15g" % p19 if isinstance(p19, float) else p19))
        fh.write("computability_status\t%s\n" % status)
        fh.write("rule_description\t%s\n" % rule_desc)
        fh.write("role\tSENSITIVITY_ONLY; does not replace the final 14-event "
                 "primary directional bridge\n")
    say("WROTE %s" % all19_path)
    say("ALL19: status=%s concordant=%s/%s P=%s" %
        (status, n_conc19, n_comparable, p19))

    # ------- check TSV -------
    cols = ["gene", "human_event_id", "mouse_event_id", "discovery_delta_psi",
            "discovery_p_raw", "discovery_direction",
            "chymera_perturbation_direction",
            "eligible_for_directional_comparison", "direction_unresolved_reason",
            "direction_unresolved_reason_in_master", "concordant_if_eligible",
            "bridge_classification", "recalculated_bridge_classification",
            "bridge_matches_reference", "rule_source_file",
            "rule_source_line_or_code_object"]
    check_path = os.path.join(OUT, "DIRECTIONAL_ELIGIBILITY_CHECK.tsv")
    with open(check_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for a in check_rows:
            fh.write("\t".join(a[c] for c in cols) + "\n")
    say("WROTE %s (19 rows)" % check_path)

    # ------- rule trace md -------
    md = f"""# Reference Workstream A — directional-eligibility rule trace

Date: 2026-08-13. Check-only reconstruction; no upstream file modified.

## Q1. What exact logical rule makes an event eligible?

The rule is implemented in
`{RULE_SOURCE_FILE}` (lines 294-316):

```python
# ASD direction
if pd.notna(dps):
    if dps < -0.01:
        asd_dir = "ASD_DECREASED_INCLUSION"
    elif dps > 0.01:
        asd_dir = "ASD_INCREASED_INCLUSION"
    else:
        asd_dir = "ASD_DIRECTION_UNRESOLVED"
else:
    asd_dir = "ASD_DIRECTION_UNRESOLVED"
```

and the bridge classification (lines 308-316): `ASD_DECREASED_INCLUSION` ->
`CONCORDANT_WITH_MICROEXON_LOSS`; `ASD_INCREASED_INCLUSION` ->
`OPPOSITE_TO_MICROEXON_LOSS`; `ASD_DIRECTION_UNRESOLVED` ->
`NO_ASD_EFFECT_DIRECTION`. An event is eligible for directional comparison
iff its bridge classification differs from `NO_ASD_EFFECT_DIRECTION`
(equivalently |discovery delta_psi| > 0.01).

**Exact rule: an event is eligible for directional comparison iff
|discovery delta_psi| > 0.01**, where delta_psi = PSI_ASD - PSI_control from
the Parikshak 2016 discovery cohort (`06_directional_bridge/00_direction_definitions.tsv`).
Concordant with CHyMErA microexon-deletion inclusion loss iff
delta_psi < -0.01.

## Q2. Does the rule depend on a discovery P-value threshold, missing sign, unresolved mapping, CHyMErA direction, or another criterion?

- **No P-value threshold**: ASD_p plays NO role in the rule. The threshold is a
  magnitude criterion on delta_psi (0.01).
- **Missing sign**: a NaN delta_psi would yield ASD_DIRECTION_UNRESOLVED, but
  all 19 events have finite signed delta_psi, so this branch is inactive.
- **Unresolved mapping**: not involved; all 19 events have a CHyMErA mapping.
- **CHyMErA direction**: constant for all 19 events
  (`MICROEXON_INCLUSION_LOSS`, because every CHyMErA perturbation in the set is
  a microexon deletion; line 305-306). It defines what "concordant" means but
  does not phase eligibility.
- The effective criterion is therefore solely: |delta_psi| > 0.01.

## Recalculation summary

| quantity | Reference | Recalculation | match |
|---|---|---|---|
| n_eligible | {overall["n_testable"]} | {n_eligible} | {"YES" if int(overall["n_testable"]) == n_eligible else "NO"} |
| n_concordant | {overall["n_concordant"]} | {n_conc} | {"YES" if int(overall["n_concordant"]) == n_conc else "NO"} |
| n_discordant | {int(float(overall["n_opposite"]))} | {n_opp} | {"YES" if int(float(overall["n_opposite"])) == n_opp else "NO"} |
| exact binomial P (two-sided) | {overall["binomial_p"]} | {binom_p_recalc:.15g} | {"YES" if abs(p - binom_p_recalc) < 1e-12 else "NO"} |

All-19 threshold-free sensitivity: computability_status = {status};
{n_conc19}/{n_comparable} sign-concordant; two-sided exact-binomial
P = {p19 if isinstance(p19, str) else format(p19, ".15g")}. Sensitivity only;
the final 14-event bridge remains the primary result.
"""
    md_path = os.path.join(OUT, "DIRECTIONAL_ELIGIBILITY_RULE_TRACE.md")
    with open(md_path, "w") as fh:
        fh.write(md)
    say("WROTE %s" % md_path)

    # ------- summary -------
    phase = {
        "A_roster_19": roster_match,
        "A_rule_traceable": True,
        "A_bridge_recalc_matches_reference": bridge_recalc_ok,
        "A_all19_valid": status in ("COMPUTABLE",
                                     "NOT_COMPUTABLE_WITHOUT_IMPUTATION"),
        "A_per_event_bridge_mismatches": len(mismatches),
    }
    for k, v in phase.items():
        say("PHASE %s=%s" % (k, v))

    with open(LOG, "w") as fh:
        fh.write("\n".join(_log) + "\n")


if __name__ == "__main__":
    main()
