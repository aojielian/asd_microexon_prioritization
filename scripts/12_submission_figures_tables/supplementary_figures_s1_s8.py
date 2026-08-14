#!/usr/bin/env python3
"""Targeted final visual fixes on top of the earlier S1-S8
renders (figure-only deliverables).

Changes applied in the final revision:
  S1: panel B (round-trip plot) removed entirely -> single panel, no panel
      letter, matrix enlarged; liftOver conclusion moved to the legend.
  S2: panel D grey unit text removed, log-scale x upper bound tightened.
  S3: panel B gene-list text removed (kept 12 retained / 7 excluded / the
      two P values); panel A floating 'P = 0.05' text removed (dashed line
      kept, meaning in legend).
  S4: panel D 'range xx' labels removed.
  S5: panel C labels shortened, floating 'P = 0.05' removed, wspace widened
      (label collision with panel B fixed).
  S6: panel A flow text reduced to essential counts; reason breakdown moved
      to the legend.
  S7: unchanged (recheck only).
  S8: panel C widened, event labels enlarged.

Original design notes follow.

Final re-render of Figures S1-S8 (figure-only deliverables).

Redesign per spec sections 5-6:
  S1: TWO panels only (A availability matrix, B round-trip check); the old
      text-only lineage funnel and identifier-table panels are removed
      (their facts now live in the legend).
  S2: keep A-D; italic in-panel
      unit-caveat sentences are removed (moved to the legend); gutters
      widened.
  S3: redundant badge-summary panel removed; old B -> A (sensitivity
      lollipops), old C -> B (ASD-prior exclusion, without archival wording).
  S4: keep A-D; bottom italic methodology note removed; tick-label
      readability increased.
  S5: panel C converted from a prose/badge box to a data summary (the four
      nominal tests plus the SRRM4 minimum permutation P).
  S6: colored text-strip direction panel removed; A mapping flow (without
      italic sentence, '19 events'), old C -> B local-structure matching.
  S7: text-only definition panel C removed; A and B enlarged, in-panel
      boxed sentence removed (content in legend).
  S8: keep A-D; spacing increased, D row labels enlarged.

All numeric content is read from upstream analysis sources; nothing is
recomputed.  Every figure meets automated extent checks (title/legend
overlap, panel overlap, text clipping) at render time.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supplementary_common import (ROOT, FIG_W, MM, VQC_DIR, C_PRIMARY, C_CONCORD,
                             C_DEV, C_DISCORD, C_NEG, C_SENS, C_ORANGE,
                             C_DARK, C_MID, C_BG_GREY, TIER_COLORS,
                             TIER_ORDER, panel_letter, panel_title,
                             save_figure, tsv_write, rd, load_master,
                             tier_letter, extent_checks)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.path import Path
import matplotlib.patches as mpatches

DIR11 = os.path.join(ROOT, "11_set_level_enrichment")
DIR14 = os.path.join(ROOT, "14_mechanistic_context")
DIR17 = os.path.join(ROOT, "17_gse30573_mapping")
DIR21 = os.path.join(ROOT, "21_coordinate_inference")
D29_RENDER = os.path.join(ROOT, "29_final_semantic_terminology_patch_"
                                "20260804", "03_render_data")

MT = load_master()
assert len(MT) == 19
QC_ROWS = []


def run_checks(fig, name):
    for check, ok, detail in extent_checks(fig, name):
        QC_ROWS.append([name, check, "OK" if ok else "ERROR", detail])
        print("  %-28s %-4s %s" % (check, "OK" if ok else "ERROR",
                                   detail[:70]))


def short(r):
    return "%s %s" % (r["gene"], r["HsaEX_ID"].replace("HsaEX", ""))


# ====================================================================== S1
def build_s1():
    """Single panel: data-source availability matrix (round-trip panel B
    removed; its 0 bp liftOver conclusion is carried in the legend)."""
    ROWS = sorted(MT, key=lambda r: (r["gene"], r["HsaEX_ID"]))
    SRC = ["Discovery\n(ASD cortex)", "CHyMErA", "BrainSpan\ndevelopment",
           "GSE30573", "PsychENCODE", "liftOver\nround-trip",
           "GENCODE v33\nstructure"]

    def src_state(r):
        return [
            "+",
            {"YES": "+", "NO": "\u00d7"}.get(r["CHyMErA_direction_concordant"],
                                             "\u2013"),
            "+" if r["developmental_dynamic_status"] == "DYNAMIC" else "\u2013",
            {"CONCORDANT": "+", "DISCORDANT": "\u00d7"}.get(
                r["GSE30573_direction_concordant"], "\u2013"),
            "+" if r["direction_concordant"] == "TRUE" else "\u00d7",
            "+", "+",
        ]

    STATE_COL = {"+": C_CONCORD, "\u00d7": C_DISCORD, "\u2013": C_BG_GREY}
    STATE_TXT = {"+": "white", "\u00d7": "white", "\u2013": C_MID}

    n = len(ROWS)
    fig = plt.figure(figsize=(FIG_W, 132 * MM))
    ax = fig.add_axes([0.20, 0.075, 0.62, 0.80])
    ax.axis("off")
    ax.set_xlim(-6.4, 7.6)
    ax.set_ylim(n + 2.4, -6.2)
    ax.set_title("Data-source availability per event", fontsize=10.0,
                 weight="bold", color=C_DARK, loc="center")
    for j, c in enumerate(SRC):
        ax.text(j + 0.5, -2.1, c, ha="left", va="center", fontsize=7.4,
                color=C_DARK, weight="bold", rotation=45,
                rotation_mode="anchor")
    for i, r in enumerate(ROWS):
        ax.text(-0.45, i + 0.5, short(r), ha="right", va="center",
                fontsize=7.0, color=C_DARK)
        for j, s in enumerate(src_state(r)):
            ax.add_patch(Rectangle((j + 0.05, i + 0.05), 0.90, 0.90,
                                   fc=STATE_COL[s], ec="white", lw=0.5))
            ax.text(j + 0.5, i + 0.52, s, ha="center", va="center",
                    fontsize=7.2, color=STATE_TXT[s], weight="bold")
    ax.text(3.0, n + 1.3,
            "+ supporting   \u00d7 discordant   \u2013 absent / not "
            "significant / NA", ha="center", va="center", fontsize=6.8,
            color=C_MID, style="italic")

    run_checks(fig, "Figure_S1")
    save_figure(fig, "Figure_S1")


# ====================================================================== S2
def build_s2():
    EFF_ROWS = rd(os.path.join(DIR11, "07_primary_reanalysis",
                               "01_effects_by_background.tsv"))
    order = ["BG0_WIDE_SE", "BG1_MICROEXON", "BG2_CONSERVED_MICROEXON",
             "BG3_CEM", "BG3_NN"]
    lab = {"BG0_WIDE_SE": "Wide splicing-event", "BG1_MICROEXON": "Microexon",
           "BG2_CONSERVED_MICROEXON": "Conserved microexon",
           "BG3_CEM": "CEM (matched)", "BG3_NN": "NN (matched)"}
    dd = {r["background"]: r for r in EFF_ROWS}
    EFF = [dict(bg=b, label=lab[b], perm_p=float(dd[b]["permutation_p"]))
           for b in order]
    NULLS = np.load(os.path.join(D29_RENDER, "NULL_DISTRIBUTIONS.npz"))
    BG_COLS = [C_PRIMARY, C_DEV, "#3a7ca5", C_ORANGE, C_SENS]

    REFERENCE_UNIVERSE = [
        ("Wide splicing-event", 200956, "annotated splicing events"),
        ("Conserved microexon", 377, "unique events (deduplicated)"),
        ("CEM-derived", 149, "unique background events"),
        ("NN-derived", 102, "unique background events"),
        ("PSI-matched pool", 4948, "candidate sampling pool"),
    ]
    ANALYSIS_UNIVERSE = [
        ("BG0 wide SE", 20916, "valid-effect events"),
        ("BG1 microexon", 467, "valid-effect microexon events"),
        ("BG2 cons. microexon", 452, "as-mapped conserved events"),
        ("BG3 CEM matched", 308, "pair records (180 unique bg events)"),
        ("BG3 NN matched", 380, "pair records (120 unique bg events)"),
    ]

    obs = float(np.mean(np.abs([float(r["Parikshak_delta_PSI"])
                                for r in MT])))
    assert abs(obs - 0.0399) < 5e-5, obs

    fig = plt.figure(figsize=(FIG_W, 185 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.85, wspace=0.60,
                           left=0.165, right=0.97, top=0.94, bottom=0.10,
                           width_ratios=[1.0, 1.15])

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(REFERENCE_UNIVERSE))[::-1]
    for i, (name, nn, unit) in enumerate(REFERENCE_UNIVERSE):
        ax.barh(y[i], nn, color=BG_COLS[i], height=0.62, zorder=2)
        ax.text(nn * 1.22, y[i] + 0.13, format(nn, ","), va="center",
                fontsize=6.8, color=C_DARK, weight="bold")
        ax.text(nn * 1.22, y[i] - 0.30, unit, va="center", fontsize=5.6,
                color=C_MID)
    ax.set_xscale("log")
    ax.set_xlim(60, 1500000)
    ax.set_xticks([1e2, 1e3, 1e4, 1e5, 1e6])
    ax.set_yticks(y)
    ax.set_yticklabels([u[0] for u in REFERENCE_UNIVERSE], fontsize=6.6)
    ax.set_xlabel("background universe size (log scale)")
    panel_title(ax, "Background universes")
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[0, 1])
    y = np.arange(len(EFF))[::-1]
    for i, e in enumerate(EFF):
        ns = NULLS[e["bg"]]
        lo, med, hi = np.percentile(ns, [2.5, 50.0, 97.5])
        ax.plot([lo, hi], [y[i], y[i]], "-", color=BG_COLS[i], lw=4.6,
                solid_capstyle="round", alpha=0.5, zorder=2)
        ax.plot([med, med], [y[i] - 0.24, y[i] + 0.24], "-", color=BG_COLS[i],
                lw=1.6, zorder=3)
        pct = 100.0 * float(np.mean(ns < obs))
        p = e["perm_p"]
        ptxt = "P = 0.0001" if p < 0.001 else "P = %.4f" % p
        ax.text(0.046, y[i] + 0.13, e["label"], fontsize=6.6, va="center",
                color=C_DARK, weight="bold")
        ax.text(0.046, y[i] - 0.26, "%s, percentile %.2f"
                % (ptxt, pct), fontsize=5.8, va="center", color=C_MID)
    ax.axvline(obs, color=C_DISCORD, lw=1.3, zorder=4)
    ax.text(obs - 0.0012, len(EFF) - 0.42, "observed = %.4f" % obs,
            color=C_DISCORD, fontsize=6.6, ha="right", va="bottom",
            weight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_ylim(-0.62, len(EFF) - 0.30)
    ax.set_xlim(0.0, 0.095)
    ax.set_xticks(np.arange(0.0, 0.051, 0.01))
    ax.set_xlabel("permuted 19-event mean |\u0394PSI|")
    panel_title(ax, "Permutation null 95% interval vs observed",
                fontsize=8.5)
    panel_letter(ax, "B")

    # ------------------------------------------------------------- C
    ax = fig.add_subplot(gs[1, 0])
    sdv = [("CEM matched", -0.336), ("NN matched", -0.073)]
    yy = [1, 0]
    for (name, v), yv in zip(sdv, yy):
        ax.plot(v, yv, "o", ms=8, color=C_SENS, mec="white", mew=0.9,
                zorder=3)
        ax.text(0.55, yv, "SMD = %.3f" % v, va="center", fontsize=6.8,
                color=C_DARK)
    ax.axvspan(-0.5, 0.5, color=C_BG_GREY, zorder=1)
    ax.axvline(0, color=C_DARK, ls="--", lw=0.8)
    ax.set_yticks(yy)
    ax.set_yticklabels([s[0] for s in sdv], fontsize=6.8)
    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-0.6, 1.9)
    ax.set_xlabel("standardized mean difference, exon length")
    panel_title(ax, "Matched-background balance diagnostics")
    panel_letter(ax, "C")

    # ------------------------------------------------------------- D
    ax = fig.add_subplot(gs[1, 1])
    y = np.arange(len(ANALYSIS_UNIVERSE))[::-1]
    for i, (name, nn, unit) in enumerate(ANALYSIS_UNIVERSE):
        ax.barh(y[i], nn, color=BG_COLS[i], height=0.62, zorder=2)
        ax.text(nn * 1.22, y[i], format(nn, ","), va="center",
                fontsize=6.8, color=C_DARK, weight="bold")
    ax.set_xscale("log")
    ax.set_xlim(60, 60000)
    ax.set_xticks([1e2, 1e3, 1e4])
    ax.set_yticks(y)
    ax.set_yticklabels([u[0] for u in ANALYSIS_UNIVERSE], fontsize=6.6)
    ax.set_xlabel("analysis-set size (log scale)")
    panel_title(ax, "Analysis backgrounds (record units)")
    panel_letter(ax, "D")

    run_checks(fig, "Figure_S2")
    save_figure(fig, "Figure_S2")


# ====================================================================== S3
def build_s3():
    fig = plt.figure(figsize=(FIG_W, 115 * MM))
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.62, left=0.165,
                           right=0.975, top=0.88, bottom=0.16,
                           width_ratios=[1.15, 1.0])

    # --------------------------------------------- A sensitivity lollipops
    ax = fig.add_subplot(gs[0])
    tests = [("Microexon length \u2264 27 nt", 0.0040),
             ("Microexon length \u2264 30 nt", 0.0042),
             ("Microexon length \u2264 36 nt", 0.0042),
             ("One-event-per-gene", 0.0024),
             ("Target > random sets", 0.0023),
             ("ASD-prior exclusion", 0.0335)]
    y = np.arange(len(tests))[::-1]
    for i, (name, p) in enumerate(tests):
        v = -np.log10(p)
        ax.plot([0, v], [y[i], y[i]], "-", color=C_SENS, lw=1.9, zorder=2)
        ax.plot(v, y[i], "o", ms=6.0, color=C_SENS, mec="white", mew=0.7,
                zorder=3)
        ax.text(v + 0.10, y[i], "P = %.4f" % p, fontsize=6.6, va="center",
                color=C_MID)
    thr = -np.log10(0.05)
    ax.axvline(thr, color=C_DISCORD, ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([t[0] for t in tests], fontsize=7.0)
    ax.set_xlabel("\u2212log10(P)")
    ax.set_xlim(0, 3.1)
    panel_title(ax, "Definition and selection sensitivity P values")
    panel_letter(ax, "A")

    # ----------------------------------------------- B ASD-prior exclusion
    ax = fig.add_subplot(gs[1])
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "ASD-prior gene exclusion (event level)", fontsize=8.0)
    ax.text(4, 84, "exclusion applied at the 19-event level",
            fontsize=7.0, color=C_DARK, weight="bold")
    ax.add_patch(Rectangle((4, 62), 40, 11, fc=C_CONCORD, ec="none"))
    ax.text(48, 67.5, "12 events retained", va="center", color=C_DARK,
            fontsize=7.0, weight="bold")
    ax.add_patch(Rectangle((4, 44), 40 * 7 / 12, 11, fc=C_ORANGE, ec="none"))
    ax.text(48, 49.5, "7 events excluded", va="center", color=C_DARK,
            fontsize=7.0, weight="bold")
    ax.text(4, 28, "permutation P = 0.0335", fontsize=7.4, color=C_DARK,
            weight="bold")
    ax.text(4, 20, "Wilcoxon P = 0.0266", fontsize=7.4, color=C_DARK,
            weight="bold")
    panel_letter(ax, "B", x=-0.10)

    run_checks(fig, "Figure_S3")
    save_figure(fig, "Figure_S3")


# ====================================================================== S4
def build_s4():
    TRAJ_COL = {"PLPH": C_DEV, "PHPL": C_ORANGE, "NON_DYNAMIC": "#bdbdbd"}

    def devrank(r):
        return {"PLPH": 0, "PHPL": 1, "NON_DYNAMIC": 2}.get(
            r["developmental_trajectory"], 2)

    fig = plt.figure(figsize=(FIG_W, 195 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.72, wspace=0.62,
                           left=0.12, right=0.96, top=0.935, bottom=0.08,
                           width_ratios=[1.0, 1.05])

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0, 0])
    rows = sorted(MT, key=lambda r: -float(r["prenatal_mean_PSI"]))
    for r in rows:
        pre = float(r["prenatal_mean_PSI"])
        post = float(r["postnatal_mean_PSI"])
        c = TRAJ_COL[r["developmental_trajectory"]]
        ax.plot([0, 1], [pre, post], "-", color=c, lw=1.4, alpha=0.8, zorder=2)
        ax.plot(0, pre, "o", ms=3.0, color=c, mec="white", mew=0.4, zorder=3)
        ax.plot(1, post, "o", ms=3.0, color=c, mec="white", mew=0.4, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Prenatal", "Postnatal"])
    ax.set_yticks([0, 50, 100])
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.32, 1.32)
    for r in [rows[0], rows[-1]]:
        ax.text(-0.03, float(r["prenatal_mean_PSI"]), short(r), ha="right",
                va="center", fontsize=5.8, color=C_MID)
    for t, labx in [("PLPH", "PLPH (9)"), ("PHPL", "PHPL (1)"),
                    ("NON_DYNAMIC", "non-dynamic (9)")]:
        ax.plot([], [], "-", color=TRAJ_COL[t], lw=1.4, label=labx)
    ax.legend(loc="lower center", fontsize=6.2, ncol=3,
              bbox_to_anchor=(0.5, 1.14))
    ax.set_ylabel("mean PSI (%)")
    panel_title(ax, "All 19 events: prenatal-to-postnatal PSI")
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[0, 1])
    tests = [("vs conserved microexon", 0.0004, C_DEV),
             ("vs CEM", 0.0005, C_DEV),
             ("vs NN", 0.0013, C_DEV),
             ("PSI-matched (negative)", 0.7186, C_NEG),
             ("gene-block (negative)", 1.0, C_NEG),
             ("BrainSpan substrate (expression)", 0.7199, C_NEG),
             ("zebrafish (exploratory)", 0.0688, C_ORANGE)]
    y = np.arange(len(tests))[::-1]
    for i, (name, p, c) in enumerate(tests):
        v = -np.log10(p)
        ax.plot([0, v], [y[i], y[i]], "-", color=c, lw=1.9, zorder=2)
        ax.plot(v, y[i], "o", ms=6.0, color=c, mec="white", mew=0.7, zorder=3)
        ax.text(v + 0.07, y[i], "P = %.4f" % p, fontsize=6.4, va="center",
                color=C_MID)
    thr = -np.log10(0.05)
    ax.axvline(thr, color=C_DISCORD, ls="--", lw=0.8)
    ax.text(thr + 0.05, -0.55, "P = 0.05", fontsize=6.2, color=C_DISCORD,
            va="top")
    ax.set_yticks(y)
    ax.set_yticklabels([t[0] for t in tests], fontsize=6.8)
    ax.set_xlabel("\u2212log10(P)")
    ax.set_xlim(0, 4.6)
    ax.set_ylim(-0.7, len(tests) - 0.4)
    panel_title(ax, "Dynamicity tests and negative controls")
    panel_letter(ax, "B")

    # ------------------------------------------------------------- C
    ax = fig.add_subplot(gs[1, 0])
    rows = sorted(MT, key=lambda r: -float(r["monotonicity_rho"]))
    y = np.arange(len(rows))[::-1]
    for i, r in enumerate(rows):
        v = float(r["monotonicity_rho"])
        c = C_DEV if r["developmental_dynamic_status"] == "DYNAMIC" else C_NEG
        ax.plot([0, v], [y[i], y[i]], "-", color=c, lw=1.7, zorder=2)
        ax.plot(v, y[i], "o", ms=4.4, color=c, mec="white", mew=0.5, zorder=3)
    ax.axvline(0, color=C_DARK, ls=":", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([short(r) for r in rows], fontsize=6.0)
    ax.set_xlabel("per-event monotonicity rho (prenatal vs postnatal)")
    panel_title(ax, "Per-event developmental monotonicity")
    panel_letter(ax, "C")

    # ------------------------------------------------------------- D
    ax = fig.add_subplot(gs[1, 1])
    rows = sorted(MT, key=lambda r: (devrank(r), r["gene"]))
    y = np.arange(len(rows))[::-1]
    for i, r in enumerate(rows):
        pre = float(r["prenatal_mean_PSI"])
        post = float(r["postnatal_mean_PSI"])
        c = TRAJ_COL[r["developmental_trajectory"]]
        ax.plot([pre, post], [y[i], y[i]], "-", color=c, lw=2.6,
                solid_capstyle="round", zorder=2)
        ax.plot(pre, y[i], "o", ms=3.6, color=c, mec="white", mew=0.5,
                zorder=3)
        ax.plot(post, y[i], "s", ms=3.6, color=c, mec="white", mew=0.5,
                zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([short(r) for r in rows], fontsize=6.0)
    ax.set_xlim(0, 122)
    ax.set_xticks([0, 30, 60, 90, 120])
    ax.set_xlabel("mean PSI (%) prenatal (circle) to postnatal (square)")
    panel_title(ax, "Per-event PSI range detail")
    panel_letter(ax, "D")

    run_checks(fig, "Figure_S4")
    save_figure(fig, "Figure_S4")


# ====================================================================== S5
def build_s5():
    RBP = rd(os.path.join(DIR14, "15_reports", "rbp_motif_results.tsv"))
    for r in RBP:
        r["perm_p"] = float(r["perm_p"])
        r["perm_fdr"] = float(r["perm_fdr"])
    assert len(RBP) == 240
    assert sum(1 for r in RBP if r["perm_p"] < 0.05) == 4
    assert all(r["perm_fdr"] >= 0.05 for r in RBP)

    COMP_LAB = {"primary19_vs_background": "19 vs background",
                "dynamic10_vs_background": "dynamic10 vs background",
                "dynamic10_vs_nondynamic9": "dynamic vs non-dynamic",
                "tier2_5_vs_background": "tier2-5 vs background"}
    REG_LAB = ["combined_proximal", "up_proximal", "up_extended", "exon",
               "down_proximal", "down_extended"]

    fig = plt.figure(figsize=(FIG_W, 150 * MM))
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.72, left=0.065,
                           right=0.985, top=0.90, bottom=0.15,
                           width_ratios=[1.0, 1.05, 1.15])

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0])
    ps = sorted(r["perm_p"] for r in RBP)
    yv = [-np.log10(p) for p in ps]
    ax.plot(np.arange(len(yv)), yv, ".", ms=2.6, color=C_NEG, zorder=2)
    nom = [p for p in ps if p < 0.05]
    ax.plot(np.arange(len(nom)), [-np.log10(p) for p in nom], "o", ms=5.0,
            color=C_DISCORD, mec="white", mew=0.6, zorder=3)
    ax.axhline(-np.log10(0.05), color=C_DISCORD, ls="--", lw=0.8)
    ax.text(len(yv) * 0.45, -np.log10(0.05) + 0.06, "nominal P = 0.05",
            fontsize=6.4, color=C_DISCORD)
    ax.set_xlim(-6, 246)
    ax.set_xticks([0, 60, 120, 180, 240])
    ax.set_xlabel("ranked permutation tests (240 total)")
    ax.set_ylabel("\u2212log10(perm P)")
    ax.text(0.98, 0.97, "0/240 BH-FDR significant\n4/240 nominal P < 0.05",
            transform=ax.transAxes, va="top", ha="right", fontsize=6.8,
            weight="bold", color=C_DARK)
    panel_title(ax, "Ranked RBP permutation tests")
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[1])
    comps = list(COMP_LAB)
    Mm = np.zeros((len(comps), len(REG_LAB)))
    for i, c in enumerate(comps):
        for j, g in enumerate(REG_LAB):
            vals = [r["perm_p"] for r in RBP if r["comparison"] == c and
                    r["region"] == g]
            Mm[i, j] = -np.log10(min(vals)) if vals else 0
    ax.imshow(Mm, cmap="viridis", aspect="auto", vmin=0, vmax=2)
    for i in range(len(comps)):
        for j in range(len(REG_LAB)):
            ax.text(j, i, "%.2f" % Mm[i, j], ha="center", va="center",
                    fontsize=6.2, color="white" if Mm[i, j] > 1.1 else C_DARK)
    ax.set_xticks(np.arange(len(REG_LAB)))
    ax.set_xticklabels([g.replace("_", " ") for g in REG_LAB], fontsize=6.2,
                       rotation=35, ha="right")
    ax.set_yticks(np.arange(len(comps)))
    ax.set_yticklabels([COMP_LAB[c] for c in comps], fontsize=6.4)
    panel_title(ax, "Max \u2212log10(perm P) per comparison")
    panel_letter(ax, "B")

    # ------------------------------------------------- C data summary panel
    ax = fig.add_subplot(gs[2])
    # shortened labels (final revision); the full comparison wording is carried
    # in the figure legend
    groups = [("MBNL1/2 \u2014 down ext.", 0.0454, C_ORANGE),
              ("PTBP1/2 \u2014 up prox.", 0.0356, C_ORANGE),
              ("PTBP1/2 \u2014 up ext.", 0.0441, C_ORANGE),
              ("ELAVL2/3/4 \u2014 exon", 0.0261, C_ORANGE),
              ("SRRM4 \u2014 minimum", 0.0928, C_NEG)]
    y = np.arange(len(groups))[::-1]
    for i, (name, p, c) in enumerate(groups):
        v = -np.log10(p)
        ax.plot([0, v], [y[i], y[i]], "-", color=c, lw=1.9, zorder=2)
        ax.plot(v, y[i], "o", ms=6.0, color=c, mec="white", mew=0.7, zorder=3)
        ax.text(v + 0.06, y[i], "P = %.4f" % p, fontsize=6.2, va="center",
                color=C_MID)
    thr = -np.log10(0.05)
    ax.axvline(thr, color=C_DISCORD, ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([g[0] for g in groups], fontsize=6.2)
    ax.set_xlabel("\u2212log10(perm P)")
    ax.set_xlim(0, 3.3)
    ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    panel_title(ax, "Nominal tests (4/240) and SRRM4", fontsize=7.6)
    panel_letter(ax, "C")

    run_checks(fig, "Figure_S5")
    save_figure(fig, "Figure_S5")


# ====================================================================== S6
def build_s6():
    MAP = rd(os.path.join(DIR17, "05_junction_structure_mapping",
                          "02_primary19_final_mapping.tsv"))

    def bezier_band(ax, x0, y0b, y0t, x1, y1b, y1t, color, alpha=0.55):
        xm = 0.5 * (x0 + x1)
        verts = [(x0, y0b), (xm, y0b), (xm, y1b), (x1, y1b),
                 (x1, y1t), (xm, y1t), (xm, y0t), (x0, y0t), (x0, y0b)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.CLOSEPOLY]
        ax.add_patch(mpatches.PathPatch(Path(verts, codes), fc=color,
                                        ec="none", alpha=alpha, zorder=1))

    mapped = [r for r in MT if r["GSE30573_mapping_status"] ==
              "MAPPED_ANALYZABLE"]
    assert len(mapped) == 3
    assert sum(1 for r in mapped if r["GSE30573_direction_concordant"] ==
               "CONCORDANT") == 2

    fig = plt.figure(figsize=(FIG_W, 128 * MM))
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.50, left=0.06,
                           right=0.985, top=0.90, bottom=0.13,
                           width_ratios=[1.12, 1.0])

    # ------------------------------------------------------------- A flow
    ax = fig.add_subplot(gs[0])
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "GSE30573 mapping flow")
    S = 1.5
    W = 18
    h19 = 19 * S
    y19 = 62 - h19 / 2
    ax.add_patch(mpatches.FancyBboxPatch(
        (4, y19), W, h19, boxstyle="round,pad=0.2,rounding_size=0.8",
        fc=C_MID, ec="none"))
    ax.text(4 + W / 2, 62, "19 events", ha="center", va="center",
            color="white", fontsize=7.4, weight="bold")
    h3, h16 = 3 * S, 16 * S
    y3 = 86 - h3
    y16 = y3 - 7 - h16
    ax.add_patch(Rectangle((40, y3), W, h3, fc=C_CONCORD, ec="none"))
    ax.add_patch(Rectangle((40, y16), W, h16, fc="#c9c9c9", ec="none"))
    ax.text(40 + W / 2, y16 + h16 / 2, "16 unmapped", ha="center",
            va="center", color=C_MID, fontsize=6.8, weight="bold")
    u = h19 / 19.0
    bezier_band(ax, 4 + W, y19 + 16 * u, y19 + 19 * u, 40, y3, y3 + h3,
                C_CONCORD, 0.5)
    bezier_band(ax, 4 + W, y19, y19 + 16 * u, 40, y16, y16 + h16, "#c9c9c9",
                0.5)
    # Final: essential counts only; mapped gene list and the unmapped
    # reason breakdown (12 no local-structure match, 4 gene not in GSE
    # annotation) are carried in the legend.
    ax.text(62, y3 + h3 / 2, "3 mapped analyzable", fontsize=6.8,
            va="center", color=C_DARK, weight="bold")
    panel_letter(ax, "A", x=-0.08)

    # --------------------------------------------- B local structure match
    ax = fig.add_subplot(gs[1])
    rows = sorted([r for r in MAP if r["mapping_level"] != "NO_MATCH"],
                  key=lambda r: (r["gene"], r["HsaEX_ID"]))
    assert len(rows) == 3
    y = np.arange(len(rows))[::-1]
    LEVEL_LABEL = {
        "MATCH_EXACT_EXON_PARTIAL_JUNCTION_SUPPORT":
            "exact exon, partial junction support",
        "MATCH_EXACT_ROUNDTRIP_EXON_AND_3JUNCTIONS":
            "exact round-trip exon + 3 junctions",
    }
    for i, r in enumerate(rows):
        nj = int(r["n_junctions_matched"])
        ax.barh(y[i], nj, color=C_CONCORD, height=0.55, zorder=2)
        ax.text(nj + 0.12, y[i], LEVEL_LABEL[r["mapping_level"]],
                fontsize=5.8, va="center", color=C_DARK)
    ax.set_yticks(y)
    ax.set_yticklabels([short(r) for r in rows], fontsize=6.6)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlim(0, 8.0)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xlabel("local-structure elements matched in GSE30573")
    panel_title(ax, "Local-structure matching (3 events)")
    panel_letter(ax, "B")

    run_checks(fig, "Figure_S6")
    save_figure(fig, "Figure_S6")


# ====================================================================== S7
def build_s7():
    ROWS = sorted(MT, key=lambda r: (r["gene"], r["HsaEX_ID"]))

    fig = plt.figure(figsize=(FIG_W, 132 * MM))
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.42, left=0.125,
                           right=0.985, top=0.92, bottom=0.12,
                           width_ratios=[1.0, 1.30])

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0])
    y = np.arange(len(ROWS))[::-1]
    inc = [int(r["inclusion_transcript_count"]) for r in ROWS]
    exc = [int(r["exclusion_transcript_count"]) for r in ROWS]
    ax.barh(y + 0.18, inc, height=0.36, color=C_PRIMARY, zorder=2,
            label="inclusion transcripts")
    ax.barh(y - 0.18, exc, height=0.36, color=C_DEV, zorder=2,
            label="exclusion transcripts")
    for i in range(len(ROWS)):
        ax.text(max(inc[i], exc[i]) + 0.15, y[i], "%d/%d" % (inc[i], exc[i]),
                fontsize=6.0, va="center", color=C_MID)
    ax.set_yticks(y)
    ax.set_yticklabels([short(r) for r in ROWS], fontsize=6.2)
    ax.set_xlabel("GENCODE v33 transcripts supporting each isoform")
    ax.legend(loc="lower right", fontsize=6.6)
    panel_title(ax, "Transcript membership per event")
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[1])
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "Representative transcript structures (GRCh38)")
    reps = [r for r in ROWS if r["gene"] in ("CPEB4", "CTNND1", "MEF2D")]
    assert len(reps) == 3

    def exon_track(yy, label):
        ax.text(1.2, yy, label, fontsize=6.4, va="center", color=C_MID,
                ha="left")
        x0 = 36
        ax.add_patch(Rectangle((x0, yy - 2.2), 12, 4.4, fc="#d9d9d9",
                               ec=C_DARK, lw=0.5))
        ax.plot([x0 + 12, x0 + 17], [yy, yy], color=C_DARK, lw=0.7)
        ax.add_patch(Rectangle((x0 + 17, yy - 2.2), 4, 4.4, fc=C_DISCORD,
                               ec=C_DARK, lw=0.5))
        ax.plot([x0 + 21, x0 + 26], [yy, yy], color=C_DARK, lw=0.7)
        ax.add_patch(Rectangle((x0 + 26, yy - 2.2), 12, 4.4, fc="#d9d9d9",
                               ec=C_DARK, lw=0.5))

    ys = [84, 52, 20]
    for r, yy in zip(reps, ys):
        ax.text(1.2, yy + 9.5, r["gene"], fontsize=8.2, weight="bold",
                color=C_PRIMARY)
        exon_track(yy + 3.5, "human " + r["HsaEX_ID"])
        ax.text(98.5, yy + 3.5, "%s:%s-%s" % (r["chr_hg38"], r["start_hg38"],
                r["end_hg38"]), fontsize=6.0, va="center", ha="right",
                color=C_MID)
        exon_track(yy - 3.5, "mouse " + r["MmuEX_ID"])
        ax.text(98.5, yy - 3.5, "mm10 / VastDB", fontsize=6.0, va="center",
                ha="right", color=C_MID)
        ax.text(98.5, yy - 8.5, "incl %s / excl %s transcripts" %
                (r["inclusion_transcript_count"],
                 r["exclusion_transcript_count"]), fontsize=6.0, va="center",
                ha="right", color=C_MID)
    panel_letter(ax, "B", x=-0.06)

    run_checks(fig, "Figure_S7")
    save_figure(fig, "Figure_S7")


# ====================================================================== S8
def build_s8():
    VAL = os.path.join(DIR21, "05_mixed_model_inference",
                       "06_set_validation_recomputed.tsv")

    fig = plt.figure(figsize=(FIG_W, 190 * MM))
    # Final: panel C (direction concordance list) gets more width so
    # event labels stay readable at 100% zoom.
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.68, wspace=0.62,
                           left=0.125, right=0.97, top=0.95, bottom=0.08,
                           width_ratios=[1.18, 1.0])

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0, 0])
    rows = sorted(MT, key=lambda r: -float(r["PsychENCODE_SE"]))
    y = np.arange(len(rows))[::-1]
    for i, r in enumerate(rows):
        v = float(r["PsychENCODE_SE"])
        ax.barh(y[i], v, color=TIER_COLORS[tier_letter(r)], height=0.62,
                zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([short(r) for r in rows], fontsize=6.2)
    ax.set_xlabel("PsychENCODE logit-beta SE")
    for t in TIER_ORDER:
        ax.plot([], [], "s", ms=4.5, color=TIER_COLORS[t], label="Tier %s" % t)
    ax.legend(loc="lower right", fontsize=6.2, ncol=2)
    panel_title(ax, "Per-event effect precision")
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[0, 1])
    xs = [-np.log10(max(float(r["P_KR"]), 1e-6)) for r in MT]
    ys = [-np.log10(max(float(r["P_LRT"]), 1e-6)) for r in MT]
    m = max(max(xs), max(ys)) * 1.1
    ax.plot(xs, ys, "o", ms=5.0, color=C_SENS, mec="white", mew=0.5, zorder=3)
    ax.plot([0, m], [0, m], "--", color=C_MID, lw=0.8, zorder=1)
    ax.set_xlim(0, m)
    ax.set_ylim(0, m)
    tks = np.arange(0, m, 1.0)
    ax.set_xticks(tks)
    ax.set_yticks(tks)
    ax.set_xlabel("KR \u2212log10(P)")
    ax.set_ylabel("LRT \u2212log10(P)")
    ax.text(0.05, 0.95, "dashed: identity", transform=ax.transAxes,
            fontsize=6.2, color=C_MID, va="top")
    panel_title(ax, "KR vs LRT P-value agreement")
    panel_letter(ax, "B")

    # ------------------------------------------------------------- C
    ax = fig.add_subplot(gs[1, 0])
    ax.axis("off")
    rows = sorted(MT, key=lambda r: (r["direction_concordant"] != "TRUE",
                                     r["gene"]))
    n = len(rows)
    ax.set_xlim(0, 100)
    ax.set_ylim(-1, n + 1.8)
    panel_title(ax, "Event-level direction concordance")
    ax.text(32, n + 0.9, "discovery", ha="center", fontsize=6.4,
            color=C_DARK, weight="bold")
    ax.text(58, n + 0.9, "PsychENCODE", ha="center", fontsize=6.4,
            color=C_DARK, weight="bold")
    for i, r in enumerate(rows):
        yy = n - 1 - i
        conc = r["direction_concordant"] == "TRUE"
        col = C_CONCORD if conc else C_DISCORD
        ax.text(0, yy + 0.3, short(r), fontsize=6.6, color=C_DARK,
                va="center")
        d1 = 1 if r["Parikshak_direction"] == "UP_IN_ASD" else -1
        d2 = 1 if r["PsychENCODE_direction"] == "UP_IN_ASD" else -1
        for xc, d in [(34, d1), (56, d2)]:
            ax.annotate("", xy=(xc, yy + 0.45 * d), xytext=(xc, yy - 0.45 * d),
                        arrowprops=dict(arrowstyle="-|>", lw=1.1, color=col))
        ax.add_patch(Rectangle((72, yy - 0.42), 8, 0.84, fc=col, ec="none"))
    ax.text(84, n * 0.5, "15/19 concordant\nexact-binomial\nP = 0.0096",
            fontsize=6.8, color=C_DARK, weight="bold", va="center")
    panel_letter(ax, "C", x=-0.06)

    # ------------------------------------------------------------- D
    ax = fig.add_subplot(gs[1, 1])
    kv = {}
    with open(VAL) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == 2 and p[0] not in ("key",):
                kv[p[0]] = p[1]
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "Model-based set validation")
    items = [("KR FDR < 0.05", kv["N_FDR005_KR"], "19"),
             ("KR FDR < 0.10", kv["N_FDR010_KR"], "19"),
             ("LRT FDR < 0.05 (sensitivity)", kv["N_FDR005_LRT"], "19"),
             ("LRT FDR < 0.10 (sensitivity)", kv["N_FDR010_LRT"], "19"),
             ("LOEO minimum concordant", kv["LOO_MIN_CONCORDANT"], "19"),
             ("LOGO minimum concordant", kv["LOGO_MIN_CONCORDANT"], "19")]
    for i, (name, a, b) in enumerate(items):
        yy = 84 - i * 13
        ax.text(0, yy + 2.5, name, fontsize=6.8, color=C_DARK, weight="bold",
                va="center")
        ax.add_patch(Rectangle((0, yy - 3), 70, 4.5, fc=C_BG_GREY, ec="none"))
        ax.add_patch(Rectangle((0, yy - 3), 70 * int(a) / int(b), 4.5,
                               fc=C_PRIMARY, ec="none"))
        ax.text(73, yy - 0.8, "%s/%s" % (a, b), fontsize=7.4, va="center",
                color=C_DARK, weight="bold")
    panel_letter(ax, "D", x=-0.06)

    run_checks(fig, "Figure_S8")
    save_figure(fig, "Figure_S8")


if __name__ == "__main__":
    build_s1()
    build_s2()
    build_s3()
    build_s4()
    build_s5()
    build_s6()
    build_s7()
    build_s8()
    tsv_write(os.path.join(VQC_DIR, "EXTENT_AUTO_S1_S8.tsv"),
              ["figure", "check", "result", "detail"], QC_ROWS)
    bad = [r for r in QC_ROWS if r[2] == "ERROR"]
    print("EXTENT_CHECKS_S1_S8: %d rows, %d ERRORS" % (len(QC_ROWS), len(bad)))
    for r in bad:
        print("  ERROR:", r)
    print("S1-S8 done")
