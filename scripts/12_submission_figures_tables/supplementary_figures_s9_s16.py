#!/usr/bin/env python3
"""Targeted final visual fixes on top of the earlier S9-S16
renders (figure-only deliverables).

Changes applied in the final revision:
  S9: panel B legend box replaced by direct gene labels at the M4 line
      endpoints; panel C title simplified to 'Direction stability across
      M0-M4' (blue meaning in legend); panel D title simplified to
      'Leave-one-region-out Tier A effects' (4/4 stability in legend).
  S10: panel A title simplified to 'Analyzability by transcript-set
      definition' (green meaning in legend).
  S11: panel D floating 'reference threshold 0.50' text removed (dashed
      line kept; explained in legend).
  S12: landscape 2 x 3 grid relayout with larger gutters; panel D title
      shortened (full wording in legend).
  S13/S14: recheck only.
  S15: unchanged (pLDDT source check OK, recorded in 02_page_notes).
  S16: panel letter A removed; UpSet enlarged vertically; top blank
      reduced; explicit Tier A-D color legend added.

Original design notes follow.

Final re-render of Figures S9-S16 (figure-only deliverables).

Redesign per spec sections 5-6:
  S9: keep A-D, larger area; 'primary model M0' / 'technical-covariate
      model M4' wording; internal provenance language never appears.
  S10: keep A-D; 'evidence tier' wording; the italic in-panel note below
      panel C is removed (content in legend).
  S11: keep A-D; panel D simplified to concise event-level metrics (the
      per-event metric text block, the bold interpretation sentence and the
      internal 'phase' label are removed; legend label 'primary model M0').
  S12: A4-landscape 2 x 3 grid; in-panel prose block below B and the
      analytic-note text in E removed (legend); 'M4'/'composition-adjusted
      M4C'/'prespecified criterion' wording only.
  S13: keep A-C; per-panel x labels and the italic probability-scale
      footnote removed; one shared x-axis label.
  S14: keep A-B; sensitivity sentence footnote removed; '(primary19)' labels
      dropped; B legend moved away from the title.
  S15: keep A-D; pLDDT-context footnote removed (caveat in legend).
  S16: UpSet enlarged vertically; internal footer removed; row labels
      enlarged.

All numeric content is read from upstream analysis sources; nothing is
recomputed.  Every figure meets automated extent checks at render time.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supplementary_common import (ROOT, FIG_W, MM, VQC_DIR, TMP, D35_SUPP_TAB,
                             D36_ADJ, D36_PROT, D36_NM, DIR33, DIR34,
                             C_PRIMARY, C_CONCORD, C_DEV, C_DISCORD, C_NEG,
                             C_SENS, C_ORANGE, C_DARK, C_MID, C_BG_GREY,
                             MICROEXON_COLOR, TIER_COLORS, TIER_ORDER,
                             LAYER_NAMES, layer_flags, panel_letter,
                             panel_title, save_figure, tsv_write, rd,
                             load_master, tier_letter, extent_checks)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

ROOT32 = os.path.join(ROOT, "32_molecular_autism_integrated_validation_and_"
                            "sensitivity_20260805")
TECH = os.path.join(ROOT32, "02_psychencode_sensitivity",
                    "TECHNICAL_MODEL_EVENT_RESULTS.tsv")
LORO_RAW = os.path.join(ROOT32, "02_psychencode_sensitivity",
                        "LEAVE_ONE_REGION_OUT_RESULTS.tsv")
STAB = os.path.join(DIR33, "03_psychencode_technical_check",
                    "TECHNICAL_STABILITY_FINAL.tsv")
LOROF = os.path.join(DIR33, "03_psychencode_technical_check",
                     "LEAVE_ONE_REGION_OUT_FINAL.tsv")
EV_MASTER = os.path.join(DIR33, "04_transcript_definition_check",
                         "TRANSCRIPT_DEFINITION_EVENT_MASTER.tsv")
SET_SUM = os.path.join(DIR33, "04_transcript_definition_check",
                       "TRANSCRIPT_DEFINITION_SET_SUMMARY.tsv")
AXIS = "ASD diagnosis coefficient for logit-transformed transcript usage"
AXIS2 = "ASD diagnosis coefficient\n(logit-transformed transcript usage)"
TIER_A = ["HsaEX0015476", "HsaEX0029786", "HsaEX0050855", "HsaEX0051138"]
GENE_OF = {"HsaEX0015476": "CLASP1", "HsaEX0029786": "HERC4",
           "HsaEX0050855": "PTK2", "HsaEX0051138": "PTPRF"}
TA_COLS = [C_PRIMARY, C_ORANGE, C_CONCORD, C_SENS]

MT = load_master()
assert len(MT) == 19
QC_ROWS = []


def run_checks(fig, name):
    for check, ok, detail in extent_checks(fig, name):
        QC_ROWS.append([name, check, "OK" if ok else "ERROR", detail])
        print("  %-28s %-4s %s" % (check, "OK" if ok else "ERROR",
                                   detail[:70]))


# ====================================================================== S9
def build_s9():
    tech = rd(TECH)
    tech_map = {(r["HsaEX_ID"], r["model"]): r for r in tech}
    EVENTS = [r["HsaEX_ID"] for r in MT]
    GENE = {r["HsaEX_ID"]: r["gene"] for r in MT}
    TIER = {r["HsaEX_ID"]: tier_letter(r) for r in MT}
    TIERA = [e for e in EVENTS if TIER[e] == "A"]
    assert TIERA == TIER_A
    MODELS = ["M0_primary", "M1_nonlinear_age", "M2_seq_batch",
              "M3_ancestry", "M4_model_matrix_covariates"]
    SHORT = dict(M0_primary="M0", M1_nonlinear_age="M1",
                 M2_seq_batch="M2", M3_ancestry="M3",
                 M4_model_matrix_covariates="M4")

    stab = {r["model_id"]: r for r in rd(STAB)}
    m4s = stab["M4_model_matrix_covariates"]
    assert abs(float(m4s["beta_pearson_r_vs_M0"]) - 0.935) < 1e-6
    assert m4s["direction_agreement_vs_M0"] == "16/19"
    assert m4s["tierA_4of4_direction"] == "4/4"
    assert int(m4s["tierA_KR_BH_FDR_lt_005"]) == 4
    loroF = rd(LOROF)
    assert len(loroF) == 11
    assert all(r["tierA_direction_stable"] == "4/4" for r in loroF)
    assert all(r["any_tierA_reversal"] == "NO" for r in loroF)
    MAX_SHIFT = max(float(r["max_tierA_abs_beta_shift"]) for r in loroF)
    assert abs(MAX_SHIFT - 0.04805) < 1e-6

    loro_map = {}
    for r in rd(LORO_RAW):
        loro_map.setdefault(r["HsaEX_ID"], []).append(r)
    REGIONS = [r["region_removed"] for r in loroF]

    fig = plt.figure(figsize=(FIG_W, 172 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.70, wspace=0.45,
                           left=0.115, right=0.985, top=0.955, bottom=0.11)

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0, 0])
    x0 = [float(tech_map[(e, "M0_primary")]["beta_ASD"])
          for e in EVENTS]
    x4 = [float(tech_map[(e, "M4_model_matrix_covariates")]["beta_ASD"])
          for e in EVENTS]
    cols = [C_DISCORD if e in TIERA else C_NEG for e in EVENTS]
    ax.scatter(x0, x4, c=cols, s=28, edgecolors="black", linewidths=0.4,
               zorder=3)
    lim = max(abs(v) for v in x0 + x4) * 1.18
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.7, zorder=1)
    for e in TIERA:
        ax.annotate(GENE[e], (float(tech_map[(e, "M0_primary")]["beta_ASD"]),
                    float(tech_map[(e, "M4_model_matrix_covariates")]["beta_ASD"])),
                    xytext=(3, 3), textcoords="offset points", fontsize=6.6)
    ax.set_xlabel("primary model M0 beta", fontsize=8.0)
    ax.set_ylabel("technical-covariate model M4 beta", fontsize=8.0)
    panel_title(ax, "M0 vs M4 effects (Pearson r = %s)"
                % m4s["beta_pearson_r_vs_M0"])
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[0, 1])
    for i, e in enumerate(TIERA):
        ys = [float(tech_map[(e, m)]["beta_ASD"]) for m in MODELS]
        ax.plot(range(5), ys, "-o", color=TA_COLS[i], ms=4.0, lw=1.1)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(5))
    ax.set_xticklabels([SHORT[m] for m in MODELS])
    ax.set_ylabel(AXIS2)
    # Final: legend box removed; the four lines are directly labelled
    # at their M4 endpoints with collision-avoided y offsets.
    ends = sorted([float(tech_map[(e, MODELS[-1])]["beta_ASD"]), i]
                  for i, e in enumerate(TIERA))
    for k in range(1, len(ends)):
        if ends[k][0] - ends[k - 1][0] < 0.014:
            ends[k][0] = ends[k - 1][0] + 0.014
    for ylab, i in ends:
        ax.text(4.10, ylab, GENE[TIERA[i]], fontsize=6.8, va="center",
                color=TA_COLS[i], weight="bold")
    ax.set_xlim(-0.35, 5.55)
    panel_title(ax, "Tier A effects across M0\u2013M4")
    panel_letter(ax, "B")

    # ------------------------------------------------------------- C
    ax = fig.add_subplot(gs[1, 0])
    matsame = []
    for e in EVENTS:
        d0 = tech_map[(e, "M0_primary")]["direction"]
        matsame.append([1 if tech_map[(e, m)]["direction"] == d0 else 0
                        for m in MODELS])
    for i, row in enumerate(matsame):
        for j, v in enumerate(row):
            ax.add_patch(Rectangle((j, i), 1, 1,
                                   facecolor=C_PRIMARY if v else C_DISCORD,
                                   edgecolor="white", lw=0.5))
    ax.set_xlim(0, 5)
    ax.set_xticks([j + 0.5 for j in range(5)])
    ax.set_xticklabels([SHORT[m] for m in MODELS])
    ax.set_yticks([i + 0.5 for i in range(19)])
    ax.set_yticklabels(["%s %s" % (GENE[e], e[-6:]) for e in EVENTS],
                       fontsize=6.0)
    ax.set_ylim(19, 0)
    # Final: simplified title; the blue/same-as-M0 meaning and the
    # 16/19 M4 agreement are explained in the legend.
    panel_title(ax, "Direction stability across M0–M4", fontsize=8.6)
    panel_letter(ax, "C")

    # ------------------------------------------------------------- D
    ax = fig.add_subplot(gs[1, 1])
    for i, e in enumerate(TIERA):
        folds = sorted(loro_map[e],
                       key=lambda r: REGIONS.index(r["region_removed"]))
        ys = [float(f["beta_ASD_loro"]) for f in folds]
        full = float(folds[0]["beta_ASD_full"])
        ax.plot(range(11), ys, "-o", color=TA_COLS[i],
                label="%s (full=%.3f)" % (GENE[e], full), ms=3.6, lw=1.0)
        ax.axhline(full, color=TA_COLS[i], ls=":", lw=0.7)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(11))
    ax.set_xticklabels(REGIONS, rotation=55, ha="right", fontsize=6.2)
    ax.set_ylabel(AXIS2)
    # Final: 4/4 stability fact moved to the legend.
    panel_title(ax, "Leave-one-region-out Tier A effects", fontsize=7.8)
    ax.legend(frameon=False, fontsize=6.4)
    panel_letter(ax, "D")

    run_checks(fig, "Figure_S9")
    save_figure(fig, "Figure_S9")


# ====================================================================== S10
def build_s10():
    tdm = {}
    for r in rd(EV_MASTER):
        tdm.setdefault(r["event_id"], {})[r["definition"]] = r
    ss = {r["definition"]: r for r in rd(SET_SUM)}
    EVENTS = [r["HsaEX_ID"] for r in MT]
    GENE = {r["HsaEX_ID"]: r["gene"] for r in MT}
    TIER = {r["HsaEX_ID"]: tier_letter(r) for r in MT}
    TIERA = [e for e in EVENTS if TIER[e] == "A"]
    DEFS = ["D0", "D1", "D2", "D3"]

    assert [int(ss[d]["analyzable_n"]) for d in DEFS] == [19, 15, 19, 18]
    assert [int(ss[d]["direction_concordant_n"]) for d in DEFS] == \
        [14, 11, 15, 13]
    assert [int(ss[d]["tierA_FDR_lt_005_n"]) for d in DEFS] == [4, 2, 0, 2]
    for e in TIERA:
        for d in DEFS:
            x = tdm[e].get(d)
            if x is not None and x["analyzable"] == "YES":
                assert x["direction"] == "DOWN_IN_ASD", (e, d)

    TIER_COLORS_LOCAL = {"A": C_DISCORD, "B": C_ORANGE, "C": C_PRIMARY,
                         "D": C_NEG}

    fig = plt.figure(figsize=(FIG_W, 168 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.70, wspace=0.45,
                           left=0.115, right=0.985, top=0.955, bottom=0.10)

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0, 0])
    for i, e in enumerate(EVENTS):
        for j, d in enumerate(DEFS):
            x = tdm[e].get(d)
            ok = x is not None and x["analyzable"] == "YES"
            ax.add_patch(Rectangle((j, i), 1, 1,
                                   facecolor=C_CONCORD if ok else C_BG_GREY,
                                   edgecolor="white", lw=0.5))
    ax.set_xlim(0, 4)
    ax.set_ylim(19, 0)
    ax.set_xticks([j + 0.5 for j in range(4)])
    ax.set_xticklabels(["%s\nn=%d" % (d, int(ss[d]["analyzable_n"]))
                        for d in DEFS])
    ax.set_yticks([i + 0.5 for i in range(19)])
    ax.set_yticklabels(["%s %s" % (GENE[e], e[-6:]) for e in EVENTS],
                       fontsize=6.0)
    # Final: green = analyzable moved to the legend.
    panel_title(ax, "Analyzability by transcript-set definition",
                fontsize=8.6)
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[0, 1])
    import random
    random.seed(35)
    for e in EVENTS:
        for j, d in enumerate(DEFS):
            x = tdm[e].get(d)
            if x is not None and x["analyzable"] == "YES":
                ax.scatter(j + random.uniform(-0.12, 0.12), float(x["beta"]),
                           c=TIER_COLORS_LOCAL[TIER[e]], s=18,
                           edgecolors="black", linewidths=0.3, zorder=3)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(4))
    ax.set_xticklabels(DEFS)
    ax.set_ylabel(AXIS2)
    panel_title(ax, "Event effects per definition")
    handles = [Line2D([], [], marker="o", ls="", color=TIER_COLORS_LOCAL[t],
                      label="Tier %s" % t, ms=4.5) for t in TIER_ORDER]
    ax.legend(handles=handles, frameon=False, fontsize=6.6)
    panel_letter(ax, "B")

    # ------------------------------------------------------------- C
    ax = fig.add_subplot(gs[1, 0])
    for i, e in enumerate(TIERA):
        xs, ys, errs = [], [], []
        for j, d in enumerate(DEFS):
            x = tdm[e].get(d)
            if x is not None and x["analyzable"] == "YES":
                lo, hi = (float(v) for v in x["CI"].strip("[]").split(","))
                b = float(x["beta"])
                xs.append(j + (i - 1.5) * 0.16)
                ys.append(b)
                errs.append((b - lo, hi - b))
        ax.errorbar(xs, ys, yerr=list(zip(*errs)), fmt="-o", ms=4.0, lw=1.1,
                    color=TA_COLS[i], label=GENE[e])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(4))
    ax.set_xticklabels(DEFS)
    ax.set_ylabel(AXIS2)
    panel_title(ax, "Tier A effects with 95% CI")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.99),
              fontsize=6.6)
    panel_letter(ax, "C")

    # ------------------------------------------------------------- D
    ax = fig.add_subplot(gs[1, 1])
    xs = range(4)
    w = 0.35
    conc = [int(ss[d]["direction_concordant_n"]) for d in DEFS]
    fdr = [int(ss[d]["tierA_FDR_lt_005_n"]) for d in DEFS]
    ax.bar([x - w / 2 for x in xs], conc, width=w, color=C_PRIMARY,
           label="direction concordant with discovery")
    ax.bar([x + w / 2 for x in xs], fdr, width=w, color=C_DISCORD,
           label="Tier A KR BH-FDR < 0.05")
    for x in xs:
        ax.text(x - w / 2, conc[x] + 0.3, str(conc[x]), ha="center",
                fontsize=7.0)
        ax.text(x + w / 2, fdr[x] + 0.3, str(fdr[x]), ha="center",
                fontsize=7.0)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(DEFS)
    ax.set_ylim(0, 20)
    ax.set_ylabel("number of events")
    panel_title(ax, "Direction concordance vs Tier A significance",
                fontsize=8.4)
    ax.legend(frameon=False, fontsize=6.4)
    panel_letter(ax, "D")

    run_checks(fig, "Figure_S10")
    save_figure(fig, "Figure_S10")


# ====================================================================== S11
def build_s11():
    LONG = rd(os.path.join(DIR34, "03_lodo", "LODO_EVENT_RESULTS_LONG.tsv"))
    TAS = rd(os.path.join(DIR34, "03_lodo", "TIER_A_LODO_SUMMARY.tsv"))
    PHASE = open(os.path.join(DIR34, "03_lodo",
                             "tier_a_lodo_summary.txt")).read()
    assert "TIER_A_LODO_OVERALL=CONFIRMED" in PHASE
    for r in TAS:
        assert float(r["direction_preservation_fraction"]) >= 0.95
        assert int(r["direction_reversal_count"]) == 0
        assert float(r["min_abs_effect_retention"]) >= 0.50
        assert float(r["max_abs_DFBETA"]) < 1.0
    M0_long = [r for r in LONG if r["model_id"] == "M0_primary"]
    for ev in TIER_A:
        assert sum(1 for r in M0_long if r["event_id"] == ev) == 80

    fig = plt.figure(figsize=(FIG_W, 198 * MM))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.62, wspace=0.48,
                           left=0.11, right=0.945, top=0.955, bottom=0.09,
                           width_ratios=[1.0, 0.85])

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0, 0])
    for i, ev in enumerate(TIER_A):
        betas = [float(r["beta"]) for r in M0_long if r["event_id"] == ev]
        full = [r for r in TAS if r["model_id"] == "M0_primary"
                and r["event_id"] == ev][0]
        fb = float(full["full_beta"])
        jx = np.random.RandomState(7).uniform(-0.18, 0.18, len(betas))
        ax.scatter(np.full(len(betas), i) + jx, betas, s=8, color=C_DEV,
                   alpha=0.75, ec="white", lw=0.3, zorder=2)
        ax.plot([i - 0.3, i + 0.3], [fb, fb], "-", color=C_PRIMARY, lw=2.6,
                zorder=3)
        ax.text(i + 0.34, fb, "full %.3f" % fb, fontsize=6.2, va="center",
                color=C_PRIMARY, weight="bold")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["%s\n%s" % (GENE_OF[e], e.replace("HsaEX", ""))
                        for e in TIER_A], fontsize=6.8)
    ax.set_ylabel("ASD diagnosis coefficient\n(logit-transformed transcript "
                  "usage)")
    panel_title(ax, "Tier A effects across all 80 donor deletions")
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[0, 1])
    donors = sorted(set(r["removed_donor"] for r in M0_long))
    Z = np.zeros((len(donors), 4))
    for j, ev in enumerate(TIER_A):
        rows = [r for r in M0_long if r["event_id"] == ev]
        full = float([r for r in TAS if r["model_id"] == "M0_primary"
                      and r["event_id"] == ev][0]["full_beta"])
        vals = np.array([float(r["beta"]) for r in
                         sorted(rows, key=lambda r: r["removed_donor"])])
        sd = vals.std(ddof=1)
        Z[:, j] = (vals - full) / sd
    im = ax.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    ax.set_xticks(range(4))
    ax.set_xticklabels([GENE_OF[e] for e in TIER_A], fontsize=6.8)
    ax.set_yticks([0, 39, 79])
    ax.set_yticklabels(["donor 1", "40", "80"], fontsize=6.4)
    panel_title(ax, "Standardized donor-deletion deviation (M0)",
                fontsize=8.4)
    cb = plt.colorbar(im, ax=ax, fraction=0.024, pad=0.01)
    cb.set_label("(beta_d \u2212 full beta) / SD", fontsize=6.2)
    cb.ax.tick_params(labelsize=6.0)
    panel_letter(ax, "B")

    # ------------------------------------------------------------- C
    ax = fig.add_subplot(gs[1, 0])
    data, labs, cols = [], [], []
    for dx, c in [("ASD", C_DISCORD), ("CTL", C_CONCORD)]:
        vals = [float(r["beta"]) for r in M0_long
                if r["removed_diagnosis"] == dx]
        data.append(vals)
        labs.append("%s-donor deletions\n(%d \u00d7 4 events)" %
                    (dx, len(vals) // 4))
        cols.append(c)
    parts = ax.violinplot(data, showmedians=True, widths=0.8)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(cols[i])
        pc.set_alpha(0.7)
    for k in ("cmins", "cmaxes", "cbars", "cmedians"):
        parts[k].set_color(C_DARK)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(labs, fontsize=6.8)
    ax.set_ylabel("donor-deletion beta")
    ax.axhline(0, color=C_MID, ls=":", lw=0.7)
    med = [float(np.median(d)) for d in data]
    ax.text(0.97, 0.97, "median beta: ASD %.3f / CTL %.3f" % (med[0], med[1]),
            transform=ax.transAxes, ha="right", va="top", fontsize=6.6,
            color=C_DARK, weight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                      alpha=0.85))
    panel_title(ax, "ASD- vs control-donor deletion effects (Tier A)")
    panel_letter(ax, "C")

    # ------------------------------------------------- D simplified metrics
    ax = fig.add_subplot(gs[1, 1])
    m0 = [r for r in TAS if r["model_id"] == "M0_primary"]
    m4 = [r for r in TAS if r["model_id"] == "M4_model_matrix_covariates"]
    y = np.arange(4)[::-1]
    for off, rows, c in [(0.14, m0, C_PRIMARY), (-0.14, m4, C_SENS)]:
        for i, ev in enumerate(TIER_A):
            r = [q for q in rows if q["event_id"] == ev][0]
            ax.plot(float(r["min_abs_effect_retention"]), y[i] + off, "o",
                    ms=6.2, color=c, mec="white", mew=0.7, zorder=3)
    # Final: floating threshold text removed (dashed line kept;
    # meaning explained in the legend).
    ax.axvline(0.5, color=C_DISCORD, ls="--", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([GENE_OF[e] for e in TIER_A], fontsize=7.0)
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylim(-0.7, 4.05)
    ax.set_xlabel("minimum absolute effect retention\nacross 80 deletions")
    ax.plot([], [], "o", color=C_PRIMARY, ms=4.5, label="primary model M0")
    ax.plot([], [], "o", color=C_SENS, ms=4.5,
            label="technical-covariate model M4")
    ax.legend(loc="lower left", fontsize=6.4)
    panel_title(ax, "Direction preservation and effect retention",
                fontsize=8.0)
    panel_letter(ax, "D")

    run_checks(fig, "Figure_S11")
    save_figure(fig, "Figure_S11")


# ====================================================================== S12
def build_s12():
    CC = os.path.join(DIR34, "04_cell_composition")
    CM = os.path.join(DIR34, "05_composition_adjusted_models")
    FRAC = rd(os.path.join(CC, "COMPOSITION_FRACTIONS_HARMONIZED.tsv"))
    MARK = rd(os.path.join(CC, "COMPOSITION_MARKER_VALIDATION.tsv"))
    LOAD = rd(os.path.join(CC, "COMPOSITION_PC_LOADINGS.tsv"))
    PCDEC = open(os.path.join(CC, "COMPOSITION_PC_DECISION.md")).read()
    M4M4C = rd(os.path.join(CM, "M4_VS_M4C_SUMMARY.tsv"))
    D0D3 = rd(os.path.join(CM, "M4C_D0_D3_EVENT_RESULTS.tsv"))
    LODO_RES = rd(os.path.join(CM, "TIER_A_M4C_LODO_RESULTS.tsv"))
    LODO_SUM = rd(os.path.join(CM, "TIER_A_M4C_LODO_SUMMARY.tsv"))
    G33 = open(os.path.join(CM, "m4c_check_values.txt")).read()

    CLASSES = ["Excitatory_neuron", "Inhibitory_neuron", "Astrocyte",
               "Oligodendrocyte", "OPC", "Microglia_immune",
               "Endothelial_mural"]
    SHORT = {"Excitatory_neuron": "Exc", "Inhibitory_neuron": "Inh",
             "Astrocyte": "Ast", "Oligodendrocyte": "Oli", "OPC": "OPC",
             "Microglia_immune": "Mic", "Endothelial_mural": "End"}
    GENE_C = {"CLASP1": C_PRIMARY, "HERC4": C_ORANGE, "PTK2": C_CONCORD,
              "PTPRF": C_SENS}

    assert "M4_VS_M4C_BETA_PEARSON_R=0.9679" in G33
    assert "TIER_A_M4C_LODO_OVERALL=CONFIRMED" in G33
    assert "TIER_A_M4C_DIRECTION_4OF4=4/4" in G33

    count_ok = 0
    inh_error_offtarget = False
    for mc in CLASSES:
        rows = [r for r in MARK if r["marker_class"] == mc]
        match = [r for r in rows if r["match_type"] == "MATCHING"][0]
        off = [float(r["rho"]) for r in rows if r["match_type"] == "OFF_TARGET"]
        rho = float(match["rho"])
        ok = rho >= 0.3 and rho > max(off)
        count_ok += ok
        if mc == "Inhibitory_neuron":
            exc_off = [float(r["rho"]) for r in rows
                       if r["match_type"] == "OFF_TARGET" and
                       r["fraction_class"] == "Excitatory_neuron"][0]
            inh_error_offtarget = (rho >= 0.3) and (rho < exc_off)
    assert count_ok == 6
    assert inh_error_offtarget

    vline = [l for l in PCDEC.splitlines()
             if l.startswith("- Variance explained per PC:")][0]
    VAR = [float(x) for x in vline.split(":", 1)[1].split(",")]
    assert len(VAR) == 6 and abs(VAR[0] - 0.730) < 1e-6
    assert abs((VAR[0] + VAR[1]) - 0.854) < 1e-3

    # Final relayout: A4-landscape canvas with a clean 2 x 3 grid and
    # generous gutters so letters/titles/axes all sit inside the page and
    # C/F are never cropped.
    fig = plt.figure(figsize=(255 * MM, 162 * MM))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=1.05, wspace=0.62,
                           left=0.06, right=0.99, top=0.93, bottom=0.10)

    # ------------------------------------------------------------- A
    ax = fig.add_subplot(gs[0, 0])
    data = [[float(r[c]) for r in FRAC] for c in CLASSES]
    cols = [C_PRIMARY, C_DISCORD, C_ORANGE, C_SENS, C_CONCORD, C_DARK, C_MID]
    bp = ax.boxplot(data, showfliers=False, widths=0.62, patch_artist=True)
    for i, (box, med) in enumerate(zip(bp["boxes"], bp["medians"])):
        box.set_facecolor(cols[i])
        box.set_alpha(0.65)
        box.set_edgecolor(C_DARK)
        box.set_linewidth(0.7)
        med.set_color(C_DARK)
        med.set_linewidth(1.4)
    for k in ("whiskers", "caps"):
        for a in bp[k]:
            a.set_color(C_DARK)
            a.set_linewidth(0.7)
    ax.set_xticks(range(1, 8))
    ax.set_xticklabels([SHORT[c] for c in CLASSES], fontsize=6.6)
    ax.set_ylabel("Estimated cell fraction")
    ax.text(0.98, 0.97, "532 cortical samples\n80 donors",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.4,
            color=C_DARK, bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                    ec="none", alpha=0.85))
    panel_title(ax, "Estimated broad cell-fraction distributions",
                fontsize=8.6)
    panel_letter(ax, "A")

    # ------------------------------------------------------------- B
    ax = fig.add_subplot(gs[0, 1])
    Z = np.zeros((7, 7))
    for i, mc in enumerate(CLASSES):
        for j, fc in enumerate(CLASSES):
            Z[i, j] = float([r for r in MARK if r["marker_class"] == mc and
                             r["fraction_class"] == fc][0]["rho"])
    im = ax.imshow(Z, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(7))
    ax.set_xticklabels([SHORT[c] for c in CLASSES], fontsize=6.4)
    ax.set_yticks(range(7))
    ax.set_yticklabels([SHORT[c] for c in CLASSES], fontsize=6.4)
    ax.set_xlabel("fraction class", fontsize=6.8)
    ax.set_ylabel("marker class", fontsize=6.8)
    for i, mc in enumerate(CLASSES):
        rows = [r for r in MARK if r["marker_class"] == mc]
        rho = float([r for r in rows if r["match_type"] == "MATCHING"][0]["rho"])
        ok = rho >= 0.3 and rho > max(
            float(r["rho"]) for r in rows if r["match_type"] == "OFF_TARGET")
        ax.text(i, i, "+" if ok else "\u00d7", fontsize=8.6, va="center",
                ha="center", color="white", weight="bold", zorder=3)
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.10)
    cb.set_label("Spearman rho", fontsize=6.4)
    panel_title(ax, "Marker-score validation (6/7 classes ok)",
                fontsize=8.6)
    panel_letter(ax, "B")

    # ------------------------------------------------------------- C
    ax = fig.add_subplot(gs[0, 2])
    bars = ax.bar(range(1, 7), VAR, color=[C_PRIMARY if i < 2 else C_MID
                                           for i in range(6)],
                  edgecolor=C_DARK, linewidth=0.6)
    cum = np.cumsum(VAR)
    ax.plot(range(1, 7), cum, "-o", color=C_ORANGE, ms=3.0, lw=1.0)
    ax.axvline(2.5, color=C_DISCORD, ls="--", lw=0.8)
    ax.text(2.62, 0.95, "k=2 retained\n(cumulative 85.4%)", fontsize=6.2,
            color=C_DISCORD, va="top")
    ax.set_xticks(range(1, 7))
    ax.set_xticklabels(["PC%d" % i for i in range(1, 7)], fontsize=6.4)
    ax.set_ylabel("variance explained", fontsize=6.8)
    panel_title(ax, "Composition-PC variance explained", fontsize=8.6)
    panel_letter(ax, "C")

    # ------------------------------------------------------------- D
    ax = fig.add_subplot(gs[1, 0])
    xb = np.array([float(r["M4_beta"]) for r in M4M4C])
    yb = np.array([float(r["M4C_beta"]) for r in M4M4C])
    dr = np.array([r["direction_retained"] == "TRUE" for r in M4M4C])
    pr = float(np.corrcoef(xb, yb)[0, 1])
    assert abs(pr - 0.9679) < 2e-3
    assert int(dr.sum()) == 18
    ax.scatter(xb[~dr], yb[~dr], s=24, color=C_MID, ec=C_DARK, lw=0.6,
               zorder=2, label="direction flipped (1/19)")
    for r in M4M4C:
        if r["event_id"] in TIER_A:
            ax.scatter(float(r["M4_beta"]), float(r["M4C_beta"]), s=32,
                       color=GENE_C[GENE_OF[r["event_id"]]], ec="white",
                       lw=0.7, zorder=4)
            ax.annotate(GENE_OF[r["event_id"]], (float(r["M4_beta"]),
                        float(r["M4C_beta"])), fontsize=6.2, xytext=(5, 5),
                        textcoords="offset points", color=C_DARK,
                        weight="bold")
    keep = dr & ~np.isin(np.array([r["event_id"] for r in M4M4C]), TIER_A)
    ax.scatter(xb[keep], yb[keep], s=18, color=C_PRIMARY, alpha=0.7,
               ec="white", lw=0.5, zorder=3, label="direction retained")
    lo = min(xb.min(), yb.min()) - 0.02
    hi = max(xb.max(), yb.max()) + 0.02
    ax.plot([lo, hi], [lo, hi], "--", color=C_DARK, lw=0.8)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("M4 beta (no composition)", fontsize=6.8)
    ax.set_ylabel("composition-adjusted M4C beta", fontsize=6.8)
    ax.text(0.03, 0.97, "Pearson r = %.3f\ndirection retained 18/19\nmedian "
            "effect retention 0.662" % pr, transform=ax.transAxes, va="top",
            fontsize=6.4, color=C_DARK,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                      alpha=0.85))
    ax.legend(fontsize=6.0, loc="lower right")
    panel_title(ax, "M4 vs composition-adjusted M4C", fontsize=8.6)
    panel_letter(ax, "D")

    # ------------------------------------------------------------- E
    ax = fig.add_subplot(gs[1, 1])
    defs = ["D0", "D1", "D2", "D3"]
    for ev in TIER_A:
        g = GENE_OF[ev]
        xs, ys, loe, hie = [], [], [], []
        for d in defs:
            rr = [q for q in D0D3 if q["HsaEX_ID"] == ev and
                  q["definition"] == d]
            if not rr:
                continue  # PTPRF not analyzable under D1
            xs.append(defs.index(d))
            ys.append(float(rr[0]["beta_ASD"]))
            loe.append(float(rr[0]["beta_ASD"]) - float(rr[0]["CI_lo"]))
            hie.append(float(rr[0]["CI_hi"]) - float(rr[0]["beta_ASD"]))
        ax.errorbar(xs, ys, yerr=[loe, hie], fmt="-o", color=GENE_C[g],
                    ms=4.0, lw=1.2, capsize=2.0, mec="white", mew=0.6,
                    label=g, zorder=3)
    ax.set_xticks(range(4))
    ax.set_xticklabels(defs, fontsize=6.8)
    ax.set_xlabel("transcript-set definition (M4C)", fontsize=6.8)
    ax.set_ylabel("ASD diagnosis coefficient\n(logit-transformed transcript "
                  "usage)", fontsize=6.8)
    ax.axhline(0, color=C_MID, ls=":", lw=0.7)
    ax.legend(fontsize=6.2, ncol=4, loc="upper right")
    panel_title(ax, "Tier A effects across transcript-set definitions",
                fontsize=8.6)
    panel_letter(ax, "E")

    # ------------------------------------------------------------- F
    ax = fig.add_subplot(gs[1, 2])
    for i, ev in enumerate(TIER_A):
        rows = [r for r in LODO_RES if r["event_id"] == ev]
        assert len(rows) == 80
        betas = [float(r["beta"]) for r in rows]
        s = [r for r in LODO_SUM if r["event_id"] == ev][0]
        assert s["TIER_A_M4C_LODO_CONFIRMED"] == "TRUE"
        full = float(s["full_beta"])
        jx = np.random.RandomState(11).uniform(-0.18, 0.18, len(betas))
        ax.scatter(np.full(len(betas), i) + jx, betas, s=8, color=C_SENS,
                   alpha=0.75, ec="white", lw=0.3, zorder=2)
        ax.plot([i - 0.3, i + 0.3], [full, full], "-", color=C_PRIMARY,
                lw=2.6, zorder=3)
        # Final: two-line annotation + tighter xlim so the event
        # labels below the panel no longer collide.
        ax.text(i + 0.34, full, "full %.3f\nretention %.2f"
                % (full, float(s["min_abs_effect_retention"])), fontsize=5.8,
                va="center", color=C_DARK)
    ax.set_xlim(-0.6, 4.7)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["%s\n%s" % (GENE_OF[e], e.replace("HsaEX", ""))
                        for e in TIER_A], fontsize=6.4)
    ax.set_ylabel("M4C ASD diagnosis coefficient\n(logit-transformed "
                  "transcript usage)", fontsize=6.8)
    panel_title(ax, "Tier A M4C donor-deletion stability", fontsize=7.8)
    panel_letter(ax, "F")

    run_checks(fig, "Figure_S12")
    save_figure(fig, "Figure_S12")


# ====================================================================== S13
def build_s13():
    FOREST_SRC = os.path.join(D36_ADJ,
                              "TIER_A_ADJUSTED_TRANSCRIPT_USAGE_M0_M4_M4C.tsv")
    rows = rd(FOREST_SRC)
    assert len(rows) == 12
    assert all(float(r["adjusted_ASD_minus_control"]) < 0 for r in rows)
    MODELS = ["M0", "M4", "M4C"]
    MODEL_LABEL = {"M0": "M0 (clinical covariates)",
                   "M4": "M4 (M0 + technical covariates)",
                   "M4C": "M4C (M4 + composition PCs)"}
    COLOR = {"M0": C_PRIMARY, "M4": C_CONCORD, "M4C": MICROEXON_COLOR}
    EVENT_ORDER = list(TIER_A)
    by_model = {m: {r["event_id"]: r for r in rows if r["model"] == m}
                for m in MODELS}

    fig, axes = plt.subplots(1, 3, figsize=(FIG_W, 105 * MM), sharey=True)
    for ax, m, letter in zip(axes, MODELS, "ABC"):
        sub = by_model[m]
        for i, ev in enumerate(EVENT_ORDER):
            r = sub[ev]
            diff = float(r["adjusted_difference_percentage_points"]) / 100.0
            lo = float(r["adjusted_difference_95CI_low"])
            hi = float(r["adjusted_difference_95CI_high"])
            y = len(EVENT_ORDER) - 1 - i
            ax.errorbar(diff, y, xerr=[[diff - lo], [hi - diff]], fmt="o",
                        color=COLOR[m], ecolor=COLOR[m], elinewidth=1.6,
                        capsize=4, capthick=1.4, markersize=6.5, zorder=3)
            ax.text(diff, y + 0.30, r["gene"], fontsize=7.6, ha="center",
                    color=C_MID, style="italic")
        ax.axvline(0.0, color="0.6", linewidth=0.9, linestyle="--", zorder=1)
        panel_title(ax, MODEL_LABEL[m])
        panel_letter(ax, letter, x=-0.18 if m == "M0" else -0.24)
        ax.set_yticks(range(len(EVENT_ORDER)))
        ax.tick_params(axis="x", labelsize=7.8)
        lo_all = min(float(sub[ev]["adjusted_difference_95CI_low"])
                     for ev in EVENT_ORDER)
        hi_all = max(float(sub[ev]["adjusted_difference_95CI_high"])
                     for ev in EVENT_ORDER)
        ax.set_xlim(lo_all - 0.012, hi_all + 0.012)
        span = (hi_all + 0.012) - (lo_all - 0.012)
        step = 0.01 if span < 0.07 else 0.02
        t0 = np.ceil((lo_all - 0.012) / step) * step
        ax.set_xticks(np.arange(t0, hi_all + 0.012 + 1e-9, step))
    axes[0].set_ylim(-0.45, 3.75)
    axes[0].set_yticklabels(["%s\n%s" % (by_model["M0"][ev]["gene"], ev)
                             for ev in EVENT_ORDER][::-1], fontsize=7.6)
    for ax in axes[1:]:
        ax.set_yticklabels([])
    fig.supxlabel("Adjusted transcript-usage difference (ASD \u2212 control)",
                  fontsize=8.5)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    run_checks(fig, "Figure_S13")
    save_figure(fig, "Figure_S13")


# ====================================================================== S14
def build_s14():
    S6_TSV = os.path.join(D35_SUPP_TAB, "Supplementary_Table_S6.tsv")
    VS_TSV = os.path.join(D36_NM, "M4C_NEURON_MERGED_VS_PRIMARY_M4C.tsv")
    NM_TSV = os.path.join(D36_NM, "M4C_NEURON_MERGED_ALL19.tsv")

    def read_s6_section1(path):
        rows, header, in_sec = [], None, False
        with open(path) as f:
            for line in f:
                if line.startswith("## section:"):
                    if rows:
                        break
                    in_sec = True
                    continue
                if not in_sec or not line.strip():
                    continue
                if header is None:
                    header = line.rstrip("\n").split("\t")
                else:
                    rows.append(dict(zip(header,
                                         line.rstrip("\n").split("\t"))))
        assert len(rows) == 19, len(rows)
        return rows

    s6 = read_s6_section1(S6_TSV)
    m4 = {r["event_id"]: (float(r["M4_beta"]), float(r["M4_SE"]),
                          r["M4_KR_FDR"]) for r in s6}
    m4c_prim = {r["event_id"]: (float(r["M4C_beta"]), float(r["M4C_SE"]),
                                r["M4C_KR_FDR"]) for r in s6}
    vs = rd(VS_TSV)
    assert len(vs) == 19
    nm = {r["event_id"]: (float(r["beta"]), float(r["SE"]),
                          float(r["KR_BH_FDR"])) for r in rd(NM_TSV)}
    C_NM = MICROEXON_COLOR

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(FIG_W, 128 * MM))
    # ------------------------------------------------------------- panel A
    non = [r for r in vs if r["is_TierA"] != "TRUE"]
    tie = [r for r in vs if r["is_TierA"] == "TRUE"]
    axa.scatter([float(r["M4_beta_reference"]) for r in non],
                [float(r["M4C_NeuronMerged_beta"]) for r in non],
                s=28, facecolors="none", edgecolors="#666666", zorder=2,
                label="non-Tier A events")
    axa.scatter([float(r["M4_beta_reference"]) for r in tie],
                [float(r["M4C_NeuronMerged_beta"]) for r in tie],
                s=44, c=C_NM, zorder=3, label="Tier A events")
    LAB_OFF = {"CLASP1": (-6, 8, "right"), "HERC4": (6, -10, "left"),
               "PTK2": (7, 6, "left"), "PTPRF": (-6, -8, "right")}
    for r in tie:
        ox, oy, ha = LAB_OFF[r["gene"]]
        axa.annotate(r["gene"], (float(r["M4_beta_reference"]),
                     float(r["M4C_NeuronMerged_beta"])),
                     xytext=(ox, oy), textcoords="offset points",
                     fontsize=7.6, ha=ha)
    xs = [float(r["M4_beta_reference"]) for r in vs]
    ys = [float(r["M4C_NeuronMerged_beta"]) for r in vs]
    lims = [min(xs + ys) * 1.15, max(xs + ys) * 1.15]
    axa.plot(lims, lims, "--", color="#999999", lw=0.8, zorder=1)

    def inrng(lo, hi):
        step = max((hi - lo) / 4.0, 1e-9)
        for cand in (0.05, 0.1, 0.2, 0.25, 0.5, 1.0):
            if cand >= step:
                step = cand
                break
        t0 = np.ceil(lo / step) * step
        return np.arange(t0, hi + 1e-9, step)
    axa.set_xticks(inrng(*lims))
    axa.set_yticks(inrng(*lims))
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    r_p = (sum((x - mx) * (y - my) for x, y in zip(xs, ys)) /
           ((sum((x - mx) ** 2 for x in xs) *
             sum((y - my) ** 2 for y in ys)) ** 0.5))
    assert 0.98 < r_p < 0.995, r_p
    axa.set_xlabel("M4 ASD diagnosis coefficient\nfor logit-transformed "
                   "transcript usage", fontsize=8.0)
    axa.set_ylabel("Neuron-merged M4C coefficient\nfor logit-transformed "
                   "transcript usage", fontsize=8.0)
    panel_title(axa, "M4 vs neuron-merged coefficients (19 events)")
    panel_letter(axa, "A")
    axa.text(0.98, 0.02, "Pearson r = %.4f" % r_p, transform=axa.transAxes,
             ha="right", va="bottom", fontsize=7.6, color=C_MID)
    axa.legend(fontsize=7.4, frameon=False)
    axa.axhline(0, color="#CCCCCC", lw=0.5, zorder=0)
    axa.axvline(0, color="#CCCCCC", lw=0.5, zorder=0)
    axa.tick_params(labelsize=7.6)
    # ------------------------------------------------------------- panel B
    events = list(TIER_A)
    xvals = []
    for i, ev in enumerate(events):
        b1, s1, f1 = m4[ev]
        b2, s2, f2 = m4c_prim[ev]
        b3, s3, f3 = nm[ev]
        for j, (b, s, f, col, lab) in enumerate([
                (b1, s1, float(f1), C_PRIMARY, "M4"),
                (b2, s2, float(f2), C_CONCORD, "M4C seven-class"),
                (b3, s3, f3, C_NM, "M4C neuron-merged")]):
            y = i + (j - 1) * 0.26
            axb.plot([b - 1.96 * s, b + 1.96 * s], [y, y], "-", color=col,
                     lw=1.7, zorder=2)
            axb.plot([b], [y], "o", color=col, ms=6.2, zorder=3,
                     label=lab if i == 0 else None)
            axb.annotate("FDR %.3f" % f, (b + 1.96 * s + 0.004, y),
                         fontsize=7.2, va="center", color=C_MID)
            xvals.extend([b - 1.96 * s, b + 1.96 * s])
    axb.set_xlim(min(xvals) - 0.015,
                 max(xvals) + (max(xvals) - min(xvals)) * 0.30)
    blo, bhi = axb.get_xlim()
    axb.set_xticks(np.arange(np.ceil(blo / 0.1) * 0.1, bhi + 1e-9, 0.1))
    axb.axvline(0, color="#999999", lw=0.8, ls="--")
    axb.set_ylim(-0.8, 3.75)
    axb.set_yticks(range(4))
    axb.set_yticklabels([TIER_A[i] for i in range(4)], fontsize=7.8)
    axb.set_xlabel("ASD diagnosis coefficient for logit-transformed\n"
                   "transcript usage with 95% CI", fontsize=8.0)
    panel_title(axb, "Tier A coefficients by composition model")
    panel_letter(axb, "B")
    axb.legend(fontsize=7.2, loc="upper right", frameon=True,
               facecolor="white", edgecolor="none", framealpha=0.92)
    axb.tick_params(labelsize=7.6)
    fig.tight_layout()
    run_checks(fig, "Figure_S14")
    save_figure(fig, "Figure_S14")


# ====================================================================== S15
def build_s15():
    GENES = [("CLASP1", "Q7Z460", "HsaEX0015476", 673, 682, 1538),
             ("HERC4", "Q5GLZ8", "HsaEX0029786", 643, 650, 1057),
             ("PTK2", "Q05397", "HsaEX0050855", 393, 393, 1052),
             ("PTPRF", "P10586", "HsaEX0051138", 772, 780, 1907)]
    FCOLS = {"Domain": "#0072B2", "Topological domain": "#56B4E9",
             "Region": "#8C8C8C", "Compositional bias": "#BFBFBF",
             "Coiled coil": "#999999", "Motif": "#F0E442",
             "Modified residue": "#009E73", "Alternative sequence": "#CC79A7",
             "Binding site": "#E69F00", "Active site": "#E69F00",
             "Site": "#E69F00"}
    SITE_COL = MICROEXON_COLOR

    def parse_plddt(acc):
        for cand in ("%s/af_%s_v6.cif" % (TMP, acc),
                     "%s/af_%s.cif" % (TMP, acc)):
            if os.path.exists(cand):
                path = cand
                break
        else:
            raise SystemExit("AlphaFold CIF cache missing for %s" % acc)
        res_b = {}
        cols = []
        in_site = False
        with open(path) as fh:
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
        return res_b

    feat_rows = rd(os.path.join(D36_PROT, "TIER_A_PROTEIN_FEATURES.tsv"))

    fig, axes = plt.subplots(4, 1, figsize=(FIG_W, 200 * MM), sharex=False)
    for ax, (gene, acc, ev, lo, hi, L), letter in zip(
            axes, GENES, "ABCD"):
        pl = parse_plddt(acc)
        assert len(pl) > 0, acc
        xsk = sorted(pl)
        ax.plot(xsk, [pl[x] for x in xsk], color="#CCCCCC", lw=0.7, zorder=1)
        ax.axhline(70, color="#AAAAAA", lw=0.5, ls=":")
        ax.axhline(50, color="#AAAAAA", lw=0.5, ls=":")
        ax.set_ylim(0, 102)
        ax.set_xlim(0, L + 1)
        ax.set_yticks([0, 25, 50, 75, 100])
        step = 250 if L < 1600 else 500
        ax.set_xticks(np.arange(0, L + 1, step))
        ax.axvspan(lo - 0.5, hi + 0.5, color=SITE_COL, alpha=0.30, zorder=2)
        ax.axvline(lo if lo == hi else (lo + hi) / 2, color=SITE_COL, lw=1.2,
                   ls="--", zorder=3)
        drawn = [r for r in feat_rows if r["gene"] == gene and
                 r["proximity"] != "distant_gt_30aa"]
        drawn.sort(key=lambda r: int(r["feature_start"]))
        box_idx = 0
        for r in drawn:
            s, e = int(r["feature_start"]), int(r["feature_end"])
            ft = r["feature_type"]
            col = FCOLS.get(ft, "#E69F00")
            if ft == "Modified residue":
                ax.plot([s], [96], marker="|", color=col, ms=10, mew=1.6,
                        zorder=4)
            else:
                y = 88 - box_idx * 6
                box_idx += 1
                ax.add_patch(mpatches.Rectangle(
                    (s, y), max(e - s, 1), 4, facecolor=col,
                    edgecolor="black" if r["proximity"] == "overlap" else col,
                    lw=0.6, alpha=0.9 if r["proximity"] == "overlap"
                    else 0.65, zorder=4))
        panel_letter(ax, letter, x=-0.055, y=1.02)
        panel_title(ax, "%s (%s, %d aa), event %s" % (gene, acc, L, ev))
        ax.set_ylabel("pLDDT", fontsize=8.0)
        ax.tick_params(labelsize=7.4)
    axes[-1].set_xlabel("residue position (UniProt canonical numbering)",
                        fontsize=8.5)
    handles = [
        mpatches.Patch(color=SITE_COL, alpha=0.4,
                       label="microexon insertion site"),
        mpatches.Patch(color=FCOLS["Domain"], label="Domain"),
        mpatches.Patch(color=FCOLS["Topological domain"],
                       label="Topological domain"),
        mpatches.Patch(color=FCOLS["Region"], label="Region (incl. disordered)"),
        mpatches.Patch(color=FCOLS["Compositional bias"],
                       label="Compositional bias"),
        Line2D([0], [0], marker="|", color=FCOLS["Modified residue"], lw=0,
               ms=10, label="Modified residue"),
        mpatches.Patch(color=FCOLS["Alternative sequence"],
                       label="Alternative sequence (isoform difference)"),
        Line2D([0], [0], color="#CCCCCC", lw=1, label="AlphaFold v6 pLDDT"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.4,
               frameon=False, bbox_to_anchor=(0.5, 0.018))
    fig.tight_layout(rect=[0, 0.085, 1, 1])
    run_checks(fig, "Figure_S15")
    save_figure(fig, "Figure_S15")


# ====================================================================== S16
def build_s16():
    flags = [tuple(layer_flags(r)) for r in MT]
    tiers = [tier_letter(r) for r in MT]
    sigs = {}
    for f, t, r in zip(flags, tiers, MT):
        sigs.setdefault(f, []).append((t, r))
    items = sorted(sigs.items(), key=lambda kv: (-len(kv[1]), -sum(kv[0]),
                                                 kv[0]))
    ncol = len(items)
    assert ncol == 14, ncol
    nl = len(LAYER_NAMES)

    # Final: single-panel figure (letter removed), UpSet matrix enlarged
    # vertically, top blank reduced, explicit tier color legend added.
    figu = plt.figure(figsize=(FIG_W, 182 * MM))
    gsu = gridspec.GridSpec(2, 1, figure=figu, height_ratios=[1, 2.6],
                            hspace=0.10, left=0.18, right=0.97, top=0.94,
                            bottom=0.04)
    ax_bar = figu.add_subplot(gsu[0])
    ax_mat = figu.add_subplot(gsu[1])
    ax_bar.axis("off")
    ax_bar.set_xlim(-7.4, ncol - 0.5)
    mx = max(len(v) for v in sigs.values())
    ax_bar.set_ylim(0, mx + 1.5)
    for j, (sig, tvs) in enumerate(items):
        yb = 0
        for t in TIER_ORDER:
            k = sum(1 for (tt, _) in tvs if tt == t)
            if k:
                ax_bar.add_patch(Rectangle((j - 0.28, yb), 0.56, k,
                                           fc=TIER_COLORS[t], ec="none",
                                           zorder=2))
                yb += k
        ax_bar.text(j, len(tvs) + 0.25, str(len(tvs)), ha="center",
                    fontsize=8.0, weight="bold", color=C_DARK)
    ax_mat.axis("off")
    ax_mat.set_xlim(-7.4, ncol - 0.5)
    ax_mat.set_ylim(-0.7, nl - 0.3)
    for i, ln in enumerate(LAYER_NAMES):
        y = nl - 1 - i
        ax_mat.text(-0.7, y, ln.replace("\n", " "), ha="right", va="center",
                    fontsize=8.0, color=C_DARK)
        ax_mat.plot([-0.5, ncol - 0.5], [y, y], color="#dcdcdc", lw=0.6,
                    zorder=1)
    for j, (sig, tvs) in enumerate(items):
        ys = [nl - 1 - i for i in range(nl) if sig[i]]
        if len(ys) > 1:
            ax_mat.plot([j, j], [min(ys), max(ys)], color=C_DARK, lw=1.5,
                        zorder=2)
        for i in range(nl):
            y = nl - 1 - i
            if sig[i]:
                ax_mat.plot(j, y, "o", ms=6.0, color=C_DARK, zorder=3)
            else:
                ax_mat.plot(j, y, "o", ms=6.0, mfc="white", mec="#bdbdbd",
                            mew=0.7, zorder=3)
    tier_handles = [mpatches.Patch(color=TIER_COLORS[t], label="Tier %s" % t)
                    for t in TIER_ORDER]
    ax_bar.legend(handles=tier_handles, loc="upper left", frameon=False,
                  fontsize=7.4, bbox_to_anchor=(0.02, 0.98))
    panel_title(ax_bar, "Evidence-layer combinations")
    run_checks(figu, "Figure_S16")
    save_figure(figu, "Figure_S16")


if __name__ == "__main__":
    build_s9()
    build_s10()
    build_s11()
    build_s12()
    build_s13()
    build_s14()
    build_s15()
    build_s16()
    tsv_write(os.path.join(VQC_DIR, "EXTENT_AUTO_S9_S16.tsv"),
              ["figure", "check", "result", "detail"], QC_ROWS)
    bad = [r for r in QC_ROWS if r[2] == "ERROR"]
    print("EXTENT_CHECKS_S9_S16: %d rows, %d ERRORS" % (len(QC_ROWS), len(bad)))
    for r in bad:
        print("  ERROR:", r)
    print("S9-S16 done")
