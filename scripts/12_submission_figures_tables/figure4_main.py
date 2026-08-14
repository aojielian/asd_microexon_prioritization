#!/usr/bin/env python3
"""Final Figure 4 final layout.

Panel A: PsychENCODE event-level effects (KR primary forest)
Panel B: Discovery vs PsychENCODE effect correspondence
Panel C: KR vs LRT sensitivity
Panel D: Set-level robustness

The former donor-accounting panel has been excluded from the main figure
because it is descriptive cohort bookkeeping rather than a data result. The
donor-flow facts (112 total, 23 overlap excluded, 9 Dup15q excluded, 80
analyzed = 38 ASD
+ 42 control; 532 cortical samples) should be stated in the Figure 4 legend and
Methods instead. No statistical result is changed.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon_main import *
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path  # noqa: F401
import matplotlib.patches as mpatches  # noqa: F401
from scipy.stats import beta as _beta

MT = load_master()
FIGDIR = FIG_DIRS[4]
MASTER_TSV = os.path.join(DIR25, "06_master_event_table",
                          "MASTER_19_EVENT_EVIDENCE_TABLE.tsv")
KV18 = {r["key"]: r["value"] for r in
        rd(os.path.join(ROOT, "18_psychencode",
                        "15_qc", "key_counts.tsv"))}


def clopper_pearson(k, n):
    """Exact binomial 95% CI (Clopper-Pearson)."""
    return _beta.ppf(0.025, k, n - k + 1), _beta.ppf(0.975, k + 1, n - k)


# QC bookkeeping
A_NODES = []     # (text, patch) for containment
A_NOTES = []     # free exclusion notes
C_ANNO = []      # 4C annotation texts
D_ANNO = []      # 4D annotation + threshold texts
E_TEXTS = []     # 4E annotation texts


# --------------------------------------------------------------------------- A
def draw_cohort(ax):
    """Large, minimal donor-accounting strip.

    The panel is deliberately not a Sankey and not a compact stack. Three
    large sequential cohort states occupy most of the panel width; exclusions
    are shown as small annotations above the arrows; the final ASD/control
    composition is shown as a single segmented bar below the 80-donor state.
    """
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "PsychENCODE donor accounting", fontsize=9.0)

    def state_box(cx, cy, w, h, number, label, fc, tc="white"):
        x = cx - w / 2
        y = cy - h / 2
        p = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.02,rounding_size=1.6",
            fc=fc, ec=C_DARK, lw=0.75, zorder=3
        )
        ax.add_patch(p)
        t1 = ax.text(cx, cy + 3.2, str(number), ha="center", va="center",
                     fontsize=12.5, weight="bold", color=tc, zorder=4)
        t2 = ax.text(cx, cy - 4.5, label, ha="center", va="center",
                     fontsize=7.3, weight="bold", color=tc, zorder=4,
                     linespacing=0.95)
        A_NODES.extend([(t1, p), (t2, p)])
        return p

    def arrow(x0, x1, y=61):
        a = FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                            mutation_scale=12, lw=1.05,
                            color=C_MID, zorder=2)
        ax.add_patch(a)

    # Main sequential accounting row: deliberately large and sparse.
    y_main = 61
    state_box(14, y_main, 24, 24, 112, "total\ndonors", C_MID)
    state_box(50, y_main, 25, 24, 89, "after overlap\nexclusion",
              C_DEV, tc=C_DARK)
    state_box(86, y_main, 24, 24, 80, "analysis\ndonors", C_PRIMARY)

    arrow(26.5, 37.0, y_main)
    arrow(62.8, 73.5, y_main)

    # Exclusion counts are annotations above the arrows, not separate boxes.
    ax.text(31.8, 74.5, "−23 overlap", ha="center", va="bottom",
            fontsize=7.2, color=C_MID, weight="bold")
    ax.text(68.2, 74.5, "−9 Dup15q", ha="center", va="bottom",
            fontsize=7.2, color=C_MID, weight="bold")

    # Final diagnostic composition: one horizontal segmented bar.
    bar_x, bar_y, bar_w, bar_h = 68.0, 26.0, 28.0, 11.0
    frac_asd = 38 / 80.0
    w_asd = bar_w * frac_asd
    w_ctl = bar_w - w_asd

    p1 = FancyBboxPatch((bar_x, bar_y), w_asd, bar_h,
                        boxstyle="round,pad=0.02,rounding_size=1.0",
                        fc=C_DISCORD, ec=C_DARK, lw=0.6, zorder=3)
    p2 = FancyBboxPatch((bar_x + w_asd, bar_y), w_ctl, bar_h,
                        boxstyle="round,pad=0.02,rounding_size=1.0",
                        fc=C_CONCORD, ec=C_DARK, lw=0.6, zorder=3)
    ax.add_patch(p1)
    ax.add_patch(p2)
    t1 = ax.text(bar_x + w_asd / 2, bar_y + bar_h / 2, "38 ASD",
                 ha="center", va="center", fontsize=7.7, color="white",
                 weight="bold", zorder=4)
    t2 = ax.text(bar_x + w_asd + w_ctl / 2, bar_y + bar_h / 2,
                 "42 control", ha="center", va="center", fontsize=7.7,
                 color="white", weight="bold", zorder=4)
    A_NODES.extend([(t1, p1), (t2, p2)])

    # Connect retained cohort to the composition bar with a simple stem.
    ax.plot([86, 86], [49, 40], color=C_MID, lw=1.0, zorder=2)
    ax.add_patch(FancyArrowPatch((86, 40), (86, 38), arrowstyle="-|>",
                                 mutation_scale=11, lw=1.0,
                                 color=C_MID, zorder=2))
    ax.text(82, 20.5, "80 donors: 38 ASD / 42 control",
            ha="center", va="top", fontsize=7.1, color=C_MID)


# --------------------------------------------------------------------------- B
def draw_forest(ax):
    rows = sorted(MT, key=lambda r: -float(r["PsychENCODE_beta"]))
    y = np.arange(len(rows))[::-1]
    los = [float(r["PsychENCODE_CI95_lower"]) for r in rows]
    his = [float(r["PsychENCODE_CI95_upper"]) for r in rows]
    ax.set_xlim(min(los) - 0.004, max(his) + 0.004)
    # Final: deterministic in-range ticks (auto ticks placed labels
    # outside the axes)
    ax.set_xticks([-1.0, -0.5, 0.0])
    for i, r in enumerate(rows):
        b = float(r["PsychENCODE_beta"])
        conc = r["direction_concordant"] == "TRUE"
        col = C_CONCORD if conc else C_DISCORD
        ax.plot([los[i], his[i]], [y[i], y[i]], "-", color=col, lw=2.2,
                solid_capstyle="round", zorder=2)
        mk = "o" if float(r["BH_FDR_KR"]) >= 0.05 else "*"
        ax.plot(b, y[i], mk, ms=6.5 if mk == "*" else 5, color=col,
                mec="white", mew=0.7, zorder=3)
    ax.axvline(0, color=C_DARK, ls="--", lw=0.9, zorder=1)
    ax.set_yticks(y)
    labs = ["%s–%s" % (r["gene"], r["HsaEX_ID"].replace("HsaEX", ""))
            for r in rows]
    ax.set_yticklabels(labs, fontsize=7.2)
    # Tier A emphasis (no new semantics: bold row labels only)
    for r, tl in zip(rows, ax.get_yticklabels()):
        if tier_letter(r) == "A":
            tl.set_fontweight("bold")
    ax.set_xlabel("ASD diagnosis coefficient for\nlogit-transformed transcript usage (95% CI)")
    ax.plot([], [], "o", color=C_CONCORD, ms=4, label="direction concordant")
    ax.plot([], [], "o", color=C_DISCORD, ms=4, label="discordant")
    ax.plot([], [], "*", color=C_DARK, ms=7, label="KR BH-FDR < 0.05")
    ax.legend(loc="upper center", fontsize=7.2, bbox_to_anchor=(0.5, -0.17),
              ncol=3, frameon=False, columnspacing=1.4, handletextpad=0.6)
    panel_title(ax, "PsychENCODE event-level effects", fontsize=8.9)


# --------------------------------------------------------------------------- C
def draw_disc_scatter(ax):
    for r in MT:
        conc = r["direction_concordant"] == "TRUE"
        col = C_CONCORD if conc else C_DISCORD
        ax.plot(float(r["Parikshak_delta_PSI"]), float(r["PsychENCODE_beta"]),
                "o", ms=4.2, color=col, mec="white", mew=0.5, zorder=3)
    OFF_C = {"HsaEX0015476": (6, -4, "left"), "HsaEX0029786": (6, -9, "left"),
             "HsaEX0050855": (6, 7, "left"), "HsaEX0051138": (-6, 7, "right")}
    for r in MT:
        if tier_letter(r) == "A":
            dx, dy, ha = OFF_C[r["HsaEX_ID"]]
            C_ANNO.append(ax.annotate(
                r["gene"],
                (float(r["Parikshak_delta_PSI"]),
                 float(r["PsychENCODE_beta"])),
                fontsize=7.3, color=C_DARK, xytext=(dx, dy), ha=ha,
                textcoords="offset points", weight="bold"))
    ax.axhline(0, color=C_MID, ls=":", lw=0.7)
    ax.axvline(0, color=C_MID, ls=":", lw=0.7)
    # Final: explicit in-range ticks (auto ticks placed out-of-range
    # labels into the panel gutters)
    ax.set_xticks([-0.10, -0.05, 0.0, 0.05])
    ax.set_yticks([-0.6, -0.4, -0.2, 0.0, 0.2])
    ax.set_xlabel("Discovery ΔPSI (ASD cortex)")
    ax.set_ylabel("ASD diagnosis coefficient\n(logit-transformed transcript usage)")
    C_ANNO.append(ax.text(0.03, 0.97, "Spearman ρ = 0.4088\nP = 0.0823",
                          transform=ax.transAxes, va="top", fontsize=7.6,
                          weight="bold", color=C_DARK))
    panel_title(ax, "Discovery vs PsychENCODE", fontsize=8.9)


# --------------------------------------------------------------------------- D
def draw_kr_lrt(ax):
    xs = [-np.log10(max(float(r["BH_FDR_KR"]), 1e-4)) for r in MT]
    ys = [-np.log10(max(float(r["BH_FDR_LRT"]), 1e-4)) for r in MT]
    ax.set_xlim(-0.25, max(xs) * 1.15)
    ax.set_ylim(-0.25, max(ys) * 1.15)
    for r, x, yv in zip(MT, xs, ys):
        ax.plot(x, yv, "o", ms=4.2, color=TIER_COLORS[tier_letter(r)],
                mec="white", mew=0.5, zorder=3)
    # Keep the 0.05 and 0.10 FDR reference lines, but move their
    # definitions to the figure legend; in-panel labels made the upper region
    # unnecessarily crowded.
    for thr in (0.05, 0.10):
        v = -np.log10(thr)
        ax.axvline(v, color=C_MID, ls="--", lw=0.7)
        ax.axhline(v, color=C_MID, ls="--", lw=0.7)
    OFF_D = {"HsaEX0015476": (6, 5, "left"), "HsaEX0029786": (6, -8, "left"),
             "HsaEX0050855": (-6, 6, "right"), "HsaEX0051138": (-6, -7,
                                                                "right")}
    for r in MT:
        if r["HsaEX_ID"] in OFF_D:
            x = -np.log10(max(float(r["BH_FDR_KR"]), 1e-4))
            yv = -np.log10(max(float(r["BH_FDR_LRT"]), 1e-4))
            dx, dy, ha = OFF_D[r["HsaEX_ID"]]
            D_ANNO.append(ax.annotate(
                r["gene"], (x, yv), fontsize=7.3, color=C_DARK,
                xytext=(dx, dy), ha=ha, textcoords="offset points",
                weight="bold"))
    # Final: explicit in-range ticks (see 4C)
    ax.set_xticks([0.0, 1.0, 2.0])
    ax.set_yticks([0.0, 1.0, 2.0])
    ax.set_xlabel("KR −log10(BH-FDR)")
    ax.set_ylabel("LRT −log10(BH-FDR)")
    # italic "LRT is sensitivity, not primary" relocated to legend (spec)
    panel_title(ax, "KR–LRT sensitivity", fontsize=8.9)


# --------------------------------------------------------------------------- E
def draw_props(ax, axt):
    """Point estimate + exact binomial 95% CI; annotations in side column."""
    panel_title(ax, "Set-level robustness", fontsize=9.0)
    items = [("All-event direction\nconcordance", 15, 19, "binom"),
             ("One-event-per-gene\nconcordance", 12, 15, "binom"),
             ("Leave-one-event-out\nminimum concordant", 14, 19, "min"),
             ("Leave-one-gene-out\nminimum concordant", 13, 19, "min")]
    y = np.arange(len(items))[::-1]
    for (name, k, n, kind), yy in zip(items, y):
        p = k / float(n)
        if kind == "binom":
            lo, hi = clopper_pearson(k, n)
            ax.plot([lo, hi], [yy, yy], "-", color=C_PRIMARY, lw=2.6,
                    solid_capstyle="round", zorder=2)
            ax.plot(p, yy, "o", ms=6.5, color=C_PRIMARY, mec="white", mew=0.8,
                    zorder=3)
            E_TEXTS.append(axt.text(0.02, yy,
                                    "%d/%d = %.2f\n95%% CI %.2f–%.2f" %
                                    (k, n, p, lo, hi),
                                    fontsize=7.2, va="center", color=C_DARK))
        else:
            ax.plot(p, yy, "D", ms=6, color=C_ORANGE, mec="white", mew=0.8,
                    zorder=3)
            E_TEXTS.append(axt.text(0.02, yy,
                                    "%d/%d = %.2f\nminimum" % (k, n, p),
                                    fontsize=7.2, va="center", color=C_DARK))
    ax.axvline(0.5, color=C_MID, ls="--", lw=0.8, zorder=1)
    ax.text(0.5, 3.62, "chance = 0.5", fontsize=7.0, color=C_MID,
            ha="center", va="bottom")
    ax.set_yticks(y)
    # Final: compact one-line row labels keep the D-E gutter wide;
    # the full item names are preserved in SOURCE_DATA and the legend.
    ax.set_yticklabels(["All-event", "Per-gene", "LOO-event", "LOO-gene"],
                       fontsize=7.2)
    ax.set_xlim(-0.04, 1.12)
    ax.set_ylim(-0.65, 3.80)
    # Final: explicit in-range ticks (see 4C)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_xlabel("direction-concordant proportion")
    axt.axis("off")
    axt.set_xticks([])
    axt.set_yticks([])
    axt.set_ylim(-0.65, 3.80)
    ax.plot([], [], "-", color=C_PRIMARY, lw=2.6,
            label="proportion + exact binomial 95% CI")
    ax.plot([], [], "D", color=C_ORANGE, ms=5,
            label="minimum robustness metric (no binomial CI)")
    ax.legend(loc="upper left", fontsize=7.0, bbox_to_anchor=(0.0, -0.24),
              ncol=1)


# ---------------------------------------------------------------------------
# Final four-panel layout: descriptive donor accounting is legend-only.
# The forest is given the full top row; the three inferential/robustness panels
# occupy the bottom row with generous gutters.
fig = plt.figure(figsize=(FIG_W, 226 * MM))
gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.12, 0.88],
                       hspace=0.58, left=0.075, right=0.985, top=0.965,
                       bottom=0.07)

# A — full-width PsychENCODE forest
axA = fig.add_subplot(gs[0])
draw_forest(axA)
panel_letter(axA, "A")

# B/C/D — lower row
gs_bot = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1],
                                          width_ratios=[0.34, 0.27, 0.39],
                                          wspace=0.72)
axB = fig.add_subplot(gs_bot[0])
draw_disc_scatter(axB)
panel_letter(axB, "B")

axC = fig.add_subplot(gs_bot[1])
draw_kr_lrt(axC)
panel_letter(axC, "C")

gsD = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_bot[2],
                                       width_ratios=[0.48, 0.52], wspace=0.08)
axD = fig.add_subplot(gsD[0])
axDt = fig.add_subplot(gsD[1])
draw_props(axD, axDt)
panel_letter(axD, "D")

# ----------------------------------------------------------- QC (computed)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
figbb = fig.get_window_extent()

# A: forest tick labels remain inside the figure and all CIs remain in x-limits
fig4a = 0
for tl in axA.get_yticklabels():
    if not containment(tl.get_window_extent(renderer=renderer), figbb, tol=1.0):
        fig4a += 1
x0, x1 = axA.get_xlim()
for r in MT:
    if float(r["PsychENCODE_CI95_lower"]) < x0 or \
       float(r["PsychENCODE_CI95_upper"]) > x1:
        fig4a += 1

# B/C: pairwise annotation overlap
def pair_overlap(texts, tol=2.0):
    n = 0
    bb = [t.get_window_extent(renderer=renderer) for t in texts]
    for i in range(len(bb)):
        for j in range(i + 1, len(bb)):
            x_ov = min(bb[i].x1, bb[j].x1) - max(bb[i].x0, bb[j].x0)
            y_ov = min(bb[i].y1, bb[j].y1) - max(bb[i].y0, bb[j].y0)
            if x_ov > tol and y_ov > tol:
                n += 1
    return n

fig4b = pair_overlap(C_ANNO)
fig4c = pair_overlap(D_ANNO)

# D: annotation texts contained in dedicated text column
fig4d = 0
tbb = axDt.get_window_extent(renderer=renderer)
for t in E_TEXTS:
    if not containment(t.get_window_extent(renderer=renderer), tbb, tol=1.0):
        fig4d += 1

# Bottom-row gutter checks
def _ax_texts(axx):
    out = []
    for t in axx.texts:
        if t.get_visible() and t.get_text().strip():
            out.append(t.get_window_extent(renderer=renderer))
    for t in list(axx.get_xticklabels()) + list(axx.get_yticklabels()):
        if t.get_visible() and t.get_text().strip():
            out.append(t.get_window_extent(renderer=renderer))
    for tt in (axx.title, axx._left_title, axx._right_title):
        if tt.get_visible() and tt.get_text().strip():
            out.append(tt.get_window_extent(renderer=renderer))
    out.append(axx.xaxis.label.get_window_extent(renderer=renderer))
    out.append(axx.yaxis.label.get_window_extent(renderer=renderer))
    return out

b_ts = _ax_texts(axB)
c_ts = _ax_texts(axC)
d_ts = _ax_texts(axD) + _ax_texts(axDt)
PX_PER_MM = fig.dpi / 25.4

def h_gap(group_l, group_r):
    g = 1e9
    for a in group_l:
        for b in group_r:
            dy = max(a.y0, b.y0) - min(a.y1, b.y1)
            if dy >= 0:
                continue
            dx = max(a.x0, b.x0) - min(a.x1, b.x1)
            g = min(g, dx)
    return g

gap_bc = h_gap(b_ts, c_ts)
gap_cd = h_gap(c_ts, d_ts)
FIG4B_C_GUTTER_OK = int(gap_bc >= 4.0 * PX_PER_MM)
FIG4C_D_GUTTER_OK = int(gap_cd >= 4.0 * PX_PER_MM)

b_bb = axB.get_window_extent(renderer=renderer)
c_bb = axC.get_window_extent(renderer=renderer)
d_bb = axD.get_window_extent(renderer=renderer)
dt_bb = axDt.get_window_extent(renderer=renderer)
fig4_bot_viol = 0
for bb in b_ts:
    if overlap_area(bb, c_bb) > 0 or overlap_area(bb, d_bb) > 0:
        fig4_bot_viol += 1
for bb in c_ts:
    if overlap_area(bb, b_bb) > 0 or overlap_area(bb, d_bb) > 0 or \
            overlap_area(bb, dt_bb) > 0:
        fig4_bot_viol += 1
for bb in d_ts:
    if overlap_area(bb, c_bb) > 0 or overlap_area(bb, b_bb) > 0:
        fig4_bot_viol += 1

checks = [("FIG4A_FOREST_CLIPPING", fig4a),
         ("FIG4B_TEXT_OVERLAP", fig4b),
         ("FIG4C_TEXT_OVERLAP", fig4c),
         ("FIG4D_TEXT_CLIPPING", fig4d),
         ("FIG4B_C_GUTTER_SUFFICIENT", FIG4B_C_GUTTER_OK),
         ("FIG4C_D_GUTTER_SUFFICIENT", FIG4C_D_GUTTER_OK),
         ("FIG4_BOTTOM_ROW_BOUNDARY_COLLISION_ZERO", fig4_bot_viol),
         ("FIG4B_C_GUTTER_PX", round(gap_bc, 1)),
         ("FIG4C_D_GUTTER_PX", round(gap_cd, 1))]
_mf4 = min_text_size_pt(fig)
tsv_write(os.path.join(QC_DIR, "figure4_layout_checks.tsv"),
          ["check", "value"], [[g, v] for g, v in checks] +
          [["MIN_TEXT_SIZE_PT", _mf4],
           ["FIG_WIDTH_MM", round(fig.get_size_inches()[0] * 25.4, 2)]])
print("FIG4 v5 checks:", checks)
print("FIG4 v5 min text size (pt):", _mf4)

# ----------------------------------------------------------- render
NAME = "Figure4_FINAL_v5"
save_final(fig, FIGDIR, NAME)

# ----------------------------------------------------------- source data TSV
rows = []

# Cohort accounting is legend-only in v5 (not a plotted panel).
kc_path = os.path.join(ROOT, "18_psychencode",
                       "15_qc", "key_counts.tsv")
n_total = int(KV18["N_DONORS_TOTAL"])
n_overlap = int(KV18["N_OVERLAP_DONORS"])
n_dup = int(KV18["N_DUP15Q_DONORS"])
n_retained = n_total - n_overlap
n_analysis = n_retained - n_dup
assert (n_total, n_retained, n_overlap, n_dup, n_analysis) == (112, 89, 23, 9, 80)
assert int(KV18["N_NONOVERLAP_ASD"]) == 38
assert int(KV18["N_NONOVERLAP_CTL"]) == 42
for label, v, loc in [
        ("112 total donors", n_total, "N_DONORS_TOTAL"),
        ("23 overlap donors excluded", n_overlap, "N_OVERLAP_DONORS"),
        ("9 Dup15q donors excluded", n_dup, "N_DUP15Q_DONORS"),
        ("80 analysis donors", n_analysis, "derived"),
        ("38 ASD donors", int(KV18["N_NONOVERLAP_ASD"]), "N_NONOVERLAP_ASD"),
        ("42 control donors", int(KV18["N_NONOVERLAP_CTL"]), "N_NONOVERLAP_CTL")]:
    rows.append(["Figure_4", "Legend", "cohort accounting", label, v,
                 kc_path, loc])
rows.append(["Figure_4", "Legend", "cohort accounting",
             "532 cortical samples", 532,
             os.path.join(DIR25, "07_check_reports",
                          "ALL_NUMBERS_FOR_MANUSCRIPT.tsv"),
             "PsychENCODE/Samples"])

# A — forest
for r in sorted(MT, key=lambda r: -float(r["PsychENCODE_beta"])):
    rows.append(["Figure_4", "A", "KR forest",
                 "%s (%s)" % (r["gene"], r["HsaEX_ID"]),
                 "beta=%.4f CI=[%.4f, %.4f] KR_FDR=%.4g concordant=%s" %
                 (float(r["PsychENCODE_beta"]),
                  float(r["PsychENCODE_CI95_lower"]),
                  float(r["PsychENCODE_CI95_upper"]),
                  float(r["BH_FDR_KR"]), r["direction_concordant"]),
                 MASTER_TSV,
                 "PsychENCODE_beta/CI95/BH_FDR_KR"])

# B — discovery/PsychENCODE correlation
rows.append(["Figure_4", "B", "effect correlation", "Spearman rho / P",
             "rho=0.4088 P=0.0823 (not significant)",
             os.path.join(DIR25, "07_check_reports",
                          "ALL_NUMBERS_FOR_MANUSCRIPT.tsv"),
             "reference value 0.4088/0.0823"])

# C — KR/LRT sensitivity
for r in MT:
    rows.append(["Figure_4", "C", "KR vs LRT",
                 "%s (%s)" % (r["gene"], r["HsaEX_ID"]),
                 "KR_FDR=%.4g LRT_FDR=%.4g" %
                 (float(r["BH_FDR_KR"]), float(r["BH_FDR_LRT"])),
                 MASTER_TSV, "BH_FDR_KR/BH_FDR_LRT"])

# D — set-level robustness
for name, k, n in [("All-event direction concordance", 15, 19),
                   ("One-event-per-gene concordance", 12, 15),
                   ("Leave-one-event-out minimum concordant", 14, 19),
                   ("Leave-one-gene-out minimum concordant", 13, 19)]:
    rows.append(["Figure_4", "D", "set-level robustness", name,
                 "%d/%d" % (k, n),
                 MASTER_TSV + " ; set-level validation", name])

tsv_write(os.path.join(FIGDIR, NAME + "_SOURCE_DATA.tsv"),
          ["figure", "panel", "series", "item", "value", "source_file",
           "source_locator"], rows)

prov = [
    ["Figure_4", "Legend", "PsychENCODE cohort accounting",
     "112 total donors; 23 discovery-overlap and 9 Dup15q donors excluded; "
     "80 analysis donors (38 ASD, 42 control); 532 cortical samples",
     "15_qc/key_counts.tsv", "counts only",
     "Descriptive cohort bookkeeping moved out of the plot and into legend/Methods"],
    ["Figure_4", "A", "PsychENCODE event-level effects",
     "ASD diagnosis coefficient for logit-transformed transcript usage, 95% CI, "
     "KR BH-FDR; Tier A labels bold", MASTER_TSV,
     "KR primary model; BH-FDR family",
     "Forest promoted to full-width top panel after removal of donor-accounting panel"],
    ["Figure_4", "B", "Discovery vs PsychENCODE",
     "discovery dPSI vs PsychENCODE coefficient; Spearman rho/P; Tier A gene labels",
     MASTER_TSV, "Spearman correlation (not significant)",
     "compact gene-only Tier A labels"],
    ["Figure_4", "C", "KR–LRT sensitivity",
     "-log10 BH-FDR per event, KR vs LRT; Tier A gene labels", MASTER_TSV,
     "KR primary / LRT sensitivity FDR families",
     "FDR threshold lines retained; textual threshold labels kept in legend"],
    ["Figure_4", "D", "Set-level robustness",
     "15/19 and 12/15 exact binomial CI; LOO minima 14/19, 13/19",
     MASTER_TSV + " ; set-level validation",
     "exact binomial CI; minima without CI",
     "dedicated annotation column retained"]
]
tsv_write(os.path.join(FIGDIR, NAME + "_PROVENANCE.tsv"),
          ["figure", "panel", "panel_title", "content", "data_source_files",
           "statistics_displayed", "corrections_vs_earlier"], prov)
print("Figure4_FINAL_v5 done")
