#!/usr/bin/env python3
"""Final Figure 2: ASD cortex enrichment and
functional-direction bridge. A multi-background forest | B permutation-null
ridgeline | C robustness/sensitivity | D CHyMErA direction alluvial |
E event-level direction pairing.

Final fixes vs the earlier version:
  - 2C/2D: titles already use the recommended short forms ("Robustness and
    sensitivity" / "CHyMErA directional bridge"); the horizontal gutter between
    the two top-row panels is widened (C right edge 0.480 -> D left edge 0.575)
    so each panel title sits clearly inside its own panel area; no subtitle
    lines remain.
All numbers identical to the analysis sources and to the Earlier
SOURCE_DATA rows. No recomputation.

Local fixes for user inspection:
  - Figure 2C: removed the in-panel P = 0.05 threshold label, centered the
    -log10(P) mini-axis title, and moved 0/1/2/3 tick labels downward.
  - Figure 2D: removed the isolated bottom P = 0.0129 annotation; the exact
    binomial P value remains available for the figure legend/source data.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import *
import matplotlib.gridspec as gridspec  # noqa: F401
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde

MT = load_master()
EFF = load_effects()
NULLS = np.load(NULL_NPZ)  # render data, read-only
MASTER_TSV = os.path.join(DIR25, "06_master_event_table",
                          "MASTER_19_EVENT_EVIDENCE_TABLE.tsv")
EFF_TSV = os.path.join(DIR11, "07_primary_reanalysis",
                       "01_effects_by_background.tsv")
FIGDIR = FIG_DIRS[2]

BG_COLS = [C_PRIMARY, C_DEV, "#3a7ca5", C_ORANGE, C_SENS]

# QC bookkeeping
D_TEXTS = []     # (text_artist, owning_patch_or_None) for axD
D_PATCHES = []   # node patches in axD
E_XTICKS = []    # endpoint label texts in axE


def bezier_band(ax, x0, y0b, y0t, x1, y1b, y1t, color, alpha=0.55):
    xm = 0.5 * (x0 + x1)
    verts = [(x0, y0b), (xm, y0b), (xm, y1b), (x1, y1b),
             (x1, y1t), (xm, y1t), (xm, y0t), (x0, y0t), (x0, y0b)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(Path(verts, codes), fc=color, ec="none",
                                    alpha=alpha, zorder=1))


# --------------------------------------------------------------------------- A
def draw_forest(ax):
    y = np.arange(len(EFF))[::-1]
    for i, e in enumerate(EFF):
        ax.plot([e["lo"], e["hi"]], [y[i], y[i]], color=BG_COLS[i], lw=2.6,
                solid_capstyle="round", zorder=2)
        ax.plot(e["effect"], y[i], "o", ms=6, color=BG_COLS[i], mec="white",
                mew=0.8, zorder=3)
        if i == 0:
            ax.text(e["effect"], y[i] - 0.28, "%.4f" % e["effect"],
                    ha="center", va="top", fontsize=7.4, color=C_DARK,
                    zorder=3)
        else:
            ax.text(e["effect"], y[i] + 0.28, "%.4f" % e["effect"],
                    ha="center", va="bottom", fontsize=7.4, color=C_DARK,
                    zorder=3)
        p = e["perm_p"]
        ptxt = "P = 0.0001" if p < 0.001 else "P = %.4f" % p
        ax.text(0.0475, y[i], ptxt, va="center", fontsize=7.4, color=C_MID)
    ax.axvline(0, color=C_DARK, ls="--", lw=0.9, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([e["label"] for e in EFF])
    ax.set_xlim(-0.006, 0.062)
    # Final: deterministic in-range ticks (auto ticks placed labels
    # outside the axes)
    ax.set_xticks([0.0, 0.02, 0.04, 0.06])
    ax.set_xlabel("Mean |ΔPSI| difference vs background\n(95% bootstrap CI)")
    panel_title(ax, "Set enrichment vs five backgrounds")


# --------------------------------------------------------------------------- B
def draw_ridgeline(ax):
    obs = np.mean(np.abs([float(r["Parikshak_delta_PSI"]) for r in MT]))
    n = len(EFF)
    for i, e in enumerate(EFF):
        ns = NULLS[e["bg"]]
        kde = gaussian_kde(ns, bw_method="scott")
        xs = np.linspace(ns.min(), ns.max() + 0.004, 400)
        dens = kde(xs)
        base = (n - 1 - i) * 1.0
        scale = 0.92 / dens.max()
        ax.fill_between(xs, base, base + dens * scale, color=BG_COLS[i],
                        alpha=0.75, zorder=2)
        ax.plot(xs, base + dens * scale, color=BG_COLS[i], lw=0.8, zorder=3)
        ax.text(0.0525, base + 0.30, e["label"], fontsize=7.4, va="center",
                color=C_DARK)
        p = e["perm_p"]
        ptxt = "P = 0.0001" if p < 0.001 else "P = %.4f" % p
        ax.text(0.0525, base + 0.04, ptxt, fontsize=7.4, va="center",
                color=C_MID)
    ax.axvline(obs, color=C_DISCORD, lw=1.3, zorder=4)
    ax.text(obs - 0.0015, 5.05, "observed 19-event mean |ΔPSI| = %.4f" % obs,
            color=C_DISCORD, fontsize=7.6, ha="right", va="bottom",
            weight="bold")
    ax.set_yticks([])
    ax.set_ylim(-0.4, 5.7)
    ax.set_xlim(0.0, 0.053)
    # Final: deterministic in-range ticks (see panel A)
    ax.set_xticks([0.0, 0.02, 0.04])
    ax.set_xlabel("Permuted 19-event mean |ΔPSI|\n(seed = 42, 10,000 permutations)")
    panel_title(ax, "Permutation null distributions")


# --------------------------------------------------------------------------- C
def draw_robustness(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    # recommended short title (spec 2/Fig 2C)
    panel_title(ax, "Robustness and sensitivity", fontsize=9.0)
    frac = [("Leave-one-gene-out", 15, 15), ("Leave-one-event-out", 19, 19)]
    for i, (name, a, b) in enumerate(frac):
        yy = 76 - i * 26
        ax.text(2, yy + 7.5, name, fontsize=7.6, color=C_DARK, weight="bold")
        w = 40.0
        ax.add_patch(Rectangle((2, yy), w, 4.5, fc=C_BG_GREY, ec="none"))
        ax.add_patch(Rectangle((2, yy), w * a / b, 4.5, fc=C_CONCORD, ec="none"))
        ax.text(44, yy + 2.2, "%d/%d stable" % (a, b), fontsize=7.4,
                va="center", color=C_DARK)
    xs0 = 60.0
    tests = [("paired CEM sensitivity", 0.0198),
             ("paired NN sensitivity", 0.0038),
             ("ASD-prior exclusion", 0.0335)]
    axh = [84, 62, 40]
    # Mini x-axis for sensitivity P values.
    ax.text(xs0 + 15, 96, "-log10(P)", fontsize=7.2, color=C_MID,
            ha="center")
    ax.plot([xs0, xs0 + 30], [30, 30], color=C_DARK, lw=0.7)
    for t in [0, 1, 2, 3]:
        xx_tick = xs0 + t * 10
        ax.plot([xx_tick, xx_tick], [30, 28.6], color=C_DARK, lw=0.6)
        # Lower the tick labels so they do not sit on the axis line.
        ax.text(xx_tick, 23.5, str(t), fontsize=7.4, ha="center",
                va="top", color=C_MID)

    # P = 0.05 reference line.  The numeric threshold itself is explained in
    # the figure legend rather than printed inside the panel.
    thr = -np.log10(0.05)
    x_thr = xs0 + thr / 3 * 30
    ax.plot([x_thr, x_thr], [34, 92], color=C_DISCORD, ls="--", lw=0.8)

    for (name, p), yy in zip(tests, axh):
        v = -np.log10(p)
        xx = xs0 + v / 3 * 30
        ax.plot([xs0, xx], [yy, yy], color=C_SENS, lw=1.6, zorder=2)
        ax.plot(xx, yy, "o", ms=5, color=C_SENS, mec="white", mew=0.7, zorder=3)
        ax.text(xs0, yy + 4.0, name, fontsize=7.2, color=C_DARK, weight="bold")
        ax.text(xx + 2, yy, "P = %.4f" % p, fontsize=7.2, va="center",
                color=C_MID)


# --------------------------------------------------------------------------- D
def draw_alluvial(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    # recommended short title (spec 2/Fig 2D)
    panel_title(ax, "CHyMErA directional bridge", fontsize=9.0)
    S = 1.05
    W = 22.0
    xs = [3.0, 39.0, 75.0]

    def reg(p, t):
        D_PATCHES.append(p)
        D_TEXTS.append((t, p))

    h19 = 19 * S
    y19 = 56 - h19 / 2
    p = FancyBboxPatch((xs[0], y19), W, h19,
                       boxstyle="round,pad=0.2,rounding_size=0.8",
                       fc=C_MID, ec="none")
    ax.add_patch(p)
    t = ax.text(xs[0] + W / 2, 56, "19\nevents", ha="center",
                va="center", color="white", fontsize=7.6, weight="bold")
    reg(p, t)
    h14, h5 = 14 * S, 9 * S  # 5-box height raised (5->9 units) so the
    # two-line '5 with |dPSI| <= 0.01' label fits inside its box
    y14 = 76 - h14
    y5 = y14 - 6 - h5
    p = Rectangle((xs[1], y14), W, h14, fc=C_PRIMARY, ec="none")
    ax.add_patch(p)
    t = ax.text(xs[1] + W / 2, y14 + h14 / 2, "14 with\n|\u0394PSI| > 0.01", ha="center",
                va="center", color="white", fontsize=7.6, weight="bold")
    reg(p, t)
    p = Rectangle((xs[1], y5), W, h5, fc="#c9c9c9", ec="none")
    ax.add_patch(p)
    t = ax.text(xs[1] + W / 2, y5 + h5 / 2, "5 with\n|\u0394PSI| \u2264 0.01", ha="center",
                va="center", color=C_MID, fontsize=7.4, weight="bold")
    reg(p, t)
    h12, h2 = 12 * S, 2 * S
    y12 = 76 - h12
    y2 = y12 - 6 - h2
    p = Rectangle((xs[2], y12), W, h12, fc=C_CONCORD, ec="none")
    ax.add_patch(p)
    t = ax.text(xs[2] + W / 2, y12 + h12 / 2, "12\nconcordant", ha="center",
                va="center", color="white", fontsize=7.6, weight="bold")
    reg(p, t)
    p = Rectangle((xs[2], y2), W, h2, fc=C_DISCORD, ec="none")
    ax.add_patch(p)
    D_PATCHES.append(p)
    # discordant box too thin for internal text: centred label in the gap,
    # clear of both neighbouring boxes (checked by the computed phase)
    t = ax.text(xs[2] + W / 2, y2 + h2 + 0.9, "2 discordant", ha="center",
                va="bottom", color=C_DISCORD, fontsize=7.4, weight="bold")
    D_TEXTS.append((t, None))
    u = h19 / 19.0
    bezier_band(ax, xs[0] + W, y19 + 5 * u, y19 + 19 * u, xs[1], y14,
                y14 + h14, C_PRIMARY, 0.5)
    bezier_band(ax, xs[0] + W, y19, y19 + 5 * u, xs[1], y5, y5 + h5,
                "#c9c9c9", 0.5)
    u14 = h14 / 14.0
    bezier_band(ax, xs[1] + W, y14 + 2 * u14, y14 + 14 * u14, xs[2], y12,
                y12 + h12, C_CONCORD, 0.5)
    bezier_band(ax, xs[1] + W, y14, y14 + 2 * u14, xs[2], y2, y2 + h2,
                C_DISCORD, 0.5)


# --------------------------------------------------------------------------- E
def draw_direction(ax):
    test = [r for r in MT if r["CHyMErA_direction_concordant"] in ("YES", "NO")]
    test.sort(key=lambda r: (r["CHyMErA_direction_concordant"] != "YES",
                             -float(r["Parikshak_delta_PSI"])))
    n = len(test)
    y = np.arange(n)[::-1]
    x_disc, x_chy = 1.0, 2.0
    for i, r in enumerate(test):
        conc = r["CHyMErA_direction_concordant"] == "YES"
        col = C_CONCORD if conc else C_DISCORD
        d_disc = 1 if r["Parikshak_direction"] == "UP_IN_ASD" else -1
        d_chy = -1  # CHyMErA perturbation: microexon inclusion loss
        ax.plot([x_disc, x_chy], [y[i] + 0.16 * d_disc, y[i] + 0.16 * d_chy],
                "-", color=col, lw=2.0, zorder=2, alpha=0.9)
        ax.plot(x_disc, y[i] + 0.16 * d_disc, "o", ms=5.5, color=col,
                mec="white", mew=0.7, zorder=3)
        ax.plot(x_chy, y[i] + 0.16 * d_chy, "s", ms=5.5, color=col, mec="white",
                mew=0.7, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(["%s–%s" % (r["gene"],
                        r["HsaEX_ID"].replace("HsaEX", "")) for r in test],
                       fontsize=8.0)
    ax.set_xticks([x_disc, x_chy])
    ax.set_xticklabels(["ASD cortex\ndirection",
                        "CHyMErA microexon-inclusion-loss\ndirection"],
                       fontsize=8.0)
    ax.set_xlim(0.55, 2.45)
    ax.set_ylim(-0.8, n - 0.2)
    # NOTE: no `loss`/`gain` side words (spec); endpoint labels only
    for tl in ax.get_xticklabels():
        E_XTICKS.append(tl)
    ax.plot([], [], "-", color=C_CONCORD, lw=2.0, label="direction concordant")
    ax.plot([], [], "-", color=C_DISCORD, lw=2.0, label="direction discordant")
    ax.legend(loc="upper right", fontsize=8.0, ncol=2,
              bbox_to_anchor=(1.0, 1.12), columnspacing=1.2)
    panel_title(ax, "Event-level direction pairing (14 events with |\u0394PSI| > 0.01)")


# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(FIG_W, 230 * MM))
# Final: top row raised and mid row lowered so the A/B two-line x-axis
# labels no longer collide with the C/D panel titles.
axA = fig.add_axes([0.075, 0.707, 0.40, 0.26])
axB = fig.add_axes([0.565, 0.707, 0.42, 0.26])
# Final: wider horizontal gutter between C and D (0.065 -> 0.095) so the
# two panel titles sit clearly within their own panel areas.
axC = fig.add_axes([0.045, 0.353, 0.435, 0.26])
axD = fig.add_axes([0.575, 0.353, 0.41, 0.26])
axE = fig.add_axes([0.085, 0.055, 0.88, 0.24])

draw_forest(axA)
panel_letter(axA, "A")
draw_ridgeline(axB)
panel_letter(axB, "B")
draw_robustness(axC)
panel_letter(axC, "C")
draw_alluvial(axD)
panel_letter(axD, "D")
draw_direction(axE)
panel_letter(axE, "E")

# ----------------------------------------------------------- QC (computed)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()

# 2D: no text may overlap a node/band patch other than its own label box,
# and no pair of texts may overlap
fig2d_violations = 0
for t, own in D_TEXTS:
    tbb = t.get_window_extent(renderer=renderer)
    for p in D_PATCHES:
        if p is own:
            continue
        pbb = p.get_window_extent(renderer=renderer)
        if overlap_area(tbb, pbb) > 0:
            fig2d_violations += 1
for i in range(len(D_TEXTS)):
    for j in range(i + 1, len(D_TEXTS)):
        a = D_TEXTS[i][0].get_window_extent(renderer=renderer)
        b = D_TEXTS[j][0].get_window_extent(renderer=renderer)
        x_ov = min(a.x1, b.x1) - max(a.x0, b.x0)
        y_ov = min(a.y1, b.y1) - max(a.y0, b.y0)
        if x_ov > 2.0 and y_ov > 2.0:
            fig2d_violations += 1

# 2E: endpoint x-tick labels must not overlap each other
fig2e_violations = 0
bbs = [t.get_window_extent(renderer=renderer) for t in E_XTICKS]
for i in range(len(bbs)):
    for j in range(i + 1, len(bbs)):
        if overlap_area(bbs[i], bbs[j]) > 0:
            fig2e_violations += 1


# Final: title-collision checks for C and D.
def _title_bb(axx):
    """Bbox of the real (non-empty) title artist. NOTE: panel_title uses
    loc='left', which sets ax._left_title rather than ax.title."""
    for tt in (axx.title, axx._left_title, axx._right_title):
        if tt.get_visible() and tt.get_text().strip():
            return tt.get_window_extent(renderer=renderer)
    return axx.title.get_window_extent(renderer=renderer)


def _all_text(axx):
    out = []
    for t in axx.texts:
        if t.get_visible() and t.get_text().strip():
            out.append(t.get_window_extent(renderer=renderer))
    for t in list(axx.get_xticklabels()) + list(axx.get_yticklabels()):
        if t.get_visible() and t.get_text().strip():
            out.append(t.get_window_extent(renderer=renderer))
    for tt in (axx.title, axx._left_title, axx._right_title,
               axx.xaxis.label, axx.yaxis.label):
        if tt.get_visible() and tt.get_text().strip():
            out.append(tt.get_window_extent(renderer=renderer))
    return out


PX_PER_MM = fig.dpi / 25.4


def _title_check(tbb, other_bbs):
    """Violation count + minimum horizontal gap between one panel's title and
    the other panel's artists. Only vertically-overlapping artists constrain
    the gutter (a title sits in the band above its own axes row)."""
    viol = 0
    gmin = 1e9
    for bb in other_bbs:
        if overlap_area(tbb, bb) > 0:
            viol += 1
            continue
        y_ov = min(tbb.y1, bb.y1) - max(tbb.y0, bb.y0)
        if y_ov <= 0:
            continue
        g = max(bb.x0, tbb.x0) - min(bb.x1, tbb.x1)
        gmin = min(gmin, g)
        if g < 6.0 * PX_PER_MM:
            viol += 1
    return viol, gmin


tc_bb = _title_bb(axC)
td_bb = _title_bb(axD)
fig2c_title_viol, fig2c_title_gap = _title_check(tc_bb, _all_text(axD))
fig2d_title_viol, fig2d_title_gap = _title_check(td_bb, _all_text(axC))

# Final check: the two-line A/B x-axis labels must not collide with the
# C/D panel titles (vertical cross-row collision seen in Earlier).
fig2_midrow_viol = 0
for axx_up, axx_lo in [(axA, axC), (axB, axD)]:
    ubb = axx_up.xaxis.label.get_window_extent(renderer=renderer)
    lbb = _title_bb(axx_lo)
    # C/D titles must stay >= 2 mm below the upper panel's x-axis label
    # (vertical gap negative => overlap)
    if (ubb.y0 - lbb.y1) < 2.0 * PX_PER_MM:
        fig2_midrow_viol += 1

checks = [("FIG2D_TEXT_OVERLAP", fig2d_violations),
         ("FIG2E_ENDPOINT_TEXT_OVERLAP", fig2e_violations),
         ("FIG2C_TITLE_COLLISION_ZERO", fig2c_title_viol),
         ("FIG2D_TITLE_COLLISION_ZERO", fig2d_title_viol),
         ("FIG2_MIDROW_TITLE_OVERLAP", fig2_midrow_viol),
         ("FIG2C_TITLE_MIN_GAP_PX", round(fig2c_title_gap, 1)),
         ("FIG2D_TITLE_MIN_GAP_PX", round(fig2d_title_gap, 1))]
_mf2 = min_text_size_pt(fig)
tsv_write(os.path.join(QC_DIR, "figure2_layout_checks.tsv"),
          ["check", "value"], [[g, v] for g, v in checks] +
          [["MIN_TEXT_SIZE_PT", _mf2],
           ["FIG_WIDTH_MM", round(fig.get_size_inches()[0] * 25.4, 2)]])
print("FIG2 checks:", checks)
print("FIG2 min text size (pt):", _mf2)

# ----------------------------------------------------------- render
NAME = "Figure2_FINAL_SUBMISSION_CLEAN"
save_final(fig, FIGDIR, NAME)

# ----------------------------------------------------------- source data TSV
# Rows identical to Earlier Figure2_FINAL_SOURCE_DATA.tsv.
rows = []
for e in EFF:
    rows.append(["Figure_2", "A", "enrichment forest", e["label"],
                 "effect=%.6g CI=[%.6g, %.6g] P=%.6g n_bg=%d" %
                 (e["effect"], e["lo"], e["hi"], e["perm_p"], e["n_bg"]),
                 EFF_TSV, e["bg"]])
obs = float(np.mean(np.abs([float(r["Parikshak_delta_PSI"]) for r in MT])))
rows.append(["Figure_2", "B", "null ridgeline", "observed 19-event mean |dPSI|",
             "%.6f" % obs, MASTER_TSV, "Parikshak_delta_PSI"])
rows.append(["Figure_2", "B", "null ridgeline", "permutation nulls",
             "seed=42, 10,000 permutations per background", NULL_NPZ,
             "render data (from the analysis backgrounds)"])
LOO_DIR = os.path.join(ROOT, "18_psychencode",
                       "11_set_level_validation")
for name, a, b, f in [("Leave-one-gene-out stable", 15, 15, "06_LOO_gene.tsv"),
                      ("Leave-one-event-out stable", 19, 19, "05_LOO_event.tsv")]:
    rows.append(["Figure_2", "C", "robustness", name, "%d/%d" % (a, b),
                 os.path.join(LOO_DIR, f), "LOO minima summaries"])
for name, p in [("Matched CEM", 0.0198),
                ("Matched NN", 0.0038),
                ("ASD-prior exclusion", 0.0335)]:
    rows.append(["Figure_2", "C", "sensitivity tests", name, "P = %.4f" % p,
                 os.path.join(DIR11, "07_primary_reanalysis"), name])
for _item, _val in [("total events", "19"),
                    ("events with |delta_PSI| > 0.01", "14"),
                    ("concordant", "12"),
                    ("discordant", "2"),
                    ("events with |delta_PSI| <= 0.01", "5"),
                    ("exact binomial P", "0.012939453125")]:
    rows.append(["Figure_2", "D", "CHyMErA bridge", _item, _val, MASTER_TSV,
                 "CHyMErA_direction_concordant"])
for r in MT:
    if r["CHyMErA_direction_concordant"] in ("YES", "NO"):
        rows.append(["Figure_2", "E", "direction pairing",
                     "%s (%s)" % (r["gene"], r["HsaEX_ID"]),
                     "discovery=%s %s" %
                     (r["Parikshak_direction"],
                      ("expected direction=microexon inclusion loss"
                       if r["CHyMErA_direction_concordant"] == "YES"
                       else "discordant with microexon inclusion loss")),
                     MASTER_TSV,
                     "Parikshak_direction + CHyMErA_direction_concordant"])
tsv_write(os.path.join(FIGDIR, NAME + "_SOURCE_DATA.tsv"),
          ["figure", "panel", "series", "item", "value", "source_file",
           "source_locator"], rows)

prov = [
    ["Figure_2", "A", "Set enrichment vs five backgrounds (forest)",
     "mean |dPSI| difference vs background, 95% bootstrap CI, permutation P",
     EFF_TSV, "permutation P (valid discovery-enrichment statistics, dir 11)",
     "font sizes/clearance raised for final width (Earlier); Final: "
     "deterministic in-range ticks"],
    ["Figure_2", "B", "Permutation null distributions (ridgeline)",
     "10,000-permutation null per background vs observed mean",
     EFF_TSV + " ; " + NULL_NPZ, "empirical permutation P",
     "font sizes raised (Earlier); Final: deterministic in-range ticks"],
    ["Figure_2", "C", "Robustness and sensitivity",
     "LOO stability 15/15 and 19/19; paired CEM/NN sensitivity and "
     "ASD-prior-exclusion Ps",
     EFF_TSV, "sensitivity P values (nominal, described as sensitivity)",
     "Final: short title retained; C-D horizontal gutter widened"],
    ["Figure_2", "D", "CHyMErA directional bridge (alluvial)",
     "19 events -> 14 eligible (|delta_PSI| > 0.01) -> "
     "12 concordant / 2 discordant",
     MASTER_TSV, "exact binomial test (compact value only)",
     "Final: short title retained; C-D horizontal gutter widened"],
    ["Figure_2", "E", "Event-level direction pairing",
     "per-event discovery direction vs CHyMErA perturbation direction",
     MASTER_TSV, "none (classification display)",
     "loss/gain side words removed; endpoint labels only (Earlier)"],
]
tsv_write(os.path.join(FIGDIR, NAME + "_PROVENANCE.tsv"),
          ["figure", "panel", "panel_title", "content", "data_source_files",
           "statistics_displayed", "corrections_vs_earlier"], prov)
print("Figure2_FINAL_v2 done")
