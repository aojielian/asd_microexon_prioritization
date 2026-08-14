#!/usr/bin/env python3
"""Final Figure 5: integrated multi-layer
prioritization of ASD-associated microexons.

Final fixes vs the earlier version:
  - B/C: wider B-C gutter and the spec-recommended B:C width ~0.44:0.56;
    taller B/C/D row on a larger canvas (no compression; spec section 4).
  - 5B: label repulsion ok (iterative display-space separation of the
    event labels) + a near-collision check (FIG5B_LABEL_CROWDING_ZERO).
  - 5C: taller rows (pitch 19 -> 20 units, cells 11 -> 12 units) and more
    header-to-cell breathing (headers raised to y=92, first row lowered).
All numbers identical to the published sources and to the Earlier
SOURCE_DATA rows (byte-identical). No recomputation. The moved old-5B UpSet
supplementary candidate is NOT rebuilt here (earlier output, untouched).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import *
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle

MT = load_master()
ROWS = sorted(MT, key=lambda r: (tier_rank(r), r["gene"], r["HsaEX_ID"]))
# Mini-ok v2: legend vertical position (figure-fraction y), tunable via env
# 0.378 chosen by checks-only scan: clears B x-ticks (+34 px) and x-axis title
# (+8.5 px), keeps >=91 px clearance to panel-D title, all overlap checks 0.
LEGEND_Y = float(os.environ.get("FIG5_LEGEND_Y", "0.378"))
FIGDIR = FIG_DIRS[5]
MASTER_TSV = os.path.join(DIR25, "06_master_event_table",
                          "MASTER_19_EVENT_EVIDENCE_TABLE.tsv")
NET_RECLASS = os.path.join(DIR33, "06_pathway_network_reclassification",
                           "NETWORK_PUBLIC_INTERPRETATION.tsv")
assert all(r["permutation_inference_status"] == "INVALID"
           for r in rd(NET_RECLASS))
FEAT_TSV = os.path.join(DIR36, "03_tierA_protein_mapping",
                        "TIER_A_PROTEIN_FEATURES.tsv")
COD_TSV = os.path.join(DIR36, "03_tierA_protein_mapping",
                       "TIER_A_CODING_CONSEQUENCE.tsv")
FEATS = rd(FEAT_TSV)
COD = rd(COD_TSV)

# short horizontal headers (definitions in legend; spec §9 suggested set)
COLS = ["CHyMErA", "Dev", "Network", "GSE-\n30573", "Psych-\nENCODE",
        "KR\nFDR<0.05", "LRT\nFDR<0.05"]
COLS_FLAT = ["CHyMErA", "Dev", "Network", "GSE30573", "PsychENCODE",
             "KR FDR<0.05", "LRT FDR<0.05"]


def cell_state(r):
    """Per-layer states, standardized symbols:
    + supporting | x discordant | en-dash absent / not significant / NA."""
    gse = "–"
    if r["GSE30573_mapping_status"] == "MAPPED_ANALYZABLE":
        gse = "+" if r["GSE30573_direction_concordant"] == "CONCORDANT" \
            else "×"
    return [
        {"YES": "+", "NO": "×"}.get(r["CHyMErA_direction_concordant"], "–"),
        "+" if r["developmental_dynamic_status"] == "DYNAMIC" else "–",
        "+" if r["network_module_or_pathway"] == "SET_LEVEL_NETWORK_MEMBERSHIP"
        else "–",
        gse,
        "+" if r["direction_concordant"] == "TRUE" else "×",
        "+" if float(r["BH_FDR_KR"]) < 0.05 else "–",
        "+" if float(r["BH_FDR_LRT"]) < 0.05 else "–",
    ]


STATE_COL = {"+": C_CONCORD, "×": C_DISCORD, "–": C_BG_GREY}
STATE_TXT = {"+": "white", "×": "white", "–": C_MID}

# QC bookkeeping
A_HEADERS = []
C_CELLS = []      # (text_artist, rect_patch)
C_HEADER_TXT = []
D_DRAWN = []      # (gene, start, end, kind, label) verified-feature register


# --------------------------------------------------------------------------- A
# per-column data widths (mm-scale): wide enough for each bold header line
W_COL = [1.20, 0.62, 1.10, 0.86, 1.10, 1.30, 1.30]
CENT, _lft = [], 0.0
for _w in W_COL:
    CENT.append(_lft + _w / 2.0)
    _lft += _w
XR = _lft + 0.15


def draw_matrix(ax):
    ax.axis("off")
    n = len(ROWS)
    ax.set_xlim(-2.45, XR)
    ax.set_ylim(n + 1.6, -2.4)
    panel_title(ax, "Multi-layer evidence matrix", fontsize=9.0)
    for j, c in enumerate(COLS):
        A_HEADERS.append(ax.text(CENT[j], -1.05, c, ha="center", va="center",
                                 fontsize=7.4, color=C_DARK, weight="bold"))
    left = 0.0
    for i, r in enumerate(ROWS):
        ax.text(-0.35, i + 0.5, compact_label(r), ha="right", va="center",
                fontsize=7.2, color=C_DARK)
        ax.text(-1.95, i + 0.5, tier_letter(r), ha="center", va="center",
                fontsize=7.4, color="white", weight="bold",
                bbox=dict(boxstyle="square,pad=0.30",
                          fc=TIER_COLORS[tier_letter(r)], ec="none"))
        left = 0.0
        for j, s in enumerate(cell_state(r)):
            gap = 0.045
            ax.add_patch(Rectangle((left + gap, i + 0.06),
                                   W_COL[j] - 2 * gap, 0.88,
                                   fc=STATE_COL[s], ec="white", lw=0.5))
            ax.text(CENT[j], i + 0.5, s, ha="center", va="center",
                    fontsize=7.2, color=STATE_TXT[s], weight="bold")
            left += W_COL[j]
    ax.text(_lft / 2.0, n + 1.05,
            "+ supporting   × discordant   – absent / not significant / NA",
            ha="center", va="center", fontsize=7.0, color=C_MID,
            style="italic")


# --------------------------------------------------------------------------- B
def draw_bubble(ax):
    for r in ROWS:
        nl = sum(layer_flags(r))
        ax.scatter(abs(float(r["Parikshak_delta_PSI"])),
                   -np.log10(max(float(r["BH_FDR_KR"]), 1e-4)),
                   s=16 + nl * 13, color=TIER_COLORS[tier_letter(r)],
                   ec="white", lw=0.6, zorder=3, alpha=0.9)
    # Compact gene-only labels reduce redundancy and prevent crowding.
    # Full event IDs remain in Figure 5A/5C and the source-data tables.
    LAB = {"HsaEX0015476": (-10, -10, "right"),
           "HsaEX0029786": (8, -14, "left"),
           "HsaEX0050855": (-4, 12, "center"),
           "HsaEX0051138": (10, 4, "left"),
           "HsaEX0038710": (10, -6, "left")}
    ax.set_xlabel("ASD cortex |ΔPSI|  (bubble size = supporting layers)")
    ax.set_ylabel("KR -log10(BH-FDR)")
    # Final: explicit in-range ticks (auto ticks may place labels
    # outside the axes)
    _xhi = ax.get_xlim()[1]
    ax.set_xticks([v for v in (0.0, 0.05, 0.10, 0.15) if v <= _xhi])
    ax.set_yticks([0.0, 0.5, 1.0, 1.5])
    # Final: labels are placed at their seed offsets and then pushed
    # apart by an iterative display-space repulsion ok.
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    labels = []
    for r in ROWS:
        if r["HsaEX_ID"] not in LAB:
            continue
        dx, dy, ha = LAB[r["HsaEX_ID"]]
        px = abs(float(r["Parikshak_delta_PSI"]))
        py = -np.log10(max(float(r["BH_FDR_KR"]), 1e-4))
        disp = ax.transData.transform((px, py))
        pos = inv.transform((disp[0] + dx, disp[1] + dy))
        labels.append(ax.text(pos[0], pos[1], r["gene"],
                              fontsize=7.5, color=C_DARK, ha=ha,
                              weight="bold" if tier_letter(r) == "A"
                              else "normal"))
    for _ in range(40):
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        bbs = [t.get_window_extent(renderer=rend) for t in labels]
        moved = False
        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                a, b = bbs[i], bbs[j]
                x_ov = min(a.x1, b.x1) - max(a.x0, b.x0)
                y_gap = max(a.y0, b.y0) - min(a.y1, b.y1)
                if x_ov > -1.0 and y_gap < 1.5:
                    push = (1.5 - y_gap) / 2.0 + 0.5
                    hi, lo = (i, j) if a.y0 >= b.y0 else (j, i)
                    for idx, sgn in ((hi, 1.0), (lo, -1.0)):
                        t = labels[idx]
                        xd, yd = t.get_position()
                        dd = ax.transData.transform((xd, yd))
                        dd[1] += sgn * push
                        np_ = inv.transform(dd)
                        t.set_position((np_[0], np_[1]))
                    moved = True
        if not moved:
            break
    th = [plt.Line2D([], [], marker="o", ms=5, color=TIER_COLORS[t], ls="",
                     label="Tier %s" % t) for t in TIER_ORDER]
    nl_vals = sorted(set(sum(layer_flags(r)) for r in MT))
    ex = [nl_vals[0], nl_vals[len(nl_vals) // 2], nl_vals[-1]]
    sh = [ax.scatter([], [], s=16 + nlv * 13, facecolors="none",
                     edgecolors=C_MID, linewidths=0.8,
                     label="%d layer%s" % (nlv, "" if nlv == 1 else "s"))
          for nlv in ex]
    # single combined one-line legend in the B/C-to-D gutter (figure
    # coordinates), so no data point is covered. Mini-ok v2 keeps the
    # same wording/symbols; only the vertical anchor is adjusted so the
    # legend clears the 0.10 x-tick label and the panel-B x-axis title.
    ax.legend(handles=th + sh, loc="center", ncol=7,
              bbox_to_anchor=(0.5, LEGEND_Y), bbox_transform=fig.transFigure,
              fontsize=7.0, columnspacing=1.0, handletextpad=0.5,
              frameon=False)
    panel_title(ax, "Integrated priority map")
    return labels


# --------------------------------------------------------------------------- C
def draw_strips(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(-2, 104)
    panel_title(ax, "Tier A evidence strips", fontsize=9.0)
    ta = [r for r in ROWS if tier_letter(r) == "A"]
    headers = ["Discovery\nΔPSI", "CHyMErA", "Dev class", "PsychENCODE\nβ",
               "KR FDR"]
    # Final: wider panel -> slightly wider cells; headers raised and
    # rows spread for breathing room (pitch 20, cell height 12).
    xs = [13, 26, 51, 74, 86]
    ws = [12, 24, 22, 11, 11]
    for j, h in enumerate(headers):
        C_HEADER_TXT.append(ax.text(xs[j] + ws[j] / 2, 92, h, ha="center",
                                    va="bottom", fontsize=7.4, color=C_DARK,
                                    weight="bold"))
    for i, r in enumerate(ta):
        y = 70 - i * 20
        ax.text(0, y + 3.2, r["gene"], fontsize=8.2, weight="bold",
                color=C_PRIMARY, va="center")
        ax.text(0, y - 3.2, r["HsaEX_ID"].replace("HsaEX", ""), fontsize=7.0,
                color=C_MID, va="center")
        vals = ["%+.3f" % float(r["Parikshak_delta_PSI"]),
                "concordant" if r["CHyMErA_direction_concordant"] == "YES"
                else "direction\nunresolved",
                r["developmental_trajectory"].replace("_", "-"),
                "%+.2f" % float(r["PsychENCODE_beta"]),
                "%.3f" % float(r["BH_FDR_KR"])]
        cols = [C_DISCORD if float(r["Parikshak_delta_PSI"]) > 0
                else C_CONCORD,
                C_CONCORD if r["CHyMErA_direction_concordant"] == "YES"
                else C_NEG,
                C_DEV if r["developmental_dynamic_status"] == "DYNAMIC"
                else C_NEG,
                C_CONCORD if r["direction_concordant"] == "TRUE"
                else C_DISCORD,
                C_PRIMARY]
        for j in range(5):
            p = Rectangle((xs[j], y - 6), ws[j], 12, fc=cols[j],
                          ec="white", lw=0.8)
            ax.add_patch(p)
            t = ax.text(xs[j] + ws[j] / 2, y, vals[j], ha="center",
                        va="center", fontsize=7.0, color="white",
                        weight="bold")
            C_CELLS.append((t, p))


# --------------------------------------------------------------------------- D
# Verified probability-scale facts only (TIER_A_PROTEIN_FEATURES.tsv canonical
# insertion_site_* rows + TIER_A_CODING_CONSEQUENCE.tsv).
D_ROWS = [
    ("CLASP1", "HsaEX0015476", (623, 733), (673, 682),
     [(662, 785, "MAPRE1/3 interaction 662–785"),
      (673, 692, "low complexity 673–692")],
     [(684, "pS684"), (688, "pS688")],
     "canonical isoform includes microexon; net +9 residues"),
    ("HERC4", "HsaEX0029786", (593, 700), (643, 650),
     [(643, 650, "alternative sequence (isoform 2) 643–650")],
     [],
     "canonical isoform includes microexon; net +8 residues"),
    ("PTK2", "HsaEX0050855", (343, 443), (393, 393),
     [(422, 680, "kinase domain 422–680")],
     [(397, "pY397"), (407, "pY407")],
     "canonical isoform excludes microexon; net +6 residues in "
     "protein-coding inclusion transcript contexts"),
    ("PTPRF", "HsaEX0051138", (722, 830), (772, 780),
     [(711, 819, "FNIII-5 domain 711–819"),
      (772, 780, "alternative sequence (isoform 2) 772–780")],
     [],
     "canonical isoform includes microexon; net +9 residues"),
]


def draw_protein(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(-66, 112)
    panel_title(ax, "Tier A local protein context", fontsize=9.0)
    ys = [100, 60, 20, -20]
    for i, (gene, ev, (w0, w1), (i0, i1), regions, ptms, note) in \
            enumerate(D_ROWS):
        y = ys[i]
        x0, x1 = 16.0, 70.0
        span = float(w1 - w0)

        def X(res):
            return x0 + (res - w0) / span * (x1 - x0)

        ax.text(0.5, y + 3, gene, fontsize=8.4, weight="bold",
                color=C_PRIMARY, va="center")
        ax.text(0.5, y - 3.5, ev.replace("HsaEX", ""), fontsize=7.0,
                color=C_MID, va="center")
        # main local-window bar
        ax.add_patch(Rectangle((x0, y - 1.6), x1 - x0, 3.2, fc="#d9d9d9",
                               ec=C_DARK, lw=0.6, zorder=2))
        # microexon insertion site / microexon-derived segment (orange)
        seg_w = max(X(i1) - X(i0), 1.6)
        seg_x = X(i0) if i1 > i0 else X(i0) - seg_w / 2
        ax.add_patch(Rectangle((seg_x, y - 2.4), seg_w, 4.8,
                               fc=MICROEXON_COLOR, ec=C_DARK, lw=0.6,
                               zorder=3))
        D_DRAWN.append((gene, i0, i1, "insertion",
                        "microexon-derived segment / insertion site"))
        # window residue ticks
        for res, xx in [(w0, x0), (w1, x1)]:
            ax.plot([xx, xx], [y - 2.6, y - 4.2], color=C_MID, lw=0.6)
            ax.text(xx, y - 4.6, str(res), fontsize=7.0, color=C_COORD,
                    ha="center", va="top")
        # verified region bars below the window, each with its label line
        for k, (r0, r1, lab) in enumerate(regions):
            yy = y - 12.5 - k * 10.5
            c0, c1 = max(r0, w0), min(r1, w1)
            ax.add_patch(Rectangle((X(c0), yy - 1.3), X(c1) - X(c0), 2.6,
                                   fc=C_DEV, ec="none", alpha=0.85, zorder=2))
            part = "" if (c0, c1) == (r0, r1) else " (part.)"
            ax.text(x0, yy - 3.4, lab + part, fontsize=7.0, color=C_DARK,
                    va="top")
            D_DRAWN.append((gene, r0, r1, "region", lab))
        # verified PTM ticks above; one combined centred label
        if ptms:
            for res, lab in ptms:
                ax.plot([X(res), X(res)], [y + 1.6, y + 4.0], color=C_SENS,
                        lw=1.2, zorder=3)
                D_DRAWN.append((gene, res, res, "ptm", lab))
            mid = sum(X(r) for r, _ in ptms) / len(ptms)
            ax.text(mid, y + 6.5, "  ".join(l for _, l in ptms),
                    fontsize=7.0, color=C_SENS, ha="center", va="center",
                    weight="bold")
    # small legend (two lines)
    ax.add_patch(Rectangle((0, -58), 4, 3, fc=MICROEXON_COLOR, ec=C_DARK,
                           lw=0.6))
    ax.text(5, -56.5, "microexon-derived segment / insertion site",
            fontsize=7.0, color=C_DARK, va="center")
    ax.add_patch(Rectangle((48, -57.4), 4, 1.8, fc=C_DEV, ec="none",
                           alpha=0.85))
    ax.text(53, -56.5, "domain / interaction / alt-seq region",
            fontsize=7.0, color=C_DARK, va="center")
    ax.plot([88, 88], [-58, -55.4], color=C_SENS, lw=1.2)
    ax.text(89.5, -56.5, "PTM", fontsize=7.0, color=C_DARK, va="center")


# ---------------------------------------------------------------------------
# Final: larger canvas (295 vs 270 mm) and a taller B/C/D row; wider
# B-C gutter with the spec B:C width ratio 0.44:0.56.
fig = plt.figure(figsize=(FIG_W, 282 * MM))
gs = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[0.35, 0.28, 0.37],
                       hspace=0.52, left=0.02, right=0.98, top=0.97,
                       bottom=0.055)
axA = fig.add_subplot(gs[0])
draw_matrix(axA)
panel_letter(axA, "A")
gs_mid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1],
                                          width_ratios=[0.45, 0.55],
                                          wspace=0.36)
axB = fig.add_subplot(gs_mid[0])
labelsB = draw_bubble(axB)
panel_letter(axB, "B")
axC = fig.add_subplot(gs_mid[1])
draw_strips(axC)
panel_letter(axC, "C")
axD = fig.add_subplot(gs[2])
draw_protein(axD)
panel_letter(axD, "D")
# Final: coding-consequence prose is moved to the figure legend.
# Keep the main figure data-focused; no text-only footnote block is drawn here.

# ----------------------------------------------------------- QC (computed)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()

# 5A: horizontal headers must not overlap
hb = [t.get_window_extent(renderer=renderer) for t in A_HEADERS]
fig5a = 0
for i in range(len(hb)):
    for j in range(i + 1, len(hb)):
        if overlap_area(hb[i], hb[j]) > 0:
            fig5a += 1

# 5C: every cell label fully inside its cell
fig5c = 0
for t, p in C_CELLS:
    if not containment(t.get_window_extent(renderer=renderer),
                       p.get_window_extent(renderer=renderer), tol=1.0):
        fig5c += 1

# 5D: every drawn annotation must exist in the the analysis tables.
feat_set = set((r["gene"], int(r["feature_start"]), int(r["feature_end"]))
               for r in FEATS)
cod_sites = set()
for r in FEATS:
    a, b = r.get("insertion_site_start"), r.get("insertion_site_end")
    if a not in (None, "", "NA") and b not in (None, "", "NA"):
        cod_sites.add((r["gene"], int(a), int(b)))
fig5d = 0
for gene, s, e, kind, lab in D_DRAWN:
    if kind in ("ptm", "region"):
        if (gene, s, e) not in feat_set:
            fig5d += 1
    elif kind == "insertion":
        if (gene, s, e) not in cod_sites and (gene, s, e) not in feat_set:
            fig5d += 1


# --------------------------------------------------------------- 37B checks
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


PX_PER_MM = fig.dpi / 25.4

# B-C gutter: minimum horizontal separation among vertically-aligned artists
b_ts, c_ts = _ax_texts(axB), _ax_texts(axC)
g_bc = 1e9
for a in b_ts:
    for b in c_ts:
        dy = max(a.y0, b.y0) - min(a.y1, b.y1)
        if dy >= 0:
            continue
        g_bc = min(g_bc, max(a.x0, b.x0) - min(a.x1, b.x1))
FIG5B_C_GUTTER_OK = int(g_bc >= 4.0 * PX_PER_MM)

# 5B label crowding: pairwise overlap OR near-collision (<1.5 px vertical
# separation while horizontally overlapping)
lb = [t.get_window_extent(renderer=renderer) for t in labelsB]
fig5b_crowd = 0
for i in range(len(lb)):
    for j in range(i + 1, len(lb)):
        a, b = lb[i], lb[j]
        x_ov = min(a.x1, b.x1) - max(a.x0, b.x0)
        y_ov = min(a.y1, b.y1) - max(a.y0, b.y0)
        if x_ov > 2.0 and y_ov > 2.0:
            fig5b_crowd += 1
        elif x_ov > 0 and (max(a.y0, b.y0) - min(a.y1, b.y1)) < 1.5:
            fig5b_crowd += 1
# labels must also stay inside the figure
figbb = fig.get_window_extent()
for t in labelsB:
    if not containment(t.get_window_extent(renderer=renderer), figbb,
                       tol=1.0):
        fig5b_crowd += 1

# 5C header breathing: vertical gap header block -> first-row cells >= 2 mm
cb = [p.get_window_extent(renderer=renderer) for _t, p in C_CELLS]
hdr_bbs = [t.get_window_extent(renderer=renderer) for t in C_HEADER_TXT]
first_row_top = max(b.y1 for b in cb[:5])
hdr_bottom = min(h.y0 for h in hdr_bbs)
hdr_gap = hdr_bottom - first_row_top
FIG5C_HEADER_OK = int(hdr_gap >= 2.0 * PX_PER_MM)

# 5C row spacing: vertical gap between consecutive row cell blocks >= 2 mm
row_tops = []
for i in range(4):
    blk = cb[i * 5:(i + 1) * 5]
    row_tops.append((min(b.y0 for b in blk), max(b.y1 for b in blk)))
row_gap = 1e9
for i in range(3):
    row_gap = min(row_gap, row_tops[i][0] - row_tops[i + 1][1])
FIG5C_ROW_OK = int(row_gap >= 2.0 * PX_PER_MM)

checks = [("FIG5A_HEADER_OVERLAP", fig5a), ("FIG5C_CELL_CLIPPING", fig5c),
         ("FIG5D_UNSUPPORTED_ANNOTATION", fig5d),
         ("FIG5B_C_GUTTER_SUFFICIENT", FIG5B_C_GUTTER_OK),
         ("FIG5B_LABEL_CROWDING_ZERO", fig5b_crowd),
         ("FIG5C_HEADER_BREATHING_OK", FIG5C_HEADER_OK),
         ("FIG5C_ROW_SPACING_OK", FIG5C_ROW_OK),
         ("FIG5B_C_GUTTER_PX", round(g_bc, 1)),
         ("FIG5C_HEADER_GAP_PX", round(hdr_gap, 1)),
         ("FIG5C_ROW_GAP_PX", round(row_gap, 1))]

# ---- mini-ok v2: legend vs panel-B tick labels / x-axis title / panel C --
def _real_title_bb(ax):
    # panel_title uses loc="left" -> live artist is _left_title, not .title
    for attr in ("_left_title", "_right_title", "title"):
        bb = getattr(ax, attr).get_window_extent(renderer=renderer)
        if bb.width > 0 and bb.height > 0:
            return bb
    return bb


leg = axB.get_legend()
leg_bb = leg.get_window_extent(renderer=renderer)
xtick_bbs = [t.get_window_extent(renderer=renderer)
             for t in axB.get_xticklabels() if t.get_visible()
             and t.get_text().strip()]
xl_bb = axB.xaxis.label.get_window_extent(renderer=renderer)
c_title_bb = _real_title_bb(axC)
c_ax_bb = axC.get_window_extent(renderer=renderer)
d_title_bb = _real_title_bb(axD)
d_ax_bb = axD.get_window_extent(renderer=renderer)
fig_bb = fig.get_window_extent()
lg_overlap_xticks = int(any(overlap_area(leg_bb, b) > 0 for b in xtick_bbs))
lg_overlap_xlabel = int(overlap_area(leg_bb, xl_bb) > 0)
lg_overlap_c_title = int(overlap_area(leg_bb, c_title_bb) > 0)
lg_overlap_c_panel = int(overlap_area(leg_bb, c_ax_bb) > 0)
lg_overlap_d_title = int(overlap_area(leg_bb, d_title_bb) > 0)
lg_overlap_d_panel = int(overlap_area(leg_bb, d_ax_bb) > 0)
lg_in_figure = int(containment(leg_bb, fig_bb, tol=1.0))
# vertical clearance (px, positive = gap; display coords, y down):
# legend bottom vs x-label top, legend top vs panel-C title bottom,
# panel-D title top vs legend bottom
lg_clear_ticks = round(min(b.y0 for b in xtick_bbs) - leg_bb.y1, 1)
lg_clear_xlabel = round(xl_bb.y0 - leg_bb.y1, 1)
lg_clear_c = round(leg_bb.y0 - c_title_bb.y1, 1)
lg_clear_d = round(d_title_bb.y0 - leg_bb.y1, 1)
leg_checks = [("LEGEND_OVERLAP_B_XTICKS", lg_overlap_xticks),
             ("LEGEND_OVERLAP_B_XLABEL", lg_overlap_xlabel),
             ("LEGEND_OVERLAP_C_TITLE", lg_overlap_c_title),
             ("LEGEND_OVERLAP_C_PANEL", lg_overlap_c_panel),
             ("LEGEND_OVERLAP_D_TITLE", lg_overlap_d_title),
             ("LEGEND_OVERLAP_D_PANEL", lg_overlap_d_panel),
             ("LEGEND_IN_FIGURE", lg_in_figure),
             ("LEGEND_CLEARANCE_TICKS_PX", lg_clear_ticks),
             ("LEGEND_CLEARANCE_XLABEL_PX", lg_clear_xlabel),
             ("LEGEND_CLEARANCE_C_PX", lg_clear_c),
             ("LEGEND_CLEARANCE_D_PX", lg_clear_d),
             ("LEGEND_Y", LEGEND_Y)]
_mf5 = min_text_size_pt(fig)
tsv_write(os.path.join(QC_DIR, "figure5_layout_checks.tsv"),
          ["check", "value"], [[g, v] for g, v in checks + leg_checks] +
          [["MIN_TEXT_SIZE_PT", _mf5],
           ["FIG_WIDTH_MM", round(fig.get_size_inches()[0] * 25.4, 2)]])
print("FIG5 checks:", checks)
print("FIG5 leg checks:", leg_checks)
print("FIG5 min text size (pt):", _mf5)

# ----------------------------------------------------------- render
# FIG5_FINAL_RENDER=0 -> checks-only ok (fast legend-position scans)
NAME = "Figure5_FINAL_SUBMISSION_CLEAN_v2"
if os.environ.get("FIG5_FINAL_RENDER", "1") != "0":
    save_final(fig, FIGDIR, NAME)

    # ---- mini-ok v2: >=200 dpi B/C region preview (visual QA) ----
    from matplotlib.transforms import Bbox as _Bbox
    _region = _Bbox.union([axB.get_window_extent(renderer=renderer),
                           axC.get_window_extent(renderer=renderer),
                           leg_bb]).expanded(1.03, 1.12)
    _inv = fig.transFigure.inverted()
    _fb = _Bbox(_inv.transform([(_region.x0, _region.y0),
                                (_region.x1, _region.y1)]))
    _sz = fig.get_size_inches()
    _bc_inches = _Bbox([(_fb.x0 * _sz[0], _fb.y0 * _sz[1]),
                        (_fb.x1 * _sz[0], _fb.y1 * _sz[1])])
    fig.savefig(os.path.join(FIGDIR, NAME + "_BC_preview_200dpi.png"),
                dpi=200, bbox_inches=_bc_inches)
    print("saved B/C 200-dpi preview")

# ------------------------------------------------- FIG5D source map TSV
maprows = []
for gene, ev, (w0, w1), (i0, i1), regions, ptms, note in D_ROWS:
    maprows.append([gene, ev, "local window", "%d-%d" % (w0, w1),
                    FEAT_TSV, "display window: insertion site +/-50 aa on "
                    "UniProt canonical sequence (uniprot_seqlen context)",
                    "DISPLAY_WINDOW"])
    maprows.append([gene, ev, "microexon-derived segment / insertion site",
                    "%d-%d" % (i0, i1), FEAT_TSV + " ; " + COD_TSV,
                    "insertion_site_start/insertion_site_end; "
                    "residue_start_inclusion/residue_end_inclusion; "
                    "site_state_on_canonical",
                    "DIRECT_DATABASE_ANNOTATION"])
    for r0, r1, lab in regions:
        maprows.append([gene, ev, "region", "%d-%d" % (r0, r1), FEAT_TSV,
                        "feature row (%s)" % lab,
                        "DIRECT_DATABASE_ANNOTATION"])
    for res, lab in ptms:
        maprows.append([gene, ev, "PTM", str(res), FEAT_TSV,
                        "modified residue row (%s)",
                        "DIRECT_DATABASE_ANNOTATION"])
    maprows.append([gene, ev, "context note", note, COD_TSV + " ; " +
                    FEAT_TSV, "coding consequence / site_state_on_canonical",
                    "DIRECT_DATABASE_ANNOTATION"])
tsv_write(os.path.join(FIGDIR, "FIG5D_TIERA_PROTEIN_CONTEXT_SOURCE_MAP.tsv"),
          ["gene", "event_id", "drawn element", "residues", "source_file",
           "source_locator", "evidence_class"], maprows)

# ----------------------------------------------------------- source data TSV
# Rows identical to the final Figure 5_FINAL_SOURCE_DATA.tsv (A matrix /
# B priority map / C Tier A strips / D protein context). The old-5B UpSet
# panel is not part of the final figure and is not emitted here.

rows = []
for r in ROWS:
    states = cell_state(r)
    rows.append(["Figure_5", "A", "evidence matrix",
                 "%s (%s) Tier %s" % (r["gene"], r["HsaEX_ID"],
                                      tier_letter(r)),
                 " ".join("%s=%s" % (c, s)
                           for c, s in zip(COLS_FLAT, states)),
                 MASTER_TSV,
                 "CHyMErA_direction_concordant/developmental_dynamic_status/"
                 "network_module_or_pathway/GSE30573_*/direction_concordant/"
                 "BH_FDR_KR/BH_FDR_LRT"])
for r in ROWS:
    rows.append(["Figure_5", "B", "priority bubble",
                 "%s (%s)" % (r["gene"], r["HsaEX_ID"]),
                 "|dPSI|=%.4f KR_FDR=%.4g layers=%d tier=%s" %
                 (abs(float(r["Parikshak_delta_PSI"])),
                  float(r["BH_FDR_KR"]), sum(layer_flags(r)),
                  tier_letter(r)), MASTER_TSV,
                 "Parikshak_delta_PSI/BH_FDR_KR/layer count"])
for r in ROWS:
    if tier_letter(r) == "A":
        rows.append(["Figure_5", "C", "Tier A strip", r["gene"],
                     "dPSI=%+.3f CHyMErA=%s dev=%s PE_beta=%+.2f KR_FDR=%.3f"
                     % (float(r["Parikshak_delta_PSI"]),
                        r["CHyMErA_direction_concordant"],
                        r["developmental_trajectory"],
                        float(r["PsychENCODE_beta"]),
                        float(r["BH_FDR_KR"])),
                     MASTER_TSV, "Tier A row"])
for gene, ev, (w0, w1), (i0, i1), regions, ptms, note in D_ROWS:
    rows.append(["Figure_5", "D", "protein context", gene,
                 "window=%d-%d site=%d-%d; %s" % (w0, w1, i0, i1, note),
                 FEAT_TSV + " ; " + COD_TSV,
                 "Probability-scale protein mapping (descriptive)"])
tsv_write(os.path.join(FIGDIR, NAME + "_SOURCE_DATA.tsv"),
          ["figure", "panel", "series", "item", "value", "source_file",
           "source_locator"], rows)

prov = [
    ["Figure_5", "A", "Multi-layer evidence matrix",
     "7 evidence layers x 19 events; +/x/- symbols; tier strip; horizontal "
     "short headers (definitions in legend)",
     MASTER_TSV, "classification display only",
     "full-width; horizontal two-line headers replace diagonal unreadable "
     "headers; compact GENE-suffix labels; centred symbols"],
    ["Figure_5", "B", "Integrated priority map",
     "|dPSI| vs KR -log10(BH-FDR); bubble size = supporting layers; Tier A "
     "labels bold + MEF2D retained",
     MASTER_TSV, "KR BH-FDR (primary family)",
     "Final final: gene-only labels for prioritized events; label repulsion; "
     "wider B-C gutter; deterministic in-range ticks"],
    ["Figure_5", "C", "Tier A evidence strips",
     "4 Tier A events; discovery dPSI / CHyMErA / dev class / PsychENCODE "
     "beta / KR FDR; PLPH uppercase; direction wording per publication "
     "vocabulary",
     MASTER_TSV, "per-event values",
     "Final: taller rows (pitch 20, cells 12) + header-to-cell "
     "breathing; rightmost clipping fixed"],
    ["Figure_5", "D", "Tier A local protein context",
     "local +/-50 aa windows with microexon-derived segment / insertion site "
     "(orange), verified domain/interaction/alt-seq regions and PTMs; no "
     "causal wording",
     FEAT_TSV + " ; " + COD_TSV,
     "none (descriptive probability-scale mapping)",
     "replaces text-only old 5E; detailed coding-consequence wording moved "
     "to the figure legend; main panel retains direct annotation only"],
]
tsv_write(os.path.join(FIGDIR, NAME + "_PROVENANCE.tsv"),
          ["figure", "panel", "panel_title", "content", "data_source_files",
           "statistics_displayed", "corrections_vs_earlier"], prov)
print("Figure5_FINAL_v2 done")
