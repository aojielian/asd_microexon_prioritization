#!/usr/bin/env python3
"""Final Figure 1 (final visual-fix revision): cross-species construction +
resource architecture. A cross-species event selection | B representative
cross-species structures | C events across host genes | D final evidence tiers.

Final fixes vs the earlier version:
  - 1B: left text block (mouse/human labels) moved further away from the exon
    graphics; arrow re-centred over the structure and far from all labels;
    GRCh38 coordinates moved fully below the human structure and centred under
    it, with a generous vertical gap from the right flanking exon; the four
    elements (left text block, arrow, structure, coordinate text) now read as
    clearly separated. Whitespace increased; no text shrunk (sizes unchanged
    except the coordinate line which keeps its Earlier size).
  - 1D: the explanatory sentence "one square per event, ordered by tier" is
    excluded from the plot area (kept in legend facts only); panel is now
    minimal and count-focused.
All numbers identical to the final master table and to the Earlier
SOURCE_DATA rows; no recomputation.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import *
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches

MT = load_master()
assert len(MT) == 19
MASTER_TSV = os.path.join(DIR25, "06_master_event_table",
                          "MASTER_19_EVENT_EVIDENCE_TABLE.tsv")
FIGDIR = FIG_DIRS[1]

# QC bookkeeping (populated while drawing; evaluated after render)
NODE_LABELS = []      # (text_artist, box_patch) for containment
ARROWS = []           # FancyArrowPatch list (Fig 1B)
COORD_TEXTS = []      # coordinate Text list (Fig 1B)
EXON_BOXES = []       # exon Rectangle list (Fig 1B)
LEFT_LABELS = []      # mouse/human row labels (Fig 1B)
GENE_LABELS = []      # gene name labels (Fig 1B)

# 1B geometry constants (axes data units, xlim/ylim = 0..100)
X0 = 42.0             # exon-track left edge (Final: 32 -> more left gutter)
TRACK_W = 42.5        # total width of the exon-intron-microexon structure


def bezier_band(ax, x0, y0b, y0t, x1, y1b, y1t, color, alpha=0.5):
    xm = 0.5 * (x0 + x1)
    verts = [(x0, y0b), (xm, y0b), (xm, y1b), (x1, y1b),
             (x1, y1t), (xm, y1t), (xm, y0t), (x0, y0t), (x0, y0b)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(Path(verts, codes), fc=color, ec="none",
                                    alpha=alpha, zorder=1))


def node(ax, x, yb, w, h, text, fc, tc="white", fs=8.2):
    """Rounded node with centred multi-line label; registers for QC."""
    p = FancyBboxPatch((x, yb), w, h,
                       boxstyle="round,pad=0.4,rounding_size=1.2",
                       fc=fc, ec=C_DARK, lw=0.6, zorder=3)
    ax.add_patch(p)
    t = ax.text(x + w / 2, yb + h / 2, text, ha="center", va="center",
                fontsize=fs, color=tc, weight="bold", zorder=4)
    NODE_LABELS.append((t, p))
    return p, t


# --------------------------------------------------------------------------- A
def draw_sankey(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "Cross-species event selection", fontsize=9.0)

    W = 20.0
    x36, x20, x19, x15 = 1.5, 27.0, 52.5, 78.0

    # 36 CHyMErA-designed events (left)
    node(ax, x36, 26, W, 48, "36\nCHyMErA-\ndesigned\nevents", C_MID, fs=8.2)

    # 20 cortex-matched (top-middle), 16 without ASD cortex match (bottom-middle)
    node(ax, x20, 52, W, 30, "20\ncortex-\nmatched\nhuman\nevents", C_DEV,
         tc=C_DARK, fs=8.2)
    node(ax, x20, 14, W, 26, "16\nno ASD\ncortex\nmatch", "#c9c9c9",
         tc=C_MID, fs=8.0)

    # 19 microexon events (top), 1 excluded stub (bottom)
    node(ax, x19, 54, W, 30, "19\nmicroexon\nevents", C_PRIMARY, fs=8.2)
    node(ax, x19, 14, W, 16, "1\nexcluded", "#c9c9c9", tc=C_MID, fs=8.0)
    excl_note = ax.text(x19 + W + 1.2, 22, "no usable\neffect estimate",
                        ha="left", va="center", fontsize=7.6, color=C_MID,
                        style="italic")
    excl_note._s37_role = "excl_note"

    # 15 host genes (right)
    node(ax, x15, 56, W, 26, "15\nhost\ngenes", C_CONCORD, fs=8.4)

    # flow bands
    bezier_band(ax, x36 + W, 26 + 16 / 36 * 48, 26 + 48, x20, 52, 52 + 30, C_DEV, 0.5)
    bezier_band(ax, x36 + W, 26, 26 + 16 / 36 * 48, x20, 14, 14 + 26, "#c9c9c9", 0.5)
    bezier_band(ax, x20 + W, 52 + 19 / 20 * 30, 52 + 30, x19, 54, 54 + 30,
                C_PRIMARY, 0.5)
    bezier_band(ax, x20 + W, 52, 52 + 19 / 20 * 30, x19, 14, 14 + 16,
                "#c9c9c9", 0.5)
    bezier_band(ax, x19 + W, 54, 54 + 30, x15, 56, 56 + 26, C_CONCORD, 0.45)


# --------------------------------------------------------------------------- B
def _exon_track(ax, y, micro_color=MICROEXON_COLOR, x0=X0):
    """Constitutive exon - intron - microexon - intron - constitutive exon."""
    boxes = []
    b1 = Rectangle((x0, y - 2.6), 13, 5.2, fc="#d9d9d9", ec=C_DARK, lw=0.6,
                   zorder=3)
    ax.add_patch(b1)
    boxes.append(b1)
    ax.plot([x0 + 13, x0 + 19], [y, y], color=C_DARK, lw=0.8, zorder=2)
    b2 = Rectangle((x0 + 19, y - 2.6), 4.5, 5.2, fc=micro_color, ec=C_DARK,
                   lw=0.6, zorder=3)
    ax.add_patch(b2)
    boxes.append(b2)
    ax.plot([x0 + 19 + 4.5, x0 + 25 + 4.5], [y, y], color=C_DARK, lw=0.8,
            zorder=2)
    b3 = Rectangle((x0 + 25 + 4.5, y - 2.6), 13, 5.2, fc="#d9d9d9", ec=C_DARK,
                   lw=0.6, zorder=3)
    ax.add_patch(b3)
    boxes.append(b3)
    return boxes


def draw_structure(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "Representative cross-species event structures", fontsize=9.0)

    # gene, mouse event, human event, GRCh38 microexon interval (final master)
    reps = [("ANK3", "MmuEX0004931", "HsaEX0004157",
             "GRCh38 chr10: 60,082,149\u201360,082,176"),
            ("CLASP1", "MmuEX0011506", "HsaEX0015476",
             "GRCh38 chr2: 121,427,403\u2013121,427,430"),
            ("PTK2", "MmuEX0037746", "HsaEX0050855",
             "GRCh38 chr8: 140,769,574\u2013140,769,592")]

    # Final: roomier vertical rhythm (gene / mouse / human / coordinate
    # rows) so that the four elements read as clearly separated.
    y0s = [96, 64, 32]
    x_mid = X0 + TRACK_W / 2.0
    for (gene, mmu, hsa, coord), y0 in zip(reps, y0s):
        y_mouse = y0 - 6
        y_human = y0 - 20
        y_coord = y0 - 30
        gt = ax.text(1.0, y0, gene, fontsize=8.8, weight="bold",
                     color=C_PRIMARY, va="center")
        GENE_LABELS.append(gt)
        # mouse row
        t_m = ax.text(1.0, y_mouse, "mouse " + mmu, fontsize=8.0, color=C_MID,
                      va="center", ha="left")
        LEFT_LABELS.append(t_m)
        EXON_BOXES.extend(_exon_track(ax, y_mouse))
        # human row
        t_h = ax.text(1.0, y_human, "human " + hsa, fontsize=8.0, color=C_MID,
                      va="center", ha="left")
        LEFT_LABELS.append(t_h)
        boxes = _exon_track(ax, y_human, micro_color=MICROEXON_COLOR)
        EXON_BOXES.extend(boxes)
        # vertical arrow between rows, centred over the structure (far from
        # the left text block)
        arr = FancyArrowPatch((x_mid, y_mouse - 3.6), (x_mid, y_human + 3.6),
                              arrowstyle="-|>", mutation_scale=14, lw=1.0,
                              color=C_MID, zorder=2)
        ax.add_patch(arr)
        ARROWS.append(arr)
        # GRCh38 coordinates fully BELOW the human structure, centred under it
        ct = ax.text(x_mid, y_coord, coord, fontsize=7.6, color=C_COORD,
                     va="center", ha="center")
        COORD_TEXTS.append(ct)


# --------------------------------------------------------------------------- C
def draw_gene_events(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "19 events across 15 host genes", fontsize=9.0)
    by_gene = {}
    for r in MT:
        by_gene.setdefault(r["gene"], []).append(r)
    genes = sorted(by_gene, key=lambda g: (-len(by_gene[g]), g))
    n = len(genes)
    row_h = 86.0 / n
    for i, g in enumerate(genes):
        evs = by_gene[g]
        y = 93 - i * row_h
        ax.text(2.0, y, g, fontsize=8.0, va="center", color=C_DARK,
                weight="bold", ha="left")
        for k, r in enumerate(evs):
            x = 34 + k * 8.5
            ax.add_patch(plt.Circle((x, y), 2.7,
                                    fc=TIER_COLORS[tier_letter(r)],
                                    ec="white", lw=0.8, zorder=3))
        ax.text(58.0, y, "%d event%s" % (len(evs), "s" if len(evs) > 1 else ""),
                fontsize=7.6, va="center", color=C_MID, ha="left")
    # tier legend (right column, dedicated area)
    for ti, t in enumerate(TIER_ORDER):
        x = 80
        y = 88 - ti * 9
        ax.add_patch(plt.Circle((x, y), 2.6, fc=TIER_COLORS[t], ec="white",
                                lw=0.8))
        ax.text(x + 5, y, "Tier %s" % t, fontsize=7.8, va="center", color=C_DARK)


# --------------------------------------------------------------------------- D
def draw_waffle(ax):
    ax.axis("off")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    panel_title(ax, "Final evidence tiers", fontsize=9.0)
    counts = {"A": 0, "B": 0, "C": 0, "D": 0}
    for r in MT:
        counts[tier_letter(r)] += 1
    order = []
    for t in TIER_ORDER:
        order += [t] * counts[t]
    n = len(order)
    cw = 92.0 / n
    for i, t in enumerate(order):
        ax.add_patch(Rectangle((4 + i * cw + 0.4, 68), cw - 0.8, 22,
                               fc=TIER_COLORS[t], ec="white", lw=0.8))
    # Final: the explanatory sentence "one square per event, ordered by
    # tier" is not drawn in the plot area (legend facts only; spec 2/Fig 1D).
    for i, t in enumerate(TIER_ORDER):
        x = 4 + (i % 2) * 50
        y = 46 - (i // 2) * 22
        ax.add_patch(Rectangle((x, y - 5), 6.5, 9, fc=TIER_COLORS[t],
                               ec="white"))
        ax.text(x + 9.5, y - 0.5, "Tier %s  n = %d" % (t, counts[t]),
                fontsize=8.2, va="center", color=C_DARK, weight="bold")
    # definitions remain in the figure legend (spec); no in-panel prose


# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(FIG_W, 220 * MM))
outer = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.16, 0.88],
                          hspace=0.36, left=0.015, right=0.985, top=0.975,
                          bottom=0.02)
gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0],
                                          width_ratios=[0.86, 1.14], wspace=0.16)
axA = fig.add_subplot(gs_top[0])
draw_sankey(axA)
panel_letter(axA, "A")
axB = fig.add_subplot(gs_top[1])
draw_structure(axB)
panel_letter(axB, "B")
gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1],
                                          width_ratios=[1.0, 1.05], wspace=0.16)
axC = fig.add_subplot(gs_bot[0])
draw_gene_events(axC)
panel_letter(axC, "C")
axD = fig.add_subplot(gs_bot[1])
draw_waffle(axD)
panel_letter(axD, "D")

# ----------------------------------------------------------- QC (computed)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()

PX_PER_MM = fig.dpi / 25.4  # display px per mm at render dpi

# 1A: every node label fully inside its node box (Earlier check retained)
fig1a_violations = 0
for t, p in NODE_LABELS:
    tbb = t.get_window_extent(renderer=renderer)
    pbb = p.get_window_extent(renderer=renderer)
    if not containment(tbb, pbb, tol=1.0):
        fig1a_violations += 1

texts_B = [(t, bb) for (t, bb) in text_boxes(fig, renderer) if t.axes is axB]

# 1B (Earlier retained): arrows must not collide with any text
fig1b_arrow_violations = 0
for arr in ARROWS:
    abb = arr.get_window_extent(renderer=renderer)
    for t, tbb in texts_B:
        if overlap_area(abb, tbb) > 0:
            fig1b_arrow_violations += 1

# 1B (Earlier retained): coordinate text must not overlap exon boxes
fig1b_coord_violations = 0
for ct in COORD_TEXTS:
    cbb = ct.get_window_extent(renderer=renderer)
    for eb in EXON_BOXES:
        ebb = eb.get_window_extent(renderer=renderer)
        if overlap_area(cbb, ebb) > 0:
            fig1b_coord_violations += 1


def bbox_gap_x(a, b):
    """Horizontal gap (px) between two bboxes (0 if horizontally overlapping)."""
    return max(0.0, max(a.x0, b.x0) - min(a.x1, b.x1))


def bbox_gap_y(a, b):
    """Vertical gap (px) between two bboxes (0 if vertically overlapping)."""
    return max(0.0, max(a.y0, b.y0) - min(a.y1, b.y1))


# Layout check: left text block -> structure horizontal spacing.
# For each mouse/human row label, measure the horizontal gap to the nearest
# exon box of the track on that same row (the leftmost exon is nearest).
left_to_struct_gaps = []
for lab in LEFT_LABELS:
    lbb = lab.get_window_extent(renderer=renderer)
    row_gaps = []
    for eb in EXON_BOXES:
        ebb = eb.get_window_extent(renderer=renderer)
        # same row only: vertical extents overlap or are within one row height
        if bbox_gap_y(lbb, ebb) < 12.0:
            row_gaps.append(bbox_gap_x(lbb, ebb))
    left_to_struct_gaps.append(min(row_gaps))
fig1b_left_gap_px = min(left_to_struct_gaps)
FIG1B_LEFT_TEXT_TO_STRUCTURE_SPACING_OK = int(fig1b_left_gap_px >= 6.0 * PX_PER_MM)

# Layout check: arrow -> nearest label spacing (any text in panel B).
arrow_label_gaps = []
for arr in ARROWS:
    abb = arr.get_window_extent(renderer=renderer)
    dmin = min(max(bbox_gap_x(abb, tbb), bbox_gap_y(abb, tbb))
               for _t, tbb in texts_B)
    arrow_label_gaps.append(dmin)
fig1b_arrow_gap_px = min(arrow_label_gaps)
FIG1B_ARROW_TO_LABEL_SPACING_OK = int(fig1b_arrow_gap_px >= 6.0 * PX_PER_MM)

# Layout check: coordinate text -> nearest exon (right flanking exon) gap.
coord_gaps = []
for ct in COORD_TEXTS:
    cbb = ct.get_window_extent(renderer=renderer)
    dmin = min(max(bbox_gap_x(cbb, eb.get_window_extent(renderer=renderer)),
                   bbox_gap_y(cbb, eb.get_window_extent(renderer=renderer)))
               for eb in EXON_BOXES)
    coord_gaps.append(dmin)
fig1b_coord_gap_px = min(coord_gaps)
FIG1B_RIGHT_COORDINATE_SPACING_OK = int(fig1b_coord_gap_px >= 4.0 * PX_PER_MM)

# Layout check (Fig 1D): explanatory sentence excluded from the plot area.
_d_texts = [t.get_text() for t in axD.texts if t.get_visible()]
FIG1D_EXPLANATORY_SENTENCE_REMOVED = int(
    not any("one square per event" in s.lower() for s in _d_texts))

checks = [("FIG1A_TEXT_OVERLAP", fig1a_violations),
         ("FIG1B_ARROW_TEXT_COLLISION", fig1b_arrow_violations),
         ("FIG1B_COORDINATE_EXON_COLLISION", fig1b_coord_violations),
         ("FIG1B_LEFT_TEXT_TO_STRUCTURE_SPACING_OK",
          FIG1B_LEFT_TEXT_TO_STRUCTURE_SPACING_OK),
         ("FIG1B_ARROW_TO_LABEL_SPACING_OK", FIG1B_ARROW_TO_LABEL_SPACING_OK),
         ("FIG1B_RIGHT_COORDINATE_SPACING_OK",
          FIG1B_RIGHT_COORDINATE_SPACING_OK),
         ("FIG1D_EXPLANATORY_SENTENCE_REMOVED",
          FIG1D_EXPLANATORY_SENTENCE_REMOVED),
         ("FIG1B_LEFT_TEXT_TO_STRUCTURE_MIN_GAP_PX",
          round(fig1b_left_gap_px, 1)),
         ("FIG1B_ARROW_TO_LABEL_MIN_GAP_PX", round(fig1b_arrow_gap_px, 1)),
         ("FIG1B_RIGHT_COORDINATE_MIN_GAP_PX", round(fig1b_coord_gap_px, 1))]
_mf1 = min_text_size_pt(fig)
tsv_write(os.path.join(QC_DIR, "figure1_layout_checks.tsv"),
          ["check", "value"], [[g, v] for g, v in checks] +
          [["MIN_TEXT_SIZE_PT", _mf1],
           ["FIG_WIDTH_MM", round(fig.get_size_inches()[0] * 25.4, 2)]])
print("FIG1 checks:", checks)
print("FIG1 min text size (pt):", _mf1)

# ----------------------------------------------------------- render
NAME = "Figure1_FINAL_SUBMISSION_CLEAN"
save_final(fig, FIGDIR, NAME)

# ----------------------------------------------------------- source data TSV
# Rows are IDENTICAL to the the released Figure1 source data.tsv rows
# (same order, same values) so NO_NUMERIC_CHANGE_FROM_EARLIER can be verified
# by a plain content diff.
by_gene = {}
for r in MT:
    by_gene.setdefault(r["gene"], []).append(r)
counts = {"A": 0, "B": 0, "C": 0, "D": 0}
for r in MT:
    counts[tier_letter(r)] += 1

rows = []
for label, value in [("CHyMErA designed events", 36), ("CHyMErA genes", 31),
                     ("CHyMErA guides", 274), ("CHyMErA cells screened", 21134),
                     ("cortex-matched human events", 20),
                     ("events without ASD cortex match", 16),
                     ("microexon events", 19),
                     ("excluded (no usable effect estimate)", 1), ("host genes", 15)]:
    rows.append(["Figure_1", "A", "selection flow", label, value, MASTER_TSV,
                 "19-event set construction summary"])
for gene, mmu, hsa, coord in [("ANK3", "MmuEX0004931", "HsaEX0004157",
                               "chr10:60082149-60082176"),
                              ("CLASP1", "MmuEX0011506", "HsaEX0015476",
                               "chr2:121427403-121427430"),
                              ("PTK2", "MmuEX0037746", "HsaEX0050855",
                               "chr8:140769574-140769592")]:
    rows.append(["Figure_1", "B", "structure", "%s mouse=%s human=%s" %
                 (gene, mmu, hsa), coord, MASTER_TSV,
                 "representative reciprocal liftOver structures"])
rows.append(["Figure_1", "B", "verification", "reciprocal liftOver verified",
             "19/19", MASTER_TSV, "round-trip verification summary (legend)"])
rows.append(["Figure_1", "B", "verification",
             "coordinate-equivalent local structures (GENCODE v33)", "19/19",
             MASTER_TSV, "coordinate lineage summary (legend)"])
rows.append(["Figure_1", "B", "color", "microexon color word",
             MICROEXON_COLOR_WORD, "this script", "MICROEXON_COLOR_WORD constant"])
for g in sorted(by_gene, key=lambda g: (-len(by_gene[g]), g)):
    rows.append(["Figure_1", "C", "events per gene", g, len(by_gene[g]),
                 MASTER_TSV, "gene column"])
for t in TIER_ORDER:
    rows.append(["Figure_1", "D", "evidence tier count", "Tier %s" % t,
                 counts[t], MASTER_TSV, "final_evidence_tier"])
tsv_write(os.path.join(FIGDIR, NAME + "_SOURCE_DATA.tsv"),
          ["figure", "panel", "series", "item", "value", "source_file",
           "source_locator"], rows)

prov = [
    ["Figure_1", "A", "Cross-species event selection",
     "36->20->19 events across 15 host genes; 16 without ASD match; 1 excluded",
     MASTER_TSV, "counts only; no statistical test",
     "labels concise + centred; bottom prose moved to legend (Earlier)"],
    ["Figure_1", "B", "Representative cross-species event structures",
     "3 exemplar mouse-human pairs; microexon in orange; GRCh38 coords below human structure",
     MASTER_TSV, "verification counts in legend only",
     "Final: left text block moved further from structure; arrow centred "
     "over structure far from labels; coordinates fully below and centred "
     "under the structure with larger vertical gap"],
    ["Figure_1", "C", "19 events across 15 host genes",
     "per-gene event counts, tier-coloured", MASTER_TSV, "counts only",
     "alignment + legend spacing improved (Earlier)"],
    ["Figure_1", "D", "Final evidence tiers (waffle)",
     "Tier A=4, B=3, C=8, D=4", MASTER_TSV, "counts only",
     "Final: explanatory sentence 'one square per event, ordered by tier' "
     "excluded from plot area (legend facts only)"],
]
tsv_write(os.path.join(FIGDIR, NAME + "_PROVENANCE.tsv"),
          ["figure", "panel", "panel_title", "content", "data_source_files",
           "statistics_displayed", "corrections_vs_earlier"], prov)
print("Figure1_FINAL_v2 done")
