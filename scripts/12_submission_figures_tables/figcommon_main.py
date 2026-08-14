#!/usr/bin/env python3
"""Final shared figure style + final-event data loading + overlap-QC utilities.

Adapted from the earlier figure-common module for the final version
visual-fix revision. All plotted numbers are read from the SAME upstream
sources (dirs 11/13/14/17/21/25/33) and the reference outputs that Earlier
used; NO statistic is recomputed and NO upstream file is modified.

Differences vs figcommon_main:
  - outputs go to dir 37B (per-figure directories 01_figure1 .. 05_figure5)
  - final deliverables named FigureN_FINAL_v2.{pdf,tiff} + journal-width/50% PNGs
  - visual QC dir is 06_visual_qc; report dir is 10_final_report
  - Earlier outputs (DIR37) are READ-ONLY inputs (source-data diff baseline)
"""
import os
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np  # noqa: F401  (kept available for figure scripts)

ROOT = os.environ.get("PROJECT_ROOT", ".")
OUT = os.path.join(ROOT, "37_main_figures")
DIR37 = os.path.join(ROOT, "37_main_figures")  # read-only
DIR25 = os.path.join(ROOT, "25_master_evidence")
DIR11 = os.path.join(ROOT, "11_set_level_enrichment")
DIR13 = os.path.join(ROOT, "13_developmental_timing_repair")
DIR14 = os.path.join(ROOT, "14_mechanistic_context")
DIR17 = os.path.join(ROOT, "17_gse30573_mapping")
DIR21 = os.path.join(ROOT, "21_coordinate_inference")
DIR29 = os.path.join(ROOT, "29_semantic_terminology")
DIR33 = os.path.join(ROOT, "33_evidence_reconciliation")
DIR36 = os.path.join(ROOT, "36_probability_scale_and_protein")

FIG_DIRS = {
    1: os.path.join(OUT, "01_figure1"),
    2: os.path.join(OUT, "02_figure2"),
    3: os.path.join(OUT, "03_figure3"),
    4: os.path.join(OUT, "04_figure4"),
    5: os.path.join(OUT, "05_figure5"),
}
QC_DIR = os.path.join(OUT, "06_visual_qc")
LEGEND_DIR = os.path.join(OUT, "07_legend_facts")
SCRIPT_DIR = os.path.join(OUT, "08_scripts")
CK_DIR = os.path.join(OUT, "09_checksums")
REP_DIR = os.path.join(OUT, "10_final_report")
NULL_NPZ = os.path.join(DIR29, "03_render_data", "NULL_DISTRIBUTIONS.npz")  # read-only
for _d in list(FIG_DIRS.values()) + [QC_DIR, LEGEND_DIR, CK_DIR, REP_DIR]:
    os.makedirs(_d, exist_ok=True)

# ---------------------------------------------------------------------- style
# Semantic, colorblind-safe palette (Okabe-Ito derived; identical to dir 29/35/37)
C_PRIMARY = "#0072B2"   # deep blue   : primary / Tier A / strong support
C_CONCORD = "#009E73"   # teal-green  : concordant direction
C_DEV     = "#56B4E9"   # blue-green  : developmental maturation
C_DISCORD = "#D55E00"   # orange-red  : discordant ; ALSO microexon fill (orange)
C_NEG     = "#999999"   # grey        : negative / absent / not direction-eligible
C_SENS    = "#CC79A7"   # purple      : sensitivity-only
C_ORANGE  = "#E69F00"   # orange      : PHPL / Tier C
C_DARK    = "#333333"
C_MID     = "#4d4d4d"
C_LIGHT   = "#f2f2f2"
C_BG_GREY = "#e8e8e8"
C_COORD   = "#8a8a8a"   # small grey for genomic coordinates (Fig 1B)

# The microexon structure is visually ORANGE (C_DISCORD = #D55E00, an orange-red).
MICROEXON_COLOR = C_DISCORD
MICROEXON_COLOR_WORD = "orange"   # caption/legend must say orange, never red

TIER_COLORS = {"A": C_PRIMARY, "B": C_DEV, "C": C_ORANGE, "D": C_NEG}
TIER_ORDER = ["A", "B", "C", "D"]

MM = 1.0 / 25.4       # mm -> inch
FIG_W = 180 * MM      # two-column journal width (180 mm)

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 8.5,             # raised from 7.5 for final-width legibility
    "axes.edgecolor": C_DARK,
    "axes.linewidth": 0.7,
    "axes.titlesize": 9.0,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.5,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 7.8,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "svg.fonttype": "none",       # keep text as editable text (QC-greppable)
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def panel_letter(ax, letter, x=-0.14, y=1.06, fontsize=10.5):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=fontsize,
            weight="bold", color=C_DARK, va="bottom", ha="left")


def panel_title(ax, text, fontsize=9.0):
    ax.set_title(text, fontsize=fontsize, weight="bold", color=C_DARK, loc="left")


def footnote(fig, text, y=0.008, fontsize=7.2):
    fig.text(0.5, y, text, ha="center", va="bottom", fontsize=fontsize,
             style="italic", color=C_MID)


def save_final(fig, directory, name):
    """Write the journal deliverables + previews for one figure.

    name should be e.g. ``Figure1_FINAL_v2``. Produces:
      - name.pdf                vector (submission)
      - name.tiff               600 dpi LZW (submission)
      - name_finalwidth.png     300 dpi, journal-width preview
      - name_50pct.png          half-size on-screen preview
      - name.svg                editable-text vector (for QC token grepping)
    """
    base = os.path.join(directory, name)
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + ".svg", bbox_inches="tight")
    tmp = base + "_tmp600.png"
    fig.savefig(tmp, dpi=600, bbox_inches="tight")
    Image.open(tmp).convert("RGB").save(base + ".tiff",
                                        compression="tiff_lzw", dpi=(600, 600))
    os.remove(tmp)
    fig.savefig(base + "_finalwidth.png", dpi=300, bbox_inches="tight")
    im = Image.open(base + "_finalwidth.png")
    im.resize((max(1, im.width // 2), max(1, im.height // 2)),
              Image.LANCZOS).save(base + "_50pct.png", dpi=(300, 300))
    plt.close(fig)
    print("saved", name, "->", directory)


def tsv_write(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    print("wrote", os.path.basename(path), "(%d rows)" % len(rows))


# ---------------------------------------------------------------------- io
def rd(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_master():
    return rd(os.path.join(DIR25, "06_master_event_table",
                           "MASTER_19_EVENT_EVIDENCE_TABLE.tsv"))


def tier_rank(r):
    t = r["final_evidence_tier"]
    return 0 if t.startswith("TIER_A") else 1 if t.startswith("TIER_B") else \
           2 if t.startswith("TIER_C") else 3


def tier_letter(r):
    return TIER_ORDER[tier_rank(r)]


def devrank(r):
    return {"PLPH": 0, "PHPL": 1, "NON_DYNAMIC": 2}.get(
        r["developmental_trajectory"], 2)


def ev_label(r, full=True):
    return "%s (%s)" % (r["gene"], r["HsaEX_ID"]) if full else \
           "%s %s" % (r["gene"], r["HsaEX_ID"].replace("HsaEX", ""))


def compact_label(r):
    """Compact event-label convention: GENE-event_suffix (e.g. CLASP1-15476)."""
    return "%s-%s" % (r["gene"], r["HsaEX_ID"].replace("HsaEX", ""))


# Evidence-layer helpers (final 19-event definitions, master table)
LAYER_NAMES = ["CHyMErA\nconcordance", "Developmental\ndynamicity",
               "Network\nmembership", "GSE30573\nmapped",
               "PsychENCODE\nconcordance", "KR FDR\n<0.05", "LRT FDR\n<0.05"]


def layer_flags(r):
    """7 evidence layers for one event, from the final master table."""
    return [
        r["CHyMErA_direction_concordant"] == "YES",
        r["developmental_dynamic_status"] == "DYNAMIC",
        r["network_module_or_pathway"] == "SET_LEVEL_NETWORK_MEMBERSHIP",
        r["GSE30573_mapping_status"] == "MAPPED_ANALYZABLE",
        r["direction_concordant"] == "TRUE",
        float(r["BH_FDR_KR"]) < 0.05,
        float(r["BH_FDR_LRT"]) < 0.05,
    ]


def load_effects():
    """Reference enrichment effects per background (dir 11)."""
    rows = rd(os.path.join(DIR11, "07_primary_reanalysis",
                           "01_effects_by_background.tsv"))
    order = ["BG0_WIDE_SE", "BG1_MICROEXON", "BG2_CONSERVED_MICROEXON",
             "BG3_CEM", "BG3_NN"]
    lab = {"BG0_WIDE_SE": "Wide splicing-event", "BG1_MICROEXON": "Microexon",
           "BG2_CONSERVED_MICROEXON": "Conserved microexon",
           "BG3_CEM": "CEM (matched)", "BG3_NN": "NN (matched)"}
    d = {r["background"]: r for r in rows}
    return [dict(bg=b, label=lab[b],
                 effect=float(d[b]["effect_mean_difference"]),
                 lo=float(d[b]["bootstrap_95CI_lower"]),
                 hi=float(d[b]["bootstrap_95CI_upper"]),
                 perm_p=float(d[b]["permutation_p"]),
                 target_mean=float(d[b]["target_mean_abs_dpsi"]),
                 bg_mean=float(d[b]["background_mean_abs_dpsi"]),
                 n_bg=int(d[b]["n_background"])) for b in order]


def load_network_edges():
    return rd(os.path.join(DIR14, "08_host_gene_network",
                           "03_observed_network_edges.tsv"))


def load_pathways():
    rows = rd(os.path.join(DIR14, "09_pathway_permutation",
                           "03_matched_permutation_results.tsv"))
    lab = {"synaptic_signaling": "Synaptic signaling",
           "neuron_projection": "Neuron projection",
           "chromatin_transcription": "Chromatin regulation",
           "cytoskeleton_organization": "Cytoskeleton",
           "cell_adhesion": "Cell adhesion",
           "axon_guidance": "Axon guidance",
           "protein_localization": "Protein localization",
           "vesicle_transport": "Vesicle trafficking",
           "calcium_signaling": "Calcium signaling",
           "ubiquitin_proteasome": "Ubiquitin/proteasome"}
    order = ["synaptic_signaling", "neuron_projection", "chromatin_transcription",
             "cytoskeleton_organization", "cell_adhesion", "axon_guidance",
             "protein_localization", "vesicle_transport", "calcium_signaling",
             "ubiquitin_proteasome"]
    d = {r["pathway"]: r for r in rows}
    return [dict(key=k, label=lab[k],
                 genes=d[k]["overlap_genes"].split(","),
                 n=int(d[k]["n_overlap"]),
                 perm_p=float(d[k]["permutation_p"])) for k in order]


def load_rbp_tests():
    """All 240 RBP motif permutation tests."""
    return rd(os.path.join(DIR14, "15_reports", "rbp_motif_results.tsv"))


# ============================================================================
# Programmatic overlap / clipping QC (hard checks are COMPUTED, never hand-set)
# ============================================================================
def _draw(fig):
    fig.canvas.draw()
    return fig.canvas.get_renderer()


def text_boxes(fig, renderer, min_len=1):
    """Display-space bounding boxes for every visible Text artist."""
    out = []
    for ax in fig.axes:
        for t in ax.texts:
            s = t.get_text()
            if s is None or len(s.strip()) < min_len:
                continue
            if not t.get_visible():
                continue
            out.append((t, t.get_window_extent(renderer=renderer)))
    return out


def patch_boxes(fig, renderer, types=None):
    """Display-space bounding boxes for patches (optionally filtered by type)."""
    out = []
    for ax in fig.axes:
        for p in ax.patches:
            if types is not None and not isinstance(p, types):
                continue
            if not p.get_visible():
                continue
            try:
                bb = p.get_window_extent(renderer=renderer)
            except Exception:
                continue
            out.append((p, bb))
    return out


def overlap_area(a, b):
    """Intersection area (px^2) of two Bbox; 0 if disjoint."""
    x0 = max(a.x0, b.x0)
    y0 = max(a.y0, b.y0)
    x1 = min(a.x1, b.x1)
    y1 = min(a.y1, b.y1)
    if x1 <= x0 or y1 <= y0:
        return 0.0
    return (x1 - x0) * (y1 - y0)


def containment(child, parent, tol=1.5):
    """True if child bbox lies inside parent bbox (within tol px)."""
    return (child.x0 >= parent.x0 - tol and child.x1 <= parent.x1 + tol and
            child.y0 >= parent.y0 - tol and child.y1 <= parent.y1 + tol)


def pairwise_text_overlap(texts, ignore_px=2.0):
    """Count pairs of text boxes whose overlap exceeds ignore_px in both dims.

    ``texts`` is a list of (artist, bbox). Adjacent labels that merely touch are
    not counted; only real area overlap above a small pixel tolerance counts.
    """
    n = 0
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = texts[i][1], texts[j][1]
            x_ov = min(a.x1, b.x1) - max(a.x0, b.x0)
            y_ov = min(a.y1, b.y1) - max(a.y0, b.y0)
            if x_ov > ignore_px and y_ov > ignore_px:
                n += 1
    return n


def figure_px(fig, renderer, pts):
    """Convert a size in points to display pixels for the given figure."""
    return pts * fig.dpi / 72.0


def min_text_size_pt(fig):
    """Smallest font size (pt) among visible Text artists in the figure."""
    sizes = []
    for ax in fig.axes:
        for t in ax.texts:
            if t.get_visible() and t.get_text() and t.get_text().strip():
                sizes.append(t.get_fontsize())
    return min(sizes) if sizes else None
