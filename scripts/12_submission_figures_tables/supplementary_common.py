#!/usr/bin/env python3
"""Shared paths, house style and QC helpers.

Final submission redesign of the 16 supplementary figures.
All numeric content is read from upstream analysis sources; NO statistic is
recomputed and NO upstream file is modified.  Supplementary tables are
the final publication tables (untouched).

Captions are rendered with a real paragraph engine (ReportLab Paragraph
with TA_JUSTIFY on the system Helvetica TTF collection); figures keep the
Earlier/38 matplotlib house style.
"""
import os
import csv
import hashlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = os.environ.get("PROJECT_ROOT", ".")
OUT = os.path.join(ROOT,
                   "38D_supplementary_figures_targeted_final_visual_fixes_"
                   "20260812")

# ------------------------------------------------------------ final inputs
DIR25 = os.path.join(ROOT, "25_master_evidence")
DIR29 = os.path.join(ROOT, "29_semantic_terminology")
DIR33 = os.path.join(ROOT,
                     "33_evidence_reconciliation"
                     "")
DIR34 = os.path.join(ROOT, "34_robustness_and_composition")
DIR35 = os.path.join(ROOT,
                     "35_final_manuscript_evidence_and_submission_assets_"
                     "20260806")
DIR36 = os.path.join(ROOT, "36_low_cost_molecular_autism_reinforcement_"
                           "20260807")
DIR37B = os.path.join(ROOT, "37_main_figures")
DIR38 = os.path.join(ROOT, "38_supplementary_package")
DIR38B = os.path.join(ROOT, "38_supplementary_tables")

D29_RENDER = os.path.join(DIR29, "03_render_data")
D35_SUPP_TAB = os.path.join(DIR35, "06_supplementary_tables")
D36_ADJ = os.path.join(DIR36, "02_adjusted_usage_effects")
D36_PROT = os.path.join(DIR36, "03_tierA_protein_mapping")
D36_NM = os.path.join(DIR36, "04_neuron_merged_composition")
D36_D2 = os.path.join(DIR36, "05_D2_transcript_check")
D36_ANC = os.path.join(DIR36, "06_ancestry_audit")
TAB_DIR = os.path.join(OUT, "03_supplementary_tables")
TQC_DIR = os.path.join(OUT, "07_table_qc")
MASTER = os.path.join(DIR25, "06_master_event_table",
                      "MASTER_19_EVENT_EVIDENCE_TABLE.tsv")
TMP = os.environ.get("SCRATCH_ROOT", "/tmp")

# ------------------------------------------------------------ output layout
# Final layout (targeted visual fixes).  Reference inputs above are shared
# and read-only; only these output dirs are written.
INP_DIR = os.path.join(OUT, "01_inputs")
NOTE_DIR = os.path.join(OUT, "02_page_notes")
FIG_DIR = os.path.join(OUT, "03_rebuilt_figures")
PAGE_DIR = os.path.join(OUT, "04_composed_pages")
PREV_DIR = os.path.join(OUT, "05_preview_png")
SCR_DIR = os.path.join(OUT, "06_scripts")
VQC_DIR = os.path.join(OUT, "07_qc")
REP_DIR = os.path.join(OUT, "08_reports")
LEG_DIR = REP_DIR
CK_DIR = os.path.join(OUT, "09_checksums")
# compat aliases used by shared helpers
AUD_DIR = INP_DIR
REG_DIR = REP_DIR
for _d in (INP_DIR, NOTE_DIR, FIG_DIR, PAGE_DIR, PREV_DIR, SCR_DIR,
           VQC_DIR, REP_DIR, CK_DIR):
    os.makedirs(_d, exist_ok=True)

# -------------------------------------------------------------------- style
# identical palette to the earlier figures (colorblind-safe Okabe-Ito derived)
C_PRIMARY = "#0072B2"
C_CONCORD = "#009E73"
C_DEV = "#56B4E9"
C_DISCORD = "#D55E00"
C_NEG = "#999999"
C_SENS = "#CC79A7"
C_ORANGE = "#E69F00"
C_DARK = "#333333"
C_MID = "#4d4d4d"
C_BG_GREY = "#e8e8e8"
MICROEXON_COLOR = C_DISCORD
TIER_COLORS = {"A": C_PRIMARY, "B": C_DEV, "C": C_ORANGE, "D": C_NEG}
TIER_ORDER = ["A", "B", "C", "D"]

MM = 1.0 / 25.4
PT2MM = 25.4 / 72.0
FIG_W = 180 * MM

# page geometry (A4 portrait); captions use 20 mm margins per spec 18-22 mm
PAGE_W_MM, PAGE_H_MM = 210.0, 297.0
ART_MARGIN_MM = 14.0
CAP_MARGIN_MM = 20.0

plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 8.5,
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
    "savefig.dpi": 300,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def panel_letter(ax, letter, x=-0.14, y=1.06, fontsize=10.5):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=fontsize,
            weight="bold", color=C_DARK, va="bottom", ha="left")


def panel_title(ax, text, fontsize=9.0):
    ax.set_title(text, fontsize=fontsize, weight="bold", color=C_DARK,
                 loc="left", pad=6.0)


def save_figure(fig, name):
    """Figure-only deliverables into 02_final_figure_only (600 dpi PNG for
    page composition plus vector PDF)."""
    base = os.path.join(FIG_DIR, name)
    fig.savefig(base + ".pdf", bbox_inches="tight")
    fig.savefig(base + "_image.png", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("saved figure", name)


def tsv_write(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    print("wrote", os.path.basename(path), "(%d rows)" % len(rows))


def rd(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def png_mm(path):
    im = Image.open(path)
    dpi = im.info.get("dpi", (300, 300))
    return im.width / float(dpi[0]) * 25.4, im.height / float(dpi[1]) * 25.4


# ----------------------------------------------- final master table helpers
def load_master():
    return rd(MASTER)


def tier_rank(r):
    t = r["final_evidence_tier"]
    return 0 if t.startswith("TIER_A") else 1 if t.startswith("TIER_B") else \
        2 if t.startswith("TIER_C") else 3


def tier_letter(r):
    return TIER_ORDER[tier_rank(r)]


LAYER_NAMES = ["CHyMErA\nconcordance", "Developmental\ndynamicity",
               "Network\nmembership", "GSE30573\nmapped",
               "PsychENCODE\nconcordance", "KR FDR\n<0.05", "LRT FDR\n<0.05"]


def layer_flags(r):
    return [
        r["CHyMErA_direction_concordant"] == "YES",
        r["developmental_dynamic_status"] == "DYNAMIC",
        r["network_module_or_pathway"] == "SET_LEVEL_NETWORK_MEMBERSHIP",
        r["GSE30573_mapping_status"] == "MAPPED_ANALYZABLE",
        r["direction_concordant"] == "TRUE",
        float(r["BH_FDR_KR"]) < 0.05,
        float(r["BH_FDR_LRT"]) < 0.05,
    ]


# ------------------------------------------------- internal-language scanner
FORBIDDEN_CAPTION_TOKENS = [
    "frozen", "source-verified", "source_verified", "registry", "gate",
    "stage", "master", "audit", "internal", "displayed", "classification",
]
METADATA_MARKERS = ["Elements:", "samples/donors:", "Statistical test:",
                    "Multiple-testing family:", "Effect scale:",
                    "Error bars/CI:", "Direction definition:",
                    "Classification:"]


def tokens(t):
    import re
    return set(x for x in re.split(r"[^A-Za-z0-9]+", t.lower()) if x)


def scan_text(text):
    """Return list of forbidden token hits (case-insensitive tokens)."""
    toks = tokens(text)
    return sorted(toks & set(FORBIDDEN_CAPTION_TOKENS))


# ------------------------------------------------- render-extent QC helpers
def ax_titles(ax):
    out = []
    for t in (ax.title, getattr(ax, "_left_title", None),
              getattr(ax, "_right_title", None)):
        if t is not None and t.get_visible() and t.get_text().strip():
            out.append(t)
    return out


def extent_checks(fig, figname, skip_outside_axes=()):
    """Automated extent checks used while rendering.

    Returns a list of (check, ok, detail) triples:
      TITLE_LEGEND_OVERLAP_ZERO, AXES_OVERLAP_ZERO, TEXT_CLIPPING_ZERO,
      TITLE_TITLE_OVERLAP_ZERO.
    skip_outside_axes: axes whose artists may sit outside the axes bbox
    (axis-off custom panels are checked against the figure bbox only).
    """
    import matplotlib.transforms as mtransforms
    fig.canvas.draw()
    problems = []
    texts = []
    titles = []
    legends = []
    axes_on = []
    for ax in fig.axes:
        if ax.get_visible() and ax.axison:
            axes_on.append(ax)
        for t in ax_titles(ax):
            titles.append((ax, t))
            texts.append(t)
        lg = ax.get_legend()
        if lg is not None and lg.get_visible():
            legends.append((ax, lg))
        if ax.axison:  # axis-off panels draw no tick labels
            for tl in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
                if tl.get_visible() and tl.get_text().strip():
                    texts.append(tl)
        for lab in (ax.xaxis.label, ax.yaxis.label):
            if lab.get_visible() and lab.get_text().strip():
                texts.append(lab)
        for tx in ax.texts:
            if tx.get_visible() and tx.get_text().strip():
                texts.append(tx)
    for t in fig.texts:
        if t.get_visible() and t.get_text().strip():
            texts.append(t)

    fb = mtransforms.Bbox([[0, 0], [fig.get_figwidth(),
                                   fig.get_figheight()]]).transformed(
        fig.dpi_scale_trans)

    def bb(artist):
        return artist.get_window_extent()

    # clipping: any text outside figure bbox
    clip = []
    for t in texts + [t for _, t in titles] :
        b = bb(t)
        if not fb.contains(*b.p0) or not fb.contains(*b.p1):
            # allow 0.5 pt tolerance for descenders at canvas edge
            if b.x0 < -0.7 or b.y0 < -0.7 or b.x1 > fb.x1 + 0.7 or \
                    b.y1 > fb.y1 + 0.7:
                clip.append(t.get_text()[:40])
    # title vs legend overlap
    tl_ov = []
    for _, lg in legends:
        lb = bb(lg)
        for ax, t in titles:
            tb = bb(t)
            if lb.overlaps(tb):
                tl_ov.append(t.get_text()[:40])
    # title vs title overlap
    tt_ov = []
    for i in range(len(titles)):
        for j in range(i + 1, len(titles)):
            bi, bj = bb(titles[i][1]), bb(titles[j][1])
            if bi.overlaps(bj):
                tt_ov.append(titles[i][1].get_text()[:30])
    # axes (panel) overlap
    ao = []
    for i in range(len(axes_on)):
        for j in range(i + 1, len(axes_on)):
            bi = axes_on[i].get_position()
            bj = axes_on[j].get_position()
            if bi.overlaps(bj):
                ao.append((i, j))
    return [
        ("TITLE_LEGEND_OVERLAP_ZERO", len(tl_ov) == 0, ";".join(tl_ov)),
        ("TITLE_TITLE_OVERLAP_ZERO", len(tt_ov) == 0, ";".join(tt_ov)),
        ("AXES_OVERLAP_ZERO", len(ao) == 0, str(ao)),
        ("TEXT_CLIPPING_ZERO", len(clip) == 0, ";".join(clip)),
    ]
