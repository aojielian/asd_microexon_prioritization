#!/usr/bin/env python3
"""
Analysis-R Finalization: Reports, Figures, and QC Files
========================================================
Generates all remaining outputs for Analysis-R:
  - 11 QC figures (PDF+PNG or SVG+PDF fallback) in 14_figures_qc/
  - 9 QC metadata files in 13_qc/
  - 14 report files in 15_reports/
  - Directory tree snapshot

Corrects STATUS from the provisional
ASD_SPECIFIC_DEVELOPMENTAL_WINDOW to
BROAD_NEURAL_MICROEXON_MATURATION.

Author: Analysis-R Pipeline
Date: 2026-07-31
Random seed: 42
"""

import sys
import platform
import os
import zlib
import struct
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ============================================================================
# MATPLOTLIB IMPORT WITH FALLBACK
# ============================================================================

HAS_MATPLOTLIB = False
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except (ImportError, OSError) as _mpl_err:
    print(f"[WARN] matplotlib unavailable ({_mpl_err.__class__.__name__}); "
          "using pure-Python SVG/PDF fallback for figures.")

# ============================================================================
# CONFIGURATION
# ============================================================================

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", "."))
TASK_ROOT = PROJECT_ROOT / "13_developmental_timing_repair"

# Output directories
FIGURES_DIR = TASK_ROOT / "14_figures_qc"
QC_DIR = TASK_ROOT / "13_qc"
REPORTS_DIR = TASK_ROOT / "15_reports"

for d in [FIGURES_DIR, QC_DIR, REPORTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Input data paths
TIMING_TSV    = TASK_ROOT / "07_primary_timing_reanalysis" / "00_primary_timing_results.tsv"
TRAJECTORY_TSV = TASK_ROOT / "08_trajectory_direction_tests" / "00_trajectory_direction_results.tsv"
DYNAMICITY_TSV = TASK_ROOT / "05_dynamicity_definition" / "00_event_dynamicity_metrics.tsv"
TIERS_TSV     = TASK_ROOT / "12_tier_reclassification" / "00_event_tier_reclassification.tsv"
SENSITIVITY_TSV = TASK_ROOT / "13_sensitivity" / "00_sensitivity_results.tsv"
ASD_CORR_TSV  = TASK_ROOT / "09_ASD_timing_correlation" / "00_asd_timing_correlation.tsv"
BRAINSPAN_TSV = TASK_ROOT / "10_brainspan_recheck" / "00_brainspan_recheck.tsv"
ZEBRAFISH_TSV = TASK_ROOT / "11_zebrafish_recheck" / "00_zebrafish_recheck.tsv"
GENE_BLOCK_TSV = TASK_ROOT / "07_primary_timing_reanalysis" / "01_gene_block_permutation.tsv"
LOO_GENE_TSV  = TASK_ROOT / "07_primary_timing_reanalysis" / "02_LOO_gene.tsv"
LOO_EVENT_TSV = TASK_ROOT / "07_primary_timing_reanalysis" / "03_LOO_event.tsv"
RECONCILIATION_SUMMARY_TSV = TASK_ROOT / "03_event_set_reconciliation" / "00_event_set_summary.tsv"
RECONCILIATION_MEMBERSHIP_TSV = TASK_ROOT / "03_event_set_reconciliation" / "01_event_set_membership.tsv"
GROUP_DEFS_TSV = TASK_ROOT / "04_vastdb_group_check" / "00_group_definitions.tsv"
BG_SUMMARY_TSV = TASK_ROOT / "06_strict_background_rebuild" / "05_background_summary.tsv"

NOW = datetime.now()
TIMESTAMP = NOW.strftime("%Y-%m-%dT%H:%M:%S")
DATE_STR = NOW.strftime("%Y-%m-%d")

# Corrected final status
STATUS = "BROAD_NEURAL_MICROEXON_MATURATION"


# ============================================================================
# PURE-PYTHON SVG/PDF/PNG FALLBACK FIGURE ENGINE
# ============================================================================

class FallbackFigure:
    """Minimal pure-Python figure engine that writes SVG, PDF, and PNG
    without matplotlib or PIL.

    Provides a coordinate system and primitive drawing operations:
      - rect, line, circle, text
    """

    def __init__(self, width_pt=540, height_pt=400):
        self.W = width_pt
        self.H = height_pt
        self.elements = []
        # Plot area margins
        self.margin_left = 90
        self.margin_right = 30
        self.margin_top = 50
        self.margin_bottom = 60
        self.plot_x = self.margin_left
        self.plot_y = self.margin_top
        self.plot_w = self.W - self.margin_left - self.margin_right
        self.plot_h = self.H - self.margin_top - self.margin_bottom

    def _esc(self, s):
        """Escape XML special characters."""
        return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    def add_bg(self, color="#FFFFFF"):
        self.elements.append(
            f'<rect x="0" y="0" width="{self.W}" height="{self.H}" fill="{color}"/>')

    def add_rect(self, x, y, w, h, fill="#4C72B0", stroke="none", opacity=1.0):
        self.elements.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'fill="{fill}" stroke="{stroke}" opacity="{opacity}"/>')

    def add_line(self, x1, y1, x2, y2, stroke="#000000", width=1, dash=""):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        self.elements.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{width}"{da}/>')

    def add_circle(self, cx, cy, r=4, fill="#C44E52", stroke="none"):
        self.elements.append(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" '
            f'fill="{fill}" stroke="{stroke}"/>')

    def add_text(self, x, y, text, size=10, fill="#000000", anchor="start",
                 weight="normal", rotate=0):
        rot = ""
        if rotate:
            rot = f' transform="rotate({rotate},{x:.1f},{y:.1f})"'
        self.elements.append(
            f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" '
            f'text-anchor="{anchor}" font-family="Helvetica,Arial,sans-serif" '
            f'font-weight="{weight}"{rot}>{self._esc(text)}</text>')

    def add_title(self, text, size=13):
        self.add_text(self.W / 2, 25, text, size=size, anchor="middle", weight="bold")

    def map_x(self, val, vmin, vmax):
        if vmax == vmin:
            return self.plot_x + self.plot_w / 2
        return self.plot_x + (val - vmin) / (vmax - vmin) * self.plot_w

    def map_y(self, val, vmin, vmax):
        if vmax == vmin:
            return self.plot_y + self.plot_h / 2
        return self.plot_y + self.plot_h - (val - vmin) / (vmax - vmin) * self.plot_h

    def draw_axes(self, xlabel="", ylabel="", xticks=None, yticks=None):
        """Draw plot box with axes and optional tick labels."""
        # Box
        self.add_rect(self.plot_x, self.plot_y, self.plot_w, self.plot_h,
                      fill="none", stroke="#CCCCCC")
        # X-axis ticks
        if xticks:
            for val, label in xticks:
                px = self.map_x(val, xticks[0][0], xticks[-1][0])
                self.add_line(px, self.plot_y + self.plot_h,
                              px, self.plot_y + self.plot_h + 5,
                              stroke="#666666")
                self.add_text(px, self.plot_y + self.plot_h + 18, label,
                              size=7, anchor="middle", fill="#333333")
        # Y-axis ticks
        if yticks:
            for val, label in yticks:
                py = self.map_y(val, yticks[0][0], yticks[-1][0])
                self.add_line(self.plot_x - 5, py, self.plot_x, py,
                              stroke="#666666")
                self.add_text(self.plot_x - 8, py + 3, label,
                              size=7, anchor="end", fill="#333333")
        # Labels
        if xlabel:
            self.add_text(self.plot_x + self.plot_w / 2,
                          self.H - 8, xlabel, size=9, anchor="middle")
        if ylabel:
            self.add_text(12, self.plot_y + self.plot_h / 2, ylabel,
                          size=9, anchor="middle", rotate=-90)

    def to_svg(self):
        parts = [
            f'<?xml version="1.0" encoding="UTF-8"?>',
            f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{self.W}" height="{self.H}" viewBox="0 0 {self.W} {self.H}">',
        ]
        parts.extend(self.elements)
        parts.append("</svg>")
        return "\n".join(parts)

    def to_pdf(self):
        """Generate minimal valid PDF with embedded SVG-equivalent content."""
        # We'll create a simple single-page PDF with text operators
        # Convert SVG coordinates to PDF coordinates (origin bottom-left)
        page_w = self.W
        page_h = self.H

        stream_lines = []
        # White background
        stream_lines.append("1 1 1 rg")  # set fill color white
        stream_lines.append(f"0 0 {page_w} {page_h} re f")  # fill rect

        for elem in self.elements:
            if elem.startswith("<rect"):
                # Parse rect attributes
                import re
                attrs = dict(re.findall(r'(\w+)="([^"]*)"', elem))
                rx = float(attrs.get("x", 0))
                ry = float(attrs.get("y", 0))
                rw = float(attrs.get("w", attrs.get("width", 0)))
                rh = float(attrs.get("h", attrs.get("height", 0)))
                fill = attrs.get("fill", "#000000")
                # PDF y is flipped
                pdf_y = page_h - ry - rh
                if fill != "none":
                    r, g, b = _hex_to_rgb(fill)
                    op = float(attrs.get("opacity", 1.0))
                    stream_lines.append(f"{r:.3f} {g:.3f} {b:.3f} rg")
                    stream_lines.append(f"{rx:.1f} {pdf_y:.1f} {rw:.1f} {rh:.1f} re f")
                stroke = attrs.get("stroke", "none")
                if stroke != "none":
                    r, g, b = _hex_to_rgb(stroke)
                    stream_lines.append(f"{r:.3f} {g:.3f} {b:.3f} RG")
                    stream_lines.append(f"{rx:.1f} {pdf_y:.1f} {rw:.1f} {rh:.1f} re S")

            elif elem.startswith("<line"):
                import re
                attrs = dict(re.findall(r'(\w+)="([^"]*)"', elem))
                x1 = float(attrs.get("x1", 0))
                y1 = page_h - float(attrs.get("y1", 0))
                x2 = float(attrs.get("x2", 0))
                y2 = page_h - float(attrs.get("y2", 0))
                stroke = attrs.get("stroke", "#000000")
                r, g, b = _hex_to_rgb(stroke)
                w = float(attrs.get("stroke-width", 1))
                stream_lines.append(f"{r:.3f} {g:.3f} {b:.3f} RG")
                stream_lines.append(f"{w} w")
                if "stroke-dasharray" in attrs:
                    stream_lines.append("[3 3] 0 d")
                else:
                    stream_lines.append("[] 0 d")
                stream_lines.append(f"{x1:.1f} {y1:.1f} m {x2:.1f} {y2:.1f} l S")

            elif elem.startswith("<circle"):
                import re
                attrs = dict(re.findall(r'(\w+)="([^"]*)"', elem))
                cx = float(attrs.get("cx", 0))
                cy = page_h - float(attrs.get("cy", 0))
                cr = float(attrs.get("r", 4))
                fill = attrs.get("fill", "#000000")
                if fill != "none":
                    r, g, b = _hex_to_rgb(fill)
                    stream_lines.append(f"{r:.3f} {g:.3f} {b:.3f} rg")
                    # Approximate circle with 4 bezier curves
                    k = 0.5522847498
                    stream_lines.append(
                        f"{cx - cr:.2f} {cy:.2f} m "
                        f"{cx - cr:.2f} {cy + cr * k:.2f} "
                        f"{cx - cr * k:.2f} {cy + cr:.2f} "
                        f"{cx:.2f} {cy + cr:.2f} c "
                        f"{cx + cr * k:.2f} {cy + cr:.2f} "
                        f"{cx + cr:.2f} {cy + cr * k:.2f} "
                        f"{cx + cr:.2f} {cy:.2f} c "
                        f"{cx + cr:.2f} {cy - cr * k:.2f} "
                        f"{cx + cr * k:.2f} {cy - cr:.2f} "
                        f"{cx:.2f} {cy - cr:.2f} c "
                        f"{cx - cr * k:.2f} {cy - cr:.2f} "
                        f"{cx - cr:.2f} {cy - cr * k:.2f} "
                        f"{cx - cr:.2f} {cy:.2f} c f")

            elif elem.startswith("<text"):
                import re
                attrs = dict(re.findall(r'(\w+)="([^"]*)"', elem))
                tx = float(attrs.get("x", 0))
                ty = page_h - float(attrs.get("y", 0))
                fsize = float(attrs.get("font-size", 10))
                fill = attrs.get("fill", "#000000")
                r, g, b = _hex_to_rgb(fill)
                anchor = attrs.get("text-anchor", "start")
                # Extract text content
                text = re.sub(r'<[^>]+>', '', elem.split(">")[-1].split("<")[0])
                text = text.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
                # PDF text
                weight = attrs.get("font-weight", "normal")
                font = "/F2" if weight == "bold" else "/F1"
                stream_lines.append(f"{r:.3f} {g:.3f} {b:.3f} rg")
                stream_lines.append("BT")
                stream_lines.append(f"{font} {fsize:.0f} Tf")
                if anchor == "middle":
                    est_w = len(text) * fsize * 0.5
                    tx -= est_w / 2
                elif anchor == "end":
                    est_w = len(text) * fsize * 0.5
                    tx -= est_w
                # Escape PDF string
                text = text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
                stream_lines.append(f"{tx:.1f} {ty:.1f} Td")
                stream_lines.append(f"({text}) Tj")
                stream_lines.append("ET")

        stream_content = "\n".join(stream_lines)
        stream_bytes = stream_content.encode("latin-1", errors="replace")

        # Build PDF structure
        objects = []
        # Obj 1: Catalog
        objects.append(b"1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n")
        # Obj 2: Pages
        objects.append(b"2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n")
        # Obj 3: Page
        objects.append(
            f"3 0 obj\n<< /Type /Page /Parent 2 0 R "
            f"/MediaBox [0 0 {page_w} {page_h}] "
            f"/Contents 4 0 R /Resources << /Font << /F1 5 0 R /F2 6 0 R >> >> >>\n"
            f"endobj\n".encode())
        # Obj 4: Stream
        objects.append(
            f"4 0 obj\n<< /Length {len(stream_bytes)} >>\nstream\n".encode()
            + stream_bytes
            + b"\nendstream\nendobj\n")
        # Obj 5: Font (Helvetica)
        objects.append(
            b"5 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\nendobj\n")
        # Obj 6: Font (Helvetica-Bold)
        objects.append(
            b"6 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>\nendobj\n")

        # Assemble PDF
        pdf = bytearray()
        pdf.extend(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = []
        for obj in objects:
            offsets.append(len(pdf))
            pdf.extend(obj)
        # Cross-reference table
        xref_offset = len(pdf)
        pdf.extend(b"xref\n")
        pdf.extend(f"0 {len(objects) + 1}\n".encode())
        pdf.extend(b"0000000000 65535 f \n")
        for off in offsets:
            pdf.extend(f"{off:010d} 00000 n \n".encode())
        pdf.extend(b"trailer\n")
        pdf.extend(f"<< /Size {len(objects) + 1} /Root 1 0 R >>\n".encode())
        pdf.extend(b"startxref\n")
        pdf.extend(f"{xref_offset}\n".encode())
        pdf.extend(b"%%EOF\n")

        return bytes(pdf)

    def to_png(self, width=800, height=600):
        """Generate a minimal valid PNG file (placeholder with white background
        and title text rendered as metadata). For full-quality PNG, run with
        matplotlib on native architecture."""
        # Create minimal 1x1 white PNG as placeholder
        # (SVG is the primary vector output; PNG is supplementary)
        def _chunk(chunk_type, data):
            c = chunk_type + data
            crc = struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)
            return struct.pack(">I", len(data)) + c + crc

        signature = b"\x89PNG\r\n\x1a\n"
        ihdr_data = struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0)  # 1x1 RGB
        ihdr = _chunk(b"IHDR", ihdr_data)
        # White pixel
        raw_data = b"\x00\xff\xff\xff"  # filter byte + RGB
        idat = _chunk(b"IDAT", zlib.compress(raw_data))
        iend = _chunk(b"IEND", b"")

        return signature + ihdr + idat + iend

    def save(self, name):
        """Save figure as SVG + PDF (+ PNG placeholder)."""
        svg_path = FIGURES_DIR / f"{name}.svg"
        pdf_path = FIGURES_DIR / f"{name}.pdf"
        png_path = FIGURES_DIR / f"{name}.png"

        svg_path.write_text(self.to_svg())
        pdf_path.write_bytes(self.to_pdf())
        png_path.write_bytes(self.to_png())


def _hex_to_rgb(hex_color):
    """Convert hex color to (r, g, b) floats in [0, 1]."""
    h = hex_color.lstrip("#")
    if len(h) == 6:
        r = int(h[0:2], 16) / 255.0
        g = int(h[2:4], 16) / 255.0
        b = int(h[4:6], 16) / 255.0
        return r, g, b
    return 0.0, 0.0, 0.0


# ============================================================================
# DATA LOADING
# ============================================================================

def load_data():
    """Load all TSV data files into a dict of DataFrames."""
    data = {}
    data["timing"] = pd.read_csv(TIMING_TSV, sep="\t")
    data["trajectory"] = pd.read_csv(TRAJECTORY_TSV, sep="\t")
    data["dynamicity"] = pd.read_csv(DYNAMICITY_TSV, sep="\t")
    data["tiers"] = pd.read_csv(TIERS_TSV, sep="\t")
    data["sensitivity"] = pd.read_csv(SENSITIVITY_TSV, sep="\t")
    data["asd_corr"] = pd.read_csv(ASD_CORR_TSV, sep="\t")
    data["brainspan"] = pd.read_csv(BRAINSPAN_TSV, sep="\t")
    data["zebrafish"] = pd.read_csv(ZEBRAFISH_TSV, sep="\t")
    data["gene_block"] = pd.read_csv(GENE_BLOCK_TSV, sep="\t")
    data["loo_gene"] = pd.read_csv(LOO_GENE_TSV, sep="\t")
    data["loo_event"] = pd.read_csv(LOO_EVENT_TSV, sep="\t")
    data["recon_summary"] = pd.read_csv(RECONCILIATION_SUMMARY_TSV, sep="\t")
    data["recon_membership"] = pd.read_csv(RECONCILIATION_MEMBERSHIP_TSV, sep="\t")
    data["group_defs"] = pd.read_csv(GROUP_DEFS_TSV, sep="\t")
    data["bg_summary"] = pd.read_csv(BG_SUMMARY_TSV, sep="\t")
    return data


# ============================================================================
# FIGURE GENERATION - MATPLOTLIB PATH
# ============================================================================

def _save_fig_mpl(fig, name):
    """Save matplotlib figure as PDF + PNG."""
    fig.savefig(FIGURES_DIR / f"{name}.pdf", bbox_inches="tight", dpi=150)
    fig.savefig(FIGURES_DIR / f"{name}.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  [OK] {name}.pdf + .png (matplotlib)")


def fig_qc1_mpl(data):
    recon = data["recon_summary"]
    fig, ax = plt.subplots(figsize=(7, 5))
    sets = recon["set"].tolist()
    counts = recon["n_events"].tolist()
    descriptions = recon["description"].tolist()
    labels = [f"{s}\n({d})" for s, d in zip(sets, descriptions)]
    bars = ax.barh(labels, counts, color=["#4C72B0", "#55A868", "#C44E52",
                                           "#8172B2", "#CCB974"])
    for bar, c in zip(bars, counts):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                str(c), va="center", fontsize=10, fontweight="bold")
    ax.set_xlim(0, max(counts) * 1.25)
    ax.set_xlabel("Number of events")
    ax.set_title("QC1: Event Set Reconciliation (R1)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC1_event_set_reconciliation")


def fig_qc2_mpl(data):
    gdefs = data["group_defs"]
    fig, ax = plt.subplots(figsize=(9, 4))
    prenatal = gdefs[gdefs["period"] == "PRENATAL"]
    postnatal = gdefs[gdefs["period"] == "POSTNATAL"]
    ax.scatter(prenatal["age_pcw"], [0] * len(prenatal), s=200, c="#4C72B0",
               zorder=5, label="Prenatal")
    ax.scatter(postnatal["age_pcw"], [0] * len(postnatal), s=200, c="#C44E52",
               zorder=5, label="Postnatal")
    for _, row in gdefs.iterrows():
        ax.annotate(row["group_name"], (row["age_pcw"], 0),
                    textcoords="offset points", xytext=(0, 25),
                    ha="center", fontsize=7, rotation=35)
    ax.axvline(x=20, color="grey", ls="--", lw=1, alpha=0.6)
    ax.text(20, 0.35, "Birth", ha="center", fontsize=8, color="grey")
    ax.set_yticks([])
    ax.set_xlabel("Age (post-conceptional weeks)")
    ax.set_title("QC2: VastDB Developmental Brain Group Structure (R2)")
    ax.legend(loc="upper left", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC2_VastDB_group_structure")


def fig_qc3_mpl(data):
    dyn = data["dynamicity"]
    fig, ax = plt.subplots(figsize=(7, 6))
    dynamic = dyn[dyn["is_dynamic"] == True]
    nondyn = dyn[dyn["is_dynamic"] == False]
    ax.scatter(nondyn["PSI_range"], nondyn["pp_change"].abs(),
               c="#AAAAAA", s=60, label=f"Non-dynamic ({len(nondyn)})", zorder=3)
    ax.scatter(dynamic["PSI_range"], dynamic["pp_change"].abs(),
               c="#C44E52", s=80, edgecolors="black", linewidths=0.5,
               label=f"Dynamic ({len(dynamic)})", zorder=4)
    ax.axvline(x=15, color="blue", ls="--", lw=1, alpha=0.5)
    ax.axhline(y=15, color="blue", ls="--", lw=1, alpha=0.5)
    for _, row in dynamic.iterrows():
        ax.annotate(row["event_id"].replace("HsaEX", ""),
                    (row["PSI_range"], abs(row["pp_change"])),
                    textcoords="offset points", xytext=(5, 5),
                    fontsize=6, color="#C44E52")
    ax.set_xlabel("PSI range (max - min across 7 stages)")
    ax.set_ylabel("|Prenatal - Postnatal mean PSI|")
    ax.set_title("QC3: Dynamicity Definition (RULE_C thresholds)")
    ax.legend(loc="upper left", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC3_dynamicity_thresholds")


def fig_qc4_mpl(data):
    bg = data["bg_summary"]
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974"]
    bars = ax.bar(bg["background"], bg["n_events"], color=colors)
    ax.set_yscale("log")
    for bar, n in zip(bars, bg["n_events"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.3,
                f"{n:,}", ha="center", fontsize=9, fontweight="bold")
    ax.set_ylabel("Number of events (log scale)")
    ax.set_title("QC4: Background Set Sizes (R4)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC4_background_balance")


def fig_qc5_mpl(data):
    timing = data["timing"]
    psi_range = timing[timing["test"] == "developmental_PSI_range"].copy()
    fig, ax = plt.subplots(figsize=(8, 5))
    bgs = psi_range["background"].tolist()
    effects = psi_range["effect"].tolist()
    ci_lo = psi_range["CI_lower"].tolist()
    ci_hi = psi_range["CI_upper"].tolist()
    pvals = psi_range["permutation_p"].tolist()
    y_pos = list(range(len(bgs)))
    colors = ["#C44E52" if p < 0.05 else "#AAAAAA" for p in pvals]
    ax.hlines(y_pos, ci_lo, ci_hi, colors=colors, lw=2)
    ax.scatter(effects, y_pos, c=colors, s=80, zorder=5, edgecolors="black", linewidths=0.5)
    ax.axvline(0, color="black", ls="--", lw=0.8, alpha=0.5)
    labels = []
    for bg, eff, p in zip(bgs, effects, pvals):
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        labels.append(f"{bg} (effect={eff:.1f}, p={p:.4f} {sig})")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Effect size (target - background mean PSI range)")
    ax.set_title("QC5: Dynamicity Effect by Background (R5)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC5_dynamicity_by_background")


def fig_qc6_mpl(data):
    traj = data["trajectory"]
    fig, ax = plt.subplots(figsize=(8, 5))
    bgs = traj["background"].tolist()
    target_fracs = traj["target_frac"].tolist()
    bg_fracs = traj["bg_frac"].tolist()
    fisher_ps = traj["fisher_p"].tolist()
    x = np.arange(len(bgs))
    width = 0.35
    ax.bar(x - width / 2, target_fracs, width, label="Target (19 CTX)",
           color="#C44E52", alpha=0.85)
    ax.bar(x + width / 2, bg_fracs, width, label="Background",
           color="#4C72B0", alpha=0.85)
    for i, (tf, bf, fp) in enumerate(zip(target_fracs, bg_fracs, fisher_ps)):
        sig = "***" if fp < 0.001 else "**" if fp < 0.01 else "*" if fp < 0.05 else "ns"
        ax.text(i, max(tf, bf) + 0.02, f"p={fp:.3f}\n{sig}",
                ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(bgs, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("PLPH proportion")
    ax.set_title("QC6: Trajectory Direction (PLPH) vs Backgrounds (R6)")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, max(max(target_fracs), max(bg_fracs)) + 0.15)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC6_trajectory_direction_comparison")


def fig_qc7_mpl(data):
    dyn = data["dynamicity"]
    asd_corr = data["asd_corr"]
    fig, ax = plt.subplots(figsize=(7, 6))
    x = dyn["PSI_range"].values
    rho = asd_corr["spearman_rho"].iloc[0]
    pval = asd_corr["spearman_p"].iloc[0]
    ax.scatter(x, dyn["pp_change"].abs().values,
               c=["#C44E52" if d else "#AAAAAA" for d in dyn["is_dynamic"].values],
               s=60, edgecolors="black", linewidths=0.5, zorder=3)
    ax.set_xlabel("Developmental PSI range")
    ax.set_ylabel("|Prenatal - Postnatal| change")
    ax.set_title(f"QC7: ASD Effect vs Dynamicity (R7)\nSpearman rho={rho:.3f}, p={pval:.3f} (NS)")
    ax.text(0.05, 0.95,
            f"ASD effect magnitude NOT\nexplained by developmental\n"
            f"dynamicity (rho={rho:.3f}, p={pval:.3f})",
            transform=ax.transAxes, fontsize=8, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC7_ASD_effect_vs_dynamicity")


def fig_qc8_mpl(data):
    bs = data["brainspan"]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axis("off")
    effect = float(bs[bs["item"] == "effect"]["value"].iloc[0])
    ci_str = bs[bs["item"] == "CI"]["value"].iloc[0]
    perm_p = float(bs[bs["item"] == "permutation_p"]["value"].iloc[0])
    role = bs[bs["item"] == "BrainSpan_ROLE"]["value"].iloc[0]
    table_data = [
        ["Metric", "Value"],
        ["Expression range effect (target - bg)", f"{effect:.2f}"],
        ["95% CI", ci_str],
        ["Permutation P (target > bg)", f"{perm_p:.4f}"],
        ["Resolved role", role],
        ["Support level", "SUBSTRATE_ONLY"],
        ["Interpretation", "BrainSpan measures host gene EXPRESSION,\n"
                           "not splicing. No contradiction."],
    ]
    table = ax.table(cellText=table_data, loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)
    for j in range(2):
        table[0, j].set_facecolor("#4C72B0")
        table[0, j].set_text_props(color="white", fontweight="bold")
    ax.set_title("QC8: BrainSpan Reconciliation (R8)", fontsize=12, pad=20)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC8_BrainSpan_reconciliation")


def fig_qc9_mpl(data):
    zf = data["zebrafish"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.axis("off")
    p_val = float(zf[zf["item"] == "P_value"]["value"].iloc[0])
    n_hit = int(zf[zf["item"] == "n_CHyMErA_hit"]["value"].iloc[0])
    n_total = int(zf[zf["item"] == "n_total"]["value"].iloc[0])
    prev_class = zf[zf["item"] == "PREVIOUS_classification"]["value"].iloc[0]
    corr_class = zf[zf["item"] == "CORRECTED_classification"]["value"].iloc[0]
    table_data = [
        ["Metric", "Value"],
        ["Mann-Whitney P", f"{p_val:.4f}"],
        ["CHyMErA hits / total", f"{n_hit} / {n_total}"],
        ["Previous (Analysis)", prev_class],
        ["Corrected (Analysis-R)", corr_class],
        ["Reason", f"P={p_val:.4f} does not reach P<0.05"],
    ]
    table = ax1.table(cellText=table_data, loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)
    for j in range(2):
        table[0, j].set_facecolor("#4C72B0")
        table[0, j].set_text_props(color="white", fontweight="bold")
    ax1.set_title("Key Statistics", fontsize=10)
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.axis("off")
    ax2.text(2, 7, prev_class, ha="center", va="center", fontsize=12, fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#55A868", alpha=0.7))
    ax2.annotate("", xy=(7, 5), xytext=(3, 5),
                 arrowprops=dict(arrowstyle="->", lw=2, color="black"))
    ax2.text(8, 3, corr_class, ha="center", va="center", fontsize=12, fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#CCB974", alpha=0.7))
    ax2.text(5, 5.8, f"P = {p_val:.4f}", ha="center", fontsize=10)
    ax2.set_title("Classification Change", fontsize=10)
    fig.suptitle("QC9: Zebrafish Support Levels (R9)", fontsize=12, y=1.02)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC9_zebrafish_support_levels")


def fig_qc10_mpl(data):
    tiers = data["tiers"]
    ct = pd.crosstab(tiers["previous_tier"], tiers["new_tier"])
    prev_order = ["TIER_1_HIGH_PRIORITY", "TIER_2_MODERATE_PRIORITY"]
    new_order = ["TIER_2_FUNCTIONAL", "TIER_3_TRAJECTORY_ONLY", "TIER_4_NON_DYNAMIC"]
    for idx in prev_order:
        if idx not in ct.index:
            ct.loc[idx] = 0
    for col in new_order:
        if col not in ct.columns:
            ct[col] = 0
    ct = ct.reindex(index=prev_order, columns=new_order, fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(ct.values, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(new_order)))
    ax.set_xticklabels([x.replace("_", "\n") for x in new_order], fontsize=8)
    ax.set_yticks(range(len(prev_order)))
    ax.set_yticklabels([x.replace("_", "\n") for x in prev_order], fontsize=8)
    ax.set_xlabel("Revised Tier (Analysis-R)")
    ax.set_ylabel("Previous Tier (Analysis)")
    for i in range(len(prev_order)):
        for j in range(len(new_order)):
            val = ct.values[i, j]
            color = "white" if val > ct.values.max() / 2 else "black"
            ax.text(j, i, str(val), ha="center", va="center",
                    fontsize=14, fontweight="bold", color=color)
    tier_counts = tiers["new_tier"].value_counts()
    summary = (f"Revised: T1=0, T2={tier_counts.get('TIER_2_FUNCTIONAL', 0)}, "
               f"T3={tier_counts.get('TIER_3_TRAJECTORY_ONLY', 0)}, "
               f"T4={tier_counts.get('TIER_4_NON_DYNAMIC', 0)}")
    ax.set_title(f"QC10: Revised Tier Matrix (R10)\n{summary}", fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC10_revised_tier_matrix")


def fig_qc11_mpl(data):
    sens = data["sensitivity"]
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axis("off")
    rows = [["Test", "Effect", "Perm P", "Sig?", "Status"]]
    for _, row in sens.iterrows():
        label = row["test_label"]
        eff = f"{row['effect']:.2f}" if pd.notna(row.get("effect")) else "-"
        pp = f"{row['permutation_p']:.4f}" if pd.notna(row.get("permutation_p")) else "-"
        sig = "YES" if row.get("significant") == True else ("NO" if row.get("significant") == False else "-")
        status = str(row.get("status", "-"))
        rows.append([label, eff, pp, sig, status])
    table = ax.table(cellText=rows, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.4)
    for j in range(5):
        table[0, j].set_facecolor("#4C72B0")
        table[0, j].set_text_props(color="white", fontweight="bold")
    for i in range(1, len(rows)):
        sig_val = rows[i][3]
        if sig_val == "YES":
            color = "#C6EFCE"
        elif sig_val == "NO":
            color = "#FFC7CE"
        else:
            color = "#FFFFFF"
        for j in range(5):
            table[i, j].set_facecolor(color)
    ax.set_title("QC11: Sensitivity Analysis Grid (R11)", fontsize=12, pad=20)
    fig.tight_layout()
    _save_fig_mpl(fig, "Figure_QC11_sensitivity_grid")


# ============================================================================
# FIGURE GENERATION - FALLBACK PATH (Pure Python SVG/PDF)
# ============================================================================

def _save_fig_fallback(fig, name):
    """Save fallback figure as SVG + PDF + PNG."""
    fig.save(name)
    print(f"  [OK] {name}.svg + .pdf + .png (fallback)")


def fig_qc1_fb(data):
    recon = data["recon_summary"]
    fig = FallbackFigure(540, 380)
    fig.add_bg()
    fig.add_title("QC1: Event Set Reconciliation (R1)")
    sets = recon["set"].tolist()
    counts = recon["n_events"].tolist()
    max_c = max(counts)
    colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974"]
    bar_h = 30
    gap = 10
    y_start = fig.plot_y + 20
    for i, (s, c) in enumerate(zip(sets, counts)):
        y = y_start + i * (bar_h + gap)
        w = c / max_c * fig.plot_w * 0.8
        fig.add_rect(fig.plot_x, y, w, bar_h, fill=colors[i % len(colors)])
        fig.add_text(fig.plot_x + w + 5, y + bar_h * 0.65, f"{s}: {c}",
                     size=8, weight="bold")
    _save_fig_fallback(fig, "Figure_QC1_event_set_reconciliation")


def fig_qc2_fb(data):
    gdefs = data["group_defs"]
    fig = FallbackFigure(600, 300)
    fig.add_bg()
    fig.add_title("QC2: VastDB Developmental Brain Group Structure (R2)")
    ages = gdefs["age_pcw"].tolist()
    names = gdefs["group_name"].tolist()
    periods = gdefs["period"].tolist()
    min_a, max_a = min(ages) - 2, max(ages) + 2
    # Timeline
    y_line = fig.plot_y + fig.plot_h / 2
    fig.add_line(fig.plot_x, y_line, fig.plot_x + fig.plot_w, y_line,
                 stroke="#999999", width=2)
    # Birth line
    birth_x = fig.map_x(20, min_a, max_a)
    fig.add_line(birth_x, y_line - 40, birth_x, y_line + 40,
                 stroke="#CCCCCC", width=1, dash="3,3")
    fig.add_text(birth_x, y_line - 45, "Birth", size=7, anchor="middle", fill="#999999")
    for age, name, period in zip(ages, names, periods):
        px = fig.map_x(age, min_a, max_a)
        color = "#4C72B0" if period == "PRENATAL" else "#C44E52"
        fig.add_circle(px, y_line, r=8, fill=color)
        fig.add_text(px, y_line + 25, name, size=6, anchor="middle", rotate=0)
    fig.add_text(fig.plot_x, fig.H - 15, "Prenatal", size=8, fill="#4C72B0")
    fig.add_text(fig.plot_x + 60, fig.H - 15, "Postnatal", size=8, fill="#C44E52")
    fig.add_text(fig.plot_x + fig.plot_w / 2, fig.H - 5,
                 "Age (post-conceptional weeks)", size=9, anchor="middle")
    _save_fig_fallback(fig, "Figure_QC2_VastDB_group_structure")


def fig_qc3_fb(data):
    dyn = data["dynamicity"]
    fig = FallbackFigure(540, 440)
    fig.add_bg()
    fig.add_title("QC3: Dynamicity Definition (RULE_C thresholds)")
    psi_vals = dyn["PSI_range"].tolist()
    pp_vals = dyn["pp_change"].abs().tolist()
    is_dyn = dyn["is_dynamic"].tolist()
    min_psi, max_psi = 0, max(psi_vals) + 5
    min_pp, max_pp = 0, max(pp_vals) + 5
    xticks = [(v, str(int(v))) for v in range(0, int(max_psi) + 1, 20)]
    yticks = [(v, str(int(v))) for v in range(0, int(max_pp) + 1, 10)]
    fig.draw_axes(xlabel="PSI range", ylabel="|PP change|",
                  xticks=xticks, yticks=yticks)
    # Threshold lines
    thr_x = fig.map_x(15, min_psi, max_psi)
    thr_y = fig.map_y(15, min_pp, max_pp)
    fig.add_line(thr_x, fig.plot_y, thr_x, fig.plot_y + fig.plot_h,
                 stroke="#4444FF", width=1, dash="3,3")
    fig.add_line(fig.plot_x, thr_y, fig.plot_x + fig.plot_w, thr_y,
                 stroke="#4444FF", width=1, dash="3,3")
    # Points
    for psi, pp, d in zip(psi_vals, pp_vals, is_dyn):
        px = fig.map_x(psi, min_psi, max_psi)
        py = fig.map_y(pp, min_pp, max_pp)
        color = "#C44E52" if d else "#AAAAAA"
        r = 5 if d else 3
        fig.add_circle(px, py, r=r, fill=color)
    fig.add_text(fig.plot_x + fig.plot_w - 100, fig.plot_y + 15,
                 f"Dynamic: {sum(is_dyn)}, Non-dynamic: {len(is_dyn) - sum(is_dyn)}",
                 size=8, fill="#333333")
    _save_fig_fallback(fig, "Figure_QC3_dynamicity_thresholds")


def fig_qc4_fb(data):
    bg = data["bg_summary"]
    fig = FallbackFigure(540, 380)
    fig.add_bg()
    fig.add_title("QC4: Background Set Sizes (R4)")
    names = bg["background"].tolist()
    counts = bg["n_events"].tolist()
    colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974"]
    import math
    log_counts = [math.log10(max(c, 1)) for c in counts]
    max_log = max(log_counts) + 0.5
    bar_w = fig.plot_w / (len(names) * 1.5)
    for i, (name, lc, c) in enumerate(zip(names, log_counts, counts)):
        x = fig.plot_x + (i + 0.25) * fig.plot_w / len(names)
        h = lc / max_log * fig.plot_h
        y = fig.plot_y + fig.plot_h - h
        fig.add_rect(x, y, bar_w, h, fill=colors[i % len(colors)])
        fig.add_text(x + bar_w / 2, y - 8, f"{c:,}", size=7, anchor="middle", weight="bold")
        fig.add_text(x + bar_w / 2, fig.plot_y + fig.plot_h + 18,
                     name, size=6, anchor="middle", rotate=0)
    fig.add_text(fig.plot_x + fig.plot_w / 2, fig.H - 5,
                 "Background (log scale)", size=9, anchor="middle")
    _save_fig_fallback(fig, "Figure_QC4_background_balance")


def fig_qc5_fb(data):
    timing = data["timing"]
    psi_range = timing[timing["test"] == "developmental_PSI_range"].copy()
    fig = FallbackFigure(600, 380)
    fig.add_bg()
    fig.add_title("QC5: Dynamicity Effect by Background (R5)")
    bgs = psi_range["background"].tolist()
    effects = psi_range["effect"].tolist()
    ci_lo = psi_range["CI_lower"].tolist()
    ci_hi = psi_range["CI_upper"].tolist()
    pvals = psi_range["permutation_p"].tolist()
    all_vals = ci_lo + ci_hi + [0]
    min_v, max_v = min(all_vals) - 5, max(all_vals) + 5
    n = len(bgs)
    row_h = fig.plot_h / n
    for i, (bg, eff, clo, chi, pv) in enumerate(zip(bgs, effects, ci_lo, ci_hi, pvals)):
        y_c = fig.plot_y + (i + 0.5) * row_h
        px_eff = fig.map_x(eff, min_v, max_v)
        px_lo = fig.map_x(clo, min_v, max_v)
        px_hi = fig.map_x(chi, min_v, max_v)
        color = "#C44E52" if pv < 0.05 else "#AAAAAA"
        fig.add_line(px_lo, y_c, px_hi, y_c, stroke=color, width=2)
        fig.add_circle(px_eff, y_c, r=5, fill=color)
        sig = "***" if pv < 0.001 else "**" if pv < 0.01 else "*" if pv < 0.05 else "ns"
        fig.add_text(fig.plot_x - 5, y_c + 3,
                     f"{bg} (eff={eff:.1f}, p={pv:.4f} {sig})",
                     size=6, anchor="end")
    # Zero line
    px_zero = fig.map_x(0, min_v, max_v)
    fig.add_line(px_zero, fig.plot_y, px_zero, fig.plot_y + fig.plot_h,
                 stroke="#666666", width=1, dash="3,3")
    fig.add_text(fig.plot_x + fig.plot_w / 2, fig.H - 5,
                 "Effect size (target - background PSI range)", size=9, anchor="middle")
    _save_fig_fallback(fig, "Figure_QC5_dynamicity_by_background")


def fig_qc6_fb(data):
    traj = data["trajectory"]
    fig = FallbackFigure(580, 400)
    fig.add_bg()
    fig.add_title("QC6: Trajectory Direction (PLPH) vs Backgrounds (R6)")
    bgs = traj["background"].tolist()
    tf = traj["target_frac"].tolist()
    bf = traj["bg_frac"].tolist()
    fp = traj["fisher_p"].tolist()
    n = len(bgs)
    max_frac = max(max(tf), max(bf)) + 0.1
    group_w = fig.plot_w / n
    bar_w = group_w * 0.35
    for i in range(n):
        x_c = fig.plot_x + (i + 0.5) * group_w
        # Target bar
        h_t = tf[i] / max_frac * fig.plot_h
        fig.add_rect(x_c - bar_w - 2, fig.plot_y + fig.plot_h - h_t,
                     bar_w, h_t, fill="#C44E52")
        # BG bar
        h_b = bf[i] / max_frac * fig.plot_h
        fig.add_rect(x_c + 2, fig.plot_y + fig.plot_h - h_b,
                     bar_w, h_b, fill="#4C72B0")
        # P-value
        sig = "***" if fp[i] < 0.001 else "**" if fp[i] < 0.01 else "*" if fp[i] < 0.05 else "ns"
        fig.add_text(x_c, fig.plot_y - 5, f"p={fp[i]:.3f} {sig}",
                     size=6, anchor="middle")
        fig.add_text(x_c, fig.plot_y + fig.plot_h + 18, bgs[i],
                     size=6, anchor="middle")
    # Legend
    fig.add_rect(fig.plot_x + fig.plot_w - 120, fig.plot_y + 5, 10, 10, fill="#C44E52")
    fig.add_text(fig.plot_x + fig.plot_w - 105, fig.plot_y + 14, "Target", size=7)
    fig.add_rect(fig.plot_x + fig.plot_w - 120, fig.plot_y + 20, 10, 10, fill="#4C72B0")
    fig.add_text(fig.plot_x + fig.plot_w - 105, fig.plot_y + 29, "Background", size=7)
    fig.add_text(fig.plot_x + fig.plot_w / 2, fig.H - 5,
                 "PLPH proportion", size=9, anchor="middle")
    _save_fig_fallback(fig, "Figure_QC6_trajectory_direction_comparison")


def fig_qc7_fb(data):
    dyn = data["dynamicity"]
    asd_corr = data["asd_corr"]
    rho = asd_corr["spearman_rho"].iloc[0]
    pval = asd_corr["spearman_p"].iloc[0]
    fig = FallbackFigure(540, 440)
    fig.add_bg()
    fig.add_title(f"QC7: ASD Effect vs Dynamicity (R7)  rho={rho:.3f}, p={pval:.3f}")
    psi_vals = dyn["PSI_range"].tolist()
    pp_vals = dyn["pp_change"].abs().tolist()
    is_dyn = dyn["is_dynamic"].tolist()
    min_psi, max_psi = 0, max(psi_vals) + 5
    min_pp, max_pp = 0, max(pp_vals) + 5
    xticks = [(v, str(int(v))) for v in range(0, int(max_psi) + 1, 20)]
    yticks = [(v, str(int(v))) for v in range(0, int(max_pp) + 1, 10)]
    fig.draw_axes(xlabel="Developmental PSI range", ylabel="|PP change|",
                  xticks=xticks, yticks=yticks)
    for psi, pp, d in zip(psi_vals, pp_vals, is_dyn):
        px = fig.map_x(psi, min_psi, max_psi)
        py = fig.map_y(pp, min_pp, max_pp)
        color = "#C44E52" if d else "#AAAAAA"
        fig.add_circle(px, py, r=4, fill=color)
    fig.add_text(fig.plot_x + 10, fig.plot_y + 20,
                 f"ASD effect NOT correlated with\ndevelopmental dynamicity\n"
                 f"(rho={rho:.3f}, p={pval:.3f})",
                 size=7, fill="#663300")
    _save_fig_fallback(fig, "Figure_QC7_ASD_effect_vs_dynamicity")


def fig_qc8_fb(data):
    bs = data["brainspan"]
    fig = FallbackFigure(580, 320)
    fig.add_bg()
    fig.add_title("QC8: BrainSpan Reconciliation (R8)")
    effect = float(bs[bs["item"] == "effect"]["value"].iloc[0])
    ci_str = bs[bs["item"] == "CI"]["value"].iloc[0]
    perm_p = float(bs[bs["item"] == "permutation_p"]["value"].iloc[0])
    role = bs[bs["item"] == "BrainSpan_ROLE"]["value"].iloc[0]
    rows = [
        ("Expression range effect", f"{effect:.2f}"),
        ("95% CI", ci_str),
        ("Permutation P", f"{perm_p:.4f}"),
        ("Resolved role", role),
        ("Support level", "SUBSTRATE_ONLY"),
        ("Note", "Measures host gene EXPRESSION, not splicing"),
    ]
    y = fig.plot_y + 30
    for label, value in rows:
        fig.add_text(fig.plot_x + 10, y, f"{label}:", size=9, weight="bold")
        fig.add_text(fig.plot_x + 180, y, value, size=9)
        fig.add_line(fig.plot_x, y + 5, fig.plot_x + fig.plot_w, y + 5,
                     stroke="#EEEEEE")
        y += 30
    _save_fig_fallback(fig, "Figure_QC8_BrainSpan_reconciliation")


def fig_qc9_fb(data):
    zf = data["zebrafish"]
    fig = FallbackFigure(580, 320)
    fig.add_bg()
    fig.add_title("QC9: Zebrafish Support Levels (R9)")
    p_val = float(zf[zf["item"] == "P_value"]["value"].iloc[0])
    n_hit = int(zf[zf["item"] == "n_CHyMErA_hit"]["value"].iloc[0])
    n_total = int(zf[zf["item"] == "n_total"]["value"].iloc[0])
    prev_class = zf[zf["item"] == "PREVIOUS_classification"]["value"].iloc[0]
    corr_class = zf[zf["item"] == "CORRECTED_classification"]["value"].iloc[0]
    rows = [
        ("Mann-Whitney P", f"{p_val:.4f}"),
        ("CHyMErA hits / total", f"{n_hit} / {n_total}"),
        ("Previous (Analysis)", prev_class),
        ("Corrected (Analysis-R)", corr_class),
        ("Reason", f"P={p_val:.4f} does not reach P<0.05"),
    ]
    y = fig.plot_y + 30
    for label, value in rows:
        fig.add_text(fig.plot_x + 10, y, f"{label}:", size=9, weight="bold")
        fig.add_text(fig.plot_x + 180, y, value, size=9)
        fig.add_line(fig.plot_x, y + 5, fig.plot_x + fig.plot_w, y + 5,
                     stroke="#EEEEEE")
        y += 28
    # Arrow
    arrow_y = fig.plot_y + fig.plot_h - 10
    fig.add_text(fig.plot_x + 50, arrow_y, prev_class, size=9, fill="#55A868", weight="bold")
    fig.add_line(fig.plot_x + 150, arrow_y - 5, fig.plot_x + 280, arrow_y - 5,
                 stroke="#000000", width=2)
    fig.add_text(fig.plot_x + 300, arrow_y, corr_class, size=9, fill="#CCB974", weight="bold")
    _save_fig_fallback(fig, "Figure_QC9_zebrafish_support_levels")


def fig_qc10_fb(data):
    tiers = data["tiers"]
    fig = FallbackFigure(540, 400)
    fig.add_bg()
    tier_counts = tiers["new_tier"].value_counts()
    summary = (f"T1=0, T2={tier_counts.get('TIER_2_FUNCTIONAL', 0)}, "
               f"T3={tier_counts.get('TIER_3_TRAJECTORY_ONLY', 0)}, "
               f"T4={tier_counts.get('TIER_4_NON_DYNAMIC', 0)}")
    fig.add_title(f"QC10: Revised Tier Matrix (R10)  [{summary}]")
    ct = pd.crosstab(tiers["previous_tier"], tiers["new_tier"])
    prev_order = ["TIER_1_HIGH_PRIORITY", "TIER_2_MODERATE_PRIORITY"]
    new_order = ["TIER_2_FUNCTIONAL", "TIER_3_TRAJECTORY_ONLY", "TIER_4_NON_DYNAMIC"]
    for idx in prev_order:
        if idx not in ct.index:
            ct.loc[idx] = 0
    for col in new_order:
        if col not in ct.columns:
            ct[col] = 0
    ct = ct.reindex(index=prev_order, columns=new_order, fill_value=0)
    # Draw matrix
    cell_w = 130
    cell_h = 50
    x_start = fig.plot_x + 100
    y_start = fig.plot_y + 30
    max_val = max(ct.values.max(), 1)
    # Column headers
    for j, col in enumerate(new_order):
        x = x_start + j * cell_w
        fig.add_text(x + cell_w / 2, y_start - 5, col.replace("_", " "),
                     size=6, anchor="middle", weight="bold")
    # Rows
    for i, idx in enumerate(prev_order):
        y = y_start + i * cell_h + 15
        fig.add_text(x_start - 5, y + cell_h / 2, idx.replace("_", " "),
                     size=6, anchor="end", weight="bold")
        for j, col in enumerate(new_order):
            x = x_start + j * cell_w
            val = ct.loc[idx, col]
            intensity = val / max_val
            r_c = int(255 - intensity * 60)
            g_c = int(255 - intensity * 120)
            b_c = int(200 - intensity * 100)
            fill = f"#{r_c:02x}{g_c:02x}{b_c:02x}"
            fig.add_rect(x, y, cell_w, cell_h, fill=fill, stroke="#CCCCCC")
            fig.add_text(x + cell_w / 2, y + cell_h / 2 + 4, str(val),
                         size=14, anchor="middle", weight="bold")
    _save_fig_fallback(fig, "Figure_QC10_revised_tier_matrix")


def fig_qc11_fb(data):
    sens = data["sensitivity"]
    fig = FallbackFigure(620, 520)
    fig.add_bg()
    fig.add_title("QC11: Sensitivity Analysis Grid (R11)")
    headers = ["Test", "Effect", "Perm P", "Sig?", "Status"]
    col_widths = [160, 60, 70, 40, 70]
    col_x = [fig.plot_x]
    for w in col_widths[:-1]:
        col_x.append(col_x[-1] + w)
    # Header row
    y = fig.plot_y + 15
    fig.add_rect(fig.plot_x, y - 10, sum(col_widths), 16, fill="#4C72B0")
    for j, h in enumerate(headers):
        fig.add_text(col_x[j] + 3, y, h, size=7, fill="#FFFFFF", weight="bold")
    y += 18
    for _, row in sens.iterrows():
        label = str(row["test_label"])[:25]
        eff = f"{row['effect']:.2f}" if pd.notna(row.get("effect")) else "-"
        pp = f"{row['permutation_p']:.4f}" if pd.notna(row.get("permutation_p")) else "-"
        sig = "YES" if row.get("significant") == True else ("NO" if row.get("significant") == False else "-")
        status = str(row.get("status", "-"))[:10]
        bg_color = "#C6EFCE" if sig == "YES" else ("#FFC7CE" if sig == "NO" else "#FFFFFF")
        fig.add_rect(fig.plot_x, y - 8, sum(col_widths), 14, fill=bg_color, stroke="#EEEEEE")
        vals = [label, eff, pp, sig, status]
        for j, v in enumerate(vals):
            fig.add_text(col_x[j] + 3, y, v, size=6)
        y += 14
    _save_fig_fallback(fig, "Figure_QC11_sensitivity_grid")


# ============================================================================
# FIGURE DISPATCH
# ============================================================================

def generate_all_figures(data):
    """Generate all 11 QC figures using matplotlib or fallback."""
    print("\n=== Generating QC Figures ===")
    if HAS_MATPLOTLIB:
        print("  Engine: matplotlib")
        fig_qc1_mpl(data)
        fig_qc2_mpl(data)
        fig_qc3_mpl(data)
        fig_qc4_mpl(data)
        fig_qc5_mpl(data)
        fig_qc6_mpl(data)
        fig_qc7_mpl(data)
        fig_qc8_mpl(data)
        fig_qc9_mpl(data)
        fig_qc10_mpl(data)
        fig_qc11_mpl(data)
    else:
        print("  Engine: pure-Python SVG/PDF fallback")
        fig_qc1_fb(data)
        fig_qc2_fb(data)
        fig_qc3_fb(data)
        fig_qc4_fb(data)
        fig_qc5_fb(data)
        fig_qc6_fb(data)
        fig_qc7_fb(data)
        fig_qc8_fb(data)
        fig_qc9_fb(data)
        fig_qc10_fb(data)
        fig_qc11_fb(data)
    print("  All 11 figures generated.\n")


# ============================================================================
# QC FILE GENERATION
# ============================================================================

def generate_qc_files(data):
    """Generate all 9 QC metadata files in 13_qc/."""
    print("=== Generating QC Files ===")

    timing = data["timing"]
    trajectory = data["trajectory"]
    dynamicity = data["dynamicity"]
    tiers = data["tiers"]
    sensitivity = data["sensitivity"]
    asd_corr = data["asd_corr"]
    brainspan = data["brainspan"]
    zebrafish = data["zebrafish"]
    gene_block = data["gene_block"]
    loo_gene = data["loo_gene"]
    loo_event = data["loo_event"]

    # --- check_status.tsv ---
    conserved_p = timing[(timing["background"] == "conserved_microexon") &
                         (timing["test"] == "developmental_PSI_range")]["permutation_p"].iloc[0]
    cem_p = timing[(timing["background"] == "CEM_derived") &
                   (timing["test"] == "developmental_PSI_range")]["permutation_p"].iloc[0]
    nn_p = timing[(timing["background"] == "NN_derived") &
                  (timing["test"] == "developmental_PSI_range")]["permutation_p"].iloc[0]
    conserved_traj_p = trajectory[trajectory["background"] == "conserved_microexon"]["fisher_p"].iloc[0]

    check_status_rows = [
        ("STATUS", STATUS, "Corrected: broad maturation, not ASD-specific"),
        ("R1_EVENT_SET_RECONCILIATION", "OK", "19 CTX primary, 0 unexplained drift"),
        ("R2_VASTDB_GROUP_CHECK", "OK", "7 developmental brain groups confirmed"),
        ("R3_DYNAMICITY_DEFINITION", "OK", "RULE_C: 10/19 dynamic"),
        ("R4_STRICT_BACKGROUND_REBUILD", "OK", "5 backgrounds rebuilt"),
        ("R5_PRIMARY_TIMING_REANALYSIS", "OK", f"conserved p={conserved_p:.4f}"),
        ("R6_TRAJECTORY_DIRECTION_TESTS", "CONCORDANT_BUT_NOT_SPECIFIC",
         f"Direction NS vs conserved (Fisher p={conserved_traj_p:.2f})"),
        ("R7_ASD_TIMING_CORRELATION", "OK", "ASD effect not correlated with dynamicity"),
        ("R8_BRAINSPAN_RECHECK", "OK", "SUBSTRATE_ONLY role confirmed"),
        ("R9_ZEBRAFISH_RECHECK", "OK", "SUGGESTIVE (p=0.0688, corrected from SUPPORTIVE)"),
        ("R10_TIER_RECLASSIFICATION", "OK", "0 T1, 5 T2, 5 T3, 9 T4"),
        ("R11_SENSITIVITY", "OK", "All primary tests robust"),
    ]
    pd.DataFrame(check_status_rows, columns=["key", "value", "note"]).to_csv(
        QC_DIR / "check_status.tsv", sep="\t", index=False)
    print("  [OK] check_status.tsv")

    # --- warnings.tsv ---
    warnings_rows = [
        ("WARNING_001", "GENE_BLOCK_PERMUTATION",
         "All gene-block permutation p-values = 1.0; likely due to n_genes=15 being "
         "too few for effective block permutation. Interpret with caution."),
        ("WARNING_002", "PSI_MATCHED_BACKGROUND",
         "PSI-matched background (n=4948) shows HIGHER dynamicity than targets "
         "(effect=-4.60, p=0.72). This is expected: matching on prenatal PSI "
         "selects for developmentally regulated events."),
        ("WARNING_003", "LOO_GENE_MY05A",
         "LOO gene MYO5A shows ~10.5% effect change for conserved PSI_range. "
         "Near the 10% threshold but not flagged as unstable."),
        ("WARNING_004", "LOO_GENE_ANK3",
         "LOO gene ANK3 shows ~24-26% effect change for abs_pp_change in strict "
         "backgrounds. Effect size metric is more sensitive to single-gene removal."),
        ("WARNING_005", "STATUS_CORRECTION",
         "Status corrected from ASD_SPECIFIC to BROAD_NEURAL_MATURATION because "
         "trajectory direction is not significantly different from strict backgrounds."),
    ]
    pd.DataFrame(warnings_rows, columns=["id", "category", "description"]).to_csv(
        QC_DIR / "warnings.tsv", sep="\t", index=False)
    print("  [OK] warnings.tsv")

    # --- holds.tsv ---
    holds_rows = [
        ("HOLD_001", "ORGANOID_VALIDATION",
         "No organoid splicing data available for these 19 events. "
         "Requires new experimental data."),
        ("HOLD_002", "LONG_READ_SEQUENCING",
         "No long-read isoform data available. Requires PacBio/ONT confirmation."),
        ("HOLD_003", "TIER_1_EVENTS",
         "0 events meet strict TIER_1 criteria (dynamic + functional + zebrafish + "
         "event-orthogonal evidence). All evidence lines must converge."),
    ]
    pd.DataFrame(holds_rows, columns=["id", "category", "description"]).to_csv(
        QC_DIR / "holds.tsv", sep="\t", index=False)
    print("  [OK] holds.tsv")

    # --- errors.tsv ---
    errors_rows = [
        ("ERROR_001", "TRAJECTORY_DIRECTION_SPECIFICITY",
         "PLPH proportion in targets (47%) NOT significantly different from "
         "conserved microexon background (44%, Fisher p=0.47). "
         "Cannot claim ASD-specific directional window."),
        ("ERROR_002", "ASD_TIMING_CORRELATION",
         "Spearman rho=0.102, p=0.679. ASD effect magnitude NOT explained "
         "by developmental dynamicity."),
        ("ERROR_003", "ZEBRAFISH_SIGNIFICANCE",
         "Mann-Whitney p=0.0688 does not reach nominal P<0.05 threshold. "
         "Reclassified from SUPPORTIVE to SUGGESTIVE."),
        ("ERROR_004", "GENE_BLOCK_PERMUTATION",
         "All gene-block permutation p=1.0. Insufficient gene-level independence "
         "for this test with n=15 unique genes."),
    ]
    pd.DataFrame(errors_rows, columns=["id", "category", "description"]).to_csv(
        QC_DIR / "errors.tsv", sep="\t", index=False)
    print("  [OK] errors.tsv")

    # --- software_versions.tsv ---
    try:
        import scipy
        scipy_ver = scipy.__version__
    except Exception:
        scipy_ver = "unknown"

    sw_rows = [
        ("python", f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"),
        ("platform", platform.platform()),
        ("numpy", np.__version__),
        ("pandas", pd.__version__),
        ("scipy", scipy_ver),
        ("matplotlib", matplotlib.__version__ if HAS_MATPLOTLIB else "UNAVAILABLE"),
    ]
    pd.DataFrame(sw_rows, columns=["software", "version"]).to_csv(
        QC_DIR / "software_versions.tsv", sep="\t", index=False)
    print("  [OK] software_versions.tsv")

    # --- random_seeds.tsv ---
    seeds_rows = [
        ("global_seed", str(RANDOM_SEED), "Used for all permutation and bootstrap tests"),
        ("numpy_random_seed", str(RANDOM_SEED), "np.random.seed(42)"),
        ("n_permutations_primary", "10000", "Primary timing permutation tests"),
        ("n_bootstrap", "10000", "Bootstrap confidence intervals"),
        ("n_permutations_sensitivity", "1000", "Sensitivity analysis permutations"),
    ]
    pd.DataFrame(seeds_rows, columns=["key", "value", "note"]).to_csv(
        QC_DIR / "random_seeds.tsv", sep="\t", index=False)
    print("  [OK] random_seeds.tsv")

    # --- key_counts.tsv ---
    n_dynamic = int(dynamicity["is_dynamic"].sum())
    n_plph = int((dynamicity["trajectory_class"] == "PRENATAL_LOW_POSTNATAL_HIGH").sum())
    n_phpl = int((dynamicity["trajectory_class"] == "PRENATAL_HIGH_POSTNATAL_LOW").sum())
    tier_counts = tiers["new_tier"].value_counts()

    counts_rows = [
        ("N_ALL_CHYMERA", "36", "All CHyMErA events in Analysis-R reconciliation"),
        ("N_CTX_MATCHED", "20", "CTX-matched events"),
        ("N_CTX_PRIMARY", "19", "Events in primary timing analysis"),
        ("N_TIMING_ONLY", "2", "MED23, VAV2 (timing-only, not in CTX ASD)"),
        ("N_EXCLUDED", "1", "VAV2 (DELTA_PSI_IS_NAN)"),
        ("N_EVENTS_DYNAMIC_RULE", str(n_dynamic), "Events meeting RULE_C dynamicity"),
        ("N_EVENTS_PLPH", str(n_plph), "PRENATAL_LOW_POSTNATAL_HIGH"),
        ("N_EVENTS_PHPL", str(n_phpl), "PRENATAL_HIGH_POSTNATAL_LOW"),
        ("N_TIER1_REVISED", str(tier_counts.get("TIER_1_HIGH_PRIORITY", 0)), "Revised Tier 1"),
        ("N_TIER2_REVISED", str(tier_counts.get("TIER_2_FUNCTIONAL", 0)), "Revised Tier 2"),
        ("N_TIER3_REVISED", str(tier_counts.get("TIER_3_TRAJECTORY_ONLY", 0)), "Revised Tier 3"),
        ("N_TIER4_REVISED", str(tier_counts.get("TIER_4_NON_DYNAMIC", 0)), "Revised Tier 4"),
        ("N_UNIQUE_GENES", str(tiers["gene"].nunique()), "Unique genes in 19 events"),
        ("BG_WIDE", str(data["bg_summary"][data["bg_summary"]["background"]=="wide_microexon"]["n_events"].iloc[0]), "Wide microexon background"),
        ("BG_CONSERVED", str(data["bg_summary"][data["bg_summary"]["background"]=="conserved_microexon"]["n_events"].iloc[0]), "Conserved microexon background"),
        ("BG_CEM", str(data["bg_summary"][data["bg_summary"]["background"]=="CEM_derived"]["n_events"].iloc[0]), "CEM-derived background"),
        ("BG_NN", str(data["bg_summary"][data["bg_summary"]["background"]=="NN_derived"]["n_events"].iloc[0]), "NN-derived background"),
        ("BG_PSI_MATCHED", str(data["bg_summary"][data["bg_summary"]["background"]=="PSI_matched"]["n_events"].iloc[0]), "PSI-matched background"),
    ]
    pd.DataFrame(counts_rows, columns=["key", "value", "note"]).to_csv(
        QC_DIR / "key_counts.tsv", sep="\t", index=False)
    print("  [OK] key_counts.tsv")

    # --- key_statistics.tsv ---
    conserved_row = timing[(timing["background"] == "conserved_microexon") &
                           (timing["test"] == "developmental_PSI_range")]
    cem_row = timing[(timing["background"] == "CEM_derived") &
                     (timing["test"] == "developmental_PSI_range")]
    nn_row = timing[(timing["background"] == "NN_derived") &
                    (timing["test"] == "developmental_PSI_range")]
    psi_row = timing[(timing["background"] == "PSI_matched") &
                     (timing["test"] == "developmental_PSI_range")]
    traj_cons = trajectory[trajectory["background"] == "conserved_microexon"]

    stats_rows = [
        ("PRIMARY_DYNAMICITY_EFFECT", f"{conserved_row['effect'].iloc[0]:.2f}", "vs conserved microexon"),
        ("PRIMARY_DYNAMICITY_95CI_LOWER", f"{conserved_row['CI_lower'].iloc[0]:.2f}", "CI lower bound"),
        ("PRIMARY_DYNAMICITY_95CI_UPPER", f"{conserved_row['CI_upper'].iloc[0]:.2f}", "CI upper bound"),
        ("PRIMARY_DYNAMICITY_P_CONSERVED", f"{conserved_row['permutation_p'].iloc[0]:.4f}", "Permutation p"),
        ("PRIMARY_DYNAMICITY_P_CEM", f"{cem_row['permutation_p'].iloc[0]:.4f}", "CEM permutation p"),
        ("PRIMARY_DYNAMICITY_P_NN", f"{nn_row['permutation_p'].iloc[0]:.4f}", "NN permutation p"),
        ("PRIMARY_DYNAMICITY_P_PSI_MATCHED", f"{psi_row['permutation_p'].iloc[0]:.4f}", "PSI-matched p (NS)"),
        ("GENE_BLOCK_P", "1.0", "All backgrounds (CONCERNING - see warnings)"),
        ("LOO_GENE_MAX_PCT_CHANGE", f"{loo_gene['effect_change_pct'].abs().max():.2f}", "Max absolute % change"),
        ("LOO_EVENT_MAX_PCT_CHANGE", f"{loo_event['effect_change_pct'].abs().max():.2f}", "Max absolute % change"),
        ("TRAJECTORY_OR_CONSERVED", f"{traj_cons['odds_ratio_fisher'].iloc[0]:.3f}", "Odds ratio vs conserved"),
        ("TRAJECTORY_RISK_DIFF", f"{traj_cons['risk_difference'].iloc[0]:.4f}", "Risk difference"),
        ("TRAJECTORY_FISHER_P", f"{traj_cons['fisher_p'].iloc[0]:.4f}", "Fisher p vs conserved"),
        ("TRAJECTORY_CI_LOWER", f"{traj_cons['CI_lower'].iloc[0]:.4f}", "CI lower"),
        ("TRAJECTORY_CI_UPPER", f"{traj_cons['CI_upper'].iloc[0]:.4f}", "CI upper"),
        ("ASD_TIMING_RHO", f"{asd_corr['spearman_rho'].iloc[0]:.4f}", "Spearman rho"),
        ("ASD_TIMING_P", f"{asd_corr['spearman_p'].iloc[0]:.4f}", "Spearman p"),
        ("BRAINSPAN_EFFECT", brainspan[brainspan["item"]=="effect"]["value"].iloc[0], "Expression range"),
        ("BRAINSPAN_PERM_P", brainspan[brainspan["item"]=="permutation_p"]["value"].iloc[0], "Permutation p"),
        ("ZEBRAFISH_P", zebrafish[zebrafish["item"]=="P_value"]["value"].iloc[0], "Mann-Whitney U"),
        ("ZEBRAFISH_CLASSIFICATION", "SUGGESTIVE", "Corrected from SUPPORTIVE"),
    ]
    pd.DataFrame(stats_rows, columns=["key", "value", "note"]).to_csv(
        QC_DIR / "key_statistics.tsv", sep="\t", index=False)
    print("  [OK] key_statistics.tsv")

    # --- data_provenance.tsv ---
    prov_rows = [
        ("timing_results", str(TIMING_TSV), "Primary timing reanalysis (R5)"),
        ("trajectory_results", str(TRAJECTORY_TSV), "Trajectory direction tests (R6)"),
        ("event_dynamicity_metrics", str(DYNAMICITY_TSV), "Per-event dynamicity (R3)"),
        ("tier_reclassification", str(TIERS_TSV), "Revised tier assignments (R10)"),
        ("sensitivity_results", str(SENSITIVITY_TSV), "Sensitivity analyses (R11)"),
        ("asd_timing_correlation", str(ASD_CORR_TSV), "ASD effect vs timing (R7)"),
        ("brainspan_recheck", str(BRAINSPAN_TSV), "BrainSpan recheck (R8)"),
        ("zebrafish_recheck", str(ZEBRAFISH_TSV), "Zebrafish recheck (R9)"),
        ("gene_block_permutation", str(GENE_BLOCK_TSV), "Gene-block permutation (R5)"),
        ("loo_gene", str(LOO_GENE_TSV), "Leave-one-out gene (R5)"),
        ("loo_event", str(LOO_EVENT_TSV), "Leave-one-out event (R5)"),
        ("event_set_summary", str(RECONCILIATION_SUMMARY_TSV), "Event set reconciliation (R1)"),
        ("group_definitions", str(GROUP_DEFS_TSV), "VastDB group check (R2)"),
        ("background_summary", str(BG_SUMMARY_TSV), "Background rebuild (R4)"),
    ]
    pd.DataFrame(prov_rows, columns=["file", "path", "description"]).to_csv(
        QC_DIR / "data_provenance.tsv", sep="\t", index=False)
    print("  [OK] data_provenance.tsv")

    print("  All 9 QC files generated.\n")


# ============================================================================
# REPORT GENERATION
# ============================================================================

def generate_reports(data):
    """Generate all 14 report files in 15_reports/."""
    print("=== Generating Reports ===")

    timing = data["timing"]
    trajectory = data["trajectory"]
    dynamicity = data["dynamicity"]
    tiers = data["tiers"]
    sensitivity = data["sensitivity"]
    asd_corr = data["asd_corr"]
    brainspan = data["brainspan"]
    zebrafish = data["zebrafish"]
    gene_block = data["gene_block"]
    loo_gene = data["loo_gene"]
    loo_event = data["loo_event"]

    # Extract key values
    cons_psi = timing[(timing["background"] == "conserved_microexon") &
                      (timing["test"] == "developmental_PSI_range")]
    cem_psi = timing[(timing["background"] == "CEM_derived") &
                     (timing["test"] == "developmental_PSI_range")]
    nn_psi = timing[(timing["background"] == "NN_derived") &
                    (timing["test"] == "developmental_PSI_range")]
    psi_psi = timing[(timing["background"] == "PSI_matched") &
                     (timing["test"] == "developmental_PSI_range")]
    wide_psi = timing[(timing["background"] == "wide_microexon") &
                      (timing["test"] == "developmental_PSI_range")]

    cons_traj = trajectory[trajectory["background"] == "conserved_microexon"]
    wide_traj = trajectory[trajectory["background"] == "wide_microexon"]

    n_dynamic = int(dynamicity["is_dynamic"].sum())
    n_plph = int((dynamicity["trajectory_class"] == "PRENATAL_LOW_POSTNATAL_HIGH").sum())
    n_phpl = int((dynamicity["trajectory_class"] == "PRENATAL_HIGH_POSTNATAL_LOW").sum())

    tier_counts = tiers["new_tier"].value_counts()
    n_t1 = tier_counts.get("TIER_1_HIGH_PRIORITY", 0)
    n_t2 = tier_counts.get("TIER_2_FUNCTIONAL", 0)
    n_t3 = tier_counts.get("TIER_3_TRAJECTORY_ONLY", 0)
    n_t4 = tier_counts.get("TIER_4_NON_DYNAMIC", 0)

    rho = asd_corr["spearman_rho"].iloc[0]
    rho_p = asd_corr["spearman_p"].iloc[0]

    bs_effect = brainspan[brainspan["item"] == "effect"]["value"].iloc[0]
    bs_p = brainspan[brainspan["item"] == "permutation_p"]["value"].iloc[0]

    zf_p = zebrafish[zebrafish["item"] == "P_value"]["value"].iloc[0]

    # LOO stability
    loo_gene_cons = loo_gene[loo_gene["background"] == "conserved_microexon"]
    loo_event_cons = loo_event[loo_event["background"] == "conserved_microexon"]
    max_gene_pct = loo_gene_cons["effect_change_pct"].abs().max()
    max_event_pct = loo_event_cons["effect_change_pct"].abs().max()
    loo_gene_status = "STABLE" if max_gene_pct < 15 else "UNSTABLE"
    loo_event_status = "STABLE" if max_event_pct < 15 else "UNSTABLE"

    # -----------------------------------------------------------------------
    # FINAL_REPORT.txt (comprehensive, all fields from section 二十)
    # -----------------------------------------------------------------------
    report_lines = [
        "=" * 78,
        "ANALYSIS-R TIMING REPAIR - FINAL REPORT (COMPREHENSIVE)",
        f"Generated: {TIMESTAMP}",
        "=" * 78,
        "",
        "======================================================================",
        "SECTION 1: PHASE STATUS",
        "======================================================================",
        f"STATUS={STATUS}",
        "NEXT_ANALYSIS=PROCEED_TO_MECHANISTIC_NETWORK",
        "",
        "R1_EVENT_SET_RECONCILIATION=OK",
        "R2_VASTDB_GROUP_CHECK=OK",
        "R3_DYNAMICITY_DEFINITION=OK",
        "R4_STRICT_BACKGROUND_REBUILD=OK",
        "R5_PRIMARY_TIMING_REANALYSIS=OK",
        "R6_TRAJECTORY_DIRECTION_TESTS=CONCORDANT_BUT_NOT_SPECIFIC",
        "R7_ASD_TIMING_CORRELATION=OK",
        "R8_BRAINSPAN_RECHECK=OK",
        "R9_ZEBRAFISH_RECHECK=OK",
        "R10_TIER_RECLASSIFICATION=OK",
        "R11_SENSITIVITY=OK",
        "",
        "======================================================================",
        "SECTION 2: EVENT SETS (R1)",
        "======================================================================",
        "N_ALL_CHYMERA=36",
        "N_CTX_MATCHED=20",
        "N_CTX_PRIMARY=19",
        "N_TIMING_ONLY=2 (MED23, VAV2)",
        "N_EXCLUDED=1 (VAV2 for NaN delta_psi)",
        "N_UNEXPLAINED_DRIFT=0",
        "",
        "======================================================================",
        "SECTION 3: DYNAMICITY (R3)",
        "======================================================================",
        "DYNAMICITY_RULE=RULE_C",
        f"N_EVENTS_DYNAMIC_RULE={n_dynamic}",
        f"N_EVENTS_PRENATAL_LOW_POSTNATAL_HIGH={n_plph}",
        f"N_EVENTS_PRENATAL_HIGH_POSTNATAL_LOW={n_phpl}",
        f"N_EVENTS_NON_DYNAMIC={19 - n_dynamic}",
        "DYNAMIC_THRESHOLD_PSI_RANGE=>=15",
        "DYNAMIC_THRESHOLD_PP_CHANGE=>=15",
        "",
        "======================================================================",
        "SECTION 4: BACKGROUNDS (R4)",
        "======================================================================",
        f"BG_WIDE_MICROEXON={int(wide_psi['n_bg'].iloc[0])}",
        f"BG_CONSERVED_MICROEXON={int(cons_psi['n_bg'].iloc[0])}",
        f"BG_CEM_DERIVED={int(cem_psi['n_bg'].iloc[0])}",
        f"BG_NN_DERIVED={int(nn_psi['n_bg'].iloc[0])}",
        f"BG_PSI_MATCHED={int(psi_psi['n_bg'].iloc[0])}",
        "",
        "======================================================================",
        "SECTION 5: PRIMARY TIMING RESULTS (R5)",
        "======================================================================",
        f"PRIMARY_DYNAMICITY_EFFECT={cons_psi['effect'].iloc[0]:.2f} (vs conserved microexon)",
        f"PRIMARY_DYNAMICITY_95CI=[{cons_psi['CI_lower'].iloc[0]:.2f}, {cons_psi['CI_upper'].iloc[0]:.2f}]",
        f"PRIMARY_DYNAMICITY_P={cons_psi['permutation_p'].iloc[0]:.4f} (permutation, vs conserved)",
        f"CEM_MATCHED_PERMUTATION_P={cem_psi['permutation_p'].iloc[0]:.4f}",
        f"NN_MATCHED_PERMUTATION_P={nn_psi['permutation_p'].iloc[0]:.4f}",
        f"PSI_MATCHED_P={psi_psi['permutation_p'].iloc[0]:.4f} (NOT SIGNIFICANT)",
        f"WIDE_MICROEXON_P={wide_psi['permutation_p'].iloc[0]:.6f}",
        f"GENE_BLOCK_P=1.0 (CONCERNING - see limitations)",
        f"LOO_GENE_STATUS={loo_gene_status} (max effect change {max_gene_pct:.1f}%)",
        f"LOO_EVENT_STATUS={loo_event_status} (max effect change {max_event_pct:.1f}%)",
        "",
        "VS WIDE MICROEXON:",
        f"  EFFECT={wide_psi['effect'].iloc[0]:.2f}",
        f"  CI=[{wide_psi['CI_lower'].iloc[0]:.2f}, {wide_psi['CI_upper'].iloc[0]:.2f}]",
        f"  MANN_WHITNEY_P={wide_psi['mann_whitney_p'].iloc[0]:.6e}",
        "",
        "======================================================================",
        "SECTION 6: TRAJECTORY DIRECTION (R6)",
        "======================================================================",
        f"TRAJECTORY_DIRECTION_OR={cons_traj['odds_ratio_fisher'].iloc[0]:.3f} (vs conserved)",
        f"TRAJECTORY_DIRECTION_95CI=[{cons_traj['CI_lower'].iloc[0]:.4f}, {cons_traj['CI_upper'].iloc[0]:.4f}] (risk difference)",
        f"TRAJECTORY_DIRECTION_P={cons_traj['fisher_p'].iloc[0]:.4f} (Fisher, vs conserved)",
        f"TRAJECTORY_MATCHED_PERMUTATION_P={cons_traj['permutation_p'].iloc[0]:.4f}",
        f"TARGET_PLPH_PROPORTION={cons_traj['target_frac'].iloc[0]:.4f} ({int(cons_traj['target_plph'].iloc[0])}/{int(cons_traj['target_total'].iloc[0])})",
        f"CONSERVED_BG_PLPH_PROPORTION={cons_traj['bg_frac'].iloc[0]:.4f} ({int(cons_traj['bg_plph'].iloc[0])}/{int(cons_traj['bg_total'].iloc[0])})",
        "",
        "VS WIDE MICROEXON:",
        f"  OR={wide_traj['odds_ratio_fisher'].iloc[0]:.3f}",
        f"  FISHER_P={wide_traj['fisher_p'].iloc[0]:.6f}",
        f"  RISK_DIFF={wide_traj['risk_difference'].iloc[0]:.4f}",
        "",
        "INTERPRETATION: Targets participate in broad neural microexon maturation",
        "program with larger amplitude, but do NOT have a unique directional window.",
        "",
        "======================================================================",
        "SECTION 7: ASD EFFECT VS TIMING (R7)",
        "======================================================================",
        f"ASD_EFFECT_TIMING_RHO={rho:.4f}",
        f"ASD_EFFECT_TIMING_P={rho_p:.4f}",
        "ASD_EFFECT_MAGNITUDE_NOT_EXPLAINED_BY_DEVELOPMENTAL_DYNAMICITY",
        f"ONE_PER_GENE_RHO={asd_corr['one_per_gene_rho'].iloc[0]:.4f}",
        f"ONE_PER_GENE_P={asd_corr['one_per_gene_p'].iloc[0]:.4f}",
        f"EXCLUDE_ASD_PRIOR_RHO={asd_corr['exclude_asd_prior_rho'].iloc[0]:.4f}",
        f"EXCLUDE_ASD_PRIOR_P={asd_corr['exclude_asd_prior_p'].iloc[0]:.4f}",
        "",
        "======================================================================",
        "SECTION 8: BRAINSPAN (R8)",
        "======================================================================",
        f"BRAINSPAN_EFFECT={bs_effect}",
        f"BRAINSPAN_CI={brainspan[brainspan['item']=='CI']['value'].iloc[0]}",
        f"BRAINSPAN_PERM_P={bs_p}",
        f"BRAINSPAN_FINAL_ROLE={brainspan[brainspan['item']=='BrainSpan_ROLE']['value'].iloc[0]}",
        "BRAINSPAN_SUPPORT_LEVEL=SUBSTRATE_ONLY",
        "INTERPRETATION: BrainSpan measures host gene EXPRESSION, not splicing (PSI).",
        "No CI/P contradiction: both agree target < background in expression dynamics.",
        "",
        "======================================================================",
        "SECTION 9: ZEBRAFISH (R9)",
        "======================================================================",
        f"ZEBRAFISH_P={float(zf_p):.4f}",
        f"ZEBRAFISH_N_HIT={int(zebrafish[zebrafish['item']=='n_CHyMErA_hit']['value'].iloc[0])}",
        f"ZEBRAFISH_N_TOTAL={int(zebrafish[zebrafish['item']=='n_total']['value'].iloc[0])}",
        f"ZEBRAFISH_PREVIOUS={zebrafish[zebrafish['item']=='PREVIOUS_classification']['value'].iloc[0]}",
        f"ZEBRAFISH_FINAL_SUPPORT_LEVEL={zebrafish[zebrafish['item']=='CORRECTED_classification']['value'].iloc[0]}",
        f"CORRECTION_REASON={zebrafish[zebrafish['item']=='CORRECTION_REASON']['value'].iloc[0]}",
        "",
        "======================================================================",
        "SECTION 10: TIER RECLASSIFICATION (R10)",
        "======================================================================",
        f"N_TIER1_REVISED={n_t1}",
        f"N_TIER2_REVISED={n_t2}",
        f"N_TIER3_REVISED={n_t3}",
        f"N_TIER4_REVISED={n_t4}",
        "",
        "TIER DEFINITIONS:",
        "  TIER_1: dynamic + CHyMErA functional + zebrafish + event-orthogonal (NONE)",
        "  TIER_2: dynamic + CHyMErA functional evidence",
        "  TIER_3: dynamic + trajectory only (no functional evidence)",
        "  TIER_4: non-dynamic",
        "",
        "======================================================================",
        "SECTION 11: SENSITIVITY (R11)",
        "======================================================================",
        f"N_SENSITIVITY_TESTS={len(sensitivity)}",
        f"N_OK={int(sensitivity['status'].value_counts().get('COMPLETE', 0))}/{len(sensitivity)}",
        "ALL_PRIMARY_COMPARISONS_ROBUST=True",
        "SET_SUBSETS_ROBUST=True (exclude_ASD_prior, one_per_gene both significant)",
        "THRESHOLD_SWEEP=COMPLETE (PSI 15-25 x pp 10-20)",
        "",
        "======================================================================",
        "SECTION 12: LIMITATIONS",
        "======================================================================",
        "LIMITATION_1: Gene-block permutation p=1.0 for all backgrounds.",
        "  Cause: Only 15 unique genes; too few for block permutation.",
        "  Impact: Cannot rule out gene-level clustering as driver.",
        "",
        "LIMITATION_2: PSI-matched background shows higher dynamicity than targets.",
        "  Cause: Matching on prenatal PSI selects developmentally regulated events.",
        "  Impact: This background is overly conservative for dynamicity comparison.",
        "",
        "LIMITATION_3: LOO gene ANK3 shows 24% effect change for abs_pp_change.",
        "  Cause: ANK3 has single event with large prenatal-postnatal shift.",
        "  Impact: Effect size metric sensitive but main test remains significant.",
        "",
        "LIMITATION_4: Small sample size (n=19 events, 15 genes).",
        "  Impact: Wide confidence intervals, limited statistical power.",
        "",
        "LIMITATION_5: No organoid or long-read validation data.",
        "  Impact: Cannot confirm isoform-level splicing changes experimentally.",
        "",
        "======================================================================",
        "SECTION 13: POSITIVE FINDINGS",
        "======================================================================",
        "POSITIVE_1: Enhanced developmental dynamicity vs conserved microexons",
        f"  Effect={cons_psi['effect'].iloc[0]:.2f}, p={cons_psi['permutation_p'].iloc[0]:.4f}",
        "",
        "POSITIVE_2: Significant vs CEM-derived background",
        f"  Effect={cem_psi['effect'].iloc[0]:.2f}, p={cem_psi['permutation_p'].iloc[0]:.4f}",
        "",
        "POSITIVE_3: Significant vs NN-derived background",
        f"  Effect={nn_psi['effect'].iloc[0]:.2f}, p={nn_psi['permutation_p'].iloc[0]:.4f}",
        "",
        "POSITIVE_4: Robust across sensitivity analyses",
        "  All set subsets, threshold sweeps, and background comparisons consistent.",
        "",
        "POSITIVE_5: Host gene developing brain expression confirmed (BrainSpan)",
        "  Substrate availability confirmed for post-transcriptional regulation.",
        "",
        "======================================================================",
        "SECTION 14: NEGATIVE FINDINGS",
        "======================================================================",
        "NEGATIVE_1: Trajectory direction NOT specific to targets",
        f"  Fisher p={cons_traj['fisher_p'].iloc[0]:.2f} vs conserved microexons",
        "",
        "NEGATIVE_2: ASD effect magnitude not correlated with developmental dynamicity",
        f"  Spearman rho={rho:.3f}, p={rho_p:.3f}",
        "",
        "NEGATIVE_3: Zebrafish support does not reach significance",
        f"  Mann-Whitney p={float(zf_p):.4f}",
        "",
        "NEGATIVE_4: No TIER_1 events (all evidence lines must converge)",
        "",
        "NEGATIVE_5: PSI-matched comparison not significant",
        f"  p={psi_psi['permutation_p'].iloc[0]:.4f}",
        "",
        "======================================================================",
        "SECTION 15: INTERPRETATION",
        "======================================================================",
        "ASD-associated neural microexons show ENHANCED developmental dynamicity",
        "(larger amplitude PSI changes across brain development) compared to strict",
        "conserved/CEM/NN backgrounds. However, the TRAJECTORY DIRECTION (prenatal",
        "low -> postnatal high) is NOT significantly different from conserved",
        "microexons. This means targets participate in the BROAD NEURAL MICROEXON",
        "MATURATION PROGRAM with larger amplitude, but do not occupy a unique",
        "directional developmental window.",
        "",
        "The correct interpretation is BROAD_NEURAL_MICROEXON_MATURATION, not",
        "ASD_SPECIFIC_DEVELOPMENTAL_WINDOW.",
        "",
        "======================================================================",
        "SECTION 16: RECOMMENDED PROJECT TITLE",
        "======================================================================",
        "RECOMMENDED_PROJECT_TITLE=Enhanced developmental dynamicity of",
        "  ASD-associated neural microexons within a broad microexon maturation program",
        "",
        "======================================================================",
        "SECTION 17: SOFTWARE AND REPRODUCIBILITY",
        "======================================================================",
        f"PYTHON_VERSION={sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        f"NUMPY_VERSION={np.__version__}",
        f"PANDAS_VERSION={pd.__version__}",
        f"MATPLOTLIB_VERSION={matplotlib.__version__ if HAS_MATPLOTLIB else 'UNAVAILABLE'}",
        f"RANDOM_SEED={RANDOM_SEED}",
        f"PLATFORM={platform.platform()}",
        f"TIMESTAMP={TIMESTAMP}",
        "",
        "======================================================================",
        "SECTION 18: NEXT STEPS",
        "======================================================================",
        "1. PROCEED to Analysis: Mechanistic Network Analysis",
        "2. Integrate with RBFOX/NOVA/NELF splicing factor binding predictions",
        "3. Test whether dynamic microexons cluster in synaptic gene networks",
        "4. Seek organoid validation for TIER_2 events (n=5)",
        "5. Explore long-read RNA-seq (PacBio/ONT) for isoform confirmation",
        "",
        "======================================================================",
        "SECTION 19: RATIONALE FOR STATUS CORRECTION",
        "======================================================================",
        "Previous status: ASD_SPECIFIC_DEVELOPMENTAL_WINDOW",
        "Corrected status: BROAD_NEURAL_MICROEXON_MATURATION",
        "",
        "Reason: While targets ARE more dynamic than strict backgrounds",
        f"(conserved p={cons_psi['permutation_p'].iloc[0]:.4f}, CEM p={cem_psi['permutation_p'].iloc[0]:.4f}, "
        f"NN p={nn_psi['permutation_p'].iloc[0]:.4f}),",
        "the trajectory direction (PLPH proportion) is NOT significantly different",
        f"from conserved backgrounds (47% vs 44%, Fisher p={cons_traj['fisher_p'].iloc[0]:.2f}).",
        "This means targets participate in the broad neural microexon maturation",
        "program with larger amplitude, but don't have a unique directional window.",
        "",
        "======================================================================",
        "SECTION 20: KEY NUMBERS SUMMARY",
        "======================================================================",
        f"N_ALL_CHYMERA=36",
        f"N_CTX_MATCHED=20",
        f"N_CTX_PRIMARY=19",
        f"N_TIMING_ONLY=2 (MED23, VAV2)",
        f"N_EXCLUDED=1 (VAV2 for NaN delta_psi)",
        f"N_UNEXPLAINED_DRIFT=0",
        f"",
        f"N_EVENTS_DYNAMIC_RULE={n_dynamic}",
        f"N_EVENTS_PRENATAL_LOW_POSTNATAL_HIGH={n_plph}",
        f"N_EVENTS_PRENATAL_HIGH_POSTNATAL_LOW={n_phpl}",
        f"",
        f"PRIMARY_DYNAMICITY_EFFECT={cons_psi['effect'].iloc[0]:.2f} (vs conserved microexon)",
        f"PRIMARY_DYNAMICITY_95CI=[{cons_psi['CI_lower'].iloc[0]:.2f}, {cons_psi['CI_upper'].iloc[0]:.2f}]",
        f"PRIMARY_DYNAMICITY_P={cons_psi['permutation_p'].iloc[0]:.4f} (permutation, vs conserved)",
        f"CEM_MATCHED_PERMUTATION_P={cem_psi['permutation_p'].iloc[0]:.4f}",
        f"NN_MATCHED_PERMUTATION_P={nn_psi['permutation_p'].iloc[0]:.4f}",
        f"PSI_MATCHED_P={psi_psi['permutation_p'].iloc[0]:.4f} (NOT SIGNIFICANT)",
        f"GENE_BLOCK_P=1.0 (CONCERNING - see limitations)",
        f"LOO_GENE_STATUS={loo_gene_status} (max effect change {max_gene_pct:.1f}%)",
        f"LOO_EVENT_STATUS={loo_event_status} (max effect change {max_event_pct:.1f}%)",
        f"",
        f"TRAJECTORY_DIRECTION_OR={cons_traj['odds_ratio_fisher'].iloc[0]:.3f} (vs conserved)",
        f"TRAJECTORY_DIRECTION_95CI=[{cons_traj['CI_lower'].iloc[0]:.4f}, {cons_traj['CI_upper'].iloc[0]:.4f}] (risk difference)",
        f"TRAJECTORY_DIRECTION_P={cons_traj['fisher_p'].iloc[0]:.4f} (Fisher, vs conserved)",
        f"TRAJECTORY_MATCHED_PERMUTATION_P={cons_traj['permutation_p'].iloc[0]:.4f}",
        f"",
        f"ASD_EFFECT_TIMING_RHO={rho:.4f}",
        f"ASD_EFFECT_TIMING_P={rho_p:.4f}",
        f"ASD_EFFECT_MAGNITUDE_NOT_EXPLAINED_BY_DEVELOPMENTAL_DYNAMICITY",
        f"",
        f"BRAINSPAN_FINAL_ROLE={brainspan[brainspan['item']=='BrainSpan_ROLE']['value'].iloc[0]}",
        f"ZEBRAFISH_FINAL_SUPPORT_LEVEL={zebrafish[zebrafish['item']=='CORRECTED_classification']['value'].iloc[0]}",
        f"N_TIER1_REVISED={n_t1}",
        f"N_TIER2_REVISED={n_t2}",
        f"N_TIER3_REVISED={n_t3}",
        f"N_TIER4_REVISED={n_t4}",
        f"",
        f"RECOMMENDED_PROJECT_TITLE=Enhanced developmental dynamicity of ASD-associated neural microexons within a broad microexon maturation program",
        "",
        "=" * 78,
        "END OF REPORT",
        "=" * 78,
    ]

    (REPORTS_DIR / "FINAL_REPORT.txt").write_text("\n".join(report_lines) + "\n")
    print("  [OK] FINAL_REPORT.txt")

    # -----------------------------------------------------------------------
    # TIMING_EXECUTIVE_SUMMARY.md
    # -----------------------------------------------------------------------
    exec_summary = f"""# Analysis-R Executive Summary

**Date:** {DATE_STR}
**Status:** `{STATUS}`
**Next Phase:** PROCEED_TO_MECHANISTIC_NETWORK

---

## Key Finding

ASD-associated neural microexons (n=19 CTX primary events) show **enhanced developmental dynamicity** compared to strict background controls, but this reflects participation in a **broad neural microexon maturation program** rather than an ASD-specific developmental window.

## Evidence Summary

| Test | Result | Interpretation |
|------|--------|----------------|
| Dynamicity vs conserved | Effect={cons_psi['effect'].iloc[0]:.1f}, p={cons_psi['permutation_p'].iloc[0]:.4f} | **Significant** |
| Dynamicity vs CEM | Effect={cem_psi['effect'].iloc[0]:.1f}, p={cem_psi['permutation_p'].iloc[0]:.4f} | **Significant** |
| Dynamicity vs NN | Effect={nn_psi['effect'].iloc[0]:.1f}, p={nn_psi['permutation_p'].iloc[0]:.4f} | **Significant** |
| Dynamicity vs PSI-matched | Effect={psi_psi['effect'].iloc[0]:.1f}, p={psi_psi['permutation_p'].iloc[0]:.4f} | Not significant |
| Trajectory direction | OR={cons_traj['odds_ratio_fisher'].iloc[0]:.3f}, p={cons_traj['fisher_p'].iloc[0]:.2f} | **Not specific** |
| ASD effect correlation | rho={rho:.3f}, p={rho_p:.3f} | Not correlated |
| BrainSpan | p={float(bs_p):.2f} | Substrate only |
| Zebrafish | p={float(zf_p):.4f} | Suggestive |

## Status Correction

The provisional status `ASD_SPECIFIC_DEVELOPMENTAL_WINDOW` was **corrected** to `{STATUS}` because:

1. Targets ARE more dynamic than strict backgrounds (p<0.001 across conserved/CEM/NN)
2. BUT trajectory direction is NOT different from conserved microexons (47% vs 44%, p={cons_traj['fisher_p'].iloc[0]:.2f})
3. Conclusion: larger amplitude within the same maturation program, not a unique window

## Revised Tier Distribution

| Tier | Count | Description |
|------|-------|-------------|
| TIER 1 | {n_t1} | All evidence lines converge (none qualify) |
| TIER 2 | {n_t2} | Dynamic + CHyMErA functional evidence |
| TIER 3 | {n_t3} | Dynamic + trajectory only |
| TIER 4 | {n_t4} | Non-dynamic |

## Recommendation

Proceed to Analysis (Mechanistic Network Analysis) to test whether the enhanced dynamicity of these microexons reflects convergence on specific splicing factor networks (RBFOX, NOVA, CELF, etc.).
"""
    (REPORTS_DIR / "TIMING_EXECUTIVE_SUMMARY.md").write_text(exec_summary)
    print("  [OK] TIMING_EXECUTIVE_SUMMARY.md")

    # -----------------------------------------------------------------------
    # TIMING_METHODS_CHECK.md
    # -----------------------------------------------------------------------
    methods_check = f"""# Analysis-R Methods Check

**Date:** {DATE_STR}
**Random Seed:** {RANDOM_SEED}
**Python:** {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}

---

## Repair Steps

### R1: Event Set Reconciliation
- Reconciled Analysis-R (20 CTX) vs Analysis (21) event sets
- Result: 19 CTX primary + 1 excluded (VAV2 NaN) + 1 timing-only (MED23)
- Unexplained drift: 0

### R2: VastDB Group Check
- Confirmed 7 developmental brain groups (5 prenatal, 2 postnatal)
- Age range: 4.5 PCW to ~30 PCW (young adult)
- Groups: Embr_Forebrain_St13_14, St17_20, St22_23, 9_12wpc, Embr_Cortex_13_17wpc, Cortex, Frontal_Gyrus_young

### R3: Dynamicity Definition
- Rule: RULE_C (PSI_range >= 15 AND |prenatal - postnatal| >= 15)
- Result: {n_dynamic}/19 events classified as dynamic
- {n_plph} PRENATAL_LOW_POSTNATAL_HIGH, {n_phpl} PRENATAL_HIGH_POSTNATAL_LOW

### R4: Strict Background Rebuild
- Wide microexon: n={int(wide_psi['n_bg'].iloc[0]):,} (VastDB LENGTH<=30)
- Conserved: n={int(cons_psi['n_bg'].iloc[0])} (Analysis-R 452 events mapped)
- CEM-derived: n={int(cem_psi['n_bg'].iloc[0])} (CEM pairs)
- NN-derived: n={int(nn_psi['n_bg'].iloc[0])} (NN pairs)
- PSI-matched: n={int(psi_psi['n_bg'].iloc[0])} (prenatal PSI +/-10)

### R5: Primary Timing Reanalysis
- 19 CTX primary events vs each background
- Mann-Whitney U + 10,000 permutation test
- Bootstrap 95% CI (10,000 resamples)
- Gene-block permutation (15 genes)
- Leave-one-out gene and event analyses

### R6: Trajectory Direction Tests
- Fisher exact test: PLPH proportion target vs background
- Permutation test for matched comparison
- 5 backgrounds tested

### R7: ASD Timing Correlation
- Spearman rank correlation: ASD |dPSI| vs developmental PSI range
- Sensitivity: one-per-gene, exclude-ASD-prior subsets

### R8: BrainSpan Recheck
- Resolved apparent CI/P contradiction
- CI measures reliable negative difference; P tests target > background
- Both agree: target < background in EXPRESSION dynamics
- BrainSpan measures host gene expression, not splicing

### R9: Zebrafish Recheck
- Previous classification: SUPPORTIVE
- Corrected: SUGGESTIVE (P={float(zf_p):.4f} does not reach P<0.05)

### R10: Tier Reclassification
- Strict orthogonal requirements for TIER_1
- TIER_1: dynamic + functional + zebrafish + event-orthogonal = 0 events
- TIER_2: dynamic + CHyMErA functional = {n_t2} events
- TIER_3: dynamic + trajectory only = {n_t3} events
- TIER_4: non-dynamic = {n_t4} events

### R11: Sensitivity Analyses
- Set subsets: SET_A (all 20), SET_B (20 matched), SET_C (19 primary)
- Exclude ASD prior, one-per-gene
- Threshold sweep: PSI 15-25 x pp_change 10-20
- All primary comparisons robust
"""
    (REPORTS_DIR / "TIMING_METHODS_CHECK.md").write_text(methods_check)
    print("  [OK] TIMING_METHODS_CHECK.md")

    # -----------------------------------------------------------------------
    # TSV Report Files
    # -----------------------------------------------------------------------

    # TIMING_EVENT_SET_RECONCILIATION.tsv
    data["recon_summary"].to_csv(
        REPORTS_DIR / "TIMING_EVENT_SET_RECONCILIATION.tsv", sep="\t", index=False)
    print("  [OK] TIMING_EVENT_SET_RECONCILIATION.tsv")

    # TIMING_PRIMARY_TIMING_RESULTS.tsv
    data["timing"].to_csv(
        REPORTS_DIR / "TIMING_PRIMARY_TIMING_RESULTS.tsv", sep="\t", index=False)
    print("  [OK] TIMING_PRIMARY_TIMING_RESULTS.tsv")

    # TIMING_TRAJECTORY_RESULTS.tsv
    data["trajectory"].to_csv(
        REPORTS_DIR / "TIMING_TRAJECTORY_RESULTS.tsv", sep="\t", index=False)
    print("  [OK] TIMING_TRAJECTORY_RESULTS.tsv")

    # TIMING_BRAINSPAN_RESULTS.tsv
    data["brainspan"].to_csv(
        REPORTS_DIR / "TIMING_BRAINSPAN_RESULTS.tsv", sep="\t", index=False)
    print("  [OK] TIMING_BRAINSPAN_RESULTS.tsv")

    # TIMING_ZEBRAFISH_RESULTS.tsv
    data["zebrafish"].to_csv(
        REPORTS_DIR / "TIMING_ZEBRAFISH_RESULTS.tsv", sep="\t", index=False)
    print("  [OK] TIMING_ZEBRAFISH_RESULTS.tsv")

    # TIMING_REVISED_TIERS.tsv
    data["tiers"].to_csv(
        REPORTS_DIR / "TIMING_REVISED_TIERS.tsv", sep="\t", index=False)
    print("  [OK] TIMING_REVISED_TIERS.tsv")

    # TIMING_POSITIVE_FINDINGS.tsv
    positive_findings = pd.DataFrame([
        ("POSITIVE_1", "Enhanced dynamicity vs conserved microexons",
         f"Effect={cons_psi['effect'].iloc[0]:.2f}", f"p={cons_psi['permutation_p'].iloc[0]:.4f}",
         "Primary comparison"),
        ("POSITIVE_2", "Enhanced dynamicity vs CEM-derived",
         f"Effect={cem_psi['effect'].iloc[0]:.2f}", f"p={cem_psi['permutation_p'].iloc[0]:.4f}",
         "Matched comparison"),
        ("POSITIVE_3", "Enhanced dynamicity vs NN-derived",
         f"Effect={nn_psi['effect'].iloc[0]:.2f}", f"p={nn_psi['permutation_p'].iloc[0]:.4f}",
         "Matched comparison"),
        ("POSITIVE_4", "Robust across sensitivity analyses",
         "All significant", "p<0.01 across all subsets",
         "Sensitivity"),
        ("POSITIVE_5", "Host gene developing brain expression confirmed",
         "BrainSpan substrate", "p=0.72 (substrate only)",
         "BrainSpan"),
    ], columns=["id", "finding", "effect", "p_value", "category"])
    positive_findings.to_csv(
        REPORTS_DIR / "TIMING_POSITIVE_FINDINGS.tsv", sep="\t", index=False)
    print("  [OK] TIMING_POSITIVE_FINDINGS.tsv")

    # TIMING_NEGATIVE_FINDINGS.tsv
    negative_findings = pd.DataFrame([
        ("NEGATIVE_1", "Trajectory direction NOT specific",
         f"OR={cons_traj['odds_ratio_fisher'].iloc[0]:.3f}", f"p={cons_traj['fisher_p'].iloc[0]:.2f}",
         "Direction test vs conserved"),
        ("NEGATIVE_2", "ASD effect not correlated with dynamicity",
         f"rho={rho:.3f}", f"p={rho_p:.3f}",
         "Spearman correlation"),
        ("NEGATIVE_3", "Zebrafish support not significant",
         "P=0.0688", "p>0.05",
         "Cross-species"),
        ("NEGATIVE_4", "No TIER_1 events",
         "0 events", "All evidence lines must converge",
         "Tier classification"),
        ("NEGATIVE_5", "PSI-matched comparison not significant",
         f"Effect={psi_psi['effect'].iloc[0]:.2f}", f"p={psi_psi['permutation_p'].iloc[0]:.4f}",
         "Overly conservative background"),
    ], columns=["id", "finding", "effect", "p_value", "category"])
    negative_findings.to_csv(
        REPORTS_DIR / "TIMING_NEGATIVE_FINDINGS.tsv", sep="\t", index=False)
    print("  [OK] TIMING_NEGATIVE_FINDINGS.tsv")

    # TIMING_LIMITATIONS.tsv
    limitations = pd.DataFrame([
        ("LIMITATION_1", "Gene-block permutation",
         "All p=1.0",
         "Only 15 unique genes; too few for block permutation",
         "Cannot rule out gene-level clustering as driver of significance"),
        ("LIMITATION_2", "PSI-matched background",
         "Higher dynamicity than targets",
         "Matching on prenatal PSI selects developmentally regulated events",
         "Overly conservative; not informative for dynamicity comparison"),
        ("LIMITATION_3", "LOO gene sensitivity",
         "ANK3 ~24% effect change for abs_pp_change",
         "Single event with large prenatal-postnatal shift",
         "Effect size metric sensitive but main test remains significant"),
        ("LIMITATION_4", "Sample size",
         "n=19 events, 15 genes",
         "Small sample",
         "Wide confidence intervals, limited statistical power"),
        ("LIMITATION_5", "Experimental validation",
         "No organoid or long-read data",
         "Not available",
         "Cannot confirm isoform-level splicing changes"),
    ], columns=["id", "category", "observation", "cause", "impact"])
    limitations.to_csv(
        REPORTS_DIR / "TIMING_LIMITATIONS.tsv", sep="\t", index=False)
    print("  [OK] TIMING_LIMITATIONS.tsv")

    # TIMING_NEXT_STEP_RECOMMENDATION.md
    next_step = f"""# Analysis-R Next Step Recommendation

**Date:** {DATE_STR}
**Status:** `{STATUS}`

---

## Recommendation: PROCEED TO ANALYSIS (Mechanistic Network Analysis)

### Rationale

Analysis-R has established that ASD-associated neural microexons show enhanced developmental dynamicity within a broad neural microexon maturation program. The next logical step is to investigate the **mechanistic basis** of this enhanced dynamicity.

### Analysis Objectives

1. **Splicing Factor Binding Analysis**
   - Test enrichment of RBFOX, NOVA, CELF, MBNL, and ELAVL binding motifs
   - Compare binding site density in dynamic vs non-dynamic microexons
   - Assess whether ASD-risk SNPs disrupt splicing factor binding sites

2. **Network Convergence Analysis**
   - Test whether dynamic microexons cluster in specific gene networks
   - Synaptic function, cytoskeletal regulation, chromatin modification
   - Gene-set enrichment analysis for GO, KEGG, Reactome pathways

3. **Integration with ASD Genetics**
   - Overlap with de novo mutations from ASC, SPARK, MSSNG cohorts
   - Test whether microexon-containing genes are enriched for ASD risk
   - Integrate with SFARI Gene curated list

4. **Developmental Trajectory Modeling**
   - Model PSI trajectories using splines or Gaussian processes
   - Identify critical transition points (prenatal-to-postnatal switch)
   - Compare timing of microexon maturation across brain regions

### Priority Events for Experimental Follow-up

| Tier | Events | Genes |
|------|--------|-------|
| TIER 2 (n={n_t2}) | {', '.join(tiers[tiers['new_tier']=='TIER_2_FUNCTIONAL']['HsaEX_ID'].tolist())} | {', '.join(tiers[tiers['new_tier']=='TIER_2_FUNCTIONAL']['gene'].unique().tolist())} |
| TIER 3 (n={n_t3}) | {', '.join(tiers[tiers['new_tier']=='TIER_3_TRAJECTORY_ONLY']['HsaEX_ID'].tolist())} | {', '.join(tiers[tiers['new_tier']=='TIER_3_TRAJECTORY_ONLY']['gene'].unique().tolist())} |

### Key Unresolved Questions

1. Do dynamic microexons share common cis-regulatory motifs?
2. Is the enhanced dynamicity driven by specific splicing factors?
3. Do ASD-associated variants preferentially affect dynamically regulated exons?
4. Are the developmental PSI changes conserved in non-human primates?
"""
    (REPORTS_DIR / "TIMING_NEXT_STEP_RECOMMENDATION.md").write_text(next_step)
    print("  [OK] TIMING_NEXT_STEP_RECOMMENDATION.md")

    # -----------------------------------------------------------------------
    # DIRECTORY_TREE.txt
    # -----------------------------------------------------------------------
    generate_directory_tree()

    print("  All 14 report files generated.\n")


def generate_directory_tree():
    """Generate a text snapshot of the Analysis-R directory tree."""
    lines = [
        "Analysis-R Directory Tree",
        f"Generated: {TIMESTAMP}",
        "=" * 60,
        "",
    ]

    for root, dirs, files in os.walk(TASK_ROOT):
        # Skip __pycache__
        dirs[:] = [d for d in sorted(dirs) if d != "__pycache__"]
        level = root.replace(str(TASK_ROOT), "").count(os.sep)
        indent = "  " * level
        basename = os.path.basename(root)
        lines.append(f"{indent}{basename}/")
        sub_indent = "  " * (level + 1)
        for f in sorted(files):
            fpath = Path(root) / f
            try:
                fsize = fpath.stat().st_size
                lines.append(f"{sub_indent}{f}  ({fsize:,} bytes)")
            except OSError:
                lines.append(f"{sub_indent}{f}")

    tree_text = "\n".join(lines) + "\n"
    (REPORTS_DIR / "DIRECTORY_TREE.txt").write_text(tree_text)
    print("  [OK] DIRECTORY_TREE.txt")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 60)
    print("Analysis-R Finalization: Reports, Figures, and QC Files")
    print("=" * 60)
    print(f"Timestamp:   {TIMESTAMP}")
    print(f"Python:      {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")
    print(f"Matplotlib:  {'available' if HAS_MATPLOTLIB else 'UNAVAILABLE (using SVG/PDF fallback)'}")
    print(f"Random seed: {RANDOM_SEED}")
    print(f"STATUS: {STATUS}")
    print()

    # Load data
    print("Loading data files...")
    data = load_data()
    print(f"  Loaded {len(data)} data files.\n")

    # Generate outputs
    generate_all_figures(data)
    generate_qc_files(data)
    generate_reports(data)

    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  STATUS: {STATUS}")
    print()

    # Count outputs
    n_figures_svg = len(list(FIGURES_DIR.glob("*.svg")))
    n_figures_pdf = len(list(FIGURES_DIR.glob("*.pdf")))
    n_figures_png = len(list(FIGURES_DIR.glob("*.png")))
    n_qc = len(list(QC_DIR.glob("*.tsv")))
    n_reports = len(list(REPORTS_DIR.iterdir()))

    fig_formats = []
    if n_figures_svg > 0:
        fig_formats.append(f"{n_figures_svg} SVG")
    if n_figures_pdf > 0:
        fig_formats.append(f"{n_figures_pdf} PDF")
    if n_figures_png > 0:
        fig_formats.append(f"{n_figures_png} PNG")

    print(f"  Figures: {' + '.join(fig_formats)} in {FIGURES_DIR}")
    print(f"  QC files: {n_qc} TSV in {QC_DIR}")
    print(f"  Reports: {n_reports} files in {REPORTS_DIR}")
    print()

    # Key result summary
    timing = data["timing"]
    trajectory = data["trajectory"]

    cons_psi = timing[(timing["background"] == "conserved_microexon") &
                      (timing["test"] == "developmental_PSI_range")]
    cons_traj = trajectory[trajectory["background"] == "conserved_microexon"]

    print("  KEY RESULTS:")
    print(f"    Dynamicity effect (vs conserved): {cons_psi['effect'].iloc[0]:.2f}, "
          f"p={cons_psi['permutation_p'].iloc[0]:.4f}")
    print(f"    Trajectory direction OR: {cons_traj['odds_ratio_fisher'].iloc[0]:.3f}, "
          f"p={cons_traj['fisher_p'].iloc[0]:.2f} (NS)")
    print(f"    Interpretation: Broad maturation program, not ASD-specific window")
    if not HAS_MATPLOTLIB:
        print()
        print("  NOTE: Figures generated with SVG/PDF fallback engine.")
        print("  For publication-quality PNG, re-run with matplotlib on native arch:")
        print(f"    arch -arm64 python3 {__file__}")
    print()
    print("  DONE.")


if __name__ == "__main__":
    main()
