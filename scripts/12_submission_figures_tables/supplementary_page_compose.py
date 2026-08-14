#!/usr/bin/env python3
"""Final A4 page composition (artwork + justified caption).

Deltas vs the earlier version: combined PDF is
Supplementary_Figures_FINAL_SUBMISSION_v2.pdf; page preview PNGs are
written to 05_preview_png; S12 remains A4 landscape.

Design notes follow.

Final A4 page composition (artwork + justified caption).

Each supplementary figure page carries the artwork on top and the
journal-style caption directly below on the same page.  Captions are set
by a real paragraph engine (ReportLab Paragraph with TA_JUSTIFY on the
system Helvetica TTF collection); no manual space insertion anywhere.

Geometry (spec section 9): A4 portrait (S12 landscape), artwork margins
14 mm, caption margins 20 mm, caption spans the usable width, artwork
scaled to fit the remaining height (cap 1.10x), whole block vertically
centred when leftover space exceeds 30 mm.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supplementary_common import (FIG_DIR, PAGE_DIR, PREV_DIR, PT2MM, PAGE_W_MM,
                             PAGE_H_MM, ART_MARGIN_MM, CAP_MARGIN_MM,
                             tsv_write)
import supplementary_legends as S38L
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle

HELV = "/System/Library/Fonts/Helvetica.ttc"
pdfmetrics.registerFont(TTFont("Helv", HELV, subfontIndex=0))
pdfmetrics.registerFont(TTFont("Helv-Bold", HELV, subfontIndex=1))
INK = HexColor("#1a1a1a")
ST_TITLE = ParagraphStyle("captitle", fontName="Helv-Bold", fontSize=9.2,
                          leading=11.0, alignment=TA_LEFT, textColor=INK)
ST_BODY = ParagraphStyle("capbody", fontName="Helv", fontSize=8.6,
                         leading=10.0, alignment=TA_JUSTIFY, textColor=INK)
GAP_MM = 7.0
TOP_MM = 12.0
BOT_MM = 18.0
SCALE_MAX = 1.10
LANDSCAPE = {12}


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def compose_page(c, n, title, body):
    land = n in LANDSCAPE
    W = PAGE_H_MM if land else PAGE_W_MM
    H = PAGE_W_MM if land else PAGE_H_MM
    png = os.path.join(FIG_DIR, "Figure_S%d_image.png" % n)
    im = Image.open(png)
    dpi = im.info.get("dpi", (600, 600))
    fw = im.width / float(dpi[0]) * 25.4
    fh = im.height / float(dpi[1]) * 25.4
    avail_w = W - 2 * ART_MARGIN_MM
    cap_w_pt = (W - 2 * CAP_MARGIN_MM) / PT2MM
    tp = Paragraph(esc(title), ST_TITLE)
    bp = Paragraph(esc(body), ST_BODY)
    _, th = tp.wrap(cap_w_pt, 1e9)
    _, bh = bp.wrap(cap_w_pt, 1e9)
    th_mm, bh_mm = th * PT2MM, bh * PT2MM
    cap_h = th_mm + 1.5 + bh_mm
    avail_art_h = H - TOP_MM - BOT_MM - GAP_MM - cap_h
    s = min(avail_w / fw, avail_art_h / fh, SCALE_MAX)
    art_w, art_h = fw * s, fh * s
    block = art_h + GAP_MM + cap_h
    leftover = H - TOP_MM - BOT_MM - block
    off = leftover / 2.0 if leftover > 30.0 else 0.0
    y_top = TOP_MM + off
    x_art = (W - art_w) / 2.0
    c.drawImage(png, x_art / PT2MM, (H - y_top - art_h) / PT2MM,
                width=art_w / PT2MM, height=art_h / PT2MM)
    y_title_bot = H - y_top - art_h - GAP_MM - th_mm
    tp.drawOn(c, CAP_MARGIN_MM / PT2MM, y_title_bot / PT2MM)
    bp.drawOn(c, CAP_MARGIN_MM / PT2MM, (y_title_bot - 1.5 - bh_mm) / PT2MM)
    usable = H - TOP_MM - BOT_MM
    return dict(n=n, land=land, W=W, H=H, art_w=art_w, art_h=art_h,
                scale=s, cap_h=cap_h, leftover=leftover,
                pct=100.0 * art_h / usable, centered=off > 0,
                y_top=y_top, x_art=x_art)


def make_preview(info, title, body):
    """Raster page preview (artwork + left-aligned caption) for quick
    visual inspection; the PDF remains the authoritative composition."""
    from PIL import ImageDraw, ImageFont
    DPI = 100.0
    Wpx = int(info["W"] / 25.4 * DPI)
    Hpx = int(info["H"] / 25.4 * DPI)
    img = Image.new("RGB", (Wpx, Hpx), "white")
    png = os.path.join(FIG_DIR, "Figure_S%d_image.png" % info["n"])
    art = Image.open(png).convert("RGB")
    aw = int(info["art_w"] / 25.4 * DPI)
    ah = int(info["art_h"] / 25.4 * DPI)
    art = art.resize((aw, ah), Image.LANCZOS)
    xa = int(info["x_art"] / 25.4 * DPI)
    ya = int(info["y_top"] / 25.4 * DPI)
    img.paste(art, (xa, ya))
    d = ImageDraw.Draw(img)
    fb = ImageFont.truetype(HELV, 12, index=1)
    fr = ImageFont.truetype(HELV, 11, index=0)
    capw = Wpx - 2 * int(CAP_MARGIN_MM / 25.4 * DPI)

    def wrap(text, font):
        out, line = [], ""
        for w in text.split():
            trial = (line + " " + w).strip()
            if font.getlength(trial) <= capw:
                line = trial
            else:
                out.append(line)
                line = w
        out.append(line)
        return out

    y = ya + ah + int(GAP_MM / 25.4 * DPI)
    for line in wrap(title, fb):
        d.text((int(CAP_MARGIN_MM / 25.4 * DPI), y), line, font=fb,
               fill=(26, 26, 26))
        y += 15
    y += 2
    for line in wrap(body, fr):
        d.text((int(CAP_MARGIN_MM / 25.4 * DPI), y), line, font=fr,
               fill=(26, 26, 26))
        y += 14
    img.save(os.path.join(PREV_DIR, "Figure_S%d_preview.png" % info["n"]),
             dpi=(DPI, DPI))


def main():
    caps = S38L.legends()
    rows = []
    comb_path = os.path.join(
        PAGE_DIR, "Supplementary_Figures_FINAL_SUBMISSION_v2.pdf")
    comb = canvas.Canvas(comb_path,
                         pagesize=(PAGE_W_MM / PT2MM, PAGE_H_MM / PT2MM))
    for n in range(1, 17):
        short, body = caps[n]
        title = "Figure S%d. %s" % (n, short)
        ps = ((PAGE_H_MM / PT2MM, PAGE_W_MM / PT2MM) if n in LANDSCAPE
              else (PAGE_W_MM / PT2MM, PAGE_H_MM / PT2MM))
        c = canvas.Canvas(os.path.join(PAGE_DIR, "Figure_S%d.pdf" % n),
                          pagesize=ps)
        info = compose_page(c, n, title, body)
        c.showPage()
        c.save()
        make_preview(info, title, body)
        comb.setPageSize(ps)
        info = compose_page(comb, n, title, body)
        comb.showPage()
        rows.append([
            n, "landscape" if info["land"] else "portrait",
            "%.1f" % info["art_w"], "%.1f" % info["art_h"],
            "%.3f" % info["scale"], "%.1f" % info["cap_h"],
            "%.1f" % info["pct"], "%.1f" % info["leftover"],
            "YES" if info["centered"] else "NO",
        ])
        print("page S%-2d %-9s art %5.1fx%5.1f mm scale %.3f cap %5.1f mm "
              "pct %.1f leftover %.1f centered %s" %
              (n, rows[-1][1], info["art_w"], info["art_h"], info["scale"],
               info["cap_h"], info["pct"], info["leftover"], rows[-1][8]))
    comb.save()
    tsv_write(os.path.join(PAGE_DIR, "PAGE_LAYOUT.tsv"),
              ["figure_n", "orientation", "artwork_w_mm", "artwork_h_mm",
               "scale", "caption_h_mm", "artwork_pct_usable_height",
               "leftover_mm", "block_centered"], rows)
    n_pages = 16
    print("PAGES_CREATED=%d/16" % n_pages)
    print("page compose done")


if __name__ == "__main__":
    main()
