#!/usr/bin/env python3
"""Final Figure 3 (final visual-fix revision): developmental maturation and
host-gene network context. A developmental PSI heatmap | B dynamic-event
slopegraph | C effect-vs-range null scatter | D host-gene network |
E gene-pathway membership map.

Final fixes vs the earlier version:
  - A/B: wider top-row gutter so panel B's y-axis, tick marks, left event
    labels and title never intrude into panel A (and A's artists never
    intrude into B).
  - D/E: bottom-row width ratio moved from 0.44:0.56 to 0.40:0.60 (more
    width for E) and a much wider D-E gutter; E's left margin expanded so the
    gene-name y-tick labels never intrude into D. Text sizes unchanged
    (margins increased first, per spec).
All numbers identical to the analysis sources and to the Earlier
SOURCE_DATA rows. No recomputation.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon_main import *
import matplotlib.gridspec as gridspec

MT = load_master()
EDGES = load_network_edges()
PATHS = load_pathways()
NET_RECLASS = os.path.join(DIR33, "06_pathway_network_reclassification",
                           "NETWORK_PUBLIC_INTERPRETATION.tsv")
PW_RECLASS = os.path.join(DIR33, "06_pathway_network_reclassification",
                          "PATHWAY_PUBLIC_INTERPRETATION.tsv")
NET_TSV = os.path.join(DIR14, "08_host_gene_network",
                       "03_observed_network_edges.tsv")
NET_METRICS = os.path.join(DIR14, "08_host_gene_network",
                           "04_network_metrics.tsv")
PW_TSV = os.path.join(DIR14, "09_pathway_permutation",
                      "03_matched_permutation_results.tsv")
MASTER_TSV = os.path.join(DIR25, "06_master_event_table",
                          "MASTER_19_EVENT_EVIDENCE_TABLE.tsv")
FIGDIR = FIG_DIRS[3]

# verify against the final reclassification before drawing
_net = rd(NET_RECLASS)
assert all(r["permutation_inference_status"] == "INVALID" for r in _net)
_full = [r for r in _net if r["edge_subset"] == "full_curated_network"][0]
assert int(_full["n_edges"]) == 25
_pw = rd(PW_RECLASS)
assert all(r["public_interpretation"] ==
           "DESCRIPTIVE_CURATED_FUNCTIONAL_ANNOTATION" for r in _pw)

TRAJ_COL = {"PLPH": C_DEV, "PHPL": C_ORANGE, "NON_DYNAMIC": C_BG_GREY}

GENES = sorted(set(r["gene"] for r in MT))
gene_tier = {}
for g in GENES:
    gene_tier[g] = min(tier_rank(r) for r in MT if r["gene"] == g)

# QC bookkeeping
D_GENE_LABELS = []   # gene-name Text artists in 3D
D_NODES = []         # node Circle artists in 3D
E_LEGEND = [None]    # 3E legend artist


def short(r):
    return "%s %s" % (r["gene"], r["HsaEX_ID"].replace("HsaEX", ""))


# --------------------------------------------------------------------------- A
def draw_heatmap(ax, axstrip):
    rows = sorted(MT, key=lambda r: (devrank(r), r["gene"], r["HsaEX_ID"]))
    M = np.array([[float(r["prenatal_mean_PSI"]),
                   float(r["postnatal_mean_PSI"])] for r in rows])
    ax.imshow(M, cmap="viridis", vmin=0, vmax=100, aspect="auto")
    for i in range(len(rows)):
        for j in range(2):
            v = M[i, j]
            ax.text(j, i, "%.0f" % v, ha="center", va="center", fontsize=7.0,
                    color="white" if v < 45 else C_DARK)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Prenatal", "Postnatal"])
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([short(r) for r in rows], fontsize=7.2)
    ax.set_title("Developmental PSI (%, BrainSpan)", fontsize=9.0,
                 weight="bold", color=C_DARK, loc="left")
    axstrip.axis("off")
    axstrip.set_ylim(len(rows) - 0.5, -0.5)
    axstrip.set_xlim(0, 1)
    for i, r in enumerate(rows):
        axstrip.add_patch(plt.Rectangle(
            (0.05, i - 0.45), 0.9, 0.9,
            fc=TRAJ_COL[r["developmental_trajectory"]], ec="white", lw=0.5))
    axstrip.set_yticks([])


# --------------------------------------------------------------------------- B
def draw_slope(ax):
    dyn = [r for r in MT if r["developmental_dynamic_status"] == "DYNAMIC"]
    dyn.sort(key=lambda r: -float(r["prenatal_mean_PSI"]))
    last = 999.0
    lab_y = {}
    for r in dyn:
        pre = float(r["prenatal_mean_PSI"])
        yy = min(pre, last - 6.0)
        lab_y[r["HsaEX_ID"]] = yy
        last = yy
    for r in dyn:
        pre = float(r["prenatal_mean_PSI"])
        post = float(r["postnatal_mean_PSI"])
        c = TRAJ_COL[r["developmental_trajectory"]]
        ax.plot([0, 1], [pre, post], "-", color=c, lw=1.7, alpha=0.85,
                zorder=2)
        ax.plot(0, pre, "o", ms=3.6, color=c, mec="white", mew=0.5, zorder=3)
        ax.plot(1, post, "o", ms=3.6, color=c, mec="white", mew=0.5, zorder=3)
        yy = lab_y[r["HsaEX_ID"]]
        if abs(yy - pre) > 0.5:
            ax.plot([-0.045, 0], [yy, pre], "-", lw=0.4, color=C_MID, zorder=1)
        ax.text(-0.06, yy, short(r), ha="right", va="center", fontsize=7.0,
                color=C_MID)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Prenatal", "Postnatal"])
    ax.set_yticks([0, 50, 100])
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=2)
    ax.set_ylim(0, 100)
    ax.set_xlim(-0.05, 1.25)
    # Compact trajectory legend: counts only; event identity is already
    # labelled directly in the slopegraph.
    ax.plot([], [], "-", color=C_DEV, lw=1.8, label="PLPH (n = 9)")
    ax.plot([], [], "-", color=C_ORANGE, lw=1.8, label="PHPL (n = 1)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98),
              fontsize=7.2, frameon=False, handlelength=2.2,
              handletextpad=0.7, borderaxespad=0.0, labelspacing=0.45)
    panel_title(ax, "Dynamic events", fontsize=8.7)


# --------------------------------------------------------------------------- C
def draw_negscatter(ax):
    cor = rd(os.path.join(DIR13, "09_ASD_timing_correlation",
                          "02_event_correlation_data.tsv"))
    kv = {r["key"]: r["value"] for r in
          rd(os.path.join(DIR13, "13_qc", "key_statistics.tsv"))}
    traj = {r["HsaEX_ID"]: r["developmental_trajectory"] for r in MT}
    for r in cor:
        c = TRAJ_COL[traj[r["event_id"]]]
        ax.plot(float(r["PSI_range"]), float(r["abs_delta_psi"]), "o", ms=4.2,
                color=c, mec="white", mew=0.5, zorder=3)
    # Final: deterministic in-range ticks (auto ticks placed labels
    # outside the axes)
    ax.set_xticks([20.0, 40.0, 60.0, 80.0])
    ax.set_yticks([0.0, 0.04, 0.08, 0.12])
    ax.set_xlabel("Developmental PSI range (%)")
    ax.set_ylabel("ASD cortex |ΔPSI|")
    # rho/P only; interpretation prose moved to the figure legend (spec)
    ax.text(0.97, 0.97, "Spearman ρ = %s\nP = %s" % (kv["ASD_TIMING_RHO"],
            kv["ASD_TIMING_P"]), transform=ax.transAxes, ha="right",
            va="top", fontsize=7.8, color=C_DARK, weight="bold")
    panel_title(ax, "Effect size vs dynamic range", fontsize=8.7)


# --------------------------------------------------------------------------- D
def draw_network(ax):
    """Host-gene network, DESCRIPTIVE ONLY (content correction retained)."""
    ax.axis("off")
    ax.set_xlim(-20, 120)
    ax.set_ylim(-24, 112)
    panel_title(ax, "Host-gene interaction network", fontsize=9.0)
    n = len(GENES)
    ang = {g: np.pi / 2 + 2 * np.pi * i / n for i, g in enumerate(GENES)}
    pos = {g: (50 + 40 * np.cos(ang[g]), 50 + 40 * np.sin(ang[g]))
           for g in GENES}
    deg = {g: 0 for g in GENES}
    for e in EDGES:
        deg[e["gene1"]] += 1
        deg[e["gene2"]] += 1
    for e in EDGES:
        x0, y0 = pos[e["gene1"]]
        x1, y1 = pos[e["gene2"]]
        if e["edge_type"] == "physical":
            ax.plot([x0, x1], [y0, y1], "-", color=C_DARK, lw=1.1, alpha=0.75,
                    zorder=1)
        else:
            ax.plot([x0, x1], [y0, y1], "--", color=C_MID, lw=0.8, alpha=0.6,
                    zorder=1)
    for g in GENES:
        x, y = pos[g]
        r = 2.6 + 0.75 * deg[g]
        ax.add_patch(plt.Circle((x, y), r,
                                fc=TIER_COLORS[TIER_ORDER[gene_tier[g]]],
                                ec="white", lw=0.9, zorder=2))
        D_NODES.append(ax.patches[-1])
        ax.text(x, y, str(deg[g]), ha="center", va="center", fontsize=7.0,
                color="white", weight="bold", zorder=3)
        # angle-aware anchoring keeps labels clear of the ring and each other
        ca, sa = np.cos(ang[g]), np.sin(ang[g])
        ha = "left" if ca > 0.35 else "right" if ca < -0.35 else "center"
        va = "bottom" if sa > 0.75 else "top" if sa < -0.75 else "center"
        # Final: centre-anchored (bottom/top) labels pushed to a larger
        # horizontal radius so adjacent bottom labels stay apart in the
        # narrower D panel.
        rad = 58.0 if ha == "center" else 50.0
        lx = 50 + rad * ca + (1.5 if ha == "left" else -1.5 if ha == "right"
                              else 0)
        ly = 50 + 50 * sa + (1.5 if va == "bottom" else -1.5 if va == "top"
                            else 0)
        D_GENE_LABELS.append(ax.text(lx, ly, g, ha=ha, va=va, fontsize=7.6,
                                     color=C_DARK, weight="bold"))
    # dedicated legend row below the network (structure counts -> legend)
    ax.plot([], [], "-", color=C_DARK, lw=1.1, label="physical interaction")
    ax.plot([], [], "--", color=C_MID, lw=0.8, label="functional association")
    ax.plot([], [], "o", ms=6, color=C_NEG, mec="white", label="size = degree")
    for t in TIER_ORDER:
        ax.plot([], [], "o", ms=4.5, color=TIER_COLORS[t], mec="white",
                label="Tier %s" % t)
    ax.legend(loc="upper center", fontsize=7.2, bbox_to_anchor=(0.5, -0.10),
              ncol=4, columnspacing=1.0, handletextpad=0.5)


# --------------------------------------------------------------------------- E
PW_LABEL2 = {"Synaptic signaling": "Synaptic\nsignaling",
             "Neuron projection": "Neuron\nprojection",
             "Chromatin regulation": "Chromatin\nregulation",
             "Cytoskeleton": "Cytoskeleton",
             "Cell adhesion": "Cell\nadhesion",
             "Axon guidance": "Axon\nguidance",
             "Protein localization": "Protein\nlocalization",
             "Vesicle trafficking": "Vesicle\ntrafficking",
             "Calcium signaling": "Calcium\nsignaling",
             "Ubiquitin/proteasome": "Ubiquitin /\nproteasome"}


def draw_pathway_heatmap(ax):
    """Gene x pathway MEMBERSHIP heatmap, DESCRIPTIVE ONLY."""
    genes = sorted(GENES, key=lambda g: (gene_tier[g], g))
    ng, npp = len(genes), len(PATHS)
    for i, g in enumerate(genes):
        for j, p in enumerate(PATHS):
            member = g in p["genes"]
            fc = TIER_COLORS[TIER_ORDER[gene_tier[g]]] if member else "#f0f0f0"
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc=fc,
                                       ec="white", lw=0.7, zorder=2))
    ax.set_xlim(-0.5, npp - 0.5)
    ax.set_ylim(ng - 0.5, -0.5)
    ax.set_xticks(range(npp))
    ax.set_xticklabels([PW_LABEL2[p["label"]].replace("\n", " ") for p in PATHS],
                       rotation=45, ha="right", rotation_mode="anchor",
                       fontsize=7.2)
    ax.tick_params(axis="x", pad=3)
    ax.set_yticks(range(ng))
    ax.set_yticklabels(genes, fontsize=7.4)
    ax.set_xticks(np.arange(npp) - 0.5, minor=True)
    ax.set_yticks(np.arange(ng) - 0.5, minor=True)
    ax.grid(which="minor", color="white", lw=0.4)
    ax.tick_params(which="minor", length=0)
    handles = [plt.Rectangle((0, 0), 1, 1, fc="#bfbfbf", ec="white",
                             label="member")]
    for t in TIER_ORDER:
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=TIER_COLORS[t],
                                     ec="white", label="Tier %s gene" % t))
    E_LEGEND[0] = ax.legend(handles=handles, loc="upper center",
                            fontsize=7.2, ncol=5, bbox_to_anchor=(0.5, -0.34),
                            columnspacing=1.0, handletextpad=0.4)
    panel_title(ax, "Curated functional memberships", fontsize=9.0)


# ---------------------------------------------------------------------------
# Final layout: taller canvas; wider A-B gutter (top row wspace 0.45 ->
# 0.62); bottom row D:E = 0.40:0.60 with a much wider gutter (0.12 -> 0.30)
# so E's gene-name labels sit fully inside the D-E gutter.
fig = plt.figure(figsize=(FIG_W, 252 * MM))
gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.0, 1.0],
                       hspace=0.5, left=0.02, right=0.98, top=0.97,
                       bottom=0.18)
gs_top = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0],
                                          width_ratios=[1.08, 0.94, 1.08],
                                          wspace=0.78)
gsA = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_top[0],
                                       width_ratios=[0.94, 0.06], wspace=0.05)
axA = fig.add_subplot(gsA[0])
axAs = fig.add_subplot(gsA[1])
draw_heatmap(axA, axAs)
panel_letter(axA, "A")
axB = fig.add_subplot(gs_top[1])
draw_slope(axB)
panel_letter(axB, "B")
axC = fig.add_subplot(gs_top[2])
draw_negscatter(axC)
panel_letter(axC, "C")
gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1],
                                          width_ratios=[0.40, 0.60],
                                          wspace=0.30)
axD = fig.add_subplot(gs_bot[0])
draw_network(axD)
panel_letter(axD, "D")
axE = fig.add_subplot(gs_bot[1])
draw_pathway_heatmap(axE)
panel_letter(axE, "E")

# ----------------------------------------------------------- QC (computed)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()

# 3D: gene labels vs gene labels and vs node circles
fig3d_violations = 0
lb = [(t, t.get_window_extent(renderer=renderer)) for t in D_GENE_LABELS]
for i in range(len(lb)):
    for j in range(i + 1, len(lb)):
        if overlap_area(lb[i][1], lb[j][1]) > 0:
            fig3d_violations += 1
    for c in D_NODES:
        if overlap_area(lb[i][1], c.get_window_extent(renderer=renderer)) > 0:
            fig3d_violations += 1

# 3E: axis labels (y gene names pairwise; x pathway labels pairwise)
yl = [t.get_window_extent(renderer=renderer) for t in axE.get_yticklabels()]
xl = [t.get_window_extent(renderer=renderer) for t in axE.get_xticklabels()]
fig3e_axis = 0
for group in (yl, xl):
    for i in range(len(group)):
        for j in range(i + 1, len(group)):
            if overlap_area(group[i], group[j]) > 2.0:
                fig3e_axis += 1

# 3E: legend must not cover any label / title
fig3e_leg = 0
legbb = E_LEGEND[0].get_window_extent(renderer=renderer)
for t in list(axE.get_yticklabels()) + list(axE.get_xticklabels()):
    if overlap_area(legbb, t.get_window_extent(renderer=renderer)) > 0:
        fig3e_leg += 1
for tt in (axE.title, axE._left_title, axE._right_title):
    if tt.get_visible() and tt.get_text().strip() and \
            overlap_area(legbb, tt.get_window_extent(renderer=renderer)) > 0:
        fig3e_leg += 1


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
    return out


a_bb = axA.get_window_extent(renderer=renderer)
as_bb = axAs.get_window_extent(renderer=renderer)
b_bb = axB.get_window_extent(renderer=renderer)
d_bb = axD.get_window_extent(renderer=renderer)
e_bb = axE.get_window_extent(renderer=renderer)

# A/B boundary: no B text inside A (or strip) bbox; no A text inside B bbox.
fig3ab_viol = 0
for bb in _ax_texts(axB):
    if overlap_area(bb, a_bb) > 0 or overlap_area(bb, as_bb) > 0:
        fig3ab_viol += 1
for bb in _ax_texts(axA):
    if overlap_area(bb, b_bb) > 0:
        fig3ab_viol += 1

# D/E boundary: no E text inside D bbox; no D text inside E bbox.
fig3de_viol = 0
for bb in _ax_texts(axE):
    if overlap_area(bb, d_bb) > 0:
        fig3de_viol += 1
for bb in _ax_texts(axD):
    if overlap_area(bb, e_bb) > 0:
        fig3de_viol += 1

# E left margin: every E gene-name label must stay >= 2 mm right of D's axes
PX_PER_MM = fig.dpi / 25.4
e_left_gap = min(bb.x0 for bb in yl) - d_bb.x1
FIG3E_LEFT_MARGIN_OK = int(e_left_gap >= 2.0 * PX_PER_MM)

checks = [("FIG3D_LABEL_OVERLAP", fig3d_violations),
         ("FIG3E_AXIS_LABEL_OVERLAP", fig3e_axis),
         ("FIG3E_LEGEND_OVERLAP", fig3e_leg),
         ("FIG3A_B_BOUNDARY_COLLISION_ZERO", fig3ab_viol),
         ("FIG3D_E_BOUNDARY_COLLISION_ZERO", fig3de_viol),
         ("FIG3E_LEFT_MARGIN_SUFFICIENT", FIG3E_LEFT_MARGIN_OK),
         ("FIG3E_LEFT_MARGIN_GAP_PX", round(e_left_gap, 1))]
_mf3 = min_text_size_pt(fig)
tsv_write(os.path.join(QC_DIR, "figure3_layout_checks.tsv"),
          ["check", "value"], [[g, v] for g, v in checks] +
          [["MIN_TEXT_SIZE_PT", _mf3],
           ["FIG_WIDTH_MM", round(fig.get_size_inches()[0] * 25.4, 2)]])
print("FIG3 checks:", checks)
print("FIG3 min text size (pt):", _mf3)

# ----------------------------------------------------------- render
NAME = "Figure3_FINAL_v2"
save_final(fig, FIGDIR, NAME)

# ----------------------------------------------------------- source data TSV
# Rows identical to Earlier Figure3_FINAL_SOURCE_DATA.tsv.
rows = []
for r in sorted(MT, key=lambda r: (devrank(r), r["gene"], r["HsaEX_ID"])):
    rows.append(["Figure_3", "A", "developmental PSI",
                 "%s (%s)" % (r["gene"], r["HsaEX_ID"]),
                 "prenatal=%.1f postnatal=%.1f trajectory=%s" %
                 (float(r["prenatal_mean_PSI"]),
                  float(r["postnatal_mean_PSI"]),
                  r["developmental_trajectory"]),
                 MASTER_TSV,
                 "prenatal_mean_PSI/postnatal_mean_PSI/"
                 "developmental_trajectory"])
_dyn = sorted((r for r in MT
               if r["developmental_dynamic_status"] == "DYNAMIC"),
              key=lambda r: -float(r["prenatal_mean_PSI"]))
assert len(_dyn) == 10
for r in _dyn:
    rows.append(["Figure_3", "B", "dynamic events slopegraph",
                 "%s (%s)" % (r["gene"], r["HsaEX_ID"]),
                 "prenatal=%.1f postnatal=%.1f trajectory=%s" %
                 (float(r["prenatal_mean_PSI"]),
                  float(r["postnatal_mean_PSI"]),
                  r["developmental_trajectory"]),
                 MASTER_TSV,
                 "developmental_dynamic_status=DYNAMIC subset"])
kv13 = {r["key"]: r["value"] for r in
        rd(os.path.join(DIR13, "13_qc", "key_statistics.tsv"))}
rows.append(["Figure_3", "C", "timing correlation", "Spearman rho / P",
             "rho=%s P=%s (not significant)" % (kv13["ASD_TIMING_RHO"],
                                                kv13["ASD_TIMING_P"]),
             os.path.join(DIR13, "13_qc", "key_statistics.tsv"),
             "ASD_TIMING_RHO / ASD_TIMING_P"])
net_metrics = rd(NET_METRICS)[0]
assert net_metrics["n_genes"] == "15" and net_metrics["n_edges"] == "25"
rows.append(["Figure_3", "D", "network structure", "nodes / edges / density",
             "%s / %s / %.3f" % (net_metrics["n_genes"],
                                 net_metrics["n_edges"],
                                 float(net_metrics["density"])),
             NET_METRICS, "n_genes / n_edges / density (legend text)"])
for e in EDGES:
    rows.append(["Figure_3", "D", "network edge",
                 "%s-%s" % (e["gene1"], e["gene2"]), e["edge_type"], NET_TSV,
                 "edge_type"])
for p in PATHS:
    rows.append(["Figure_3", "E", "pathway membership", p["label"],
                 "overlap=%d/15 genes=%s" % (p["n"], ",".join(p["genes"])),
                 PW_TSV, p["key"]])
tsv_write(os.path.join(FIGDIR, NAME + "_SOURCE_DATA.tsv"),
          ["figure", "panel", "series", "item", "value", "source_file",
           "source_locator"], rows)

prov = [
    ["Figure_3", "A", "Developmental PSI heatmap (BrainSpan)",
     "prenatal/postnatal mean PSI per event + trajectory strip",
     MASTER_TSV, "none (descriptive means)",
     "fonts raised for final width (Earlier)"],
    ["Figure_3", "B", "Dynamic events slopegraph",
     "10 dynamic events (9 PLPH, 1 PHPL) prenatal->postnatal PSI",
     MASTER_TSV, "none (descriptive)",
     "Final: wider A-B gutter; B y-axis/labels clear of A"],
    ["Figure_3", "C", "Effect size vs developmental dynamic range",
     "per-event |dPSI| vs PSI range; Spearman rho/P only",
     "09_ASD_timing_correlation + 13_qc/key_statistics.tsv",
     "Spearman correlation (non-significant negative control)",
     "italic interpretation prose removed (Earlier); Final: "
     "deterministic in-range ticks"],
    ["Figure_3", "D", "Host-gene interaction network",
     "15 nodes / 25 edges / density 0.238, tier-coloured (counts in legend)",
     NET_TSV + " ; " + NET_METRICS + " ; " + NET_RECLASS,
     "NONE: network permutation P not shown (direction-inference status "
     "reclassified in the evidence reconciliation)",
     "Final: D-E gutter widened; D:E width 0.40:0.60"],
    ["Figure_3", "E", "Curated functional memberships",
     "gene membership in 10 curated pathway sets (observed overlaps)",
     PW_TSV + " ; " + PW_RECLASS,
     "NONE: matched-permutation P not shown (pathway memberships are "
     "descriptive-only)",
     "Final: wider panel + left margin; gene labels never intrude into D"],
]
tsv_write(os.path.join(FIGDIR, NAME + "_PROVENANCE.tsv"),
          ["figure", "panel", "panel_title", "content", "data_source_files",
           "statistics_displayed", "corrections_vs_earlier"], prov)
print("Figure3_FINAL_v2 done")
