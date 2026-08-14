#!/usr/bin/env python3
"""Final journal-style supplementary figure legends.

Rewrites all 16 legends in publication style (spec sections 3, 6, 8):
no internal metadata blocks (Elements:/Statistical test:/...), no internal
language, essential
statistics and moved panel facts integrated in prose, 80-180 words each.
Writes 04_final_legends/SUPPLEMENTARY_FIGURE_LEGENDS_FINAL.md.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supplementary_common import (LEG_DIR, D36_ADJ, D36_NM, D36_PROT, DIR34,
                             rd, load_master, layer_flags, scan_text,
                             METADATA_MARKERS)

LEGEND_MD = os.path.join(LEG_DIR, "SUPPLEMENTARY_FIGURE_LEGENDS_FINAL.md")


def legends():
    d = {}

    # Final: S1 is now a single-panel figure (round-trip plot removed;
    # the liftOver conclusion is carried here in the legend).
    d[1] = (
        "Data-source availability and coordinate cross-checks for the "
        "19-event analysis set.",
        "Availability per event for the seven evidence layers: discovery "
        "ASD cortex, CHyMErA, BrainSpan developmental, GSE30573, "
        "PsychENCODE, reciprocal liftOver and GENCODE v33 local structure; "
        "'+' denotes supporting or available, '×' discordant and "
        "'–' absent or not applicable. In the reciprocal liftOver "
        "round-trip check (hg38–hg19–hg38) all 19 events returned "
        "identical round-trip coordinates: maximum absolute start and end "
        "differences were 0 bp. Of 36 designed events, 20 were "
        "cortex-matched human events and 19 entered the analysis set; "
        "57/57 liftOver chain steps completed and all 19 events show "
        "coordinate-equivalent local structures in GENCODE v33. Full human "
        "and mouse identifiers with GRCh38 coordinates are provided in the "
        "Supplementary Table.")

    d[2] = (
        "Background universes, permutation nulls and matching diagnostics.",
        "(A) Background universes in unique-event or annotation-catalog "
        "units: wide splicing-event universe (200,956 events), conserved "
        "microexon (377), CEM-derived (149), NN-derived (102) and "
        "PSI-matched sampling pool (4,948). (B) Permutation null 95% "
        "intervals with medians for each analysis background against the "
        "observed 19-event mean |ΔPSI| = 0.0399 (seed = 42, 10,000 "
        "permutations), annotated with the permutation P and the observed "
        "percentile. (C) Matched-background balance: standardized mean "
        "difference of exon length, CEM −0.336 and NN −0.073 (shaded "
        "region |SMD| < 0.5). (D) Analysis backgrounds in post-filter "
        "analysis-record units: wide SE 20,916; microexon 467; conserved "
        "microexon 452; CEM 308 pair records (180 unique); NN 380 pair "
        "records (120 unique). Panels A and D use different statistical "
        "units and are not directly comparable; the microexon-only "
        "background exists only at the analysis level. Enrichment remained "
        "significant across all five stringency levels.")

    d[3] = (
        "Definition, selection and prior-exclusion sensitivity analyses.",
        "(A) Sensitivity P values on the −log10 scale with the P = 0.05 "
        "reference: microexon length thresholds ≤27 nt (P = 0.0040), "
        "≤30 nt (P = 0.0042), ≤36 nt (P = 0.0042), one-event-per-gene "
        "deduplication (P = 0.0024), target versus random same-size "
        "19-event sets drawn from the conserved-microexon background "
        "(empirical one-sided P = 0.0023; 23/10,001 random sets reached or "
        "exceeded the observed statistic) and ASD-prior gene exclusion "
        "(P = 0.0335). (B) ASD-prior exclusion applied at the 19-event "
        "level: 12 events retained and 7 excluded (permutation "
        "P = 0.0335; Wilcoxon P = 0.0266); the excluded events map to five "
        "genes (ANK3, CAMTA1, CTNND1, MEF2A, PTK2) of the nine-gene "
        "ASD-prior list; NRG1, RIMS2, UNC13A and UNC13B carry no event in "
        "the analysis set. Mapping level, "
        "background stringency and deduplication were stable and no "
        "single-gene driver was present (15/15); sensitivity values are "
        "reported unadjusted.")

    d[4] = (
        "Developmental timing analyses for all 19 events.",
        "(A) Slopegraph of prenatal-to-postnatal mean PSI for all 19 "
        "events, colored by trajectory class (PLPH 9, PHPL 1, non-dynamic "
        "9). (B) Dynamicity tests against conserved microexon "
        "(P = 0.0004), CEM (P = 0.0005) and NN (P = 0.0013) backgrounds, "
        "together with negative and exploratory controls: PSI-matched "
        "P = 0.7186, gene-block P = 1.0, a BrainSpan expression-level "
        "substrate test (P = 0.7199; not a splicing test) and an "
        "exploratory zebrafish Mann–Whitney test (P = 0.0688, suggestive "
        "only). (C) Per-event developmental monotonicity (Spearman rho), "
        "colored by dynamic status. (D) Per-event developmental PSI range "
        "(prenatal circle to postnatal square). PLPH denotes "
        "prenatal-low/postnatal-high and PHPL prenatal-high/postnatal-low "
        "trajectories; dynamicity tests are reported unadjusted as "
        "developmental context.")

    d[5] = (
        "RBP motif permutation landscape across all 240 tests.",
        "(A) Ranked −log10 permutation P for all 240 motif-enrichment "
        "tests (4 comparisons × 6 regions × 10 motifs across 26 curated "
        "RBP genes); 4/240 tests were nominal (P < 0.05) and none survived "
        "BH-FDR. (B) Maximum −log10 permutation P per comparison and "
        "region. (C) The four nominal tests grouped by RBP/motif: MBNL1/2 "
        "(dynamic versus non-dynamic, down-extended, P = 0.0454), PTBP1/2 "
        "(19-event set versus background, up-proximal P = 0.0356 and "
        "up-extended P = 0.0441) and ELAVL2/3/4 (tier2-5 versus background, "
        "exon, P = 0.0261); none remained significant after BH-FDR, and "
        "SRRM4 was not nominally significant (minimum permutation "
        "P = 0.0928). The dashed vertical line marks the nominal P = 0.05 "
        "reference; the shortened panel labels denote: down ext., dynamic "
        "versus non-dynamic down-extended; up prox. and up ext., 19-event "
        "set versus background up-proximal and up-extended; exon, tier2-5 "
        "versus background exon; minimum, minimum over all 240 tests. "
        "Motif enrichment is exploratory context; all results "
        "are BH-FDR non-significant.")

    d[6] = (
        "GSE30573 cross-cohort mapping and local-structure matching.",
        "(A) Mapping flow: 3 of 19 events were mapped and analyzable "
        "(ANK3, FBXO25, HERC4); 16 were unmapped (12 without a "
        "local-structure match and 4 whose gene was absent from the GSE "
        "annotation). (B) Local event-structure matching for the three "
        "mapped events; the 16 unmapped events are omitted. Among mapped "
        "events, 2 of 3 were directionally concordant with discovery "
        "(FBXO25, HERC4) and 1 discordant (ANK3; exact-binomial P = 0.5). "
        "GSE30573 provides limited cross-cohort directional context only "
        "and is not an independent validation dataset.")

    d[7] = (
        "Transcript membership and representative transcript structures.",
        "(A) GENCODE v33 inclusion- and exclusion-transcript counts "
        "supporting each isoform for all 19 events. (B) Representative "
        "transcript structures for CPEB4 HsaEX0016976 "
        "(chr5:173943025-173943049), CTNND1 HsaEX0017872 "
        "(chr11:57806460-57806478) and MEF2D HsaEX0038710 "
        "(chr1:156476493-156476514) with mouse orthologous events; "
        "microexons are shown in red, and all 19 events show "
        "coordinate-equivalent local structures in GENCODE v33. "
        "PsychENCODE effects are transcript-usage logit betas quantified "
        "as RSEM effective-length-normalized ratios of inclusion versus "
        "exclusion transcript sets, i.e. a transcript-usage measure rather "
        "than junction PSI.")

    d[8] = (
        "PsychENCODE diagnostics and set-level robustness.",
        "(A) Per-event PsychENCODE logit-beta standard errors, colored by "
        "evidence tier. (B) Kenward–Roger versus likelihood-ratio-test "
        "P-value agreement with the identity line. (C) Per-event direction "
        "concordance between discovery and PsychENCODE; 15/19 events were "
        "concordant (exact-binomial P = 0.0096). (D) Model-based set "
        "validation: KR BH-FDR < 0.05, 4/19; KR BH-FDR < 0.10, 7/19; LRT "
        "BH-FDR < 0.05, 6/19 and < 0.10, 10/19 (sensitivity); "
        "leave-one-event-out minimum concordant 14/19; leave-one-gene-out "
        "minimum concordant 13/19. Analyses used 532 cortical samples from "
        "80 donors (38 ASD / 42 control); BH-FDR was applied across the 19 "
        "event-level Kenward–Roger tests.")

    d[9] = (
        "Technical-covariate and regional stability of PsychENCODE "
        "effects.",
        "(A) Per-event effects under the primary model M0 versus the "
        "technical-covariate model M4 for all 19 events (Pearson "
        "r = 0.935; direction agreement 16/19; Tier A events in color). "
        "(B) Tier A effects under models M0–M4; all four directions are "
        "stable under every model and remain KR BH-FDR < 0.05 under M4. "
        "(C) Direction-stability heatmap (19 events × M0–M4); blue denotes "
        "the M0 direction and red a flip; M4 retains 16/19, and the three "
        "flips involve non-Tier-A events. (D) Leave-one-region-out effects "
        "for the four Tier A events across 11 region-exclusion folds; "
        "dotted lines mark the full-data effects; direction was stable in "
        "4/4 events in all folds (no reversal; maximum |beta "
        "shift| = 0.04805). Coefficients are ASD diagnosis coefficients "
        "for logit-transformed transcript usage (532 cortical samples from "
        "80 donors).")

    d[10] = (
        "Transcript-set definition sensitivity (D0–D3).",
        "(A) Event analyzability under definitions D0 (reference sets), "
        "D1 (strict-local), D2 (deterministic representative pair) and D3 "
        "(protein-coding); green denotes analyzable and grey "
        "non-analyzable event\u2013definition pairs (analyzable counts 19, "
        "15, 19 and 18). (B) "
        "Per-event effects across definitions, colored by evidence tier; "
        "directions are broadly stable while magnitudes vary by "
        "definition. (C) Tier A effects with 95% confidence intervals "
        "under each analyzable definition; direction remains down in ASD "
        "for all four events under every analyzable definition, while KR "
        "BH-FDR < 0.05 counts are 4, 2, 0 and 2 for D0–D3 (significance is "
        "definition-sensitive; PTPRF is not analyzable under D1). (D) "
        "Direction-concordant events (14, 11, 15, 13) versus Tier A "
        "significant events (4, 2, 0, 2) per definition. Coefficients are "
        "ASD diagnosis coefficients for logit-transformed transcript "
        "usage.")

    d[11] = (
        "Leave-one-donor-out influence analysis for Tier A events.",
        "(A) Tier A full-model effects with all 80 single-donor-deletion "
        "estimates per event under the primary model M0. (B) Standardized "
        "donor-deletion deviation heatmap (80 donors × 4 events; "
        "(beta_d − full beta)/within-event SD). (C) Donor-deletion beta "
        "distributions after deleting ASD donors (152 estimates: 38 "
        "donors × 4 events) versus control donors (168: 42 × 4), with "
        "medians. (D) Minimum absolute effect retention across the 80 "
        "deletions for M0 and M4 with the 0.50 reference threshold: "
        "per-event direction preservation was 1.0 (no reversal), minimum "
        "retention 0.755–0.885 and maximum |DFBETA| ≤ 0.859. All eight "
        "Tier A model rows met the prespecified criteria (direction "
        "preservation ≥ 0.95, minimum retention ≥ 0.50, no "
        "|DFBETA| ≥ 1.0); no single donor drove the Tier A effects. "
        "Analyses used 532 cortical samples from 80 donors (38 ASD / 42 "
        "control).")

    d[12] = (
        "Cell-composition estimates and composition-adjusted (M4C) "
        "sensitivity.",
        "(A) Estimated broad cell-fraction distributions across seven "
        "classes (532 cortical samples, 80 donors). (B) Marker-score "
        "validation matrix (7 × 7 Spearman rho; '+' where the class meets "
        "the prespecified criterion, matching rho ≥ 0.3 and greater than "
        "all off-target rhos); six of seven classes met the criterion; the "
        "inhibitory-neuron estimate did not exceed its excitatory-neuron "
        "off-target correlation. (C) Composition-PC variance explained "
        "(PC1–PC6: 0.730, 0.124, 0.087, 0.033, 0.017, 0.008; k = 2 "
        "retained, cumulative 85.4%) with PC1/PC2 loadings. (D) M4 versus "
        "composition-adjusted M4C coefficients (Pearson r = 0.9679, "
        "Spearman rho = 0.9702, direction retained 18/19, median effect "
        "retention 0.662). (E) Tier A effects across transcript-set "
        "definitions under M4C with 95% CI (PTPRF not analyzable under "
        "D1). (F) Tier A M4C donor-deletion stability (80 deletions per "
        "event; all four rows met the prespecified criteria).")

    d[13] = (
        "Tier A model-adjusted transcript-usage differences on the "
        "probability scale.",
        "(A–C) Model-adjusted transcript-usage differences (ASD − "
        "control) for the four Tier A events on the probability scale "
        "under (A) model M0 (clinical covariates), (B) model M4 (M0 plus "
        "technical covariates) and (C) model M4C (M4 plus "
        "cell-composition principal components). Point estimates are "
        "marginal standardized fixed-effect predictions over the "
        "532-sample covariate distribution; 95% confidence intervals were "
        "obtained by parametric simulation of the fixed-effect "
        "variance–covariance matrix (1,000 draws, seed = 42). All 12 "
        "model-by-event estimates are negative, matching the direction of "
        "the logit-scale diagnosis coefficients: ASD is associated with "
        "lower adjusted usage of the inclusion transcript set under every "
        "model. Event-level significance is reported with each model "
        "elsewhere; this figure displays effect magnitudes and intervals.")

    d[14] = (
        "Neuron-merged cell-composition sensitivity.",
        "(A) Per-event ASD diagnosis coefficients under the neuron-merged "
        "composition model versus model M4 for all 19 events (dashed "
        "line = identity; Pearson r = 0.9856). The neuron-merged model "
        "merges excitatory and inhibitory neuron fractions before PCA "
        "(six broad classes; k = 3 composition PCs, 91.5% cumulative "
        "variance) and otherwise matches M4C. (B) Tier A coefficients with "
        "95% CI under M4, seven-class M4C and the neuron-merged model, "
        "with BH-FDR annotated per estimate. All four Tier A directions "
        "are preserved; CLASP1, PTK2 and PTPRF remain BH-FDR < 0.05 "
        "(0.0335, 0.0148 and 0.0026) while HERC4 is sub-threshold (0.0586), "
        "mirroring the seven-class pattern. The seven-class M4C results "
        "are the primary analysis; the neuron-merged model is shown as a "
        "sensitivity analysis.")

    d[15] = (
        "Tier A microexon insertion sites in protein context.",
        "(A–D) AlphaFold v6 pLDDT profiles (grey line) with UniProt "
        "features near each Tier A insertion site (red band): (A) CLASP1 "
        "(Q7Z460), residues 673–682 within a disordered MAPRE1/MAPRE3 "
        "interaction region (662–785); (B) HERC4 (Q5GLZ8), residues "
        "643–650; (C) PTK2 (Q05397), insertion at the junction before "
        "E393, four residues from the regulatory autophosphorylation site "
        "Y397; (D) PTPRF (P10586), residues 772–780 within fibronectin "
        "type-III repeat 5 (711–819). All four microexons are in-frame "
        "(27, 24, 18 and 27 nt; 9, 8, 6 and 9 added residues). Colored "
        "boxes and ticks mark UniProt features within 30 residues of the "
        "site; the complete feature inventory is provided in Supplementary "
        "Table S8. pLDDT is a computational confidence measure shown as "
        "structural context only and does not establish functional or "
        "binding effects.")

    # Final: single-panel figure (letter removed); tier colors stated
    # explicitly rather than 'same as main figures'.
    d[16] = (
        "Evidence-layer combinations across the 19-event set.",
        "Bars show event counts for each observed combination of the "
        "seven evidence layers, stacked by evidence tier (14 combinations "
        "among 19 events). The layers are CHyMErA direction concordance, "
        "developmental dynamicity, host-gene network membership, GSE30573 "
        "mappability, PsychENCODE direction concordance, Kenward–Roger "
        "BH-FDR < 0.05 and LRT BH-FDR < 0.05. Each column represents one "
        "observed combination; the matrix rows mark layer membership, and "
        "columns are ordered by event count and then by number of layers. "
        "Tier colors match the main figures: Tier A blue, Tier B light "
        "blue, Tier C orange and Tier D grey.")
    return d


def check():
    d = legends()
    for n in range(1, 17):
        title, body = d[n]
        words = len((title + " " + body).split())
        assert 80 <= words <= 185, (n, words)
        hit = scan_text(title + " " + body)
        assert not hit, (n, hit)
        for m in METADATA_MARKERS:
            assert m not in title + body, (n, m)
        # Final: S1 and S16 are single-panel figures without a panel
        # letter, so their legends carry no "(A)" reference.
        assert "(A)" in body or n in (1, 16), n
    # spot-assert key transcribed values against the source tables
    rows = rd(os.path.join(D36_ADJ,
                           "TIER_A_ADJUSTED_TRANSCRIPT_USAGE_M0_M4_M4C.tsv"))
    assert len(rows) == 12 and all(
        float(r["adjusted_ASD_minus_control"]) < 0 for r in rows)
    vs = rd(os.path.join(D36_NM, "M4C_NEURON_MERGED_VS_PRIMARY_M4C.tsv"))
    assert len(vs) == 19
    MT = load_master()
    assert len(MT) == 19
    assert len({tuple(layer_flags(r)) for r in MT}) == 14
    print("LEGEND_CHECKS_OK (word counts, terminology, metadata)")
    return d


def write(d):
    L = ["# Supplementary Figure Legends — final submission version",
         "",
         "Journal-style legends for the 16 supplementary figures of the",
         "Molecular Autism submission. Each legend states what the figure",
         "shows, describes the panels, and carries only essential",
         "methodological and statistical context.",
         "",
         "Shared conventions: 19-event / 15-gene microexon analysis set",
         "(GRCh38/hg38); colorblind-safe semantic palette; '+' supporting,",
         "'×' discordant, '–' absent / not significant / NA; discovery",
         "enrichment permutations use seed = 42 with 10,000 permutations;",
         "PsychENCODE effects are ASD diagnosis coefficients for",
         "logit-transformed transcript usage.",
         ""]
    for n in range(1, 17):
        title, body = d[n]
        L += ["## Figure S%d" % n, "", "**Figure S%d. %s** %s" % (n, title,
                                                                 body), ""]
    with open(LEGEND_MD, "w") as f:
        f.write("\n".join(L))
    print("wrote", os.path.basename(LEGEND_MD))
    return {n: (len(d[n][1].split()) + len(d[n][0].split()))
            for n in d}


if __name__ == "__main__":
    d = check()
    wc = write(d)
    for n in range(1, 17):
        print("Figure S%d: %d words" % (n, wc[n]))
