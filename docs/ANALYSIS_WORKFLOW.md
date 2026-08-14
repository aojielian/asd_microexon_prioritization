# Analysis workflow

End-to-end order of the ASD microexon analysis. Each numbered step maps to
one `scripts/NN_*` directory. Scripts resolve all inputs through the
environment variables documented in `config/paths_template.yaml` and never
require the author's private directory layout.

## Step 1 — Event definition and cross-species mapping (`01_candidate_mapping`)

- Filter CHyMErA perturbation-defined microexon events (mouse).
- Map mouse events to human orthologous microexons using VastDB orthology.
- Convert coordinates to GRCh38 with reciprocal liftOver verification.
- Check GENCODE v33 local structure at each mapped event.
- Construct the final set of 19 events and the master event table
  (also released as Table S1 and `metadata/event_dictionary.tsv`).

Scripts: `core_mapping.py`, `coordinate_matching.py`, `coordinate_inference.R`,
`transcript_membership.R`, `master_event_table.py`.

## Step 2 — Autism set-level enrichment (`02_set_level_enrichment`)

- Extract Parikshak (GSE64018) delta-PSI values.
- Build the wide / microexon / conserved / CEM / NN backgrounds.
- Run 10,000-permutation enrichment with bootstrap confidence intervals.
- Apply BH and Bonferroni multiplicity correction; run one-event-per-gene
  and length-threshold sensitivities; run random same-size set comparison
  and matched-background sensitivity.

Scripts: `backgrounds_and_enrichment.py`, `event_reconciliation.py`,
`multiplicity_correction.py`, `sensitivity_and_negative_controls.py`.

## Step 3 — CHyMErA directional concordance (`03_directional_concordance`)

- Apply the |delta-PSI| > 0.01 testability rule to CHyMErA perturbation
  results (testability is based on effect magnitude, not P-value).
- Compute primary directional concordance (12/14 eligible events) and the
  all-event sensitivity (16/19).
- Perform one-to-one perturbation mapping with guide-sharing and multi-exon
  mapping checks.

Scripts: `chymera_directional_analysis.py`, `chymera_direction_verification.py`,
`directional_eligibility_rule.py`.

## Step 4 — BrainSpan developmental timing (`04_brainspan_development`)

- Extract developmental PSI from BrainSpan.
- Classify each event as prenatal-low/postnatal-high (PLPH),
  prenatal-high/postnatal-low (PHPL), or non-dynamic.
- Test dynamicity with background sensitivities; compute monotonicity and
  PSI range; correlate developmental range with the ASD effect.

Scripts: `developmental_timing_analysis.py`, `developmental_reports.py`.

## Step 5 — Network and RBP context (`05_network_and_rbp_context`)

- Host-gene network context (STRING/GO/Reactome-based curated network).
- RBP motif analysis at microexon-flanking regions (240-test corrected
  battery).
- GSE30573 limited cross-cohort directional context (small-n; reported as
  context only, not as independent validation).

Scripts: `mechanistic_context_analysis.py`, `context_reports.py`,
`gse30573_limited_context.py`.

## Step 6 — PsychENCODE event-level models (`06_psychencode_models`)

- De-duplicate donors with exact subject-ID overlap filtering
  (final analysis set: 80 donors, 38 ASD / 42 controls; 532 cortical
  samples).
- Fit mixed models M0 (primary) through M4 (technical covariates) for
  transcript usage at each event; Kenward-Roger inference with LRT
  sensitivity and BH correction.
- One-event-per-gene concordance, diagnosis-by-region omnibus, and Tier A
  diagnosis-by-sex / diagnosis-by-age sensitivity.

Scripts: `donor_deduplication.R`, `mixed_model_analysis.R`,
`model_reproduction.R`, `region_interaction_omnibus.R`,
`tierA_moderator_sensitivity.R`.

## Step 7 — Transcript-set definition sensitivity (`07_transcript_definition_sensitivity`)

- Construct transcript-set definitions D0-D3.
- Select deterministic D2 representative transcript pairs (highest median
  support; documented tie-breaking).
- Re-run event-level models under each definition; verify transcript support.

Scripts: `d2_transcript_selection.R`, `psychencode_sensitivity.R`.

## Step 8 — Donor robustness: LODO (`08_lodo`)

- Leave-one-donor-out fits for M0 and M4 (and M4C) across all 19 events.
- Tier A influence metrics: direction retention, minimum retention,
  maximum |DFBETA|.

Scripts: `lodo.R`, `m4c_lodo.R`.

## Step 9 — Cell composition (`09_cell_composition`)

- Build gene-level TPM from PsychENCODE gene TPMs.
- Deconvolve 532 cortical samples with non-negative least squares (NNLS)
  against a local brainSCOPE adult human brain reference, aggregated from
  24 subtypes to 7 broad classes (log2-scale pseudo-bulk linearized; NNLS
  used because the reference is not raw counts).
- Compute composition principal components; run the neuron-merged
  sensitivity (6-class reference, k = 3 PCs).

Scripts: `gene_tpm.R`, `marker_pc.R`, `deconvolution_nnls.py`,
`neuron_merged_sensitivity.R`.

## Step 10 — Tier A probability-scale effects (`10_probability_scale`)

- Model-adjusted transcript-usage effects for the four Tier A events on the
  probability scale (marginal standardized fixed-effect predictions with
  parametric bootstrap confidence intervals).

Scripts: `adjusted_usage_probability_scale.R`.

## Step 11 — Protein context (`11_protein_context`)

- Determine in-frame status of Tier A microexon insertions on canonical
  transcripts.
- Map insertion sites to UniProt features and AlphaFold pLDDT context.

Scripts: `protein_source_retrieval.py`, `protein_context_mapping.py`.

## Step 12 — Submission figures and tables (`12_submission_figures_tables`)

- Render main figures 1-5 and supplementary figures S1-S16 (PDF, SVG,
  TIFF, PNG).
- Compose supplementary figure pages and the combined PDF.
- Build the supplementary tables workbook (S1-S12) with submission-clean,
  reader-facing vocabulary.
- Build the publication-exact Main Table 1 (Primary Tier A microexon
  events and model stability) and Main Table 2 (Tier A robustness,
  probability-scale effects and protein context) source files
  (`final_main_tables_1_2.py`, 4 rows each).

Scripts: `figcommon.py`, `figcommon_main.py`, `figure1_main.py` ..
`figure5_main.py`, `supplementary_figures_s1_s8.py`,
`supplementary_figures_s9_s16.py`, `supplementary_legends.py`,
`supplementary_page_compose.py`, `prepare_supplementary_tables.py`,
`final_main_tables_1_2.py`, `supplementary_common.py`,
`finalize_table_s1.py`.
