# Reproducibility

This document explains how to run the analysis code and what a faithful
reproduction requires.

## Environment

- R >= 4.5.1 with `lme4` 1.1.38, `pbkrtest` 0.5.5 and `MASS` 7.3.65
  (see `environment/R_sessionInfo.txt` and `environment/software_versions.tsv`).
- Python 3.13.2 with the pinned package versions in
  `environment/python_environment.txt` (matplotlib, numpy, pandas, scipy,
  pyliftover, openpyxl, Pillow, reportlab).
- External executables used in upstream steps: UCSC `liftOver` (coordinate
  conversion) and `bedtools` (event definition). Exact versions are not
  recorded.

## Inputs

All scripts read inputs through the neutral environment variables documented
in `config/paths_template.yaml`:

- `PROJECT_ROOT` — pipeline workspace root containing the numbered analysis
  directories.
- `REFERENCE_ROOT` — reference resources (GENCODE v33 GTF, VastDB tables,
  liftOver chains, brainSCOPE reference profiles, ...).
- `DATA_ROOT` / `OUTPUT_ROOT` / `SCRATCH_ROOT` / `LIFTOVER_PATH` — optional
  overrides (see the template).

The upstream public/restricted data sources are listed in
`metadata/source_dataset_manifest.tsv` with accessions. Source data are not
redistributed; scripts expect the input formats produced by the preceding
steps of the pipeline (for example, junction PSI tables from CHyMErA,
gene-level TPM from PsychENCODE).

## Recommended order

Run the steps in the order documented in `docs/ANALYSIS_WORKFLOW.md`
(1 candidate mapping -> 2 set-level enrichment -> 3 directional
concordance -> 4 BrainSpan -> 5 context -> 6 PsychENCODE models -> 7
transcript definitions -> 8 LODO -> 9 composition -> 10 probability scale
-> 11 protein context -> 12 figures and tables). Each step reads the
outputs of the previous steps.

## What a faithful reproduction produces

- The 19-event master table with cross-species identifiers, GRCh38
  coordinates and per-event evidence (`metadata/event_dictionary.tsv` is a
  reader-facing projection; Table S1 is the full evidence table).
- Set-level enrichment values (10,000 permutations), multiplicity-corrected
  results and sensitivity tables.
- CHyMErA directional concordance (12/14 eligible; 16/19 all-event).
- BrainSpan developmental classification (9 PLPH / 1 PHPL / 9 non-dynamic;
  10 dynamic).
- PsychENCODE model results (Kenward-Roger P, BH FDR) with donor overlap
  filtering (80 donors; 532 cortical samples).
- Transcript-definition sensitivity (D0-D3), LODO robustness, cell
  composition fractions and M4C composition-adjusted results.
- Tier A probability-scale effects and protein-context tables.
- Figures 1-5, Supplementary Figures S1-S16 and Supplementary Tables
  S1-S12.

## Notes

- Exact numeric reproduction of stochastic analyses (10,000-permutation
  enrichment, bootstrap CIs) requires identical random seeds; seed values
  are fixed inside the corresponding scripts.
- The final publication values were verified in a separate verification
  step (upstream numeric-integrity checks); the released derived data in the
  companion Zenodo package are the authoritative values.
- Where full reproduction depends on restricted data (PsychENCODE),
  apply for access through the original resource and provide inputs in the
  documented formats.
