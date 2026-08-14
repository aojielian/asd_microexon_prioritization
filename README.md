# ASD-associated microexon events: perturbation-guided prioritization analysis code

## Purpose

This repository contains the analysis code for perturbation-guided
prioritization of autism-associated microexon splicing events. Starting
from CHyMErA perturbation-defined microexon events in mouse, the pipeline
maps events to the human cortex, tests set-level enrichment for autism,
establishes developmental and directional context, models event-level
diagnosis effects in the independent PsychENCODE resource, adds
robustness, cell-composition and transcript-definition sensitivity
analyses, and finally prioritizes events through protein-context evidence.

## Study overview

```
CHyMErA perturbation-defined events
-> human cortical event mapping
-> autism set-level enrichment
-> developmental and directional context
-> independent PsychENCODE mixed-model analysis
-> robustness / composition / transcript-definition sensitivity
-> protein-context prioritization
```

## Repository organization

```
asd_microexon_prioritization/
├── README.md                  This file
├── LICENSE                    License
├── CITATION.cff               Citation metadata
├── environment/               Captured R/Python/tool versions
├── config/                    paths_template.yaml + documentation
├── scripts/                   12 numbered analysis families + utils/
│   ├── 01_candidate_mapping/          CHyMErA event filtering, mouse-human mapping, coordinates, 19-event master table
│   ├── 02_set_level_enrichment/       backgrounds, 10,000-permutation enrichment, multiplicity correction
│   ├── 03_directional_concordance/    CHyMErA |delta-PSI| rule, 12/14 and 16/19 directional concordance
│   ├── 04_brainspan_development/      developmental PSI, PLPH/PHPL/non-dynamic classification
│   ├── 05_network_and_rbp_context/    host-gene network, RBP motif battery, GSE30573 limited context
│   ├── 06_psychencode_models/         donor de-duplication, M0-M4 mixed models, Kenward-Roger inference
│   ├── 07_transcript_definition_sensitivity/  D0-D3 definitions, deterministic D2 pairs
│   ├── 08_lodo/                       leave-one-donor-out (M0/M4/M4C), Tier A influence metrics
│   ├── 09_cell_composition/           NNLS deconvolution of 532 cortical samples, composition PCs
│   ├── 10_probability_scale/          Tier A probability-scale effects
│   ├── 11_protein_context/            in-frame status, UniProt features, AlphaFold pLDDT
│   ├── 12_submission_figures_tables/  figures 1-5, S1-S16, tables S1-S12 assembly
│   └── utils/                         shared path helpers
├── metadata/                  event_dictionary.tsv, analysis_dictionary.tsv, source_dataset_manifest.tsv
├── examples/                  quick-start notes
└── docs/                      ANALYSIS_WORKFLOW.md, REPRODUCIBILITY.md, DATA_AVAILABILITY.md, FILE_CROSSWALK.md
```

## Data sources

Public/external sources (see `metadata/source_dataset_manifest.tsv` for
accessions and access notes; source data are not redistributed):

- CHyMErA / GSE291610 / GSE291608 — perturbation-defined microexon events
- Parikshak / GSE64018 — autism transcriptome delta-PSI values
- BrainSpan — developmental transcriptome
- VastDB — microexon annotation, orthology, PSI tables
- GENCODE v33 — reference annotation (GRCh38)
- GSE30573 — limited cross-cohort directional context only (not an
  independent validation dataset)
- PsychENCODE processed resource — event-level mixed-model analysis
- adult human brain reference (brainSCOPE) — cell-composition analysis
- UniProt; AlphaFold — protein-context annotation

## Reproduction

Scripts resolve all paths from the neutral environment variables documented
in `config/paths_template.yaml` (`PROJECT_ROOT`, `DATA_ROOT`, `OUTPUT_ROOT`,
`REFERENCE_ROOT`, `SCRATCH_ROOT`, `LIFTOVER_PATH`) and never require the
author's private directory structure.

Run the 12 families in order (details in `docs/ANALYSIS_WORKFLOW.md`):

1. `scripts/01_candidate_mapping/` — event definition and mapping
2. `scripts/02_set_level_enrichment/` — enrichment and multiplicity
3. `scripts/03_directional_concordance/` — CHyMErA directional rule
4. `scripts/04_brainspan_development/` — developmental timing
5. `scripts/05_network_and_rbp_context/` — context analyses
6. `scripts/06_psychencode_models/` — PsychENCODE models
7. `scripts/07_transcript_definition_sensitivity/` — transcript definitions
8. `scripts/08_lodo/` — leave-one-donor-out robustness
9. `scripts/09_cell_composition/` — composition deconvolution
10. `scripts/10_probability_scale/` — probability-scale effects
11. `scripts/11_protein_context/` — protein context
12. `scripts/12_submission_figures_tables/` — figures and tables

Full reproduction depends on the external/restricted data listed above;
each script documents the expected input format. See
`docs/REPRODUCIBILITY.md` and `environment/` for the software environment.

## Final outputs

- **Figures 1-5** — produced by `scripts/12_submission_figures_tables/`
  (`figure1_main.py` .. `figure5_main.py`).
- **Supplementary Figures S1-S16** — produced by
  `supplementary_figures_s1_s8.py`, `supplementary_figures_s9_s16.py`,
  `supplementary_page_compose.py` (with `supplementary_legends.py`).
- **Main Table 1 (Primary Tier A microexon events and model stability)
  and Main Table 2 (Tier A robustness, probability-scale effects and
  protein context)** — produced by
  `scripts/12_submission_figures_tables/final_main_tables_1_2.py`
  (publication-exact 4-row source data; also released in the Zenodo
  package under `13_tables/`).
- **Supplementary Tables S1-S12** — produced by
  `prepare_supplementary_tables.py` / `finalize_table_s1.py`
  (with `supplementary_common.py`).

## Data availability

Derived analysis data underlying every figure and table are released as a
companion Zenodo data package:

- **Package**: `asd_microexon_prioritization_analysis_data`
- **DOI**: `10.5281/zenodo.21928557`

See `docs/DATA_AVAILABILITY.md` for the package contents and access notes
for the original source data.

## Citation

Please cite the manuscript and this code repository using the metadata in
`CITATION.cff`.
