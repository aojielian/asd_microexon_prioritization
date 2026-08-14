# File crosswalk

Map between manuscript deliverables (figures, tables, claims) and the
analysis outputs / scripts that produce them.

## Main figures

| Deliverable | Producing analysis family | Key outputs | Zenodo package location |
| --- | --- | --- | --- |
| Figure 1 (cross-species construction of the perturbation-defined microexon set) | 01_candidate_mapping | Master event table; mapping QC | 01_event_master/ and 12_figure_source_data/ |
| Figure 2 (autism cortical enrichment and CHyMErA directional comparison) | 02_set_level_enrichment | Enrichment and multiplicity tables | 02_set_level_enrichment/, 03_chymera_direction/, and 12_figure_source_data/ |
| Figure 3 (developmental timing and host-gene biological context) | 04_brainspan_development | Trajectory classification and statistics | 04_brainspan_development/ and 12_figure_source_data/ |
| Figure 4 (cross-cohort PsychENCODE support for the perturbation-defined event set) | 06_psychencode_models | Model results and concordance tables | 05_psychencode_models/ and 12_figure_source_data/ |
| Figure 5 (integrated prioritization and local protein context of Tier A microexons) | Multiple (see panel provenance in the legend facts) | Summary tables | 09_composition_adjustment/, 10_probability_scale/, 11_protein_context/, and 12_figure_source_data/ |

## Supplementary figures

| Deliverable | Producing analysis family |
| --- | --- |
| S1 (data-source availability and coordinate cross-checks for the 19-event analysis set) | 01_candidate_mapping |
| S2 (background universes, permutation nulls and matching diagnostics) | 02_set_level_enrichment |
| S3 (definition, selection and prior-exclusion sensitivity analyses) | 01_candidate_mapping |
| S4 (developmental timing analyses for all 19 events) | 04_brainspan_development |
| S5 (RBP motif permutation landscape across all 240 tests) | 05_network_and_rbp_context, 03_directional_concordance |
| S6 (GSE30573 cross-cohort mapping and local-structure matching) | 05_network_and_rbp_context |
| S7 (transcript membership and representative transcript structures) | 05_network_and_rbp_context |
| S8 (PsychENCODE diagnostics and set-level robustness) | 06_psychencode_models |
| S9 (technical-covariate and regional stability of PsychENCODE effects) | 06_psychencode_models, 08_lodo |
| S10 (transcript-set definition sensitivity, D0-D3) | 07_transcript_definition_sensitivity |
| S11 (leave-one-donor-out influence analysis for Tier A events) | 08_lodo |
| S12 (cell-composition estimates and composition-adjusted (M4C) sensitivity) | 09_cell_composition |
| S13 (Tier A model-adjusted transcript-usage differences on the probability scale) | 10_probability_scale |
| S14 (neuron-merged cell-composition sensitivity) | 09_cell_composition |
| S15 (Tier A microexon insertion sites in protein context) | 11_protein_context |
| S16 (evidence-layer combinations across the 19-event set) | 02_set_level_enrichment, 06_psychencode_models |

## Supplementary tables

| Deliverable | Title (abridged) | Producing analysis family |
| --- | --- | --- |
| Table S1 | Microexon event set and cross-species mapping: complete evidence for the 19 events | 01_candidate_mapping; 12_submission_figures_tables (finalize_table_s1.py) |
| Table S2 | PsychENCODE diagnosis coefficients under M0-M4 | 06_psychencode_models |
| Table S3 | PsychENCODE diagnosis coefficients under transcript-set definitions D0-D3 | 07_transcript_definition_sensitivity |
| Table S4 | Leave-one-donor-out (LODO) influence summary (M0, M4) | 08_lodo |
| Table S5 | Cell-composition reference coverage, aggregate fraction summaries, marker validation, and composition-PC loadings/variance. | 09_cell_composition |
| Table S6 | Composition-adjusted (M4C) event-level results | 08_lodo, 09_cell_composition |
| Table S7 | Tier A probability-scale effects (M0/M4/M4C) | 10_probability_scale |
| Table S8 | Tier A coding consequences and protein annotation | 11_protein_context |
| Table S9 | Neuron-merged composition sensitivity | 09_cell_composition |
| Table S10 | Deterministic D2 representative transcript pairs | 07_transcript_definition_sensitivity |
| Table S11 | Reported ancestry categories and model encoding | 06_psychencode_models |
| Table S12 | Additional sensitivity analyses of directional concordance, multiplicity, regional and demographic heterogeneity, donor de-duplication, transcript support, and perturbation-to-human event mapping. | 02/03/06 families |

## Main tables

| Deliverable | Title | Producing analysis family |
| --- | --- | --- |
| Table 1 | Primary Tier A microexon events and model stability (4 rows) | 12_submission_figures_tables (final_main_tables_1_2.py), from 01_candidate_mapping master table, 10_probability_scale adjusted usage, 08_lodo LODO summary |
| Table 2 | Tier A robustness, probability-scale effects and protein context (4 rows) | 12_submission_figures_tables (final_main_tables_1_2.py), from 08_lodo LODO summary, 10_probability_scale adjusted usage, 09_cell_composition neuron-merged summary, 11_protein_context |

## Script-to-output map

Full per-script input/output contracts are listed in
`metadata/analysis_dictionary.tsv` (12 families, one row each) and in each
script's docstring. File-level provenance and the mapping from publication
items to released scripts and source data are provided in
PUBLICATION_CROSSWALK.tsv and ZENODO_TO_SCRIPT_MAP.tsv.

## Notes

- The final publication figures/tables were assembled from these outputs
  by family `12_submission_figures_tables` with submission-clean,
  reader-facing vocabulary; the numeric content is identical to the
  upstream derived tables (verified by numeric-integrity checks at
  packaging time).
- The companion Zenodo package mirrors the file-level locations in the
  first table's last column; see its `DATA_DICTIONARY.md` for column
  definitions.
