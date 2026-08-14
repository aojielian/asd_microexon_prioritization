# Data availability

## Derived analysis data (Zenodo)

All derived analysis data needed to verify and reuse the paper are released
as a companion Zenodo data package:

- **Package**: `asd_microexon_prioritization_analysis_data`
- **DOI**: to be updated at release (placeholder)

The Zenodo package contains the derived data underlying every figure and
table of the paper (event master data, enrichment and sensitivity tables,
developmental timing, PsychENCODE model results, transcript-definition
sensitivity, LODO robustness, cell-composition estimates, composition-
adjusted (M4C) results, probability-scale effects, protein-context tables,
reproducibility checks and the source-to-output crosswalk). See the
package README and `DATA_DICTIONARY.md` for the full file list and column
definitions.

This repository (`asd_microexon_prioritization`, GitHub) contains the analysis code;
the Zenodo package contains the derived data. Both are versioned together
at release.

## Public/restricted source data

Source data are **not redistributed** in either release; they remain
available from their original repositories (see
`metadata/source_dataset_manifest.tsv` for accessions and access notes):

- CHyMErA perturbation atlas (GSE291610, GSE291608)
- Parikshak et al. autism transcriptome (GSE64018)
- BrainSpan developmental transcriptome
- VastDB (microexon annotation, orthology, PSI tables)
- GENCODE v33 (GRCh38)
- GSE30573 fetal brain transcriptome
- PsychENCODE processed resource (restricted; apply through the portal)
- brainSCOPE adult human brain cell-type reference
- UniProt; AlphaFold protein structure database

## Citation of the data package

See `CITATION.cff` in this repository (and the `CITATION` section of the
Zenodo record) for the recommended citation.
