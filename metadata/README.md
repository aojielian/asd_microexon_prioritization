# metadata/

Reader-facing metadata tables for the ASD microexon analysis.

- `event_dictionary.tsv` — the 19 microexon events: cross-species
  identifiers (HsaEX_ID, MmuEX_ID), gene, GRCh38 coordinates, main evidence
  tier, CHyMErA perturbation direction, developmental dynamicity and
  trajectory, GSE30573 mapping, and PsychENCODE concordance/significance.
- `analysis_dictionary.tsv` — the 12 script families: purpose, inputs,
  scripts, outputs, and manuscript destination (figure/table) of each
  analysis family.
- `source_dataset_manifest.tsv` — public/external data sources used in the
  study with accession/DOI identifiers, their role, and redistribution
  status. Source data are not redistributed; derived analysis data are
  provided in the companion Zenodo data package.

All tables are tab-separated with a header row; lines beginning with `#`
are comments. Values use reader-facing vocabulary (no internal pipeline
enums).
