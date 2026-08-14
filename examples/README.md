# examples/

Run-through examples for the analysis code.

This directory is reserved for small worked examples (how to run a script
against small inputs, expected output shapes, and how to interpret a result
file). The full pipeline order and per-script input/output contracts are
documented in `docs/ANALYSIS_WORKFLOW.md`, `docs/REPRODUCIBILITY.md` and
`docs/FILE_CROSSWALK.md`.

Quick start (see also `config/README.md`):

```bash
export PROJECT_ROOT=/path/to/your/analysis/workspace
export REFERENCE_ROOT=/path/to/reference/resources

# Python scripts (activate a Python 3.13 environment with the packages
# listed in environment/python_environment.txt)
python scripts/01_candidate_mapping/core_mapping.py

# R scripts (R >= 4.5 with lme4, pbkrtest, MASS)
Rscript scripts/06_psychencode_models/donor_deduplication.R
```

Because full reproduction depends on external/restricted data (CHyMErA,
Parikshak, BrainSpan, VastDB, GENCODE v33, PsychENCODE; see
`metadata/source_dataset_manifest.tsv`), scripts read their inputs through
the neutral environment variables documented in `config/paths_template.yaml`
and expect the input formats produced by the upstream steps of the pipeline.
