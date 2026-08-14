# environment/

Captured software environment used to run the analysis code in this
repository.

- `R_sessionInfo.txt` — R session information at packaging time.
- `python_environment.txt` — Python version and pinned package versions
  (pip freeze of the analysis virtual environment) plus the core packages
  used by the analysis scripts.
- `software_versions.tsv` — per-tool version table for R, Python, external
  executables (liftOver, bedtools) and the key analysis libraries.

Version honesty policy: only versions that were actually captured are
reported. Where the exact original version could not be recovered the value
is written as `not recorded`; no versions are invented.
