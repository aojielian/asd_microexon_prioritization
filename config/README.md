# config/

Runtime configuration for the analysis code.

- `paths_template.yaml` — documents the environment variables used by all
  scripts (`PROJECT_ROOT`, `DATA_ROOT`, `OUTPUT_ROOT`, `REFERENCE_ROOT`,
  `SCRATCH_ROOT`, `LIFTOVER_PATH`) and the expected directory layout.

Scripts never require the author's private directory structure; they resolve
every file from these neutral environment variables. To run the pipeline on
your own machine:

```bash
export PROJECT_ROOT=/path/to/your/analysis/workspace
export REFERENCE_ROOT=/path/to/reference/resources
export DATA_ROOT=/path/to/raw/data            # optional
export LIFTOVER_PATH=/path/to/liftOver        # optional; needed for coordinate conversion
```

Do not commit machine-specific paths; keep a copy of this template outside
the repository.
