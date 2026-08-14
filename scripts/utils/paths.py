"""Shared path helpers for the ASD microexon analysis code.

All public scripts must resolve file locations from the neutral environment
variables documented in ``config/paths_template.yaml`` and must never require
the author's private directory layout.

Environment variables:

- ``PROJECT_ROOT``  -- pipeline workspace root (the directory containing the
  numbered analysis directories, e.g. ``11_set_level_enrichment``,
  ``34_robustness_and_composition``, ``41_submission_figures_and_tables``).
- ``DATA_ROOT``     -- optional raw-data root (used when raw inputs live
  outside the workspace).
- ``OUTPUT_ROOT``   -- optional override for generated outputs.
- ``REFERENCE_ROOT``-- reference-resource root (GENCODE, VastDB, liftOver
  chains, brainSCOPE, ...).
- ``SCRATCH_ROOT``  -- scratch directory for large intermediates.
- ``LIFTOVER_PATH`` -- path to the UCSC liftOver executable.

All helpers return ``os.PathLike`` objects; missing variables fall back to the
current working directory.
"""
import os


def project_root() -> str:
    """Pipeline workspace root (see config/paths_template.yaml)."""
    return os.environ.get("PROJECT_ROOT", ".")


def data_root() -> str:
    """Optional raw-data root (see config/paths_template.yaml)."""
    return os.environ.get("DATA_ROOT", project_root())


def output_root() -> str:
    """Optional output override (see config/paths_template.yaml)."""
    return os.environ.get("OUTPUT_ROOT", project_root())


def reference_root() -> str:
    """Reference-resource root (see config/paths_template.yaml)."""
    return os.environ.get("REFERENCE_ROOT", data_root())


def scratch_root() -> str:
    """Scratch directory for large intermediates."""
    return os.environ.get("SCRATCH_ROOT", os.path.join(project_root(), "scratch"))


def liftover_path() -> str:
    """Path to the UCSC liftOver executable (empty if not configured)."""
    return os.environ.get("LIFTOVER_PATH", "")
