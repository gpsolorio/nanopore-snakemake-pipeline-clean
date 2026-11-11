##############################################################################
# common.smk — Shared globals, paths, helpers, and shell utilities
# Included by nanopore_adaptive_pipeline.smk before all rule includes.
##############################################################################

import os
import csv
import shutil
import sys
from pathlib import Path
from glob import glob

import pandas as pd

base = Path(workflow.basedir).resolve()


# ============================================================================
# Directory Layout
# ============================================================================

WORKFLOW_ROOT = (base.parent).resolve()
REPO_ROOT     = (base.parent.parent).resolve()

ENVS_DIR       = (WORKFLOW_ROOT / "envs").resolve()
SCRIPTS_DIR    = (WORKFLOW_ROOT / "scripts").resolve()
CONTAINERS_DIR = (WORKFLOW_ROOT / "containers").resolve()
TOOLS_DIR      = (WORKFLOW_ROOT / "tools").resolve()

sys.path.insert(0, str(SCRIPTS_DIR))


# ============================================================================
# Default Paths (overridable via config.paths)
# ============================================================================

_DEFAULT_RESULTS   = (REPO_ROOT / "results").resolve()
_DEFAULT_LOGS      = (REPO_ROOT / "logs").resolve()
_DEFAULT_SLURM     = (REPO_ROOT / "slurm_logs").resolve()
_DEFAULT_RESOURCES = (REPO_ROOT / "resources").resolve()


def as_path(value, fallback):
    if not value:
        return Path(fallback).resolve()
    return Path(value).expanduser().resolve()


def _resolve_repo_path(value):
    """
    Resolve a config-supplied path that may be either absolute (e.g. a
    reference genome living on a shared mount) or relative
    (e.g. a BED file checked into the repo under resources/). Absolute
    values are returned as-is; relative values are resolved against
    REPO_ROOT. Returns None for falsy input so callers can pass optional
    config keys straight through without a separate guard.
    """
    if not value:
        return None
    p = Path(value).expanduser()
    return (p if p.is_absolute() else (REPO_ROOT / p)).resolve()


# ============================================================================
# Config Unpacking  (must come before BARCODES resolution)
# ============================================================================

_paths   = config.get("paths", {})
_inputs  = _paths.get("inputs", {})
_outputs = _paths.get("outputs", {})
_tpaths  = _paths.get("tools", {})

EXPERIMENT_CFG = config.get("experiment", {})
STUDY = EXPERIMENT_CFG.get("study", "default")
TOOLS_CFG      = config.get("tools", {})
RESOURCE_CFG   = config.get("resources", {})
ANALYSIS_CFG   = config.get("analysis", {})

RUN_DORADO_CFG = TOOLS_CFG.get("dorado", TOOLS_CFG.get("run_dorado", {}))
CLAIR3_CFG     = TOOLS_CFG.get("clair3",  TOOLS_CFG.get("run_clair3", {}))


# ============================================================================
# Genome Mode  (drives BED / reference path selection throughout)
# ============================================================================

GENOME_MODE = config.get("genome", EXPERIMENT_CFG.get("genome", "hg38")).lower()

_ANNOTATIONS_CFG = config.get("annotations", {})

def _get_annotation(key):
    """
    Return the BED path for `key` matching the current genome mode.
    Raises clearly if the path is missing from config or doesn't exist on disk.
    """
    mode_cfg = _ANNOTATIONS_CFG.get(GENOME_MODE)
    if not mode_cfg or key not in mode_cfg:
        raise ValueError(
            f"annotations.{GENOME_MODE}.{key} is not set in config.\n"
            f"Add it under annotations: in your config YAML.\n"
            f"Current genome mode: {GENOME_MODE}"
        )
    path = _resolve_repo_path(mode_cfg[key])
    if not path.exists():
        raise FileNotFoundError(
            f"Annotation BED for '{key}' (genome={GENOME_MODE}) not found: {path}\n"
            + (
                "Run scripts/make_mini_beds.py to generate mini BEDs from hg38 BEDs."
                if GENOME_MODE == "mini" else
                "Check annotations path in config."
            )
        )
    return path.resolve()


# ============================================================================
# Prebuilt Alignment Mode  (independent of GENOME_MODE)
#
# GENOME_MODE selects which reference/annotation FILES to use (mini vs
# hg38) -- it says nothing about whether reads still need to be aligned.
# Some "mini" runs are static test fixtures with alignments already baked
# in (e.g. the AR-locus GIAB smoke tests) and should skip
# basecall_pod5/merge_bams/sort_index_bam entirely, pulling BAMs straight
# from a curated fixture directory. Other "mini" runs (e.g. a smaller,
# faster production reference) still need the real alignment pipeline to
# run over real POD5/BAM inputs. These are orthogonal, so they're two
# separate config flags -- conflating them into a single GENOME_MODE check
# silently sends real samples down the fixture-only path and produces
# confusing MissingInputException errors for any sample that was never
# manually pre-aligned into the fixture directory.
# ============================================================================

USE_PREBUILT_ALIGNMENTS = bool(EXPERIMENT_CFG.get("use_prebuilt_alignments", False))


# ============================================================================
# DMR Grouping Config  (drives the metadata TSV header + DSS group selection)
#
# group_column / group_values are entirely config-driven so the same
# dss_call_dmrs rule / run_dss.R script work for any two-level split --
# diagnosis (e.g. ASD vs TD), sex, family, batch, etc. Nothing here is
# hardcoded to a particular cohort.
# ============================================================================

_DMR_CFG = config.get("dmr", {})
DMR_GROUP_COLUMN = _DMR_CFG.get("group_column", "diagnosis")
DMR_GROUP_VALUES = _DMR_CFG.get("group_values", ["ASD", "TD"])

# Optional: restrict DMR *calling* (DSS/dmrseq) to one role (e.g. "PRI"),
# independent of which samples get processed (BARCODES). Lets a family
# sheet's `roles:` include every relative -- so they all get basecalled /
# aligned / variant-called / piled-up -- while the actual ASD-vs-TD
# statistical comparison only ever sees probands, not their relatives
# (who inherit their family's diagnosis label, not their own). Empty/unset
# means no role restriction -- every sample matching group_column/
# group_values is compared, same as before this option existed.
DMR_RESTRICT_ROLE = _DMR_CFG.get("restrict_role", "")

# Explicit per-sample group labels from config, for datasets with no
# family-style diagnosis sheet (e.g. public GIAB samples grouped by sex).
# Keys/values are just strings supplied in experiment.sample_groups --
# nothing hardcoded here. Still supported as an override/escape hatch for
# datasets with no sheet at all; when a sample_sheet IS configured (see
# sample_ids below), prefer deriving labels from the sheet instead of
# hand-copying them here, so there's one source of truth.
_EXPLICIT_SAMPLE_GROUPS = EXPERIMENT_CFG.get("sample_groups", {})


# ============================================================================
# Sample / Barcode Resolution
# ============================================================================

# Every role column that can appear in a family-style sample sheet. Used
# as the default column set to scan when experiment.sample_sheet.sample_ids
# is given without an explicit `roles:` override -- see sample_ids handling
# in _parse_sample_sheet() below for why this needs to be the full set
# rather than the usual ["PRI"] default.
ALL_FAMILY_ROLES = ["PRI", "MTH", "FTH", "MGM", "MGF", "PGM", "PGF", "SIB1"]


def _norm_id(value) -> str:
    """
    Normalize a sheet cell to a plain ID string for comparison/storage.

    openpyxl/pandas reads whole-number Excel cells (Family ID, sample IDs
    typed as numbers rather than text) back as float64, so a family ID
    that's "1001" in the spreadsheet UI becomes the Python float
    1001.0, and a naive str(x) turns that into "1001.0" -- which then
    fails to match a config family_ids entry of "1001" (isin() compares
    strings) and would show up wrong in SAMPLE_FAMILY / the metadata TSV.
    This strips a trailing ".0" for any value that's an integer-valued
    float, while leaving genuine non-integer or already-string values
    (e.g. a family ID that really does contain letters) untouched.
    """
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _resolve_barcodes_legacy(exp_cfg):
    """Original range/list mode — used when no sample_sheet is configured."""
    bc_cfg = exp_cfg.get("barcodes", {})
    mode   = bc_cfg.get("mode", "list")
    if mode == "list":
        values = bc_cfg.get("values", [])
        if not values:
            raise ValueError("experiment.barcodes.mode is 'list' but no values provided")
        return values
    if mode == "range":
        start   = bc_cfg["start"]
        end     = bc_cfg["end"]
        prefix  = bc_cfg.get("prefix", "barcode")
        padding = bc_cfg.get("padding", 2)
        return [f"{prefix}{str(i).zfill(padding)}" for i in range(start, end + 1)]
    raise ValueError(f"experiment.barcodes.mode must be 'range' or 'list', got '{mode}'")


def _parse_sample_sheet(exp_cfg):
    """
    Parse a sample sheet (CSV or Excel) -> (barcodes list, sample_id->diagnosis
    dict, sample_id->role dict, sample_id->family_id dict).
    Pure data parsing only — no filesystem writes here.
    Falls back to legacy range/list mode if no sample_sheet is configured.

    Two modes:
      - id_col set:   simple CSV/Excel with one column of sample IDs + diagnosis
                       column. No role concept -- every sample maps to role="".
                       No family concept either -- family_id="" for all.
      - id_col unset: family Excel sheet with sample IDs spread across role
                       columns (PRI, MTH, FTH, MGM, MGF, PGM, PGF, SIB1, ...).
                       The role dict records which column each sample_id came
                       from -- e.g. so a DMR-calling step can be restricted to
                       role="PRI" (probands) even when `roles:` includes every
                       family member so they all still get processed through
                       basecalling/alignment/variant-calling/pileup. Without
                       this, "process everyone" and "compare everyone" were
                       the same knob -- there was no way to do one without
                       the other, and doing both for a family sheet means
                       comparing a family's overall diagnosis label (which
                       every relative in that row shares) as if it were each
                       individual's own diagnosis, which it is not.
                       Family ID (sheet column 0) is recorded per-sample
                       alongside role, so downstream rules/results can be
                       organized or filtered by family without re-parsing
                       the sheet themselves.

    Optional `sample_ids:` allowlist (family-sheet mode only):
      family_ids + roles are inherently uniform filters -- every included
      family contributes the SAME set of role columns. That falls apart
      for a curated sample set that's asymmetric per family (e.g. take
      PRI+MTH from one family, but only PGM+PGF from another, because
      that's simply which relatives were sequenced/relevant for that
      family) -- no combination of family_ids/roles can express that
      without also pulling in each family's PRI/MTH/etc. that you didn't
      want. `sample_ids` sidesteps this: it's a flat list of exact sample
      IDs to keep, independent of which role or family they came from.
      When set, `roles` defaults to ALL_FAMILY_ROLES (rather than the
      usual ["PRI"]) so every column is scanned and no matter which role
      an allowlisted ID happens to sit under, it's found -- the caller is
      choosing samples by ID, not by role, so the role scan itself
      shouldn't be the thing silently excluding a wanted ID. diagnosis,
      role, and family_id are still read off the sheet per-sample as
      normal; only the final membership list is overridden by the
      allowlist. This also means family_ids becomes unnecessary (though
      harmless) once sample_ids is set, since sample_ids is authoritative
      regardless of family.
    """
    ss_cfg = exp_cfg.get("sample_sheet", {})
    path   = ss_cfg.get("path")

    if not path:
        return _resolve_barcodes_legacy(exp_cfg), {}, {}, {}

    path = _resolve_repo_path(path)

    if not path.exists():
        raise FileNotFoundError(
            "\n\n"
            "============================================================\n"
            "  Sample sheet not found\n"
            "============================================================\n"
            f"  The pipeline is looking for this file:\n"
            f"    {path}\n"
            "\n"
            "  This file doesn't exist yet, so the pipeline doesn't know\n"
            "  which samples to run or what their diagnoses are.\n"
            "\n"
            "  How to fix this — pick whichever applies to you:\n"
            "\n"
            "  1. If you have the sample sheet, place (or copy) it at the\n"
            "     path shown above.\n"
            "\n"
            "  2. If your sample sheet lives somewhere else, update the\n"
            "     path in your config file:\n"
            "\n"
            f"       experiment:\n"
            f"         sample_sheet:\n"
            f"           path: \"path/to/your/sample_sheet.xlsx\"\n"
            "\n"
            "     (Edit the config file you passed with --configfile.)\n"
            "\n"
            "  3. If you're just running the quick mini test and don't\n"
            "     have real sample data, check that your configs/mini_data/*.yml config\n"
            "     points at the mini test's own sample sheet, not the\n"
            "     full study's sample sheet.\n"
            "============================================================\n"
        )

    sheet  = ss_cfg.get("sheet", 0)
    dx_col = ss_cfg.get("diagnosis_col", "Diagnosis")
    id_col = ss_cfg.get("id_col", None)
    sample_ids_filter = ss_cfg.get("sample_ids", None)
    roles  = ss_cfg.get("roles", ALL_FAMILY_ROLES if sample_ids_filter else ["PRI"])
    skip_empty_dx    = ss_cfg.get("skip_empty_diagnosis", False)
    filter_diagnoses = ss_cfg.get("diagnoses_to_run", None)
    family_ids       = ss_cfg.get("family_ids", None)

    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path, header=0)
    else:
        df = pd.read_excel(path, sheet_name=sheet, header=0)

    # Drop duplicate columns (family-sheet Excel has two "Family ID" columns)
    df = df.loc[:, ~df.columns.duplicated()]

    # For Excel family sheets, drop non-data rows by requiring numeric first column
    if path.suffix.lower() != ".csv":
        df = df[pd.to_numeric(df.iloc[:, 0], errors="coerce").notna()].copy()

    # Optional allowlist by Family ID (column 0) -- for a multi-family pilot
    # sheet where only some families have actually been sequenced, this lets
    # `roles:` safely include every relative (MTH, FTH, MGM, ...) without
    # pulling in every family in the sheet, most of which have no BAM data
    # yet and would otherwise blow up BARCODES with sample IDs
    # get_bam_paths() can never resolve. Family ID is assumed to be column 0
    # (true for the family-style sheet this was built for); id_col-mode
    # sheets without a meaningful column 0 shouldn't set this. Not needed
    # (though harmless) when sample_ids is set -- see sample_ids docstring
    # note above.
    if family_ids is not None:
        _wanted_fam = {_norm_id(x) for x in family_ids}
        _before = len(df)
        df = df[df.iloc[:, 0].apply(_norm_id).isin(_wanted_fam)]
        print(f"[sample_sheet] family_ids filter: kept {len(df)}/{_before} rows "
              f"(families: {sorted(_wanted_fam)})")

    df[dx_col] = df[dx_col].fillna("").astype(str).str.strip()

    if skip_empty_dx:
        df = df[df[dx_col] != ""]
    if filter_diagnoses:
        df = df[df[dx_col].isin(filter_diagnoses)]

    barcodes         = []
    sample_diagnosis = {}
    sample_role      = {}
    sample_family    = {}

    if id_col:
        # Simple mode: one column holds sample IDs, another holds diagnosis.
        # No role column in this mode -- every sample gets role="".
        # No family concept in this mode either -- family_id="" for all.
        for _, row in df.iterrows():
            sample_id = str(row.get(id_col, "")).strip()
            if not sample_id or sample_id.lower() in ("nan", "", "none"):
                continue
            dx = str(row.get(dx_col, "")).strip()
            if sample_id in sample_diagnosis:
                continue
            barcodes.append(sample_id)
            sample_diagnosis[sample_id] = dx
            sample_role[sample_id] = ""
            sample_family[sample_id] = ""
    else:
        # Family sheet mode: sample IDs spread across role columns (PRI, MTH, etc.)
        # Family ID is column 0 (same column family_ids filtering above uses).
        for _, row in df.iterrows():
            dx  = row.get(dx_col, "")
            fam = _norm_id(row.iloc[0])
            for role in roles:
                sample_id = str(row.get(role, "")).strip()
                if not sample_id or sample_id.lower() in ("nan", "", "none"):
                    continue
                if not sample_id.isdigit():
                    continue
                if sample_id in sample_diagnosis:
                    continue
                barcodes.append(sample_id)
                sample_diagnosis[sample_id] = dx
                sample_role[sample_id] = role
                sample_family[sample_id] = fam

    # sample_ids allowlist: narrow the parsed set down to exactly these IDs,
    # regardless of which role/family they came from. Order follows the
    # allowlist itself (not sheet row order) so the config's own ordering
    # is what shows up in BARCODES/logs/the metadata TSV. Any requested ID
    # not found under the scanned `roles` is a hard error rather than a
    # silent drop -- silently proceeding with fewer samples than the
    # config explicitly asked for is the "zero samples remain" class of surprise;
    # better to fail loudly with the missing ID named.
    if sample_ids_filter:
        _wanted_ids = [str(x).strip() for x in sample_ids_filter]
        _missing_ids = [sid for sid in _wanted_ids if sid not in sample_diagnosis]
        if _missing_ids:
            raise ValueError(
                f"experiment.sample_sheet.sample_ids requested ID(s) not found "
                f"under roles {roles} in {path}: {_missing_ids}\n"
                f"If the ID exists under a role not currently scanned, add it "
                f"to experiment.sample_sheet.roles explicitly."
            )
        barcodes = _wanted_ids
        sample_diagnosis = {sid: sample_diagnosis[sid] for sid in barcodes}
        sample_role      = {sid: sample_role[sid]      for sid in barcodes}
        sample_family    = {sid: sample_family[sid]    for sid in barcodes}
        print(f"[sample_sheet] sample_ids filter: kept exactly the "
              f"{len(barcodes)} requested ID(s)")

    if not barcodes:
        raise ValueError(
            f"No sample IDs found in sample sheet {path}. "
            f"Check that 'id_col' or 'roles' values match column headers exactly.\n"
            f"  id_col={id_col}, roles={roles}\n"
            f"  Available columns: {df.columns.tolist()}"
        )

    print(f"[sample_sheet] {len(barcodes)} samples loaded from {path}")
    for sid in barcodes:
        dx = sample_diagnosis.get(sid, "")
        role_suffix = f"  [{sample_role[sid]}]" if sample_role.get(sid) else ""
        fam_suffix  = f"  (family {sample_family[sid]})" if sample_family.get(sid) else ""
        print(f"    {sid}  ->  {dx or '(no diagnosis)'}{role_suffix}{fam_suffix}")

    return barcodes, sample_diagnosis, sample_role, sample_family


if "mini_data" in config:
    BARCODES = [str(x) for x in config["mini_data"]["samples"].keys()]
    _, _all_diagnosis, _all_role, _all_family = _parse_sample_sheet(EXPERIMENT_CFG)
    SAMPLE_DIAGNOSIS = {b: _all_diagnosis.get(b, "") for b in BARCODES}
    SAMPLE_ROLE      = {b: _all_role.get(b, "")      for b in BARCODES}
    SAMPLE_FAMILY    = {b: _all_family.get(b, "")    for b in BARCODES}

    _missing_dx = [b for b in BARCODES if not SAMPLE_DIAGNOSIS[b]]
    if _missing_dx:
        print(f"[mini_data] WARNING: no diagnosis found for: {_missing_dx}")

    print(f"[mini_data] {len(BARCODES)} samples loaded from config")
    for _b in BARCODES:
        print(f"    {_b}  ->  {SAMPLE_DIAGNOSIS[_b] or '(no diagnosis)'}")
else:
    BARCODES, SAMPLE_DIAGNOSIS, SAMPLE_ROLE, SAMPLE_FAMILY = _parse_sample_sheet(EXPERIMENT_CFG)

if _EXPLICIT_SAMPLE_GROUPS:
    SAMPLE_DIAGNOSIS = {
        b: _EXPLICIT_SAMPLE_GROUPS.get(b, SAMPLE_DIAGNOSIS.get(b, ""))
        for b in BARCODES
    }
    print(f"[sample_groups] {len(_EXPLICIT_SAMPLE_GROUPS)} explicit group label(s) applied "
          f"(column='{DMR_GROUP_COLUMN}')")


# ============================================================================
# Resolved Paths  (must come after config unpacking, before TSV write)
# ============================================================================

RESULTS_ROOT = as_path(_outputs.get("results"), _DEFAULT_RESULTS)

STUDY_ROOT = RESULTS_ROOT / STUDY
STUDY_ROOT.mkdir(parents=True, exist_ok=True)

LOGS_ROOT       = as_path(_outputs.get("logs"),       _DEFAULT_LOGS)
SLURM_LOGS_ROOT = as_path(_outputs.get("slurm_logs"), _DEFAULT_SLURM)

REFERENCE_GENOME  = _get_annotation("reference")
print(f"[genome] mode        = {GENOME_MODE}")
print(f"[genome] reference   = {REFERENCE_GENOME}")
print(f"[genome] use_prebuilt_alignments = {USE_PREBUILT_ALIGNMENTS}")


def mini_aligned_root():
    """
    Root directory containing per-barcode pre-aligned fixture BAMs
    (<root>/<barcode>_AR/minimap2.bam), used when USE_PREBUILT_ALIGNMENTS
    is true. Moved up here (from its old spot under "Path Helper
    Functions") because _has_indexed_bam_data(), just below, needs to
    check this location instead of BAM_INPUT_DIR whenever
    USE_PREBUILT_ALIGNMENTS is set -- fixture-based mini runs never
    populate BAM_INPUT_DIR at all, so checking only that path there
    excluded every sample and emptied BARCODES.
    """
    return as_path(
        EXPERIMENT_CFG.get("mini_aligned_dir"),
        RESULTS_ROOT / "mini-test" / "aligned",
    )

BAM_INPUT_DIR     = as_path(_inputs.get("bam")       or config.get("bam_input_dir"),    _DEFAULT_RESOURCES / "bam")
POD5_INPUT_DIR    = as_path(_inputs.get("pod5")       or config.get("pod5_input_dir"),   _DEFAULT_RESOURCES / "pod5")


# ============================================================================
# Drop sample-sheet entries with no BAM data indexed yet
#
# A multi-role family sample sheet (experiment.sample_sheet.roles including
# relatives) lists every sample ID in those columns across the WHOLE sheet,
# regardless of whether that specific person has actually been sequenced/
# delivered/indexed yet -- a planning spreadsheet is naturally ahead of the
# data. Without this check, BARCODES would include IDs with nothing under
# BAM_INPUT_DIR at all, and get_bam_paths() raising ValueError for the
# first such sample hard-fails DAG construction for the ENTIRE pipeline,
# not just that one sample.
#
# Only ever looks inside BAM_INPUT_DIR/<barcode>/ -- never walks the rest
# of BAM_INPUT_DIR -- so one absent barcode can't trigger a scan of every
# OTHER sample's data tree. Within the barcode dir it tries the flat
# <barcode>/*.bam layout first, then falls back to an early-exit recursive
# walk (follows symlinks), because on the HPC bam_index each <barcode> is a
# symlink to a run dir that nests its BAMs one level deeper, e.g.
#   bam_index/sample_006 -> .../barcode07_sample_006_5mC-5hmC/  (BAMs live in here)
# Samples with no BAMs at either depth are treated as not yet available,
# so scripts/index_ont_bam_dirs.sh (with --verify) stays the one place new samples
# get added.
# ============================================================================

def _iter_barcode_bams(barcode_dir):
    """Yield every *.bam under barcode_dir at any depth, following symlinks.

    os.walk(followlinks=True) is used instead of Path.glob("**/*.bam")
    because pathlib's "**" does not descend into symlinked subdirectories
    on Python < 3.13. Lazy, so callers that only need existence can stop
    at the first hit.
    """
    for root, _dirs, names in os.walk(str(barcode_dir), followlinks=True):
        for name in names:
            if name.endswith(".bam"):
                yield os.path.join(root, name)


def _has_indexed_bam_data(barcode):
    if USE_PREBUILT_ALIGNMENTS:
        fixture_bam = mini_aligned_root() / f"{barcode}_AR" / "minimap2.bam"
        return fixture_bam.exists()

    # BAM_INPUT_DIR may itself be a single BAM file (see get_bam_paths()).
    if Path(BAM_INPUT_DIR).is_file() and Path(BAM_INPUT_DIR).suffix == ".bam":
        return True

    direct = Path(BAM_INPUT_DIR) / barcode
    if not direct.is_dir():
        return False

    # Fast path: BAMs directly under bam_index/<barcode>/
    if any(direct.glob("*.bam")):
        return True

    # BAMs nested below the barcode dir (HPC symlink-to-run-dir layout).
    return next(_iter_barcode_bams(direct), None) is not None


# The directory _has_indexed_bam_data() actually checked, purely for the
# warning message below -- BAM_INPUT_DIR is never touched in
# USE_PREBUILT_ALIGNMENTS mode, so a warning that names it there points
# at the wrong place and sends anyone debugging a "zero samples remain"
# failure down the wrong path (a dangling symlink under
# mini_aligned_root() is a common cause).
_BAM_CHECK_DIR = mini_aligned_root() if USE_PREBUILT_ALIGNMENTS else BAM_INPUT_DIR

_missing_data = sorted(b for b in BARCODES if not _has_indexed_bam_data(b))
if _missing_data:
    print(f"[sample_sheet] WARNING: {len(_missing_data)} sample(s) from the sheet have no "
          f"BAM data under {_BAM_CHECK_DIR} yet -- excluding from this run:")
    for _b in _missing_data:
        print(f"    {_b}  (role={SAMPLE_ROLE.get(_b, '') or '?'}, "
              f"family={SAMPLE_FAMILY.get(_b, '') or '?'}, "
              f"{DMR_GROUP_COLUMN}={SAMPLE_DIAGNOSIS.get(_b, '') or '(none)'})")
    _missing_set = set(_missing_data)
    BARCODES         = [b for b in BARCODES if b not in _missing_set]
    SAMPLE_DIAGNOSIS = {b: v for b, v in SAMPLE_DIAGNOSIS.items() if b not in _missing_set}
    SAMPLE_ROLE      = {b: v for b, v in SAMPLE_ROLE.items()      if b not in _missing_set}
    SAMPLE_FAMILY    = {b: v for b, v in SAMPLE_FAMILY.items()    if b not in _missing_set}

if not BARCODES:
    raise ValueError(
        f"After excluding sample-sheet entries with no indexed BAM data, zero samples "
        f"remain. Check {_BAM_CHECK_DIR} and your experiment.sample_sheet config."
    )

print(f"[sample_sheet] {len(BARCODES)} sample(s) with indexed BAM data will be processed.")
DORADO_MODELS_DIR = as_path(_tpaths.get("dorado_models"), TOOLS_DIR / "dorado")

OUTPUT_BASE_DIR = STUDY_ROOT


# ============================================================================
# Apptainer / Singularity Binary Resolution
#
# Repo-local install (workflow/tools/apptainer/bin/apptainer) is preferred so the pipeline doesn't depend on
# a personal ~/.local/bin or a cluster module being loaded. Falls back to
# whatever's on PATH (a module-provided singularity/apptainer, or a
# personal install) if no repo-local copy exists.
# ============================================================================

def _resolve_sif_engine_bin(engine: str) -> str:
    repo_local = TOOLS_DIR / engine / "bin" / engine
    if repo_local.exists():
        return str(repo_local.resolve())
    on_path = shutil.which(engine)
    if on_path:
        return on_path
    raise FileNotFoundError(
        f"'{engine}' not found at {repo_local} or on PATH.\n"
        f"Install it (see docs/no_sudo_install.md) or run 'module load {engine}'\n"
        f"if your cluster provides one."
    )


# ============================================================================
# Tool Paths
# ============================================================================

CLAIR3_DIR = (TOOLS_DIR / "clair3").resolve()

PREFER_POD5_INPUTS = bool(
    EXPERIMENT_CFG.get("prefer_pod5_inputs", config.get("prefer_pod5_inputs", False))
)


def _resolve_pypy3_bin():
    for p in sorted(CLAIR3_DIR.glob("pypy3.11-v7.3.20-*/bin/pypy3")):
        if p.exists():
            return p.resolve()
    return (CLAIR3_DIR / "pypy3.11-v7.3.20-macos_arm64" / "bin" / "pypy3").resolve()

PYPY3_BIN = _resolve_pypy3_bin()


# ============================================================================
# Validation
# ============================================================================

_REQUIRED_CONFIG_KEYS   = ["experiment", "tools", "resources"]
_VALID_EXPERIMENT_TYPES = ["adaptive", "nonadaptive"]

# Container engines that consume local .sif files / docker:// refs via the
# same CLI surface (exec, pull, --bind, --nv). Apptainer is the actively
# maintained fork of what used to be distributed as Singularity/SingularityCE;
# some newer cluster images ship only the `apptainer` binary with no
# `singularity` compatibility shim, so both names need to be treated as
# first-class options rather than hardcoding "singularity" everywhere.
SIF_ENGINES = ("singularity", "apptainer")


def _require_keys(d, keys, prefix=""):
    for key in keys:
        full = f"{prefix}.{key}" if prefix else key
        if key not in d:
            raise ValueError(f"Missing required config key: {full}")


def _ensure_exists(path, label, setup_script, expect_nonempty_dir=False, extra_hint=""):
    if not path.exists():
        msg = f"{label} not found at {path}\nPlease run: {setup_script}"
        if extra_hint:
            msg += f"\n{extra_hint}"
        raise FileNotFoundError(msg)
    if expect_nonempty_dir and not any(path.iterdir()):
        msg = f"{label} is empty at {path}\nPlease run: {setup_script}"
        if extra_hint:
            msg += f"\n{extra_hint}"
        raise FileNotFoundError(msg)
    print(f"[validate] {label}: {path}")


def _validate_config(cfg):
    _require_keys(cfg, _REQUIRED_CONFIG_KEYS)
    if "paths" in cfg:
        _require_keys(cfg["paths"], ["inputs", "outputs"], prefix="paths")
        _require_keys(cfg["paths"]["inputs"],  ["bam", "pod5"], prefix="paths.inputs")
        _require_keys(cfg["paths"]["outputs"], ["results", "logs", "slurm_logs"], prefix="paths.outputs")

    tools = cfg["tools"]
    if "dorado" not in tools and "run_dorado" not in tools:
        raise ValueError("Missing required config key: tools.dorado")
    if "clair3" not in tools and "run_clair3" not in tools:
        raise ValueError("Missing required config key: tools.clair3")

    res = cfg["resources"]
    for key in ("basecall_pod5", "merge_bams"):
        if key not in res:
            raise ValueError(f"Missing required config key: resources.{key}")
    if "call_variants" not in res and "clair3" not in res:
        raise ValueError("Missing required config key: resources.clair3")

    exp_type = cfg["experiment"]["type"]
    if exp_type not in _VALID_EXPERIMENT_TYPES:
        raise ValueError(f"Invalid experiment.type '{exp_type}'. Must be one of {_VALID_EXPERIMENT_TYPES}")


def _validate_tools():
    print("\n" + "=" * 78)
    print("VALIDATING TOOL INSTALLATIONS")
    print("=" * 78 + "\n")

    engine = config.get("cluster", {}).get("container_engine", "docker")
    valid_engines = ("docker",) + SIF_ENGINES
    if engine not in valid_engines:
        raise ValueError(
            f"cluster.container_engine must be one of {valid_engines}, got '{engine}'"
        )

    # Only docker gets a plain shutil.which() check here -- it's never
    # vendored into the repo, so PATH is the only place it can live. SIF
    # engines (singularity/apptainer) are resolved later via
    # _resolve_sif_engine_bin()/CONTAINER_ENGINE_BIN, which checks the
    # repo-local workflow/tools/<engine>/bin/<engine> location first and
    # only falls back to PATH -- checking shutil.which() here as well
    # would print a misleading "(None)" for a valid repo-local install
    # that simply isn't on PATH.
    if engine == "docker" and not shutil.which("docker"):
        raise FileNotFoundError(
            "Container engine 'docker' not found in PATH.\n"
            f"Install docker, or set cluster.container_engine to one of {SIF_ENGINES}."
        )

    print("[validate] Clair3: using containerized run_clair3.sh -- "
          "skipping host source-install checks")
    print("\n" + "=" * 78 + "\n")


_validate_config(config)
_validate_tools()

EXPERIMENT_TYPE = config["experiment"]["type"]


# ============================================================================
# Cluster Profile
# ============================================================================

_cluster_profile  = config.get("cluster", {}).get("profile", "cpu")
_container_engine = config.get("cluster", {}).get("container_engine", "docker")
_CONTAINER_ENGINE = _container_engine  # logical name -- "docker" / "singularity" / "apptainer"

# Absolute path to the actual binary to invoke. For docker this is just
# "docker" (assumed on PATH, since Docker isn't something we vendor into
# the repo); for SIF engines it's resolved repo-local-first, falling back
# to PATH (e.g. a cluster module) if no repo-local copy exists.
CONTAINER_ENGINE_BIN = (
    _resolve_sif_engine_bin(_CONTAINER_ENGINE)
    if _CONTAINER_ENGINE in SIF_ENGINES
    else _CONTAINER_ENGINE
)
print(f"[validate] Container engine: {_CONTAINER_ENGINE} ({CONTAINER_ENGINE_BIN})")

_clair3_res_cfg = RESOURCE_CFG.get("call_variants", {})
_clair3_profile = _clair3_res_cfg.get(_cluster_profile, _clair3_res_cfg.get("cpu", {}))

RUN_CLAIR3_USE_GPU = bool(_clair3_profile.get("use_gpu", False))


# ============================================================================
# Clair3 Container Images
# ============================================================================

_CLAIR3_IMAGE_CPU        = "hkubal/clair3:v2.0.2"
_CLAIR3_IMAGE_CPU_SIF    = str(TOOLS_DIR / "singularity" / "clair3_v2.0.2.sif")
_CLAIR3_IMAGE_GPU_DOCKER = "hkubal/clair3:v2.0.2_gpu"
_CLAIR3_IMAGE_GPU_SIF    = str(TOOLS_DIR / "singularity" / "clair3_v2.0.2_gpu.sif")

if RUN_CLAIR3_USE_GPU:
    _CLAIR3_CONTAINER = _CLAIR3_IMAGE_GPU_SIF if _CONTAINER_ENGINE in SIF_ENGINES else _CLAIR3_IMAGE_GPU_DOCKER
else:
    _CLAIR3_CONTAINER = _CLAIR3_IMAGE_CPU_SIF if _CONTAINER_ENGINE in SIF_ENGINES else _CLAIR3_IMAGE_CPU

# When running under a SIF-based engine (singularity or apptainer), always
# use a pre-built local .sif instead of a bare docker:// tag -- passing a
# docker:// reference to `<engine> exec` makes it re-pull *and* re-convert
# the image to SIF format on every single invocation (observed costing
# ~6-7 minutes of wall time per call_variants run for the CPU image alone,
# with zero actual variant calling happening during that time). Building
# the .sif once and reusing it turns that overhead into effectively
# nothing on every run after the first. This used to only be done for the
# GPU image; the CPU path is the far more commonly hit one, so it needed
# the same treatment. The .sif format itself is shared between Singularity
# and Apptainer, so the same cached file works for either engine.
if _CONTAINER_ENGINE in SIF_ENGINES:
    _sif_path, _sif_docker_tag = (
        (_CLAIR3_IMAGE_GPU_SIF, _CLAIR3_IMAGE_GPU_DOCKER) if RUN_CLAIR3_USE_GPU
        else (_CLAIR3_IMAGE_CPU_SIF, _CLAIR3_IMAGE_CPU)
    )
    _ensure_exists(
        Path(_sif_path),
        f"Clair3 {'GPU' if RUN_CLAIR3_USE_GPU else 'CPU'} {_CONTAINER_ENGINE} image",
        f"{_CONTAINER_ENGINE} pull {_sif_path} docker://{_sif_docker_tag}",
    )


def singularity_ref(image: str) -> str:
    """
    Build the reference `<engine> exec`/`<engine> pull` expects (works for
    both singularity and apptainer, which share the same CLI/SIF format).
    Local .sif paths are used as-is; anything else is treated as a
    Docker Hub tag and gets a docker:// prefix.
    """
    return image if str(image).endswith(".sif") else f"docker://{image}"


# ============================================================================
# Container Bind Targets  (portability fix)
#
# Any host path handed to a container (Clair3, modkit, plotting scripts,
# etc.) needs its containing directory explicitly bound in --
# singularity/apptainer only auto-mount $HOME/$TMP/cwd by default, and
# `docker run` mounts nothing unless told to. The previous version of this
# pipeline hardcoded a fixed list of four "blessed" roots (workflow root,
# results dir, resources dir, examples dir) into each rule's shell block --
# any config path living outside those four (a reference genome, GTF, BED,
# or BAM index living on a different mount) was silently invisible inside
# the container, producing a misleading "file not found" even though the
# file plainly exists on the host.
#
# Instead of hardcoding roots, binds are derived from wherever the
# *configured* paths actually resolve to on this cluster. A reference
# genome, annotation BED, GTF, or BAM/POD5 index can each live on a
# completely different filesystem/mount and still get mounted correctly --
# porting this pipeline to a new HPC or a new dataset needs zero edits
# here, only a new config file with different (possibly non-overlapping)
# paths.
# ============================================================================

def container_binds(*paths):
    """
    Return a deduplicated, minimal set of directories to --bind so that
    every path in `paths` is visible inside a container. "Minimal" means a
    directory is dropped if it's already covered by a parent directory
    also in the set, so overlapping repo dirs don't produce redundant
    --bind flags. Falsy entries (None, "") are ignored so callers can pass
    optional config paths straight through.
    """
    dirs = set()
    for p in paths:
        if not p:
            continue
        p = Path(p).resolve()
        dirs.add(str(p if p.is_dir() else p.parent))
    minimal = []
    for d in sorted(dirs):
        if not any(d != other and d.startswith(other.rstrip("/") + "/") for other in dirs):
            minimal.append(d)
    return minimal


# Config-driven paths that might live anywhere, gathered generically so a
# new dataset/cluster's config just works without touching this file.
_mode_annotation_paths = [
    _resolve_repo_path(v) for v in _ANNOTATIONS_CFG.get(GENOME_MODE, {}).values()
]
_gtf_path = _resolve_repo_path(ANALYSIS_CFG.get("gtf_file"))

CONTAINER_BIND_TARGETS = container_binds(
    WORKFLOW_ROOT,
    RESULTS_ROOT,
    REPO_ROOT / "resources",
    REPO_ROOT / "examples",
    BAM_INPUT_DIR,
    POD5_INPUT_DIR,
    DORADO_MODELS_DIR,
    _gtf_path,
    *_mode_annotation_paths,
)

# Pre-rendered --bind / -v flag strings, ready to interpolate straight
# into a rule's `shell:` block via params. Built once here rather than
# per-rule so every containerized rule stays in sync automatically -- add
# a new config path to CONTAINER_BIND_TARGETS above and every rule that
# uses these picks it up, instead of needing a `--bind` added by hand in
# N different places.
SINGULARITY_BIND_ARGS = " ".join(f"--bind {d}:{d}" for d in CONTAINER_BIND_TARGETS)
DOCKER_BIND_ARGS       = " ".join(f"-v {d}:{d}"     for d in CONTAINER_BIND_TARGETS)

print(f"[validate] Container bind targets ({len(CONTAINER_BIND_TARGETS)}):")
for _d in CONTAINER_BIND_TARGETS:
    print(f"    {_d}")


# ============================================================================
# Device Detection (for native tools e.g. Dorado)
# ============================================================================

import platform as _platform
import subprocess as _subprocess


def _detect_device() -> str:
    """Return 'metal', 'cuda', or 'cpu' based on the current host."""
    if _platform.system() == "Darwin" and _platform.machine() == "arm64":
        return "metal"
    try:
        if _subprocess.run(["nvidia-smi"], capture_output=True).returncode == 0:
            return "cuda"
    except FileNotFoundError:
        pass
    return "cpu"


RUN_DEVICE      = _detect_device()
RUN_DEVICE_FLAG = f"--device {RUN_DEVICE}" if RUN_DEVICE != "cpu" else ""


# ============================================================================
# PLINK2 Binary Resolution
#
# plink2 has no osx-arm64 build in bioconda. Rather than forcing the whole
# conda ecosystem (including R/Bioconductor envs) into osx-64/Rosetta just
# for this one binary, plink2 is resolved as an external executable:
#   - tools.plink2_path in config, if set (e.g. a native arm64 binary
#     downloaded from https://www.cog-genomics.org/plink/2.0/)
#   - otherwise, whatever "plink2" resolves to on PATH (works fine on
#     Linux/HPC and Intel Mac, where the conda-installed plink2 is native)
# ============================================================================

def _resolve_plink2():
    configured = TOOLS_CFG.get("plink2_path")
    if configured:
        path = Path(configured).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(
                f"tools.plink2_path is set to {path} but that file doesn't exist.\n"
                f"Download a native plink2 binary for your platform from "
                f"https://www.cog-genomics.org/plink/2.0/ and update the path."
            )
        return str(path)

    if _platform.system() == "Darwin" and _platform.machine() == "arm64":
        raise ValueError(
            "Running on Apple Silicon (arm64) but tools.plink2_path is not set.\n"
            "plink2 has no osx-arm64 build in bioconda, so it can't be installed "
            "into the conda env on this machine.\n"
            "Download the native arm64 binary from "
            "https://www.cog-genomics.org/plink/2.0/ and set in your config:\n\n"
            "  tools:\n"
            "    plink2_path: /path/to/plink2\n"
        )

    on_path = shutil.which("plink2")
    if on_path:
        return on_path
    return "plink2"  # let conda-env PATH resolve it at rule execution time

PLINK2_BIN = _resolve_plink2()
print(f"[validate] plink2: {PLINK2_BIN}")


IS_APPLE_SILICON = (_platform.system() == "Darwin" and _platform.machine() == "arm64")

def plink2_env_yml():
    """
    plink2 has no osx-arm64 bioconda build. On Apple Silicon, use the env
    that omits plink2 itself (PLINK2_BIN resolves to an external native
    binary via tools.plink2_path instead, so this env's solve stays
    native and never needs Rosetta). Everywhere else, use the env that
    installs plink2 directly from bioconda.
    """
    return str(ENVS_DIR / ("plink2-env-arm64.yml" if IS_APPLE_SILICON else "plink2-env.yml"))

# ============================================================================
# Memory Budget
# ============================================================================

def _get_total_ram_gb() -> float:
    """Return total system RAM in GB, with a safe fallback."""
    try:
        import psutil
        return psutil.virtual_memory().total / 1024 ** 3
    except ImportError:
        pass
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) / 1024 / 1024
    except FileNotFoundError:
        pass
    try:
        out = _subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True)
        if out.returncode == 0:
            return int(out.stdout.strip()) / 1024 ** 3
    except FileNotFoundError:
        pass
    return 16.0


_TOTAL_RAM_GB    = _get_total_ram_gb()
_RAM_HEADROOM    = 0.85
RUN_LOCAL_MEM_GB = _TOTAL_RAM_GB * _RAM_HEADROOM

print(f"[sys] RAM: {_TOTAL_RAM_GB:.1f} GB | budget: {RUN_LOCAL_MEM_GB:.1f} GB")


# ============================================================================
# Path Helper Functions
# ============================================================================

def sample(barcode, subdir, *parts):
    """Per-sample result path: results/{STUDY}/samples/{barcode}/{subdir}/..."""
    return os.path.join(str(OUTPUT_BASE_DIR), "samples", str(barcode), subdir, *parts)

def cohort(*parts):
    """Study-level output path: results/{STUDY}/..."""
    return os.path.join(str(OUTPUT_BASE_DIR), *parts)

def family_dir(family_id, *parts):
    """
    Family-level result path: results/{STUDY}/families/{family_id}/...
    Mirrors sample()/cohort() -- a place for outputs that operate on a
    whole family at once (e.g. a per-family trio/pedigree summary) rather
    than a single sample or the whole cohort. family_id should be one of
    the values in SAMPLE_FAMILY; "" (unknown/no-family samples, e.g.
    id_col-mode sheets) is a valid value and lands under families/_none/.
    """
    fam = str(family_id).strip() or "_none"
    return os.path.join(str(OUTPUT_BASE_DIR), "families", fam, *parts)

def sample_log(barcode, name):
    """Per-sample log path: logs/samples/{barcode}/{name}"""
    return os.path.join(str(LOGS_ROOT), "samples", str(barcode), name)

def cohort_log(name):
    """Cohort-level log path: logs/cohort/{name}"""
    return os.path.join(str(LOGS_ROOT), "cohort", name)


# ============================================================================
# Write sample metadata TSV  (needs both RESULTS_ROOT and BARCODES)
# ============================================================================

_META_TSV = STUDY_ROOT / "sample_metadata.tsv"
_META_TSV.parent.mkdir(parents=True, exist_ok=True)

with open(_META_TSV, "w", newline="") as _f:
    _w = csv.writer(_f, delimiter="\t")
    _w.writerow(["sample_id", DMR_GROUP_COLUMN, "role", "family_id", "bed_path"])

    for _sid in BARCODES:
        _w.writerow([
            _sid,
            SAMPLE_DIAGNOSIS.get(_sid, ""),
            SAMPLE_ROLE.get(_sid, ""),
            SAMPLE_FAMILY.get(_sid, ""),
            sample(_sid, "methylation", f"{_sid}_modkit.bed"),
        ])

SAMPLE_METADATA_TSV = str(_META_TSV)
print(f"[sample_sheet] metadata TSV written -> {_META_TSV} (group column: '{DMR_GROUP_COLUMN}')")


# ============================================================================
# Input Resolution Functions
# ============================================================================

def get_bam_paths(wc):
    barcode_dir = Path(BAM_INPUT_DIR) / wc.barcode

    # BAMs directly under bam_index/<barcode>/
    files = sorted(glob(os.path.join(str(barcode_dir), "*.bam")))
    if files:
        return files

    # BAMs nested below bam_index/<barcode>/ (symlink targets included).
    # Scoped to this barcode's directory only -- never scans the whole
    # BAM_INPUT_DIR tree.
    if barcode_dir.is_dir():
        files = sorted(_iter_barcode_bams(barcode_dir))
        pass_only = [f for f in files if "/pass/" in f]
        if pass_only:
            return pass_only
        if files:
            return files

    # Special case: BAM_INPUT_DIR itself is a single BAM file.
    bam_input = Path(BAM_INPUT_DIR)
    if bam_input.is_file() and bam_input.suffix == ".bam":
        return [str(bam_input.resolve())]

    raise ValueError(f"No BAM files found for {wc.barcode} in {BAM_INPUT_DIR}")


def get_pod5_paths(wc):
    files = glob(os.path.join(str(POD5_INPUT_DIR), wc.barcode, "*.pod5"))
    files = [p for p in files if wc.barcode in os.path.basename(p)]
    if files:
        return sorted(files)

    all_pod5 = glob(os.path.join(str(POD5_INPUT_DIR), "**", "*.pod5"), recursive=True)
    if not all_pod5:
        raise ValueError(f"No POD5 files found for {wc.barcode} in {POD5_INPUT_DIR}")

    matched = [p for p in all_pod5 if wc.barcode in os.path.basename(p)]
    if matched:
        return sorted(matched)

    raise ValueError(f"No POD5 files found for {wc.barcode} in {POD5_INPUT_DIR}")


# ============================================================================
# Diagnosis Group Helpers  (used by DMR rules)
# ============================================================================

def get_barcodes_by_diagnosis(diagnosis):
    """Return barcodes whose diagnosis matches the given string."""
    return [b for b in BARCODES if SAMPLE_DIAGNOSIS.get(b, "") == diagnosis]


def get_diagnosis_param(group):
    """Comma-separated barcode list for a diagnosis group, for passing to R."""
    return ",".join(get_barcodes_by_diagnosis(group))


# ============================================================================
# Family Helpers  (mirror the diagnosis helpers above, keyed on family_id
# instead of diagnosis -- for rules that need to operate per-family, e.g.
# a trio/pedigree summary, rather than per-sample or per-diagnosis-group)
# ============================================================================

def get_barcodes_by_family(family_id):
    """Return barcodes belonging to the given family_id, in BARCODES order."""
    fam = str(family_id).strip()
    return [b for b in BARCODES if SAMPLE_FAMILY.get(b, "") == fam]


def get_family_param(family_id):
    """Comma-separated barcode list for a family, for passing to R/scripts."""
    return ",".join(get_barcodes_by_family(family_id))


def get_all_family_ids():
    """
    All distinct non-empty family IDs present in the current BARCODES set,
    in first-seen order. Empty/unknown family_id ("") is excluded -- callers
    that want to iterate "one output per family" (e.g. a family_dir() target
    per pedigree) shouldn't get a bogus family for id_col-mode sheets or
    legacy range/list barcode configs, where family_id is always "".
    """
    seen = []
    for b in BARCODES:
        fam = SAMPLE_FAMILY.get(b, "")
        if fam and fam not in seen:
            seen.append(fam)
    return seen


# ============================================================================
# Resource Helper
# ============================================================================

def res(rule_name, key, default=None):
    if rule_name.startswith("_"):
        return default
    r = config["resources"].get(rule_name, {})
    if key in r:
        return r[key]
    if default is not None:
        return default
    raise KeyError(f"Missing resource key '{key}' for rule '{rule_name}'")


# ============================================================================
# Shell Utilities
# ============================================================================

_SHELL_UTILS = (SCRIPTS_DIR / "utils" / "shell_utils.sh").resolve()

SHELL_PRELUDE = (
    "set -euo pipefail\n"
    "exec > {log} 2>&1\n"
    f"source {_SHELL_UTILS}\n"
    'shell_start "{rule}"\n'
    "echo ''\n"
)

SHELL_EPILOG = (
    "echo ''\n"
    'shell_end "{rule}"\n'
)

SHELL_PRELUDE_ORDERED = (
    "set -euo pipefail\n"
    f"source {_SHELL_UTILS}\n"
    'shell_start "{rule}" >> {log} 2>&1\n'
    "echo '' >> {log}\n"
)

SHELL_EPILOG_ORDERED = (
    "echo '' >> {log}\n"
    'shell_end "{rule}" >> {log} 2>&1\n'
)

_PTY_WRAPPER = (SCRIPTS_DIR / "utils" / "run_with_pty.py").resolve()
