# Workflow Scripts

Helper scripts called by Snakemake rules.

Scripts in this directory are **not** intended to be run manually during
normal pipeline execution — they are invoked by rules under
`workflow/snakemake/rules/`.

## Layout

```
workflow/scripts/
├── methylation/     Methylation-specific preprocessing utilities
├── methylqtl/        Scripts used for TensorQTL and genotype preprocessing
├── mini_genomes/       Scripts used to generate mini reference datasets
├── utils/            Shared utilities and helper files
└── visualization/     Plotting and visualization wrappers
```

## Directories

### `methylation/`

Preprocessing utilities specific to methylation data.

| File | Purpose |
|---|---|
| `prepare_methylation_bed.py` | Prepares methylation BED input for downstream rules |

### `methylqtl/`

Scripts supporting genotype preprocessing and TensorQTL analysis.

| File | Purpose |
|---|---|
| `format_pca_covariates.py` | Formats PCA covariates for TensorQTL |
| `make_plink_aux_files.py` | Generates PLINK2 auxiliary files |
| `run_tensorqtl.py` | Runs the TensorQTL analysis |

### `mini_genomes/`

Scripts used to generate mini reference resources for smoke tests.

### `utils/`

Shared helpers and config files used across multiple rule modules.

| File | Purpose |
|---|---|
| `run_with_pty.py` | PTY wrapper used by rules that need a pseudo-terminal |
| `chr_name_conv.txt` | Chromosome name conversion table |
| `shell_utils.sh` | Shared shell helper functions |

### `visualization/`

Plotting and visualization wrappers.

## Adding a new script

Place a new script in the directory corresponding to the pipeline stage
that uses it:

- Methylation preprocessing → `methylation/`
- methylQTL / genotype prep → `methylqtl/`
- Mini test-data generation → `mini_genomes/`
- Shared across multiple stages → `utils/`
- Plotting / reporting → `visualization/`

## Environment handling

Rules activate the correct conda environment themselves via
`conda: str(ENVS_DIR / "<env>.yml")` — scripts should not activate or
`conda run` into an environment on their own. If you find a script that
does (e.g. an old `conda run -n <env> ...` wrapper), it's leftover from
before this was moved into Snakemake's own environment handling and should
be removed rather than kept as a parallel path.