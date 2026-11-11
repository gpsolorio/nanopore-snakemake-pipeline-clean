# Setup Scripts

This directory contains one-time setup utilities for installing external tools
and conda environments required by the nanopore pipeline.

Supported platforms: **Linux** (x86_64, with or without NVIDIA GPU),
**macOS on Apple Silicon** (arm64: M1/M2/M3/M4), and **macOS on Intel**
(x86_64).

## Quick reference: what to run on each platform

| Step                        | Linux                     | macOS (Apple Silicon / arm64) | macOS (Intel / x86_64) |
|-----------------------------|---------------------------|--------------------------------|--------------------------|
| Clair3                      | `setup_clair3_linux.sh`   | `setup_clair3_macos.sh`        | `setup_clair3_macos.sh`  |
| Dorado                      | `setup_dorado.sh`         | `setup_dorado.sh`               | `setup_dorado.sh`        |
| dss-env / wgcna-env / pheatmap-env | none needed — solve natively | `setup_apple_silicon_envs.sh` | none needed — solve natively |
| plink2                      | resolved via bioconda, no script needed | installed by `setup_apple_silicon_envs.sh` (native arm64 binary) | resolved via bioconda, no script needed |

Everything else (all other conda envs) is created automatically by
Snakemake's normal `--use-conda` flow on every platform — no manual step
required.

---

## Linux

### Clair3 (recommended)

Use the Singularity/Apptainer-based installer:

```bash
./setup/setup_clair3_linux.sh
```

This script:

- Detects whether an NVIDIA GPU is available.
- Downloads the appropriate Clair3 CPU and/or GPU Singularity images.
- Verifies the container installation.
- Places the images under `workflow/tools/singularity/`.

### Conda environments

No manual steps needed. `bioconductor-bsseq` / `bioconductor-dss` (and
everything else) have native `linux-64` builds, so
`snakemake --use-conda` solves and creates every environment on its own.

### plink2

Resolved automatically from bioconda (`plink2-env.yml`) — no separate
binary install required, unlike macOS.

---

## macOS

### Clair3

Apple Silicon and Intel Macs do not support the Linux Singularity
containers used by the pipeline. Instead, install Clair3 from source —
this script works on **both** arm64 and x86_64 Macs:

```bash
./setup/setup_clair3_macos.sh
```

This script:

- Creates a dedicated conda environment.
- Installs Clair3 from source.
- Installs PyTorch.
- Builds the required C components.
- Installs PyPy.
- Optionally downloads pretrained models.

### Conda environments: `dss-env`, `wgcna-env`, `pheatmap-env`

`bioconductor-bsseq` and `bioconductor-dss` have no native `osx-arm64`
bioconda build (confirmed absent at every version; present for `osx-64`
and `linux-64`). This affects three environments that depend on them:
`dss-env`, `wgcna-env`, and `pheatmap-env`.

**On Intel Macs (x86_64):** no action needed. These environments solve
*natively* as `osx-64` — the same architecture the bioconda packages are
built for. Just run the pipeline normally with `--use-conda`.

**On Apple Silicon (arm64):** Snakemake's normal `--use-conda` auto-create
will always fail for these three, because it solves natively for
`osx-arm64`, where the packages don't exist. Run this once per fresh
clone (or per `.snakemake/` wipe) before your first real pipeline run:

```bash
./setup/setup_apple_silicon_envs.sh \
    -s workflow/snakemake/nanopore_adaptive_pipeline.smk \
    -c configs/<your_config>.yml
```

This script:

- Exits immediately (no-op) if not run on Darwin/arm64.
- Iteratively runs Snakemake's `--conda-create-envs-only` pass, letting
  every ordinary env solve natively.
- For `dss-env`, `wgcna-env`, and `pheatmap-env`, catches the expected
  arm64 solve failure and rebuilds each one manually under Rosetta
  (`CONDA_SUBDIR=osx-64`) at the **exact hash-named path Snakemake itself
  computed** for it. Once built there, Snakemake sees the env as already
  present and never tries to solve it natively again.
- Distinguishes *expected* failures (the three envs above) from
  *unexpected* ones — a genuine error in some other env will stop the
  script with full diagnostic output rather than being silently retried.
- Installs a native `arm64` `plink2` binary (plink2 has no `osx-arm64`
  bioconda package at all, so this is handled outside conda entirely).
- Is idempotent — safe to re-run; it skips anything already built.

> **Note: `dss-env`/`wgcna-env` manual-prefix scripts are retired.**
> This directory previously also had `setup_dss_env_macos_arm64.sh` and
> `setup_wgcna_env_macos_arm64.sh`, which built those two envs at a fixed,
> hand-chosen prefix (`.snakemake/conda/dss-env-manual`,
> `.../wgcna-env-manual`) intended to be picked up via `DSS_ENV` /
> `WGCNA_ENV` variables in `common.smk`. A check of `common.smk` confirms
> neither variable is actually defined or referenced there — the same
> file that *does* centralize equivalent platform resolution for other
> tools (`PLINK2_BIN`, `CONTAINER_ENGINE_BIN`) — so this mechanism was
> never wired into the pipeline. Both scripts (and any leftover
> `.snakemake/conda/*-manual` directories) have been removed.
> `setup_apple_silicon_envs.sh` is the sole mechanism for all three envs
> (`dss-env`, `wgcna-env`, `pheatmap-env`), via Snakemake's own hash-path
> resolution.

### plink2

Installed automatically by `setup_apple_silicon_envs.sh` on Apple
Silicon (native arm64 binary, downloaded outside conda). On Intel Macs,
`plink2-env.yml` solves natively via bioconda like any other env — no
separate script needed.

---

## Dorado

Install Oxford Nanopore's Dorado basecaller (same script on all
platforms):

```bash
./setup/setup_dorado.sh
```

Follow the script prompts to install Dorado and configure your
environment.

---

## Logs

Setup logs are written to:

```
setup/logs/
```

`setup_apple_silicon_envs.sh` additionally writes detailed per-round logs
to a temporary directory printed at the end of its run (e.g.
`/tmp/apple_silicon_env_setup.XXXXXX/`) — useful for diagnosing anything
beyond the expected arm64 bsseq/dss gap.