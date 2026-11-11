# Nanopore Adaptive-Sampling Methylation Pipeline

A Snakemake pipeline that takes Oxford Nanopore adaptive-sampling data from raw
signal to haplotype-aware methylation results: basecalling, alignment, variant
calling and phasing, per-CpG methylation, differentially methylated region (DMR)
calling, and methylation QTL mapping. It runs on Slurm clusters (CPU or GPU)
and on a laptop, including Apple Silicon, using a small test reference.

## Pipeline overview

```
POD5 ──► Dorado basecalling + minimap2 alignment
      ──► merge / sort / index BAM (streamed)
      ──► Clair3 variant calling + WhatsHap phasing
      ──► modkit per-CpG 5mC pileup ──► methylartist plots, CpG-island histograms, XCI summaries
      ──► DMR calling (DSS or dmrseq) ──► smoothed betas, WGCNA modules,
                                           DMR annotation, pathway enrichment, heatmaps
      ──► methylQTL: merge VCFs ► PLINK2 ► PCA covariates ► TensorQTL
```

## Design highlights

- **Config-driven.** Samples, references, resources and tool choices come from a YAML
  config and a sample sheet (CSV or Excel). Nothing about a particular cohort is hard-coded.
- **Container-agnostic.** Docker, Apptainer or Singularity, selected with one config key.
  Container bind mounts are derived automatically from the paths in the config.
- **Runs anywhere.** Device detection picks Metal, CUDA or CPU for Dorado. Slurm and local
  profiles are in `profiles/`.
- **Disk-friendly.** BAM merge and sort are streamed, so no large unsorted intermediate is written.
- **Fails early with clear messages.** Config validation, tool checks, and samples with no
  data are reported before the DAG is built.
- **Reproducible.** Every wrapper launch writes a run manifest (git commit, config, tool versions).
- **Pluggable analysis.** DSS or dmrseq for DMR calling. A grouping column such as diagnosis,
  sex or batch is configurable, and calling can be restricted to one family role.
- **Two modes.** A small "mini" genome of a single locus for fast tests, and full hg38 for real runs.

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/nanopore-snakemake-pipeline.git
cd nanopore-snakemake-pipeline
```

### 2. Install Miniforge

This pipeline uses Conda to manage its software dependencies. We recommend
[Miniforge](https://github.com/conda-forge/miniforge). Follow the instructions for your
operating system, then restart your shell so the `conda` command is available.

### 3. Create the Snakemake environment

```bash
# Linux / HPC
conda env create -f workflow/envs/smk-env-linux.yml
conda activate "$(grep '^name:' workflow/envs/smk-env-linux.yml | awk '{print $2}')"

# macOS: use workflow/envs/smk-env-macos.yml instead
```

This environment provides Snakemake and what is needed to launch the workflow. All
tool-specific environments are created automatically by Snakemake on the first run.

#### Apple Silicon (M1/M2/M3/M4)

Three R/Bioconductor packages (bsseq, DSS, and their WGCNA/pheatmap dependents) have no
native osx-arm64 conda build. Before your first run:

```bash
./setup/setup_apple_silicon_envs.sh
```

This builds the three affected environments under Rosetta (osx-64) at the paths Snakemake
expects, so `--use-conda` treats them as already built. Run it once per fresh clone, or any
time you delete `.snakemake/conda/`.

You also need a native arm64 `plink2` binary, since no conda package exists. Download it
from <https://www.cog-genomics.org/plink/2.0/> and set `tools.plink2_path` in your config.

Intel Macs and Linux/HPC do not need this step.

### 4. Install a container engine

Clair3 and some other tools run in containers. Install Docker, Apptainer or Singularity,
or load one on your cluster (`module load apptainer`), and set `cluster.container_engine`
in your config.

With Apptainer or Singularity, pull the Clair3 image once (use `clair3_v2.0.2_gpu` for GPU runs):

```bash
mkdir -p workflow/tools/singularity
apptainer pull workflow/tools/singularity/clair3_v2.0.2.sif docker://hkubal/clair3:v2.0.2
```

Without root access, see [`docs/no_sudo_install.md`](docs/no_sudo_install.md).

### 5. Install Dorado

```bash
setup/setup_dorado.sh
```

See [`setup/README.md`](setup/README.md) for platform-specific details.

### 6. Prepare the test data

The mini configs use public GIAB-derived fixtures. See
[`docs/mini_test_data_giab_2023.md`](docs/mini_test_data_giab_2023.md) for how to obtain them.

## Run the pipeline

First do a dry run to check that Snakemake resolves the workflow (example: GPU workstation):

```bash
snakemake all \
    -s workflow/snakemake/nanopore_adaptive_pipeline.smk \
    --configfile configs/mini_data/mini_data_public_2023_workstation_gpu.yml \
    --cores 6 \
    --use-conda \
    --conda-frontend conda \
    -n -p
```

`-n` does a dry run and shows what would execute. `-p` prints the shell commands. To run
for real, remove `-n`. The first execution can take several minutes while Snakemake builds
the tool-specific environments from `workflow/envs/`. See
[`workflow/envs/README.md`](workflow/envs/README.md) for details, including Apple Silicon notes.

The wrappers `scripts/run_mac.sh` and `scripts/run_hpc.sh` add preflight checks and write a
run manifest:

```bash
bash scripts/run_mac.sh --configfile configs/mini_data/mini_data_public_2023_mac_local.yml
bash scripts/run_hpc.sh --configfile configs/full_data/TEMPLATE_full_data.yml
```

Other mini configs cover HPC (A100 or CPU-only) and workstation GPUs.

### Execution modes and targets

| Setting | Values | Meaning |
| --- | --- | --- |
| `experiment.genome` | `mini`, `hg38` | Mini test reference vs. full human reference |
| `targets` | `full`, `modkit_pileup`, `plots` | Complete pipeline, stop after methylation pileup, or plots only |
| `dmr.caller` | `dss`, `dmrseq` | DSS suits small cohorts, dmrseq larger ones |

## Running on your own data

1. Copy `configs/full_data/TEMPLATE_full_data.yml` and set your input paths, reference and sample sheet.
2. Point `paths.inputs.bam` at a directory containing one folder per sample ID. The helper
   `scripts/index_ont_bam_dirs.sh` can build this index from your sequencing run folders.
3. Set `cluster.container_engine` and the `resources:` block for your scheduler.
4. Launch with `scripts/run_hpc.sh --configfile your_config.yml`.

## Repository layout

| Path | Contents |
| --- | --- |
| `workflow/snakemake/` | Entry point, shared helpers (`common.smk`), and rules |
| `workflow/scripts/` | Helper scripts, including mini-genome builders |
| `workflow/envs/` | Conda environments, one per tool group |
| `analysis/dmr/` | DMR calling and downstream R analysis |
| `configs/` | Config templates and mini-test configs |
| `profiles/` | Snakemake profiles for local and Slurm runs |
| `resources/` | Reference BEDs and the mini single-locus genome |
| `scripts/` | Launch wrappers and helpers |
| `setup/` | Tool installers (Dorado, Clair3) and Apple Silicon environment setup |
| `docs/` | Setup notes and test-data documentation |
| `examples/` | Example sample manifests with placeholder IDs |

## Data and privacy

This repository contains **no participant or patient data**. Test data are derived from
public GIAB samples. Real sample manifests and cohort paths are kept outside version control.

## Author

<Your Name>, PhD candidate in bioinformatics, UC Davis.

## License

MIT