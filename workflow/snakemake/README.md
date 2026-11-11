# Snakemake Workflow Structure

This directory contains the Snakemake workflow definition for the nanopore
adaptive sequencing pipeline.

The top-level entry point is:

```
workflow/snakemake/nanopore_adaptive_pipeline.smk
```

This file defines the overall workflow and imports the individual rule
files located under `rules/`.

## Workflow overview

The pipeline DAG follows this general structure:

```
POD5 data
   │
   ▼
Basecalling (Dorado)
   │
   ▼
Read alignment
   │
   ▼
Variant calling (Clair3)
   │
   ▼
Haplotype phasing
   │
   ├──────────────────────┐
   ▼                       ▼
Methylation           Variant processing
(Modkit)               (VCF/PGEN)
   │                       │
   ▼                       ▼
DMR analysis           methylQTL analysis
(DSS/DSSseq)            (TensorQTL)
   │
   ▼
Visualization
(methylartist, CpG histograms)
```

## Directory structure

```
workflow/snakemake/
├── common.smk
├── nanopore_adaptive_pipeline.smk
└── rules/
    ├── alignment/
    ├── basecalling/
    ├── dmr/
    ├── methylation/
    ├── methylqtl/
    ├── mini_genomes/
    ├── plots/
    └── variants/
```

## Top-level files

### `nanopore_adaptive_pipeline.smk`

Main workflow entry point. Responsibilities:

- Loads pipeline configuration
- Defines sample metadata
- Includes rule modules
- Connects the full analysis DAG

Run the pipeline using:

```bash
snakemake all \
    -s workflow/snakemake/nanopore_adaptive_pipeline.smk \
    --configfile <config.yml>
```

### `common.smk`

Shared workflow utilities. Contains:

- Common paths
- Configuration parsing
- Validation helpers
- Shared resource definitions
- Tool/container paths

Rule files should use values defined here instead of hard-coding paths.

## Rule modules

### `rules/basecalling/` — Oxford Nanopore basecalling

| File | Purpose |
|---|---|
| `basecall_pod5.smk` | Runs Dorado basecalling on POD5 files |

### `rules/alignment/` — Read alignment and BAM preparation

| File | Purpose |
|---|---|
| `sort_index_bam.smk` | Sorts and indexes BAM files |
| `merge_bams.smk` | Combines BAM files when required |

### `rules/variants/` — Variant calling workflow

| File | Purpose |
|---|---|
| `call_variants.smk` | Runs Clair3 variant calling and phasing |

### `rules/methylation/` — DNA methylation processing

| File | Purpose |
|---|---|
| `modkit_pileup.smk` | Generates methylation pileups using Modkit |

### `rules/dmr/` — Differential methylation analysis

| File | Purpose |
|---|---|
| `dss_smooth_betas.smk` | Smooths methylation beta values |
| `dss_call_dmrs.smk` | Calls DMRs using DSS |
| `dmrseq_call_dmrs.smk` | Alternative DMR calling using DMRseq |
| `annotate_dmrs.smk` | Annotates DMR regions |

### `rules/methylqtl/` — Methylation quantitative trait locus analysis

| File | Purpose |
|---|---|
| `compute_pca_covariates.smk` | Generates PCA covariates |
| `merge_vcfs.smk` | Combines variant files |
| `reheader_vcf.smk` | Adjusts VCF sample headers |
| `vcf_to_pgen.smk` | Converts VCF to PLINK2 format |
| `run_tensorqtl.smk` | Runs TensorQTL analysis |

### `rules/plots/` — Visualization and reporting

| File | Purpose |
|---|---|
| `methylartist.smk` | Generates methylation visualizations |
| `plot_cpg_hists.smk` | Creates CpG methylation distribution plots |

### `rules/mini_genomes/` — Resources for generating small test datasets

| File | Purpose |
|---|---|
| `mini_data.smk` | Builds mini reference resources used for smoke tests |

## Adding or modifying rules

When adding a new analysis step:

1. Create a new rule file under the appropriate directory.
2. Add the file to `nanopore_adaptive_pipeline.smk` using `include:`.
3. Add required software dependencies to `workflow/envs/`.
4. Update configuration options under `configs/`.
5. Test using a dry run:

```bash
snakemake all \
    -s workflow/snakemake/nanopore_adaptive_pipeline.smk \
    --configfile configs/mini_data.yml \
    --use-conda \
    -n \
    -p
```

This also gives future contributors a map of where a change belongs:

- Changing Clair3 → `rules/variants`
- Changing methylation calling → `rules/methylation`
- Changing DMR statistics → `rules/dmr`
- Adding a new tool → add an env YAML + rule module

## Development notes

- Avoid hard-coded paths; use configuration values from YAML files.
- Keep tool-specific dependencies in `workflow/envs/`.
- Use the shared helpers in `common.smk` when possible.
- Test changes against the mini genome configurations before running
  production-scale datasets.