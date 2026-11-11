# Analysis

Downstream analyses built on top of pipeline outputs from `workflow/`.
Each subdirectory is a self-contained analysis with its own scripts,
configuration, and SLURM submission scripts.

## Structure

```
analysis/
└── dmr/       # Differential methylation region calling and annotation
```

Planned:
- `cnv/`       # Copy number variation
- `sv/`        # Structural variant characterization

## How to add a new analysis

```
analysis/
└── my_analysis/
    ├── README.md              # what it does, inputs, outputs, tools
    ├── config/                # R/Python config, sample sheets
    ├── data/                  # static reference files (gittracked, small)
    ├── 01_first_step/         # scripts in order
    ├── 02_second_step/
    ├── submit/                # or colocate submit scripts with step dirs
    └── results/               # gitignored outputs
```

## Relationship to workflow/

Pipeline outputs (BAMs, CpG pileups, variant calls) are produced by
`workflow/snakemake/` and written to `results/`. The analyses here treat
those outputs as inputs — they are not run by Snakemake.
