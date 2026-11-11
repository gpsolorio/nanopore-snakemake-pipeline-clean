# workflow/envs

Conda environments used by Snakemake rules via `conda: str(ENVS_DIR / "<env>.yml")`.

## Organization principle

**One environment per distinct dependency footprint, not one per pipeline
"stage."** Every rule pointing at the same `.yml` file shares one solved
conda environment — if that env fails to solve, every rule using it fails
together. Keeping heavy or fragile package trees (WGCNA, annotatr/
AnnotationHub) in their own file means a solve failure or dependency
change in one never takes down an unrelated rule. Some duplication across
files (e.g. `bsseq` appearing in three different envs) is the accepted
cost of that isolation — don't try to eliminate it by merging files back
together.

## Current envs

| File | Used by (rule) | Scope |
|---|---|---|
| `dss-env.yml` | `dss_call_dmrs`, `dss_smooth_betas` | DSS DML testing, BSmooth beta extraction. **osx-arm64 note:** see "Apple Silicon" section below — this env cannot be auto-solved natively on Apple Silicon and requires a one-time manual setup step. |
| `wgcna-env.yml` | `wgcna_modules` | WGCNA co-methylation module detection. **osx-arm64 note:** neither `bioconductor-wgcna` nor `bioconductor-bsseq` has an arm64 bioconda build — see "Apple Silicon" section below. Same manual fixed-prefix setup as `dss-env`. |
| `pheatmap-env.yml` | `plot_dmr_heatmap` | DMR z-score heatmap plotting |
| `dmr-annotation-env.yml` | `annotate_dmr_context` | GTF gene overlap + CpG/genic context annotation. `bioconductor-annotatr` confirmed available for `osx-arm64` through r45 — solves natively on Apple Silicon, no special handling needed. (`bioconductor-org.hs.eg.db`, also in this env, hasn't been explicitly checked for arm64 but is a plain annotation-data package, not compiled code, so is expected to be fine.) |
| `enrichment-env.yml` | `pathway_enrichment` | EnrichR pathway enrichment, for `pathway_enrichment.R`. Deliberately excludes annotatr/AnnotationHub — solves natively everywhere. |
| `bcftools-env.yml` | VCF/BCF rules | bcftools |
| `plink2-env.yml` | `vcf_to_pgen`, `compute_pca_covariates` | Python/pandas support only — **plink2 itself is not installed via conda.** It has no `osx-arm64` bioconda build; rather than forcing this (and every other) env into `osx-64`/Rosetta, plink2 is resolved as an external binary at runtime via `tools.plink2_path` in config (see `PLINK2_BIN` in `common.smk`). On Linux/Intel Mac it falls back to whatever `plink2` resolves to on `PATH`. |
| `tensorqtl-env.yml` | `run_tensorqtl` | tensorqtl, torch |
| `methylartist-env.yml` | `methylartist` | methylartist plotting |
| `plotting-env.yml` | `plot_cpg_hists` | general plotting |
| `ngs-utils-env.yml` | misc | samtools/general NGS tools |
| `smk-env.yml` | *(the shell env Snakemake itself runs in, not a rule `conda:` target)* | snakemake, mamba |

## Before adding a package to an existing env

Ask: does this package share a *purpose* with what's already in this
file, or does it just happen to run at a nearby pipeline stage? If it's
the latter — e.g. a new R package needed only by one new script — give it
its own env file rather than adding it to the nearest existing one. This
is exactly the mistake that caused `annotate_dmr_context` to fail with
"no package called 'tidyr'": that dependency belonged in its own
`dmr-annotation-env.yml`, not folded into whatever env happened to be
active.

## Apple Silicon

Three unrelated compatibility gaps exist on `osx-arm64`, handled two
different ways. Don't conflate them — none of them require forcing the
*whole* conda ecosystem into `osx-64`/Rosetta, which causes its own
problems (mixed-architecture R/Bioconductor postlink failures, e.g.
`GO.db`'s `Biobase.so` load erroring with "incompatible architecture").

**plink2** — no `osx-arm64` bioconda build, but a native arm64 binary is
available directly from
[cog-genomics.org/plink/2.0](https://www.cog-genomics.org/plink/2.0/)
("macOS M1" build). `plink2-env.yml` no longer installs plink2 via conda
at all; instead, download the native binary and point your config at it:

```yaml
tools:
  plink2_path: /path/to/plink2
```

`common.smk` raises a clear error at startup on Apple Silicon if
`tools.plink2_path` isn't set, rather than silently falling through.
Linux and Intel Mac need no config change — the conda-installed `plink2`
on `PATH` is used automatically.

**bsseq / DSS** — `bioconductor-bsseq` and `bioconductor-dss` have no
`osx-arm64` bioconda build at *any* version (confirmed via
`conda search --platform osx-arm64`; both are available for `osx-64` and
`linux-64` through current Bioconductor releases). Unlike plink2, there's
no standalone native binary to swap in — these are compiled R packages
deep in the `dss-env.yml` dependency tree. The fix is to pre-create this
one env, forced to `osx-64` (Rosetta), at a **fixed prefix** that
Snakemake reuses directly instead of trying to auto-solve the env
natively:

```bash
setup/setup_dss_env_macos_arm64.sh
```

(Equivalent to running `CONDA_SUBDIR=osx-64 conda env create -p
.snakemake/conda/dss-env-manual -f workflow/envs/dss-env.yml` followed by
`conda config --set subdir osx-64 --prefix .snakemake/conda/dss-env-manual`
manually, plus a verification step and logging — see the script.)

`common.smk`'s `DSS_ENV` resolves to this fixed prefix automatically on
Apple Silicon (and raises a clear error telling you to run the above if
it doesn't exist yet), while pointing at plain `dss-env.yml` unchanged on
Linux and Intel Mac, where it solves natively with no special handling.

**WGCNA** — `bioconductor-wgcna` also has no `osx-arm64` bioconda build,
independent of the `bsseq`/`DSS` gap above. (`wgcna-env.yml` duplicates
`bioconductor-bsseq` intentionally, since `wgcna_methylation.R` needs
`DSS::makeBSseqData` and `bsseq::getMeth`/`BSmooth` for building the input
BSseq object — see the comment in `wgcna-env.yml` — so this env has two
independent arm64 gaps stacked in the same file, not just one.) Same
fixed-prefix Rosetta treatment as `dss-env`:

```bash
setup/setup_wgcna_env_macos_arm64.sh
```

`common.smk`'s `WGCNA_ENV` resolves to this fixed prefix automatically on
Apple Silicon, exactly like `DSS_ENV`.

---

This keeps every other environment (`pheatmap-env`, `dmr-annotation-env`,
`enrichment-env`, `tensorqtl-env`, etc.) solving natively as `osx-arm64`
— only `dss-env` and `wgcna-env` run under Rosetta, and only because
bioconda genuinely has no arm64 alternative for `bsseq`/`DSS`/`WGCNA` yet.
If bioconda ever ships arm64 builds for these, the manual-prefix scripts
and the `is_apple_silicon` branches in `common.smk` can simply be removed
— nothing else in the pipeline depends on the workaround existing.