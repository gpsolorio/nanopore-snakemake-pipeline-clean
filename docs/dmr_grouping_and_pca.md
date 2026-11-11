# Config-driven DMR grouping & PCA sizing

This doc covers two generalizations made to unblock the public/mini GIAB
smoke test, both of which also apply to (and now benefit) the the private cohort
production pipeline.

## 1. DMR grouping is no longer hardcoded to diagnosis/ASD/TD

`run_dss.R` previously hardcoded a `diagnosis` column and `ASD`/`TD` values.
That meant any cohort without the private cohort's diagnosis sheet (e.g. public GIAB
reference samples) failed with `No samples with ASD or TD diagnosis found`.

Grouping is now fully config-driven via two keys under `dmr:`:

```yaml
dmr:
  group_column: diagnosis   # or: sex, family, batch, anything with 2 levels
  group_values: [ASD, TD]   # must be exactly 2 values
```

- `group_column` is the metadata TSV column `run_dss.R` reads.
- `group_values` are the two labels it splits samples into.

The metadata TSV itself (`sample_metadata.tsv`, written by `common.smk`) now
uses `group_column` as its header instead of a fixed `"diagnosis"` string.

### Where the group labels come from

- **the private cohort** (`configs/nanopore_config.yml`, `configs/mini_data.yml`): labels
  still come from the real diagnosis sheet (`experiment.sample_sheet`), same
  as before — only the column name is now explicit in config rather than
  implicit in code.
- **Public/GIAB** (`configs/mini_data_public.yml`): there's no diagnosis
  sheet at all, so labels are supplied directly via a new
  `experiment.sample_groups` map:

  ```yaml
  experiment:
    sample_groups:
      hg002: male
      hg003: male
      hg004: female
  ```

  This is a generic sample_id -> label map — nothing in the code expects
  "male"/"female" specifically. Any two-level split (family, batch, etc.)
  works the same way for a future cohort.

### ⚠️ Sign convention — read this before interpreting `diff.methy`

`run_dss.R` calls:

```r
DMLtest(BSobj, group1 = group_values[1], group2 = group_values[2], ...)
```

So **`diff.methy` in `DMR_results.csv` is `group_values[1]` relative to
`group_values[2]`.**

- the private cohort config uses `group_values: [ASD, TD]` → `diff.methy` is **ASD
  relative to TD** (positive = higher methylation in ASD).
- Public config uses `group_values: [male, female]` → `diff.methy` is
  **male relative to female**.

If you need the opposite direction for a given cohort, **swap the order of
`group_values` in that cohort's YAML config** — do not change the order of
arguments in `run_dss.R` itself, since that would silently flip the
convention for every other config that doesn't expect it.

> Note: the original hardcoded version of this script called
> `DMLtest(group1 = td_samples, group2 = asd_samples)`, i.e. TD-relative-to-
> ASD. The the private cohort config here intentionally uses `[ASD, TD]` (ASD-relative-
> to-TD) instead, since that's the more natural reading direction ("compare
> ASD to TD") for anyone consuming `DMR_results.csv` going forward. If any
> downstream analysis script assumed the old TD-relative-to-ASD sign, it
> will need updating — the magnitude of `diff.methy` is unchanged, only the
> sign.

## 2. PCA covariates scale with actual sample count

`compute_pca_covariates.smk` previously always requested a fixed
`--pca {n_pcs}` from `tensorqtl.n_pcs` (default 5), regardless of how many
samples were actually in the cohort. This fails outright below 6 samples
(`Error: Too few samples to compute 5 PCs`) since you can never compute more
than `N_samples - 1` PCs.

`n_pcs` and the `approx`/exact algorithm choice are now both resolved from
the real sample count in the `.psam` file at runtime:

- `n_pcs = min(tensorqtl.n_pcs, N_samples - 1)`
- `approx` is only added once `N_samples >= tensorqtl.pca_approx_min_n`
  (default 50) — `approx` mode is a randomized algorithm meant for large-N
  biobank data and isn't appropriate (or necessary) for small cohorts.

For the 3-sample public smoke test this resolves to `--pca 2` with no
`approx` flag, instead of erroring. `format_pca_covariates.py` needed no
changes since it already builds its PC column list from `--n-pcs` rather
than assuming a fixed 5.

Both `tensorqtl.n_pcs` and `tensorqtl.pca_approx_min_n` remain configurable
per-cohort in each config's `tensorqtl:` block.
