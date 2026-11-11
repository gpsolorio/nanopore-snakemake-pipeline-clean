"""
run_tensorqtl.py — cis-mQTL mapping script for Snakemake

Called by rule run_tensorqtl via shell: directive.
CLI args are passed explicitly from the rule shell block.

Stage 1 — permutation pass  → one row per DMR, q-values via Storey-Tibshirani FDR
Stage 2 — nominal pass      → all SNP-DMR pairs for FDR-significant DMRs only
Stage 3 — validation check  → top variant from Stage 1 must appear in Stage 2

Q-values computed via tensorqtl.post.calculate_qvalues (R's qvalue package via rpy2),
matching the original analysis script exactly.
"""

import argparse
import glob
import os
import shutil
import sys

# R_LIBS_USER is intentionally NOT set here. rpy2 (required by
# post.calculate_qvalues) picks up R and its library path from the active
# conda env (tensorqtl-env.yml sets R_LIBS_USER="" so R only sees the env's
# own library). Setting a hardcoded user-specific path here would override
# that and risk pointing at a library built against a different R ABI on
# whatever node this job lands on.

import numpy as np
import pandas as pd
import torch

import tensorqtl
from tensorqtl import pgen, cis, post
import tensorqtl.core as _tqtl_core

# Workaround for a tensorqtl packaging bug: core.py does a bare `import rfunc`
# inside a try/except when probing for R/rpy2/qvalue support. That bare import
# almost never succeeds from a normal site-packages install (rfunc.py is a
# package-relative module, not top-level), so core.has_rpy2 stays False and
# `rfunc` is never bound in core's namespace -- even when R, rpy2, and the
# qvalue package are all correctly installed and working (verified manually).
# tensorqtl.post.calculate_qvalues() then fails with NameError. Fix: import
# rfunc the correct way ourselves and inject it where tensorqtl expects it.
from tensorqtl import rfunc as _rfunc
_tqtl_core.rfunc = _rfunc
_tqtl_core.has_rpy2 = True
post.rfunc = _rfunc


# ── CLI args ──────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--plink-prefix",  required=True)
    p.add_argument("--bed",           required=True)
    p.add_argument("--covariates",    required=True)
    p.add_argument("--out-dir",       required=True)
    p.add_argument("--prefix",        required=True)
    p.add_argument("--cis-window",    type=int,   default=1_000_000)
    p.add_argument("--fdr",           type=float, default=0.05)
    p.add_argument("--seed",          type=int,   default=123456)
    p.add_argument("--perm-out",      required=True)
    p.add_argument("--nom-out",       required=True)
    return p.parse_args()

args = parse_args()

plink_prefix     = args.plink_prefix
expression_bed   = args.bed
covariates_file  = args.covariates
out_dir          = args.out_dir
prefix           = args.prefix
cis_window       = args.cis_window
fdr_threshold    = args.fdr
seed             = args.seed
perm_output_path = args.perm_out
nom_output_path  = args.nom_out


# ── GPU check ─────────────────────────────────────────────────────────────
print(f"PyTorch version: {torch.__version__}")
print(f"CUDA available:  {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU:             {torch.cuda.get_device_name(0)}")

# ── Output directory setup ────────────────────────────────────────────────
if os.path.exists(out_dir):
    print(f"Cleaning up old output directory: {out_dir}")
    shutil.rmtree(out_dir)
os.makedirs(out_dir)

# ── 1. Load data ──────────────────────────────────────────────────────────
print("\n── Loading phenotypes and covariates ──")
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
print(f"  Phenotypes loaded : {phenotype_df.shape[0]} DMRs × {phenotype_df.shape[1]} samples")

covariates_df = pd.read_csv(covariates_file, sep="\t", index_col=0).T
print(f"  Covariates loaded : {covariates_df.shape[1]} covariates × {covariates_df.shape[0]} samples")

# ── Align sample order to phenotype matrix ────────────────────────────────
if not covariates_df.index.equals(phenotype_df.columns):
    print("  WARNING: sample order mismatch — reindexing covariates to phenotype columns")
    print(f"    Phenotype samples : {phenotype_df.columns.tolist()}")
    print(f"    Covariate samples : {covariates_df.index.tolist()}")
    covariates_df = covariates_df.loc[phenotype_df.columns]

print(f"  Covariates loaded : {covariates_df.shape[1]} covariates × {covariates_df.shape[0]} samples")

print("\n── Loading genotypes ──")
pgr = pgen.PgenReader(plink_prefix)
genotype_df = pgr.load_genotypes()
variant_df  = pgr.variant_df
print(f"  Variants loaded   : {genotype_df.shape[0]}")

# ── 2. Stage 1: Permutation pass ──────────────────────────────────────────
print(f"\n── Stage 1: Permutation pass (window={cis_window:,} bp, seed={seed}) ──")
cis_df = cis.map_cis(
    genotype_df, variant_df, phenotype_df, phenotype_pos_df,
    covariates_df=covariates_df, window=cis_window, seed=seed,
)

print("Calculating q-values...")
post.calculate_qvalues(cis_df, fdr=fdr_threshold, qvalue_lambda=0.85)
cis_df.to_csv(perm_output_path, sep="\t")
print(f"Stage 1 results saved → {perm_output_path}")

n_sig = (cis_df["qval"] < fdr_threshold).sum()
print(f"  Total DMRs tested : {len(cis_df)}")
print(f"  Significant DMRs  : {n_sig} (FDR < {fdr_threshold})")

# ── 3. Stage 2: Nominal pass ──────────────────────────────────────────────
if n_sig == 0:
    print("\nNo DMRs passed FDR threshold — nominal pass skipped.")
    pd.DataFrame().to_csv(nom_output_path, sep="\t", index=False)
else:
    print(f"\n── Stage 2: Nominal pass ({n_sig} significant DMRs) ──")
    significant_dmrs     = cis_df[cis_df["qval"] < fdr_threshold].index
    phenotype_df_sig     = phenotype_df.loc[significant_dmrs]
    phenotype_pos_df_sig = phenotype_pos_df.loc[significant_dmrs]

    cis.map_nominal(
        genotype_df, variant_df, phenotype_df_sig, phenotype_pos_df_sig,
        prefix, covariates_df=covariates_df, window=cis_window, output_dir=out_dir,
    )

    parquet_files = sorted(glob.glob(f"{out_dir}/{prefix}.cis_qtl_pairs.*.parquet"))
    print(f"  Parquet files generated: {len(parquet_files)}")

    if not parquet_files:
        print("WARNING: No parquet files found.")
        pd.DataFrame().to_csv(nom_output_path, sep="\t", index=False)
    else:
        all_pairs = pd.concat([pd.read_parquet(f) for f in parquet_files])
        all_pairs.to_csv(nom_output_path, sep="\t", index=False)
        print(f"  Total SNP-DMR pairs : {len(all_pairs)}")
        print(f"Nominal results saved → {nom_output_path}")

        print("\n── Validation ──")
        test_pheno         = significant_dmrs[0]
        top_variant_stage1 = cis_df.loc[test_pheno, "variant_id"]
        match = all_pairs[
            (all_pairs["phenotype_id"] == test_pheno) &
            (all_pairs["variant_id"]   == top_variant_stage1)
        ]
        if not match.empty:
            print(f"[PASSED] '{top_variant_stage1}' found in nominal results.")
        else:
            print(f"[FAILED] '{top_variant_stage1}' MISSING — check chr naming or window size.")

print("\n── tensorQTL run complete ──")