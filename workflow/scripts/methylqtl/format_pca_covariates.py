"""
format_pca_covariates.py — reformat plink2 eigenvec → tensorQTL covariates TSV
"""

import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--eigenvec", required=True)
    p.add_argument("--n-pcs",    type=int, default=5)
    p.add_argument("--meta",     default="")
    p.add_argument("--metadata", default="",
                   help="sample_metadata.tsv from pipeline (sample_id + diagnosis)")
    p.add_argument("--out",      required=True)
    return p.parse_args()


def main():
    args = parse_args()

    pcs = pd.read_csv(args.eigenvec, sep="\t")
    pcs = pcs.rename(columns={"#IID": "sample_id"}).set_index("sample_id")
    pc_cols = [f"PC{i}" for i in range(1, args.n_pcs + 1)]
    pcs = pcs[pc_cols]

    # ── Validate against metadata TSV (no renaming needed — IDs already match) ──
    if args.metadata:
        meta_df = pd.read_csv(args.metadata, sep="\t", dtype=str)
        known_ids = set(meta_df["sample_id"].str.strip())
        unrecognized = [x for x in pcs.index if str(x) not in known_ids]
        if unrecognized:
            raise ValueError(
                f"Sample IDs in eigenvec not found in metadata TSV: {unrecognized}\n"
                f"Known IDs: {sorted(known_ids)}"
            )
        print(f"[PCA] All {len(pcs)} sample IDs validated against metadata TSV")

    cov = pcs.T

    if args.meta:
        meta = pd.read_csv(args.meta, index_col=0)
        meta = meta.loc[meta.index.isin(cov.columns)]
        cov = pd.concat([meta.T, cov])

    cov.to_csv(args.out, sep="\t")
    print(f"[PCA] Covariates written → {args.out}  (shape: {cov.shape})")


if __name__ == "__main__":
    main()