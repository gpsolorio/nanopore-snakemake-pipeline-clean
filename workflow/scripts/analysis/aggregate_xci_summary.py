#!/usr/bin/env python3
"""
Combine per-sample XCI summary TSVs into one population-level workbook,
flag abnormal Xi/Xa (and X?/Xa, X?/Xi) ratios per gene, and plot each
gene's ratio distribution across the cohort.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RATIO_COLS = ["ratio_Xi_Xa", "ratio_Xq_Xa", "ratio_Xq_Xi", "haplotype_Xi_Xa_ratio"]


def flag_outliers_iqr(series, k=1.5):
    """
    Tukey-fence outlier flag. Used instead of a mean/std z-score because
    the cohort here is small (pilot-scale) and these ratios are naturally
    right-skewed (Xa dominates most genes, so Xi/Xa is usually small with
    an occasional large value) -- IQR fences don't assume normality and
    aren't dragged around by the very outliers you're trying to flag.
    """
    q1, q3 = series.quantile(0.25), series.quantile(0.75)
    iqr = q3 - q1
    if iqr == 0 or pd.isna(iqr):
        return pd.Series(False, index=series.index)
    lower, upper = q1 - k * iqr, q3 + k * iqr
    return (series < lower) | (series > upper)


def plot_population_ratio(df_gene, ratio_col, gene, out_path, diagnosis_col):
    df_gene = df_gene.dropna(subset=[ratio_col])
    if df_gene.empty:
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    diagnoses = sorted(d for d in df_gene[diagnosis_col].dropna().unique())
    x_positions = {dx: i for i, dx in enumerate(diagnoses) or ["(none)"]}
    colors = {"ASD": "#d55e00", "TD": "#0072b2"}
    rng = np.random.default_rng(42)

    for dx in diagnoses:
        sub = df_gene[df_gene[diagnosis_col] == dx]
        jitter = rng.uniform(-0.12, 0.12, size=len(sub))
        x = np.full(len(sub), x_positions[dx]) + jitter
        flag_col = f"flag_{ratio_col}"
        edgecolors = np.where(sub[flag_col], "red", "black")
        linewidths = np.where(sub[flag_col], 2.0, 0.5)
        ax.scatter(x, sub[ratio_col], c=colors.get(dx, "#999999"),
            edgecolor=edgecolors, linewidth=linewidths, s=70, alpha=0.85, zorder=3)
        for _, row in sub[sub[flag_col]].iterrows():
            ax.annotate(row["sample"], (x_positions[dx], row[ratio_col]),
                fontsize=8, xytext=(6, 0), textcoords="offset points")

    ax.set_xticks(list(x_positions.values()))
    ax.set_xticklabels(list(x_positions.keys()))
    ax.set_ylabel(ratio_col.replace("_", " "))
    ax.set_title(f"{gene}: {ratio_col.replace('_', ' ')}")
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate per-sample XCI summaries into a cohort workbook.")
    p.add_argument("--xci-summary-tsvs", nargs="+", required=True)
    p.add_argument("--sample-metadata", required=True,
        help="TSV with sample_id / <diagnosis_col> / role / bed_path (common.smk's SAMPLE_METADATA_TSV)")
    p.add_argument("--diagnosis-col", default="Diagnosis")
    p.add_argument("--output-xlsx", required=True)
    p.add_argument("--output-plots-dir", required=True)
    p.add_argument("--iqr-k", type=float, default=1.5)
    return p.parse_args()


def main():
    args = parse_args()

    dfs = [pd.read_csv(p, sep="\t") for p in args.xci_summary_tsvs]
    df = pd.concat(dfs, ignore_index=True)

    meta = pd.read_csv(args.sample_metadata, sep="\t")
    meta = meta.rename(columns={"sample_id": "sample"})
    meta["sample"] = meta["sample"].astype(str)
    df["sample"] = df["sample"].astype(str)
    df = df.merge(meta[["sample", args.diagnosis_col, "role"]], on="sample", how="left")

    # Flag abnormal ratios per gene (a sample's ratio is only "abnormal"
    # relative to its own gene's distribution -- different genes have very
    # different baseline Xi/Xa balance, so pooling across genes first would
    # flag genuinely normal genes just for being globally unusual).
    for ratio_col in RATIO_COLS:
        flag_col = f"flag_{ratio_col}"
        df[flag_col] = False
        for region, idx in df.groupby("region").groups.items():
            df.loc[idx, flag_col] = flag_outliers_iqr(df.loc[idx, ratio_col], k=args.iqr_k)

    Path(args.output_xlsx).parent.mkdir(parents=True, exist_ok=True)
    plots_dir = Path(args.output_plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    flag_cols = [f"flag_{c}" for c in RATIO_COLS]
    flagged = df[df[flag_cols].any(axis=1)].sort_values(["region", "sample"])

    pop_stats = (
        df.groupby("region")[RATIO_COLS]
        .agg(["median", lambda s: s.quantile(0.25), lambda s: s.quantile(0.75)])
    )
    pop_stats.columns = [f"{c}_{stat}" for c, stat in
        zip([c for c in RATIO_COLS for _ in range(3)], ["median", "q1", "q3"] * len(RATIO_COLS))]
    pop_stats = pop_stats.reset_index()

    with pd.ExcelWriter(args.output_xlsx, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="per_sample_per_gene", index=False)
        flagged.to_excel(writer, sheet_name="flagged_abnormal", index=False)
        pop_stats.to_excel(writer, sheet_name="population_stats", index=False)
    print(f"Wrote workbook: {args.output_xlsx}")
    print(f"  {len(df)} sample x gene rows, {len(flagged)} flagged as abnormal on >=1 ratio")

    for region, df_gene in df.groupby("region"):
        for ratio_col in RATIO_COLS:
            if df_gene[ratio_col].notna().sum() < 2:
                continue
            out_path = plots_dir / f"{region}_{ratio_col}.png"
            plot_population_ratio(df_gene, ratio_col, region, out_path, args.diagnosis_col)

    print(f"Wrote population plots to: {plots_dir}")


if __name__ == "__main__":
    main()