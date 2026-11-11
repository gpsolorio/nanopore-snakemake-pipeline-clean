#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from sklearn.mixture import GaussianMixture


def load_bed(bed_path):
    """
    Load BED file with target regions.

    Column 4 is "{GENE}_{transcript_accession}" (e.g. "PLCXD1_NM_001370371").
    Many rows share identical (chrom, start, end) and differ only in which
    transcript's promoter annotation produced that window -- those are
    true duplicates and get collapsed to one region. Other rows share a
    gene symbol but have genuinely different (start, end) -- e.g. PLCXD1
    has 3 distinct promoter windows from 3 transcripts with different
    TSS's -- and are kept as separate regions rather than merged, since
    merging would misrepresent three real ~3kb promoters as one
    inaccurate ~10kb span.

    Region name is set to "{gene}_{start}_{end}" so every plot title and
    output filename shows the gene symbol plus its exact coordinates --
    this also makes the name naturally unique per distinct interval
    without needing a separate dedup counter.
    """
    unique_intervals = {}  # (chrom, start, end) -> gene symbol
    with open(bed_path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 4:
                raise ValueError("BED must have at least 4 columns: chr start end name")
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            gene = row[3].split("_", 1)[0]
            unique_intervals[(chrom, start, end)] = gene

    if not unique_intervals:
        raise ValueError("No valid regions found in BED file")

    regions = [
        (chrom, start, end, f"{gene}_{start}_{end}")
        for (chrom, start, end), gene in unique_intervals.items()
    ]
    regions.sort(key=lambda r: (r[0], r[1]))
    return regions


def load_parental_map(parental_map_path):
    """Load optional parental origin mapping."""
    if parental_map_path is None:
        return {}

    df = pd.read_csv(parental_map_path)
    required = {"sample", "hp1_origin", "hp2_origin"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Parental map missing required columns: {sorted(missing)}")

    mapping = {}
    for _, row in df.iterrows():
        mapping[str(row["sample"])] = {
            1: str(row["hp1_origin"]),
            2: str(row["hp2_origin"]),
        }
    return mapping


def sample_name_from_bam_path(bam_path):
    """Extract sample name from BAM file path."""
    name = bam_path.name
    if name.endswith("_merged.bam"):
        return name[: -len("_merged.bam")]
    if name.endswith(".bam"):
        return name[: -len(".bam")]
    return name


def hp_label(sample, hp_value, parental_map):
    """Resolve haplotype label from sample and HP tag."""
    if hp_value in (1, 2):
        if sample in parental_map:
            return parental_map[sample].get(hp_value, "Unknown")
        return f"HP{hp_value}"
    return "Unknown"


def per_read_methylation(bam_path, regions, parental_map, min_mapq, min_cpg_per_read, mod_code="m"):
    """Extract methylation data with one row per read."""
    sample = sample_name_from_bam_path(bam_path)
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    results = []
    read_count = 0

    for chrom, start, end, region_name in regions:
        try:
            iterator = bam.fetch(chrom, start, end)
        except ValueError:
            continue

        for read in iterator:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            hp_value = read.get_tag("HP") if read.has_tag("HP") else None
            origin = hp_label(sample, hp_value, parental_map)

            modified = read.modified_bases
            if not modified:
                continue

            aligned_pairs = read.get_aligned_pairs(matches_only=True)
            if not aligned_pairs:
                continue
            qpos_to_rpos = dict(aligned_pairs)

            read_cpg_probs = []
            for key, vals in modified.items():
                base, _strand, this_mod_code = key
                if base != "C":
                    continue
                if this_mod_code != mod_code:
                    continue

                for qpos, qual in vals:
                    rpos = qpos_to_rpos.get(qpos)
                    if rpos is None:
                        continue
                    if start <= rpos < end:
                        read_cpg_probs.append(qual / 255.0 * 100.0)

            if len(read_cpg_probs) >= min_cpg_per_read:
                read_count += 1
                results.append(
                    {
                        "sample": sample,
                        "region": region_name,
                        "origin": origin,
                        "read_name": read.query_name,
                        "mean_methylation": sum(read_cpg_probs) / len(read_cpg_probs),
                        "n_cpg_calls": len(read_cpg_probs),
                        "mod_code": mod_code,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "read_start": read.reference_start,
                        "read_end": read.reference_end,
                        "mapq": read.mapping_quality,
                    }
                )

    bam.close()
    print(f"  Found {read_count} reads with modification '{mod_code}' data")
    return results


def fit_two_component_gmm(values, gmm_params):
    """Fit a 2-component Gaussian mixture model."""
    X = np.array(values).reshape(-1, 1)

    gmm = GaussianMixture(
        n_components=gmm_params["n_components"],
        covariance_type=gmm_params["covariance_type"],
        max_iter=gmm_params["max_iter"],
        random_state=gmm_params["random_state"],
    )
    gmm.fit(X)

    means = gmm.means_.flatten()
    means.sort()
    return means, gmm


def gmm_half_max_cutoffs(gmm, x_min=0.0, x_max=100.0):
    """Compute cutoff lines from Gaussian half-max boundaries."""
    means = gmm.means_.flatten()

    if gmm.covariance_type == "full":
        variances = np.array([c[0, 0] for c in gmm.covariances_])
    elif gmm.covariance_type == "diag":
        variances = np.array([c[0] for c in gmm.covariances_])
    elif gmm.covariance_type == "spherical":
        variances = np.array(gmm.covariances_)
    elif gmm.covariance_type == "tied":
        variances = np.array([gmm.covariances_[0, 0]] * len(means))
    else:
        raise ValueError(f"Unsupported covariance_type: {gmm.covariance_type}")

    stds = np.sqrt(np.maximum(variances, 1e-12))

    order = np.argsort(means)
    low_idx, high_idx = order[0], order[1]

    low_mean = means[low_idx]
    high_mean = means[high_idx]

    halfmax_scale = np.sqrt(2.0 * np.log(2.0))
    low_cut = low_mean + halfmax_scale * stds[low_idx]
    high_cut = high_mean - halfmax_scale * stds[high_idx]

    low_cut = float(np.clip(low_cut, x_min, x_max))
    high_cut = float(np.clip(high_cut, x_min, x_max))

    if low_cut > high_cut:
        mid = (low_cut + high_cut) / 2.0
        low_cut = mid
        high_cut = mid

    return {
        "low_mean": float(low_mean),
        "high_mean": float(high_mean),
        "low_cut": low_cut,
        "high_cut": high_cut,
    }


def calculate_section_percentages(values, low_cut, high_cut):
    """Calculate percentage of values in three sections."""
    values = np.array(values)
    total = len(values)

    low_count = np.sum(values < low_cut)
    mid_count = np.sum((values >= low_cut) & (values <= high_cut))
    high_count = np.sum(values > high_cut)

    low_pct = (low_count / total * 100) if total > 0 else 0
    mid_pct = (mid_count / total * 100) if total > 0 else 0
    high_pct = (high_count / total * 100) if total > 0 else 0

    return {
        "low_pct": low_pct,
        "mid_pct": mid_pct,
        "high_pct": high_pct,
        "low_count": low_count,
        "mid_count": mid_count,
        "high_count": high_count,
    }


def _safe_ratio(numerator, denominator):
    """NaN-safe division for ratio columns -- avoids inf/ZeroDivisionError
    when a bin has 0% of reads (common for Xq, which is often near-empty)."""
    if denominator is None or pd.isna(denominator) or denominator == 0:
        return np.nan
    if numerator is None or pd.isna(numerator):
        return np.nan
    return numerator / denominator


def compute_xci_stats(df, gmm_params):
    """
    One GMM fit per (sample, region), reused by both the plot renderer and
    the summary table. Previously this was fit twice per region (grid
    subplot loop + individual-plot loop) and its output -- the Xa/X?/Xi
    percentages -- was only ever drawn as plot annotation text and never
    saved anywhere. This is now the single source of truth for both.

    Also computes haplotype-resolved Xi/Xa when a region has reads from
    exactly two non-"Unknown" origins (HP1/HP2, or Maternal/Paternal if
    a parental map was supplied) -- the higher-methylation haplotype is
    treated as Xi, the lower as Xa. This is a read-count-independent
    signal (mean methylation per haplotype) that complements the
    bin-percentage ratios, which is what you want when phasing coverage
    is thin (like the BCORL1 example: 27+16 phased reads vs 197 Unknown --
    the haplotype ratio is still meaningful even though most reads can't
    be split into HP1/HP2 for the binned percentages).
    """
    rows = []
    for (sample_name, region), region_data in df.groupby(["sample", "region"]):
        values = region_data["mean_methylation"].values
        n_reads = len(values)

        row = {
            "sample": sample_name, "region": region, "n_reads": n_reads,
            "low_cut": np.nan, "high_cut": np.nan,
            "pct_Xa": np.nan, "pct_Xq": np.nan, "pct_Xi": np.nan,
            "n_Xa": np.nan, "n_Xq": np.nan, "n_Xi": np.nan,
        }

        if n_reads > gmm_params["min_reads_for_gmm"]:
            try:
                _, gmm = fit_two_component_gmm(values, gmm_params)
                cutoffs = gmm_half_max_cutoffs(gmm, x_min=0.0, x_max=100.0)
                pcts = calculate_section_percentages(values, cutoffs["low_cut"], cutoffs["high_cut"])
                row.update({
                    "low_cut": cutoffs["low_cut"], "high_cut": cutoffs["high_cut"],
                    "pct_Xa": pcts["low_pct"], "pct_Xq": pcts["mid_pct"], "pct_Xi": pcts["high_pct"],
                    "n_Xa": pcts["low_count"], "n_Xq": pcts["mid_count"], "n_Xi": pcts["high_count"],
                })
            except Exception as exc:
                print(f"Warning: GMM fitting failed for {sample_name}/{region}: {exc}")

        row["ratio_Xi_Xa"] = _safe_ratio(row["pct_Xi"], row["pct_Xa"])
        row["ratio_Xq_Xa"] = _safe_ratio(row["pct_Xq"], row["pct_Xa"])
        row["ratio_Xq_Xi"] = _safe_ratio(row["pct_Xq"], row["pct_Xi"])

        haplotype_origins = sorted(o for o in region_data["origin"].unique() if o != "Unknown")
        row["n_haplotypes_present"] = len(haplotype_origins)
        row["haplotype_higher"] = ""
        row["haplotype_lower"] = ""
        row["haplotype_Xi_mean_pct"] = np.nan
        row["haplotype_Xa_mean_pct"] = np.nan
        row["haplotype_Xi_Xa_ratio"] = np.nan

        if len(haplotype_origins) == 2:
            means_by_origin = {
                o: region_data.loc[region_data["origin"] == o, "mean_methylation"].mean()
                for o in haplotype_origins
            }
            hi = max(means_by_origin, key=means_by_origin.get)
            lo = min(means_by_origin, key=means_by_origin.get)
            row["haplotype_higher"] = hi
            row["haplotype_lower"] = lo
            row["haplotype_Xi_mean_pct"] = means_by_origin[hi]
            row["haplotype_Xa_mean_pct"] = means_by_origin[lo]
            row["haplotype_Xi_Xa_ratio"] = _safe_ratio(means_by_origin[hi], means_by_origin[lo])

        rows.append(row)

    return pd.DataFrame(rows)


def sanitize_region_name(region_name):
    """Create a filesystem-safe region identifier for output filenames."""
    safe = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in str(region_name))
    return safe.strip("_") or "region"


def title_to_filename(title):
    """Convert a plot title into a clean filename stem."""
    safe = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in str(title))
    safe = safe.replace("_5mC", "_5mc")
    return "_".join(part for part in safe.split("_") if part)


def make_histogram_plots(df, output_png, bins, separate_origins, viz, xci_stats, title_suffix=""):
    """Create both a combined grid plot and one histogram image per CpG island."""
    if df.empty:
        raise ValueError("No methylation data to plot")

    samples_in_df = df["sample"].unique()
    if len(samples_in_df) != 1:
        raise ValueError(
            "make_histogram_plots expects reads from exactly one sample -- "
            "mixing samples into one region's GMM fit would corrupt the "
            f"Xa/Xi classification. Got: {sorted(samples_in_df)}"
        )
    stats_lookup = xci_stats[xci_stats["sample"] == samples_in_df[0]].set_index("region")

    regions = sorted(df["region"].unique())
    n_regions = len(regions)
    n_cols = min(3, n_regions)
    n_rows = int(np.ceil(n_regions / n_cols))

    palette = {
        "Maternal": "#d55e00", "Paternal": "#0072b2", "Unknown": "#e6c200",
        "HP1": "#d55e00", "HP2": "#0072b2",
    }

    output_png.parent.mkdir(parents=True, exist_ok=True)
    written_paths = []

    def _draw(ax, region, region_data, plot_title):
        if separate_origins:
            for origin in sorted(region_data["origin"].unique()):
                origin_data = region_data[region_data["origin"] == origin]
                ax.hist(origin_data["mean_methylation"], bins=bins, alpha=0.6,
                    label=f"{origin} (n={len(origin_data)})",
                    color=palette.get(origin, "#999999"), edgecolor="black",
                    linewidth=viz["hist_edgewidth"])
            ax.legend(fontsize=viz["legend_fontsize"])
        else:
            ax.hist(region_data["mean_methylation"], bins=bins, alpha=viz["hist_alpha"],
                color=viz["hist_color"], edgecolor="black", linewidth=viz["hist_edgewidth"])
            ax.text(0.95, 0.95, f"n={len(region_data)}", transform=ax.transAxes,
                fontsize=viz["annotation_fontsize"], va="top", ha="right",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

        ax.set_xlabel("Mean methylation per read (%)", fontsize=viz["label_fontsize"])
        ax.set_ylabel("Number of reads", fontsize=viz["label_fontsize"])
        ax.set_title(plot_title, fontsize=viz["title_fontsize"], fontweight="bold")
        ax.set_xlim(0, 100)
        ax.tick_params(axis="both", labelsize=viz["tick_fontsize"])
        ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)

        median_val = region_data["mean_methylation"].median()
        ax.axvline(median_val, color="red", linestyle="--", linewidth=viz["median_linewidth"], alpha=viz["median_alpha"])

        stats = stats_lookup.loc[region]
        if pd.notna(stats["low_cut"]):
            low_cut, high_cut = stats["low_cut"], stats["high_cut"]
            ax.axvline(low_cut, color="green", linestyle="--", linewidth=viz["gmm_linewidth"], alpha=viz["gmm_alpha"])
            ax.axvline(high_cut, color="purple", linestyle="--", linewidth=viz["gmm_linewidth"], alpha=viz["gmm_alpha"])

            y_max = ax.get_ylim()[1]
            xa_x = max(low_cut * 0.5, 2.0)
            xq_x = (low_cut + high_cut) / 2.0
            xi_x = high_cut + (100.0 - high_cut) * 0.5

            for x, pct, label, color in (
                (xa_x, stats["pct_Xa"], "Xa", "lightblue"),
                (xq_x, stats["pct_Xq"], "X?", "lightyellow"),
                (xi_x, stats["pct_Xi"], "Xi", "lightgreen"),
            ):
                ax.text(x, y_max * (0.35 if label == "X?" else 0.5), f"{label}: {pct:.1f}%",
                    fontsize=viz["annotation_fontsize"], ha="center", va="center",
                    bbox=dict(boxstyle="round", facecolor=color, alpha=0.6))

    # 1) Combined grid plot
    grid_fig, grid_axes = plt.subplots(n_rows, n_cols,
        figsize=(viz["subplot_width"] * n_cols, viz["subplot_height"] * n_rows))
    grid_axes = np.array(grid_axes).reshape(-1)

    for idx, region in enumerate(regions):
        _draw(grid_axes[idx], region, df[df["region"] == region], f"{region}{title_suffix}")
    for idx in range(n_regions, len(grid_axes)):
        grid_axes[idx].axis("off")

    grid_fig.tight_layout()
    grid_fig.savefig(output_png, dpi=viz["dpi"], bbox_inches="tight")
    plt.close(grid_fig)
    written_paths.append(output_png)

    # 2) Individual region plots
    for region in regions:
        fig, ax = plt.subplots(figsize=(viz["subplot_width"], viz["subplot_height"]))
        plot_title = f"{region}{title_suffix}"
        _draw(ax, region, df[df["region"] == region], plot_title)

        region_output = output_png.parent / f"{title_to_filename(plot_title)}{output_png.suffix}"
        plt.savefig(region_output, dpi=viz["dpi"], bbox_inches="tight")
        plt.close(fig)
        written_paths.append(region_output)

    for path in written_paths:
        print(f"Wrote histogram plot: {path}")


def make_summary_stats(df):
    """Generate summary statistics for each region."""
    summary = []
    for region in sorted(df["region"].unique()):
        region_data = df[df["region"] == region]
        stats = {
            "region": region,
            "n_reads": len(region_data),
            "mean_methylation": region_data["mean_methylation"].mean(),
            "median_methylation": region_data["mean_methylation"].median(),
            "std_methylation": region_data["mean_methylation"].std(),
            "min_methylation": region_data["mean_methylation"].min(),
            "max_methylation": region_data["mean_methylation"].max(),
        }
        summary.append(stats)
    return pd.DataFrame(summary)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate CpG island methylation histogram analysis plots."
    )

    parser.add_argument(
        "--bam-dir", required=True,
        help="Directory containing BAM files (e.g. resources/... or workflow/data/bam)",
    )
    parser.add_argument(
        "--bed-file", required=True,
        help="CpG island BED file, e.g. resources/beds/{mini,hg38}_CpG_island_coords.bed",
    )
    parser.add_argument("--parental-map-file", default=None)
    parser.add_argument(
        "--output-prefix", required=True,
        help="Prefix for output plot files, e.g. results/MY_STUDY/cohort/plots/cpg_hist",
    )

    parser.add_argument("--bam-pattern", default="*_merged.bam")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-cpg-per-read", type=int, default=1)
    parser.add_argument("--bins", type=int, default=20)
    parser.add_argument("--separate-origins", action="store_true")

    parser.add_argument("--subplot-width", type=float, default=8.5)
    parser.add_argument("--subplot-height", type=float, default=6.0)
    parser.add_argument("--title-fontsize", type=int, default=18)
    parser.add_argument("--label-fontsize", type=int, default=18)
    parser.add_argument("--tick-fontsize", type=int, default=14)
    parser.add_argument("--annotation-fontsize", type=int, default=14)
    parser.add_argument("--legend-fontsize", type=int, default=14)

    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main():
    args = parse_args()

    viz_settings = {
        "subplot_width": args.subplot_width,
        "subplot_height": args.subplot_height,
        "title_fontsize": args.title_fontsize,
        "label_fontsize": args.label_fontsize,
        "tick_fontsize": args.tick_fontsize,
        "annotation_fontsize": args.annotation_fontsize,
        "legend_fontsize": args.legend_fontsize,
        "median_linewidth": 1.5,
        "gmm_linewidth": 2,
        "median_alpha": 0.7,
        "gmm_alpha": 0.8,
        "hist_alpha": 0.7,
        "hist_edgewidth": 0.5,
        "hist_color": "#0072b2",
        "dpi": args.dpi,
    }

    gmm_settings = {
        "n_components": 2,
        "covariance_type": "full",
        "max_iter": 100,
        "random_state": 42,
        "min_reads_for_gmm": 20,
    }

    bam_dir_path = Path(args.bam_dir)
    bed_path = Path(args.bed_file)
    output_prefix_path = Path(args.output_prefix)

    if not bam_dir_path.exists():
        raise FileNotFoundError(f"BAM directory not found: {bam_dir_path}")
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    bams = sorted(bam_dir_path.glob(args.bam_pattern))
    if not bams:
        raise FileNotFoundError(
            f"No BAM files found in {bam_dir_path} matching pattern '{args.bam_pattern}'"
        )

    regions = load_bed(bed_path)
    parental_map = load_parental_map(args.parental_map_file)

    print(f"Found {len(bams)} BAM files")
    print(f"Found {len(regions)} target regions")

    all_rows_5mc = []
    all_rows_5hmc = []

    for bam_path in bams:
        print(f"Processing {bam_path.name}")

        rows_5mc = per_read_methylation(
            bam_path=bam_path,
            regions=regions,
            parental_map=parental_map,
            min_mapq=args.min_mapq,
            min_cpg_per_read=args.min_cpg_per_read,
            mod_code="m",
        )
        all_rows_5mc.extend(rows_5mc)

        rows_5hmc = per_read_methylation(
            bam_path=bam_path,
            regions=regions,
            parental_map=parental_map,
            min_mapq=args.min_mapq,
            min_cpg_per_read=args.min_cpg_per_read,
            mod_code="h",
        )
        all_rows_5hmc.extend(rows_5hmc)

    df_5mc = pd.DataFrame(all_rows_5mc)
    df_5hmc = pd.DataFrame(all_rows_5hmc)

    print(f"Total 5mC reads extracted: {len(df_5mc)}")
    print(f"Total 5hmC reads extracted: {len(df_5hmc)}")

    if not df_5mc.empty:
        mc_dir = output_prefix_path / "5mC"
        tsvs_dir_5mc = mc_dir / "tsvs"
        tsvs_dir_5mc.mkdir(parents=True, exist_ok=True)

        output_tsv_5mc = tsvs_dir_5mc / "histogram_5mC.tsv"
        df_5mc.to_csv(output_tsv_5mc, sep="\t", index=False)
        print(f"Wrote 5mC table: {output_tsv_5mc}")

        summary_5mc = make_summary_stats(df_5mc)
        summary_tsv_5mc = tsvs_dir_5mc / "histogram_5mC_summary.tsv"
        summary_5mc.to_csv(summary_tsv_5mc, sep="\t", index=False)
        print(f"Wrote 5mC summary: {summary_tsv_5mc}")

        xci_stats_5mc = compute_xci_stats(df_5mc, gmm_settings)
        xci_summary_tsv_5mc = tsvs_dir_5mc / "xci_summary_5mC.tsv"
        xci_stats_5mc.to_csv(xci_summary_tsv_5mc, sep="\t", index=False)
        print(f"Wrote 5mC XCI summary: {xci_summary_tsv_5mc}")

        output_png_5mc = mc_dir / "histogram_5mC.png"
        make_histogram_plots(
            df=df_5mc, output_png=output_png_5mc, bins=args.bins,
            separate_origins=args.separate_origins, viz=viz_settings,
            xci_stats=xci_stats_5mc, title_suffix=" (5mC)",
        )

    if not df_5hmc.empty:
        hmc_dir = output_prefix_path / "5hmC"
        tsvs_dir_5hmc = hmc_dir / "tsvs"
        tsvs_dir_5hmc.mkdir(parents=True, exist_ok=True)

        output_tsv_5hmc = tsvs_dir_5hmc / "histogram_5hmC.tsv"
        df_5hmc.to_csv(output_tsv_5hmc, sep="\t", index=False)
        print(f"Wrote 5hmC table: {output_tsv_5hmc}")

        summary_5hmc = make_summary_stats(df_5hmc)
        summary_tsv_5hmc = tsvs_dir_5hmc / "histogram_5hmC_summary.tsv"
        summary_5hmc.to_csv(summary_tsv_5hmc, sep="\t", index=False)
        print(f"Wrote 5hmC summary: {summary_tsv_5hmc}")

        viz_5hmc = dict(viz_settings)
        viz_5hmc["hist_color"] = "#e69f00"

        output_png_5hmc = hmc_dir / "histogram_5hmC.png"
        
        xci_stats_5hmc = compute_xci_stats(df_5hmc, gmm_settings)
        xci_summary_tsv_5hmc = tsvs_dir_5hmc / "xci_summary_5hmC.tsv"
        xci_stats_5hmc.to_csv(xci_summary_tsv_5hmc, sep="\t", index=False)
        print(f"Wrote 5hmC XCI summary: {xci_summary_tsv_5hmc}")

        output_png_5hmc = hmc_dir / "histogram_5hmC.png"
        make_histogram_plots(
            df=df_5hmc,
            output_png=output_png_5hmc,
            bins=args.bins,
            separate_origins=args.separate_origins,
            viz=viz_5hmc,
            xci_stats=xci_stats_5hmc,
            title_suffix=" (5hmC)",
        )

if __name__ == "__main__":
    main()