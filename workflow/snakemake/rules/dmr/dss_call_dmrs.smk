##############################################################################
# DSS DMR Calling — DML test + callDMR across all barcodes
#
# Input:  all per-barcode modkit BEDs (collected via expand)
# Output: DMR_results.csv  — group-level DMR calls (chr, start, end, stats)
#
# This rule is not wildcarded — it aggregates all barcodes in one R job.
# The R script reads dataset config from analysis/dmr/config/ directly.
#
# Grouping is config-driven (see common.smk: DMR_GROUP_COLUMN /
# DMR_GROUP_VALUES, sourced from config's dmr.group_column / dmr.group_values)
# rather than hardcoded to the private cohort's diagnosis/ASD/TD -- this lets the same
# rule run against any cohort (diagnosis, sex, family, batch, etc.).
##############################################################################

_DSS_SCRIPT = (REPO_ROOT / "analysis" / "dmr" / "02_call_dmrs" / "run_dss.R").resolve()

DMR_RESULTS_CSV = cohort("dmrs", "dss_results", "DMR_results.csv")


rule dss_call_dmrs:
    input:
        beds = expand(
            sample("{barcode}", "methylation", "{barcode}_modkit.bed"),
            barcode=BARCODES
        )
    output:
        dmrs = cohort("dmrs", "dss_results", "DMR_results.csv")
    log:
        cohort_log("dss_call_dmrs.log")
    conda:
        str(ENVS_DIR / "dss-env.yml")
    params:
        script       = str(_DSS_SCRIPT),
        repo_root    = str(REPO_ROOT),
        metadata     = SAMPLE_METADATA_TSV,
        group_col    = DMR_GROUP_COLUMN,
        group_vals   = ",".join(DMR_GROUP_VALUES),
        restrict_arg = (f"--restrict_role {DMR_RESTRICT_ROLE}" if DMR_RESTRICT_ROLE else ""),
    threads:
        res("dss_call_dmrs", "threads")
    resources:
        mem_gb    = res("dss_call_dmrs", "mem_gb"),
        time      = res("dss_call_dmrs", "time"),
        account   = res("dss_call_dmrs", "account"),
        partition = res("dss_call_dmrs", "partition")
    shell:
        SHELL_PRELUDE + """
        echo "Input BEDs:  {input.beds}"
        echo "Output DMRs: {output.dmrs}"
        echo "Group column: {params.group_col}  |  Group values: {params.group_vals}"
        echo "Role restriction: {params.restrict_arg}"
        echo ""

        mkdir -p $(dirname {output.dmrs})
        
        "$CONDA_PREFIX/bin/Rscript" {params.script} \\
            --outfile   {output.dmrs} \\
            --repo_root {params.repo_root} \\
            --metadata  {params.metadata} \\
            --group_column {params.group_col} \\
            --group_values "{params.group_vals}" \\
            {params.restrict_arg}
        """ + SHELL_EPILOG
