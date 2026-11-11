##############################################################################
# dmrseq DMR Calling
##############################################################################

_DMRSEQ_SCRIPT = (REPO_ROOT / "analysis" / "dmr" / "02_call_dmrs" / "run_dmrseq_pipeline.R").resolve()

DMR_RESULTS_CSV = cohort("dmrs", "dmrseq_results", "DMR_results.csv")


rule dmrseq_call_dmrs:
    input:
        beds = expand(
            sample("{barcode}", "methylation", "{barcode}_modkit.bed"),
            barcode=BARCODES
        )
    output:
        dmrs = cohort("dmrs", "dmrseq_results", "DMR_results.csv")
    log:
        cohort_log("dmrseq_call_dmrs.log")
    conda:
        str(ENVS_DIR / "dmrseq-env.yml")
    params:
        script       = str(_DMRSEQ_SCRIPT),
        repo_root    = str(REPO_ROOT),
        metadata     = SAMPLE_METADATA_TSV,
        group_col    = DMR_GROUP_COLUMN,
        group_vals   = ",".join(DMR_GROUP_VALUES),
        restrict_arg = (f"--restrict_role {DMR_RESTRICT_ROLE}" if DMR_RESTRICT_ROLE else ""),
    threads:
        res("dmrseq_call_dmrs", "threads")
    resources:
        mem_gb    = res("dmrseq_call_dmrs", "mem_gb"),
        time      = res("dmrseq_call_dmrs", "time"),
        account   = res("dmrseq_call_dmrs", "account"),
        partition = res("dmrseq_call_dmrs", "partition")
    shell:
        SHELL_PRELUDE + """
        Rscript {params.script} \\
            --outfile   {output.dmrs} \\
            --repo_root {params.repo_root} \\
            --metadata  {params.metadata} \\
            --group_column {params.group_col} \\
            --group_values "{params.group_vals}" \\
            {params.restrict_arg}
        """ + SHELL_EPILOG