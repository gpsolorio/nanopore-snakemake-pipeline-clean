##############################################################################
# WGCNA Co-methylation Module Detection
##############################################################################

WGCNA_SCRIPT = (REPO_ROOT / "analysis" / "dmr" / "03_integrate_signals" / "wgcna_methylation.R").resolve()

# Hardcodes the dss_call_dmrs output path. If dmr.caller == "dmrseq" in
# this config, dmrseq_call_dmrs.smk's actual output path/filename needs
# checking -- it may not be "dmrs/dss_results/DMR_results.csv". Fix this
# conditional once you confirm dmrseq_call_dmrs.smk's real output.
_DMRS_CSV = (
    cohort("dmrs", "dss_results", "DMR_results.csv")
    if DMR_CALLER == "dss"
    else cohort("dmrs", "dmrseq_results", "DMR_results.csv")  # UNVERIFIED path
)


rule wgcna_modules:
    input:
        metadata  = SAMPLE_METADATA_TSV,
        dmrs_file = _DMRS_CSV
    output:
        modules_rds = cohort("dmrs", "wgcna", "DMR_Modules.rds"),
        bsseq_rds   = cohort("dmrs", "wgcna", "bsseq_smoothed.rds"),
        modules_tsv = cohort("dmrs", "wgcna", "wgcna_modules.tsv")
    log:
        cohort_log("wgcna_modules.log")
    params:
        repo_root    = str(REPO_ROOT),
        outdir       = cohort("dmrs", "wgcna"),
        group_col    = DMR_GROUP_COLUMN,
        group_labels = ",".join(DMR_GROUP_VALUES),
        min_cluster_size = (
            f"--min_cluster_size {config['dmr']['wgcna_min_cluster_size']}"
            if config.get("dmr", {}).get("wgcna_min_cluster_size")
            else ""
        ),
        script = str(WGCNA_SCRIPT)
    threads:
        res("wgcna_modules", "threads")
    resources:
        mem_gb    = res("wgcna_modules", "mem_gb"),
        time      = res("wgcna_modules", "time"),
        account   = res("wgcna_modules", "account"),
        partition = res("wgcna_modules", "partition")
    conda:
        str(ENVS_DIR / "wgcna-env.yml")
    shell:
        SHELL_PRELUDE + """
        mkdir -p {params.outdir}

        "$CONDA_PREFIX/bin/Rscript" {params.script} \
            --repo_root {params.repo_root} \
            --metadata {input.metadata} \
            --dmr_file {input.dmrs_file} \
            --outdir {params.outdir} \
            --group_col {params.group_col} \
            --group_labels {params.group_labels} \
            {params.min_cluster_size}
        """ + SHELL_EPILOG
