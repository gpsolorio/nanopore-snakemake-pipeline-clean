##############################################################################
# DSS BSmooth + Individual Beta Extraction
#
# Input:  DMR_results.csv from dss_call_dmrs
#         all per-barcode modkit BEDs (to rebuild BSseq object)
# Output: DMR_individual_smoothed_methylation.bed.gz + .tbi — per-sample
#         smoothed beta values per DMR, formatted for tensorQTL
##############################################################################

_BSMOOTH_SCRIPT = (REPO_ROOT / "analysis" / "dmr" / "03_individual_betas" / "run_bsmooth.R").resolve()


rule dss_smooth_betas:
    input:
        dmrs = cohort("dmrs", "dss_results", "DMR_results.csv"),
        beds = expand(
            sample("{barcode}", "methylation", "{barcode}_modkit.bed"),
            barcode=BARCODES
        ),
    output:
        txt     = cohort("dmrs", "smoothed", "DMR_individual_smoothed_methylation.txt"),
        bed     = cohort("dmrs", "smoothed", "DMR_individual_smoothed_methylation.bed"),
        bed_gz  = cohort("dmrs", "smoothed", "DMR_individual_smoothed_methylation.bed.gz"),
        bed_tbi = cohort("dmrs", "smoothed", "DMR_individual_smoothed_methylation.bed.gz.tbi"),
    log:
        cohort_log("dss_smooth_betas.log")
    conda:
        str(ENVS_DIR / "dss-env.yml")
    params:
        script           = str(_BSMOOTH_SCRIPT),
        repo_root        = str(REPO_ROOT),
        dataset          = config.get("experiment", {}).get("dataset", "full"),
        metadata         = SAMPLE_METADATA_TSV,
        group_col        = config.get("experiment", {}).get("group_col", "diagnosis"),
        group_labels_arg = (
            f"--group_labels {config['experiment']['group_labels']}"
            if config.get("experiment", {}).get("group_labels", "")
            else ""
        ),
    threads:
        res("dss_smooth_betas", "threads")
    resources:
        mem_gb    = res("dss_smooth_betas", "mem_gb"),
        time      = res("dss_smooth_betas", "time"),
        account   = res("dss_smooth_betas", "account"),
        partition = res("dss_smooth_betas", "partition")
    shell:
        SHELL_PRELUDE + """

        mkdir -p $(dirname {output.txt})
        
        "$CONDA_PREFIX/bin/Rscript" {params.script} \\
            --dmr_file     {input.dmrs} \\
            --out_txt      {output.txt} \\
            --out_bed      {output.bed} \\
            --threads      {threads} \\
            --repo_root    {params.repo_root} \\
            --metadata     {params.metadata} \\
            --group_col    {params.group_col} \\
            {params.group_labels_arg}

        echo ""
        echo "Sorting, bgzipping, and indexing BED for tensorQTL..."
        
        head -1 {output.bed} > {output.bed}.sorted.tmp
        tail -n +2 {output.bed} | sort -k1,1 -k2,2n >> {output.bed}.sorted.tmp
        mv {output.bed}.sorted.tmp {output.bed}
        
        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed {output.bed_gz}
        """ + SHELL_EPILOG
