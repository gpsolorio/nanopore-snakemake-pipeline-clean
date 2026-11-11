##############################################################################
# DMR Heatmap
##############################################################################

DMR_HEATMAP_SCRIPT = (REPO_ROOT / "analysis" / "dmr" / "05_visualize" / "dmr_heatmap.R").resolve()


rule plot_dmr_heatmap:
    input:
        bsseq_file = rules.wgcna_modules.output.bsseq_rds,
        dmrs_file  = _DMRS_CSV
    output:
        heatmap = cohort("summaries", "dmr_heatmap", "heatmap.pdf")
    log:
        cohort_log("plot_dmr_heatmap.log")
    params:
        test_covariate = DMR_GROUP_COLUMN,
        script = str(DMR_HEATMAP_SCRIPT)
    threads:
        res("plot_dmr_heatmap", "threads")
    resources:
        mem_gb    = res("plot_dmr_heatmap", "mem_gb"),
        time      = res("plot_dmr_heatmap", "time"),
        account   = res("plot_dmr_heatmap", "account"),
        partition = res("plot_dmr_heatmap", "partition")
    conda:
        str(ENVS_DIR / "pheatmap-env.yml")
    shell:
        SHELL_PRELUDE + """
        mkdir -p $(dirname {output.heatmap})

        "$CONDA_PREFIX/bin/Rscript" {params.script} \
            --bsseq-file {input.bsseq_file} \
            --dmrs-file {input.dmrs_file} \
            --test-covariate {params.test_covariate} \
            --output-file {output.heatmap}
        """ + SHELL_EPILOG
