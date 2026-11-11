##############################################################################
# CpG Island Methylation Histogram Plots
##############################################################################

_CPG_HIST_SCRIPT = (WORKFLOW_ROOT / "scripts" / "visualization" / "histogram-per-CpG-island-plot.py").resolve()


rule plot_cpg_hists:
    input:
        bam = sample("{barcode}", "variants", "{barcode}.bam"),
        bai = sample("{barcode}", "variants", "{barcode}.bam.bai"),
    output:
        outdir = directory(sample("{barcode}", "histograms"))
    log:
        sample_log("{barcode}", "cpg_hists.log")
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")
    params:
        bam_dir     = lambda wc, input: str(Path(input.bam).parent),
        bam_pattern = "{barcode}.bam",
        # bed_file    = str(_CPG_ISLAND_BED),
        bed_file    = str(_METHYLARTIST_BED),   # was: str(_CPG_ISLAND_BED)
        script      = str(_CPG_HIST_SCRIPT),
    threads:
        res("plot_cpg_hists", "threads")
    resources:
        mem_gb    = res("plot_cpg_hists", "mem_gb"),
        time      = res("plot_cpg_hists", "time"),
        account   = res("plot_cpg_hists", "account"),
        partition = res("plot_cpg_hists", "partition"),
        gpus      = res("plot_cpg_hists", "gpus", 0)
    shell:
        SHELL_PRELUDE + """

        mkdir -p {output.outdir}

        python {params.script} \\
            --bam-dir {params.bam_dir} \\
            --bam-pattern {params.bam_pattern:q} \\
            --bed-file {params.bed_file} \\
            --separate-origins \\
            --output-prefix {output.outdir}
        """ + SHELL_EPILOG