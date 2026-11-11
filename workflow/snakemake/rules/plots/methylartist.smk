##############################################################################
# Reheader phased BAMs — add SM:<barcode> tag for methylartist sample tracking
#
# Clair3's docker branch omits --sample_name, so @RG headers have no SM: tag.
# methylartist uses the filename stem as the sample label when SM: is absent,
# causing all samples to collapse into one track when BAMs share a filename.
#
# This rule reheaders each BAM and names the output after the barcode so that
# both the SM: tag and the filename correctly identify the sample.
#
# Once call_variants is rerun with --sample_name fixed in the docker branch,
# this rule can be removed and plot_methylation pointed back at phased_output.bam.
##############################################################################

rule reheader_phased_bams:
    input:
        bam = sample("{barcode}", "variants", "phased_output.bam"),
        bai = sample("{barcode}", "variants", "phased_output.bam.bai"),
    output:
        bam = sample("{barcode}", "variants", "{barcode}.bam"),
        bai = sample("{barcode}", "variants", "{barcode}.bam.bai"),
    log:
        sample_log("{barcode}", "reheader_phased_bams.log")
    params:
        sample_name = lambda wc: wc.barcode
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")
    shell:
        SHELL_PRELUDE + """
        samtools addreplacerg \\
            -w \\
            -r "ID:{params.sample_name}\\tSM:{params.sample_name}" \\
            -m overwrite_all \\
            -O bam \\
            -o {output.bam} \\
            {input.bam}
        samtools index {output.bam}

        echo "[check] Verifying reheadered BAM integrity..."
        samtools quickcheck {output.bam}

        n_rg=$(samtools view -H {output.bam} | grep -c '^@RG' || true)
        if [ "$n_rg" -ne 1 ]; then
            echo "[FATAL] Expected exactly 1 @RG line in {output.bam}, found $n_rg" >&2
            exit 1
        fi
        echo "[check] OK: single @RG line (ID:{params.sample_name}, SM:{params.sample_name}) confirmed"
        """ + SHELL_EPILOG

##############################################################################
# X-Chromosome Gene Methylartist Plots
##############################################################################

_METHYLARTIST_BED    = _get_annotation("asd_genes")
_CPG_ISLAND_BED      = _get_annotation("cpg_islands")
_METHYLARTIST_SCRIPT = (
    str(SCRIPTS_DIR / "visualization" / "run_methylartist_mini.sh")
    if GENOME_MODE == "mini"
    else str(SCRIPTS_DIR / "visualization" / "run_methylartist.sh")
)


print(f"[genome] asd_genes   = {_METHYLARTIST_BED}")
print(f"[genome] cpg_islands = {_CPG_ISLAND_BED}")


rule plot_methylation:
    input:
        bam = expand(sample("{barcode}", "variants", "{barcode}.bam"),     barcode=BARCODES),
        bai = expand(sample("{barcode}", "variants", "{barcode}.bam.bai"), barcode=BARCODES),
    output:
        outdir = directory(cohort("summaries", "methylartist"))
    log:
        cohort_log("methylartist.log")
    params:
        bed_file  = str(_METHYLARTIST_BED),
        script    = str(_METHYLARTIST_SCRIPT),
        reference = str(REFERENCE_GENOME),
        bam_list  = lambda wc, input: ",".join(input.bam)
    threads:
        res("plot_methylation", "threads")
    resources:
        mem_gb    = res("plot_methylation", "mem_gb"),
        time      = res("plot_methylation", "time"),
        account   = res("plot_methylation", "account"),
        partition = res("plot_methylation", "partition"),
        gpus      = res("plot_methylation", "gpus", 0)
    conda:
        str(ENVS_DIR / "methylartist-env.yml")
    shell:
        SHELL_PRELUDE + """

        mkdir -p {output.outdir}
        
        bash {params.script} \\
            --bam {params.bam_list:q} \\
            --bed-file {params.bed_file} \\
            --reference {params.reference} \\
            --threads {threads} \\
            --output-dir {output.outdir}
        """ + SHELL_EPILOG