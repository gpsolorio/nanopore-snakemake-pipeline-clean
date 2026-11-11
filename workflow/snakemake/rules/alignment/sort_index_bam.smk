##############################################################################
# Sort + Index BAM
##############################################################################

rule sort_index_bam:
    input:
        sample("{barcode}", "bam", "{barcode}_merged_unsorted.bam")
    output:
        bam = sample("{barcode}", "bam", "{barcode}_merged_sorted.bam"),
        bai = sample("{barcode}", "bam", "{barcode}_merged_sorted.bam.bai")
    log:
        sample_log("{barcode}", "sort_index_bam.log")
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")
    threads:
        res("sort_index_bam", "threads")
    resources:
        mem_gb    = res("sort_index_bam", "mem_gb"),
        time      = res("sort_index_bam", "time"),
        account   = res("sort_index_bam", "account"),
        partition = res("sort_index_bam", "partition"),
        gpus      = res("sort_index_bam", "gpus", 0)
    shell:
        SHELL_PRELUDE + """
        samtools sort {input} -o {output.bam}
        samtools index {output.bam}
        """ + SHELL_EPILOG
