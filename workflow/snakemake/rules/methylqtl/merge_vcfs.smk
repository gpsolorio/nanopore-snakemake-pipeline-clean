##############################################################################
# Merge per-barcode phased VCFs → one multi-sample VCF
#
# Input:  one phased_merge_output.vcf.gz per barcode (from call_variants)
# Output: all_samples_merged.vcf.gz
##############################################################################

def _merge_vcfs_input_vcf(barcode):
    if CLAIR3_CFG.get("enable_phasing", True):
        return sample(barcode, "variants", "phased_merge_output.vcf.gz")
    return sample(barcode, "variants", "merge_output.vcf.gz")


rule merge_vcfs:
    input:
        vcfs = [_merge_vcfs_input_vcf(b) for b in BARCODES],
    output:
        vcf_list = cohort("methylqtl", "vcfs", "gvcf_list.txt"),
        vcf      = cohort("methylqtl", "vcfs", "all_samples_merged.vcf.gz"),
        tbi      = cohort("methylqtl", "vcfs", "all_samples_merged.vcf.gz.tbi")
    log:
        cohort_log("merge_vcfs.log")
    params:
        ref        = str(REFERENCE_GENOME),
        chr_rename = str(SCRIPTS_DIR / "utils" / "chr_name_conv.txt"),
    threads:
        res("merge_vcfs", "threads")
    resources:
        mem_gb    = res("merge_vcfs", "mem_gb"),
        time      = res("merge_vcfs", "time"),
        account   = res("merge_vcfs", "account"),
        partition = res("merge_vcfs", "partition"),
    conda:
        str(ENVS_DIR / "bcftools-env.yml")
    shell:
        SHELL_PRELUDE + """
        printf '%s\\n' {input.vcfs} > {output.vcf_list}
        
        bcftools merge \\
            --gvcf {params.ref} \\
            -l {output.vcf_list} \\
            -Oz -o {output.vcf_list}.unrenamed.vcf.gz
        
        bcftools annotate \\
            --rename-chrs {params.chr_rename} \\
            -Oz -o {output.vcf} \\
            {output.vcf_list}.unrenamed.vcf.gz
        
        bcftools index -t {output.vcf}
        rm {output.vcf_list}.unrenamed.vcf.gz
        """ + SHELL_EPILOG