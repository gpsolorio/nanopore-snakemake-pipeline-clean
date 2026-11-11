##############################################################################
# Rename VCF sample IDs: strip everything after the first underscore so that
# e.g. "barcode01_JLYZ0004_AS_5mCG-5hmCG_sorted_NHIP_geno" → "barcode01"
#
# OPTIONAL — only runs if the merged VCF sample names contain underscores
# (i.e. are longer than bare barcode IDs). Controlled by:
#   tensorqtl.rename_samples: true   # default: false
#
# This aligns VCF sample names with the barcode IDs used everywhere else
# in the pipeline without needing a hand-crafted sample_map.txt file.
##############################################################################

if config["tensorqtl"].get("rename_samples", False):

    rule reheader_vcf:
        input:
            vcf = cohort("methylqtl", "vcfs", "all_samples_merged.vcf.gz"),
        output:
            vcf = cohort("methylqtl", "vcfs", "all_samples_final.vcf.gz"),
            tbi = cohort("methylqtl", "vcfs", "all_samples_final.vcf.gz.tbi"),
        log:
            cohort_log("reheader_vcf.log")
        threads:
            res("reheader_vcf", "threads")
        resources:
            mem_gb    = res("reheader_vcf", "mem_gb"),
            time      = res("reheader_vcf", "time"),
            account   = res("reheader_vcf", "account"),
            partition = res("reheader_vcf", "partition"),
        conda:
            str(ENVS_DIR / "bcftools-env.yml")
        shell:
            SHELL_PRELUDE + """
            bcftools query -l {input.vcf} \\
                | sed 's/_.*//' \\
                > {output.vcf}.sample_map.txt
            
            echo "Sample name mapping:"
            
            paste \\
                <(bcftools query -l {input.vcf}) \\
                {output.vcf}.sample_map.txt
            
            bcftools reheader \\
                -s {output.vcf}.sample_map.txt \\
                {input.vcf} \\
                -o {output.vcf}
            
            bcftools index -t {output.vcf}
            rm {output.vcf}.sample_map.txt
            """ + SHELL_EPILOG