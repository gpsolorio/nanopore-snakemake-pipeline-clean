##############################################################################
# Convert merged VCF → PLINK2 pgen format + genotype QC
#
# Uses --make-pgen (not --make-bed) because tensorQTL's PgenReader
# requires the .pgen/.pvar/.psam format.
#
# Input VCF is the reheadered one if tensorqtl.rename_samples is true,
# otherwise the raw merged VCF is used directly.
#
# QC filters (matching Yujin's run):
#   --maf 0.01        minor allele frequency
#   --geno 0.0        no missing genotypes
#   --vcf-min-dp 5    minimum read depth
#   --vcf-min-gq 20   minimum genotype quality
##############################################################################

def _vcf_to_pgen_input(wc):
    if config["tensorqtl"].get("rename_samples", False):
        return cohort("methylqtl", "vcfs", "all_samples_final.vcf.gz")
    return cohort("methylqtl", "vcfs", "all_samples_merged.vcf.gz")


rule vcf_to_pgen:
    input:
        vcf      = _vcf_to_pgen_input,
        metadata = SAMPLE_METADATA_TSV,
    output:
        pgen = cohort("methylqtl", "pgen", "genotypes_qc.pgen"),
        pvar = cohort("methylqtl", "pgen", "genotypes_qc.pvar"),
        psam = cohort("methylqtl", "pgen", "genotypes_qc.psam"),
    log:
        cohort_log("vcf_to_pgen.log")
    params:
        plink2 = PLINK2_BIN,
        out_prefix  = lambda wc: cohort("methylqtl", "pgen", "genotypes_qc"),
        psam_file   = cohort("methylqtl", "pgen", "import.psam"),
        rename_file = cohort("methylqtl", "pgen", "id_rename.txt"),
        scripts_dir = str(SCRIPTS_DIR),
        maf         = config["tensorqtl"].get("maf_filter",  0.01),
        geno        = config["tensorqtl"].get("geno_filter", 0.1),
        genome_mode = GENOME_MODE,
    threads:
        res("vcf_to_pgen", "threads")
    resources:
        mem_gb    = res("vcf_to_pgen", "mem_gb"),
        time      = res("vcf_to_pgen", "time"),
        account   = res("vcf_to_pgen", "account"),
        partition = res("vcf_to_pgen", "partition")
    conda:
        plink2_env_yml()
    shell:
        SHELL_PRELUDE + """
        mkdir -p $(dirname {params.out_prefix})
        
        $CONDA_PREFIX/bin/python3 {params.scripts_dir}/methylqtl/make_plink_aux_files.py \\
            --metadata {input.metadata} \\
            --psam-file {params.psam_file} \\
            --rename-file {params.rename_file}
        
        {params.plink2} \\
            --vcf {input.vcf} \\
            --split-par b38 \\
            --psam {params.psam_file} \\
            --update-ids {params.rename_file} \\
            --make-pgen \\
            --max-alleles 2 \\
            --maf {params.maf} \\
            --allow-extra-chr \\
            --geno {params.geno} \\
            --out {params.out_prefix} \\
            --threads {threads}
        
        if [ "{params.genome_mode}" = "hg38" ]; then
            awk 'BEGIN{{OFS="\\t"}} /^#/{{print; next}} {{$1="chr"$1; print}}' \\
                {params.out_prefix}.pvar > {params.out_prefix}.pvar.tmp
            mv {params.out_prefix}.pvar.tmp {params.out_prefix}.pvar
        fi
        """ + SHELL_EPILOG