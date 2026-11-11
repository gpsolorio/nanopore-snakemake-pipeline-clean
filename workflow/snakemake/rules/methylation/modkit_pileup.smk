##############################################################################
# modkit pileup — CpG methylation pileup from aligned BAM
#
# Input:  sorted, haplotagged BAM + index from sort_index_bam
# Output: per-barcode modkit BED (5mC counts, CpG-only, strands combined)
#
# One job per barcode. Output BEDs are consumed by dss_call_dmrs.
##############################################################################

def _modkit_pileup_bam(wc):
    if GENOME_MODE == "mini":
        return str(mini_aligned_root() / f"{wc.barcode}_AR" / "minimap2.bam")
    return sample(wc.barcode, "bam", f"{wc.barcode}_merged_sorted.bam")

def _modkit_pileup_bai(wc):
    if GENOME_MODE == "mini":
        return str(mini_aligned_root() / f"{wc.barcode}_AR" / "minimap2.bam.bai")
    return sample(wc.barcode, "bam", f"{wc.barcode}_merged_sorted.bam.bai")


rule modkit_pileup:
    input:
        bam = _modkit_pileup_bam,
        bai = _modkit_pileup_bai,
        ref = str(REFERENCE_GENOME)
    output:
        bed = sample("{barcode}", "methylation", "{barcode}_modkit.bed")
    log:
        sample_log("{barcode}", "modkit_pileup.log")
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")
    threads:
        res("modkit_pileup", "threads")
    resources:
        mem_gb    = res("modkit_pileup", "mem_gb"),
        time      = res("modkit_pileup", "time"),
        account   = res("modkit_pileup", "account"),
        partition = res("modkit_pileup", "partition")
    shell:
        SHELL_PRELUDE + """
        echo "Input BAM:  {input.bam}"
        echo "Reference:  {input.ref}"
        echo "Output BED: {output.bed}"
        echo "Threads:    {threads}"
        echo ""
        mkdir -p $(dirname {output.bed})
        
        modkit pileup \\
            {input.bam} \\
            {output.bed} \\
            --cpg \\
            --combine-strands \\
            --modified-bases 5mC \\
            --ref {input.ref} \\
            --threads {threads}
        
        LINE_COUNT=$(wc -l < {output.bed})
        echo "[CHECK] modkit produced $LINE_COUNT lines in {output.bed}"
        if [ "$LINE_COUNT" -eq 0 ]; then
            echo "[ERROR] Empty BED file — modkit may have failed"
            exit 1
        fi
        """ + SHELL_EPILOG