##############################################################################
# Merge + Sort + Index BAM  (single streaming rule)
##############################################################################

# Shared-storage root for the bamlist file specifically. This MUST be on a
# filesystem visible to both the Snakemake controller process (which writes
# it, on the submission host) and the compute node that runs the job (which
# reads it) -- unlike the batch merge working directory below, which is
# created and consumed entirely within the same shell script on the same
# node and can safely live on node-local /tmp.
BAMLIST_SHARED_DIR = config.get("scratch_dir") or os.environ.get("TMPDIR", "/tmp")


def get_merge_bams_input(wc):
    pod5_matches = glob(os.path.join(str(POD5_INPUT_DIR), wc.barcode, "*.pod5"))
    if not pod5_matches:
        all_pod5 = glob(os.path.join(str(POD5_INPUT_DIR), "*.pod5"))
        pod5_matches = [p for p in all_pod5 if wc.barcode in os.path.basename(p)]

    if pod5_matches:
        ckpt = checkpoints.basecall_pod5.get(barcode=wc.barcode)
        bam_files = glob(
            os.path.join(ckpt.output.outdir, "**", "bam_pass", "*.bam"),
            recursive=True
        )
        if not bam_files:
            raise ValueError(
                f"basecall_pod5 completed for {wc.barcode} but produced no BAMs "
                f"in {ckpt.output.outdir}. Check the Dorado log."
            )
        return sorted(bam_files)

    return get_bam_paths(wc)


def _write_bamlist(wc, input):
    """
    Side-effecting params function: writes the full input BAM list to a file
    on SHARED storage (BAMLIST_SHARED_DIR), not resources.tmpdir. This keeps
    the rendered shell string small regardless of input file count -- with
    ~63k input BAMs for some samples, embedding the list directly into the
    shell command text blew past exec()/ARG_MAX at the singularity/bash
    layer. Writing here in Python and returning a path sidesteps that.

    Deliberately NOT written under resources.tmpdir: this function runs on
    the Snakemake controller host (login2) during job preparation, before
    the job is ever submitted -- but the job itself executes on a different
    compute node. If tmpdir is node-local (e.g. /tmp), a file written here
    would land on login2's /tmp and be invisible to the compute node's own
    /tmp, causing the job to fail immediately with "No such file or
    directory". BAMLIST_SHARED_DIR must stay on network storage so the file
    written by the controller is visible to whichever node runs the job.
    """
    os.makedirs(BAMLIST_SHARED_DIR, exist_ok=True)
    path = os.path.join(BAMLIST_SHARED_DIR, f"{wc.barcode}_bamlist.txt")
    with open(path, "w") as f:
        f.write("\n".join(input.bams) + "\n")
    return path


rule merge_bams:
    input:
        bams = get_merge_bams_input
    output:
        bam = sample("{barcode}", "bam", "{barcode}_merged_sorted.bam"),
        bai = sample("{barcode}", "bam", "{barcode}_merged_sorted.bam.bai")
    log:
        sample_log("{barcode}", "merge_bams.log")
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")
    threads:
        res("merge_bams", "threads")
    resources:
        mem_gb    = res("merge_bams", "mem_gb"),
        time      = res("merge_bams", "time"),
        account   = res("merge_bams", "account"),
        partition = res("merge_bams", "partition"),
        gpus      = res("merge_bams", "gpus", 0)
        # NOTE: resources.tmpdir (from profile default-resources, currently
        # /tmp) is used below ONLY for the batch merge working directory --
        # that directory is created and consumed entirely within this shell
        # script, on this same compute node, so node-local /tmp is safe and
        # fast for it. The bamlist itself is NOT under resources.tmpdir --
        # see _write_bamlist docstring above.
    params:
        outdir     = lambda wc: sample(wc.barcode, "bam"),
        bamlist    = _write_bamlist,
        n_input    = lambda wc, input: len(input.bams),
        first_bam  = lambda wc, input: input.bams[0],
        batch_size = 1000
    shell:
        SHELL_PRELUDE + """
        mkdir -p {params.outdir}

        n=$(wc -l < {params.bamlist})
        if [ "$n" -eq 1 ]; then
            samtools sort -@ {threads} -T {resources.tmpdir}/{wildcards.barcode}_sort \
                -o {output.bam} {params.first_bam}
        else
            batch_dir={resources.tmpdir}/{wildcards.barcode}_batches
            rm -rf $batch_dir
            mkdir -p $batch_dir

            split -l {params.batch_size} {params.bamlist} $batch_dir/batch_

            total_batches=$(find "$batch_dir" -maxdepth 1 -name 'batch_*' | wc -l)
            i=0
            for batch in "$batch_dir"/batch_*; do
                samtools merge -b "$batch" -@ {threads} -o "$batch_dir/merged_${{i}}.bam"
                i=$((i+1))
                printf '\r[%d/%d] Merging BAM batches...' "$i" "$total_batches"
            done
            printf '\n'

            printf '[final] Sorting merged BAM...\n'
            samtools merge -@ {threads} - "$batch_dir"/merged_*.bam | \
                samtools sort -@ {threads} -T {resources.tmpdir}/{wildcards.barcode}_sort \
                    -o {output.bam} -

            rm -rf $batch_dir
        fi
        printf '[final] Indexing BAM...\n'
        samtools index -@ {threads} {output.bam}
        printf '[done] %s\n' "{wildcards.barcode}"
        """ + SHELL_EPILOG
