"""
Reformat dmrseq beta output → tensorQTL phenotype BED.

Input TSV (snakemake.input.betas):
    chr  start  end  sample1  sample2  ...

Output BED (bgzipped + tabixed):
    #chr  start  end  phenotype_id  sample1  sample2  ...
"""
import pandas as pd, subprocess, sys

df = pd.read_csv(snakemake.input.betas, sep="\t")

# Build a unique phenotype_id from coordinates
df.insert(3, "phenotype_id",
          df["chr"].astype(str) + ":" +
          df["start"].astype(str) + "-" +
          df["end"].astype(str))

# Rename chr col to #chr (tensorQTL convention)
df.rename(columns={"chr": "#chr"}, inplace=True)

# Sort by chromosome then start position
df.sort_values(["#chr", "start"], inplace=True)

tmp = snakemake.output.bed.replace(".gz", "")
df.to_csv(tmp, sep="\t", index=False)

subprocess.run(["bgzip", tmp], check=True)
subprocess.run(["tabix", "-p", "bed", snakemake.output.bed], check=True)