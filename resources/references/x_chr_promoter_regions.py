#!/usr/bin/env python3
import requests

# 1. API endpoint configuration
# UCSC's REST API (NOT the marketing site at ucsc.edu).
# 'ncbiRefSeq' pulls curated and predicted RefSeq annotations.
API_URL = "https://api.genome.ucsc.edu/getData/track"
PARAMS = {
    "genome": "hg38",   # NOTE: this API uses "genome", not "db"
    "track": "ncbiRefSeq",
    "chrom": "chrX",
}

# 2. Define your promoter offsets relative to the TSS
UPSTREAM_BP = 2000
DOWNSTREAM_BP = 1000   # matches hg38_chrX_promoters_up2000_down1000.bed

print("Fetching gene annotations from UCSC API...")
try:
    response = requests.get(API_URL, params=PARAMS, timeout=60)
    response.raise_for_status()
    data = response.json()
except Exception as e:
    print(f"Error fetching data: {e}")
    raise SystemExit(1)

# 3. Process the gene data and calculate promoter regions
# The API returns a dict where the track name maps to a list of records.
track_data = data.get("ncbiRefSeq", [])

if not track_data:
    print(f"Warning: no records returned. Top-level keys were: {list(data.keys())}")

output_lines = []
for record in track_data:
    chrom = record.get("chrom")
    strand = record.get("strand")
    accession = record.get("name")    # e.g. NM_001370370
    symbol = record.get("name2")      # e.g. PLCXD1

    # txStart and txEnd are zero-based internal coordinates
    tx_start = int(record.get("txStart"))
    tx_end = int(record.get("txEnd"))

    # Calculate boundaries depending on the strand direction
    # (upstream/downstream are relative to transcription direction)
    if strand == "+":
        tss = tx_start
        promoter_start = tss - UPSTREAM_BP
        promoter_end = tss + DOWNSTREAM_BP
    elif strand == "-":
        tss = tx_end
        promoter_start = tss - DOWNSTREAM_BP
        promoter_end = tss + UPSTREAM_BP
    else:
        continue  # skip unstranded records, if any

    if promoter_start < 0:
        promoter_start = 0

    name_field = f"{symbol}_{accession}"
    output_lines.append(f"{chrom}\t{promoter_start}\t{promoter_end}\t{name_field}\t0\t{strand}")

# Sort by start coordinate, matching the style of the reference file
output_lines.sort(key=lambda line: int(line.split("\t")[1]))

out_path = "test.bed"
with open(out_path, "w") as fh:
    fh.write("\n".join(output_lines) + "\n")

print(f"Wrote {len(output_lines)} promoter regions to {out_path}")