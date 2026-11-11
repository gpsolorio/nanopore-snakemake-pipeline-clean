import argparse
import csv


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--metadata",    required=True, help="sample_metadata.tsv from pipeline")
    p.add_argument("--psam-file",   required=True)
    p.add_argument("--rename-file", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.metadata) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # PSAM using sample_id directly (matches VCF IDs)
    with open(args.psam_file, "w") as o:
        o.write("#IID\tSEX\n")
        for row in rows:
            # metadata TSV has no sex column — write 0 (unknown)
            # add a sex_col to the TSV later if needed
            o.write(f"{row['sample_id']}\t0\n")

    # Rename file is a no-op now (old ID == new ID) but plink2 still needs it
    with open(args.rename_file, "w") as o:
        o.write("#OLD_FID\tOLD_IID\tNEW_FID\tNEW_IID\n")
        for row in rows:
            sid = row["sample_id"]
            o.write(f"0\t{sid}\t0\t{sid}\n")

    print(f"psam file   → {args.psam_file}")
    print(f"rename file → {args.rename_file}")


if __name__ == "__main__":
    main()