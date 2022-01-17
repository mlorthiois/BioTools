#!/usr/bin/env python3
# ==========================================================================
# Matthias Lorthiois, last update on 11-2021.

# This script aims to convert seqname (1st column) of GTF files from a db
# id to other
#
# usage example:
# convert_seqname.py -f ncbi -t ensembl -c assembly_report.txt -i my.gtf
#
# Where:
#   - assembly_report.txt is found on NCBI FTP for each assembly
#   - my.gtf is the gtf you want to convert
#
# Output:
#   - converted gtf will be writtent to stdout
# ==========================================================================


def ensembl_seqname(sequence_role, assigned_molecule, genbank):
    if sequence_role == "assembled-molecule":
        return assigned_molecule
    else:
        return genbank


def check_ucsc_style_name(ucsc_name, seqrole, molecule, genbank):
    if ucsc_name != "na":
        return ucsc_name

    if seqrole == "assembled-molecule":
        return f"chr{molecule}"
    else:
        return f"chrUn_{genbank.replace('.', 'v')}"


def parse_config_line(line: str):
    # seqname, seqrole, assigned-molecule, location/type, genbank, relationship, refseq, assembly-uniq, length, ucsc
    l = line.rstrip().split("\t")

    return {
        "ucsc": check_ucsc_style_name(l[9], l[1], l[2], l[4]),
        "ncbi": l[6],
        "ensembl": ensembl_seqname(l[1], l[2], l[4]),
    }


def parse_config_file(file, from_db: str):
    config = {}

    for line in file:
        if line.startswith("#"):
            continue
        l = parse_config_line(line)
        key = l[from_db]
        config[key] = l

    return config


# ==========================================================================
def convert_gtf_record(config, line: str, db_to):
    l = line.rstrip().split("\t")
    seqname = l[0]
    l[0] = config[seqname][db_to]
    return "\t".join(l)


def convert_gtf(file, config, db_to):
    for line in file:
        if line.startswith("#"):
            print(line.rstrip())
            continue
        print(convert_gtf_record(config, line, db_to))


# ==========================================================================
if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Convert gtf seqname from ucsc/ncbi/ensembl to ucsc/ncbi/ensembl"
    )
    parser.add_argument(
        "-f",
        "--current",
        choices=["ucsc", "ncbi", "ensembl"],
        type=str,
        help="seqname style currently in your gtf",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--to",
        choices=["ucsc", "ncbi", "ensembl"],
        type=str,
        help="seqname style you want.",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Path to your GTF file. Use stdin by default",
        type=argparse.FileType("r"),
        default=(None if sys.stdin.isatty() else sys.stdin),
    )

    parser.add_argument(
        "-c",
        "--config",
        help="Path to ncbi assembly report file from ftp.ncbi.nlm.nih.gov/genomes/.../..._assembly_report.txt",
        type=argparse.FileType("r"),
        required=True,
    )

    args = parser.parse_args()
    config = parse_config_file(args.config, args.current)
    convert_gtf(args.input, config, args.to)
