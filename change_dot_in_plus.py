from GTF import GTF
import sys
import argparse

parser = argparse.ArgumentParser(description="Utility tools for your GTF files.")
parser.add_argument(
    "-i",
    "--input-file",
    help="Path to your GTF file. Use stdin by default",
    type=argparse.FileType("r"),
    default=(None if sys.stdin.isatty() else sys.stdin),
)
args = parser.parse_args()

for record in GTF.parse_by_line(args.input_file):
    if record.strand == ".":
        record.strand = "+"
    print(record)
