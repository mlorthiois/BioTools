#!/usr/bin/env python3


class GTFRecord:
    def __init__(self, line):
        # Remove comment
        line_splitted = line.split("#")[0].rstrip().split("\t")
        if len(line_splitted) != 9:
            sys.exit(f"{line}\nUnable to parse this line, maybe badly formatted")
        self.seqname = line_splitted[0]
        self.source = line_splitted[1]
        self.feature = line_splitted[2]
        self.start = int(line_splitted[3])
        self.end = int(line_splitted[4])
        self.score = line_splitted[5]
        self.strand = line_splitted[6]
        self.frame = line_splitted[7]
        self._parse_attributes_as_dict(line_splitted[8])

    def __contains__(self, attribute):
        return attribute in self.attributes

    def __str__(self):
        begin = "\t".join(
            [
                self.seqname,
                self.source,
                self.feature,
                str(self.start),
                str(self.end),
                self.score,
                self.strand,
                self.frame,
            ]
        )
        end = " ".join(
            [f'{attribute} "{value}";' for attribute, value in self.attributes.items()]
        )
        return begin + "\t" + end

    def __len__(self):
        return abs(self.end - self.start)

    def __getitem__(self, item):
        return self.attributes[item]

    def __setitem__(self, key, value):
        self.attributes[key] = value

    def _parse_attributes_as_dict(self, att):
        self.attributes = dict(
            x.split("|")
            for x in att.replace("; ", ";")
            .replace(' "', "|")
            .replace('";', ";")
            .split(";")
            # x != "" for last empty case in split (finish by ;)
            if x != ""
        )


class GTFRecordWithChildren(GTFRecord):
    def __init__(self, line):
        super().__init__(line)
        self.children = []

    def add_child(self, child, check_position=False):
        if check_position:
            if child.start < self.start:
                self.start = child.start
            if child.end > self.end:
                self.end = child.end
        self.children.append(child)

    def format_to_gtf(self):
        gtf_seq = str(self) + "\n"
        for child in self.children:
            if isinstance(child, GTFRecordWithChildren):
                gtf_seq += child.format_to_gtf() + "\n"
            elif isinstance(child, GTFRecord):
                gtf_seq += str(child) + "\n"
        return gtf_seq.rstrip()


class Gene(GTFRecordWithChildren):
    @property
    def transcripts(self):
        return self.children

    @property
    def exons(self):
        return [exon for transcript in self.transcripts for exon in transcript.exons]


class Transcript(GTFRecordWithChildren):
    @property
    def exons(self):
        return self.children


class GTF:
    @staticmethod
    def parse(fd, feature=None, strand=None, attributes=None, by_line=False):
        if fd.readline() == "":
            sys.exit("There is nothing to parse...")
        # Pass attributes as {"gene_id":"ENSG001", "transcript_biotype":"lncRNA"}
        if by_line:
            for line in fd:
                if line.startswith("#"):
                    continue
                record = GTFRecord(line)
                if (feature is not None) and (record.feature != feature):
                    continue

                if strand is not None and record.strand != strand:
                    continue

                attributes_check = True
                for attribute, value in record.attributes.items():
                    if attribute not in record or record[attribute] != value:
                        attributes_check = False
                        break

                if attributes_check:
                    yield record
                else:
                    continue
        else:
            gene = None
            for line in fd:
                if line.startswith("#"):
                    continue

                record = GTFRecord(line)
                if record.feature == "gene" and gene is not None:
                    yield gene
                if record.feature == "gene":
                    gene = Gene(line)
                elif record.feature == "transcript":
                    transcript = Transcript(line)
                    gene.add_child(transcript)
                else:
                    transcript.add_child(record)

            yield gene

    @staticmethod
    def reconstruct_full_gtf(file):
        gene = None
        transcript = None
        for record in GTF.parse(file, by_line=True):
            if gene is None or record["gene_id"] != gene["gene_id"]:
                if gene is not None:
                    gene.add_child(transcript)
                    yield gene
                gene = Gene(str(record))
                gene.feature = "gene"
                for attribute in list(gene.attributes.keys()):
                    if "exon" in attribute or "transcript" in attribute:
                        del gene.attributes[attribute]

            if (
                transcript is None
                or record["transcript_id"] != transcript["transcript_id"]
            ):
                if transcript is not None and gene["gene_id"] == transcript["gene_id"]:
                    gene.add_child(transcript)
                transcript = Transcript(str(record))
                transcript.feature = "transcript"
                for attribute in list(transcript.attributes.keys()):
                    if "exon" in attribute:
                        del transcript.attributes[attribute]

            if record.start < gene.start:
                gene.start = record.start
            if record.start < transcript.start:
                transcript.start = record.start

            if record.end > gene.end:
                gene.end = record.end
            if record.end > transcript.end:
                transcript.end = record.end

            transcript.add_child(record)
        gene.add_child(transcript)
        yield gene

    @staticmethod
    def stats(file):
        exons = 0
        transcripts = set()
        genes = set()
        for exon in GTF.parse(file, feature="exon", by_line=True):
            exons += 1
            genes.add(exon["gene_id"])
            transcripts.add(exon["transcript_id"])
        return len(genes), len(transcripts), exons


##################################################
if __name__ == "__main__":
    import sys
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Utility tools for your GTF files.")
    parser.add_argument(
        "mode",
        choices=["stats", "format"],
        type=str,
        help="Basic stats about your file | Format a gtf to with exon lines to gene and transcript levels",
    )
    parser.add_argument(
        "-i",
        "--input-file",
        help="Path to your GTF file. Use stdin by default",
        type=argparse.FileType("r"),
        default=(None if sys.stdin.isatty() else sys.stdin),
    )

    args = parser.parse_args()

    if args.input_file is None:
        print(
            "\033[91mPlease specify your GTF file or use stdin... See below for usage:\n\x1b[0m"
        )

        sys.exit(parser.print_help())

    if args.mode == "format":
        for gene in GTF.reconstruct_full_gtf(args.input_file):
            print(gene.format_to_gtf())

    elif args.mode == "stats":
        genes, transcripts, exons = GTF.stats(args.input_file)
        if args.input_file.name != "<stdin>":
            print(f"FILE: {os.path.abspath(args.input_file.name)}")
        print(f"# genes:\t{genes}\n# transcripts:\t{transcripts}\n# exons:\t{exons}")
