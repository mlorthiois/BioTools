#! /usr/bin/python3
import os


class GTFRecord:
    def __init__(self, line):
        line = line.rstrip().split("\t")
        self.seqname = line[0]
        self.source = line[1]
        self.feature = line[2]
        self.start = int(line[3])
        self.end = int(line[4])
        self.score = line[5]
        self.strand = line[6]
        self.frame = line[7]
        self._parse_attributes_as_dict(line[8])

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

    def add_child(self, child):
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
        return [transcript for transcript in self.children]

    @property
    def exons(self):
        return [exon for transcript in self.transcripts for exon in transcript.exons]


class Transcript(GTFRecordWithChildren):
    @property
    def exons(self):
        return self.children


class GTF:
    @staticmethod
    def parse(file, feature=None, strand=None, attributes=None, by_line=False):
        # Pass attributes as {"gene_id":"ENSG001", "transcript_biotype":"lncRNA"}
        with open(file) as fd:
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
    def reconstruct_full_gtf():
        pass

    @staticmethod
    def stats(file):
        exons = 0
        transcripts = set()
        genes = set()
        for exon in GTF.parse(file, feature="exon", by_line=True):
            exons += 1
            genes.add(exon["gene_id"])
            transcripts.add(exon["transcript_id"])

        return f"FILE: {os.path.abspath(file)}\n# genes:\t{len(genes)}\n# transcripts:\t{len(transcripts)}\n# exons:\t{exons}"


##################################################
if __name__ == "__main__":
    import sys

    print(GTF.stats(sys.argv[1]))

    # for gene in GTF.parse(sys.argv[1]):
    #     print(gene.format_to_gtf())
