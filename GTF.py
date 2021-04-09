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


class GTF:
    @staticmethod
    def parse(file, feature=None, strand=None, attributes=None):
        # Pass attributes as {"gene_id":"ENSG001", "transcript_biotype":"lncRNA"}
        with open(file) as fd:
            for line in fd:
                record = GTFRecord(line)
                if (feature is not None) and (record.feature != feature):
                    continue

                if strand is not None and record.strand != strand:
                    continue

                attributes_check = True
                for attribute, value in attributes.items():
                    if attribute not in record or record[attribute] != value:
                        attributes_check = False
                        break

                if attributes_check:
                    yield record
                else:
                    continue


if __name__ == "__main__":
    for line in GTF.parse(
        "short_jeq.gtf",
        attributes={"gene_id": "XLOC_000003", "transcript_id": "TCONS_00000004"},
    ):
        print(line)
