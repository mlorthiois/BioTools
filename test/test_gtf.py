from ..GTF import *


class TestGtfRecord:
    line = '1\tCufflinks\texon\t1\t76\t.\t+\t.\tgene_id "XLOC_000004"; transcript_id "TCONS_00000026"; exon_number "1"; gene_name "CTDP1";'

    # Test Record
    def test_init(self):
        record = GTFRecord(self.line)
        assert record.seqname == "1"
        assert record.source == "Cufflinks"
        assert record.feature == "exon"
        assert record.start == 1
        assert record.end == 76
        assert record.score == "."
        assert record.strand == "+"
        assert record.frame == "."

        # Test _parse_attributes_as_dict
        assert len(record.attributes) == 4
        assert record.attributes["gene_id"] == "XLOC_000004"
        assert record.attributes["transcript_id"] == "TCONS_00000026"
        assert record.attributes["exon_number"] == "1"
        assert record.attributes["gene_name"] == "CTDP1"

    # Test contains
    def test_all_specific_method(self):
        record = GTFRecord(self.line)
        # __contains__
        assert ("gene_id" in record) == True

        # __str__
        assert (str(record)) == self.line

        # __len__
        assert len(record) == 75

        # __setitem__
        record["gene_id"] = "new_gene_id"
        assert record["gene_id"] == "new_gene_id"

        # __delitem__
        del record["gene_name"]
        assert ("gene_name" not in record) == True

    def test_GeneClass(self):
        with open("test/short.CanFam3.gtf") as fd:
            gtf = GTF.parse(fd)

        gene = [gene for gene in gtf.values()][0]
        assert len(gene.transcripts) == 2
        assert (
            len(gene.transcripts[0].children) == 4
        )  # Test with 1 CDS not counted as exon
        assert len(gene.exons) == 5

    def test_remove_attributes(self):
        record = GTFRecord(self.line)
        assert len(record.attributes) == 4
        record.remove_attributes(["gene_name"])
        assert len(record.attributes) == 3

    def test_filter_attributes(self):
        record = GTFRecord(self.line)
        assert len(record.attributes) == 4
        record.filter_attributes(["gene_name", "transcript_id"])
        assert len(record.attributes) == 2
        assert record.attributes["gene_name"] == "CTDP1"
        assert record.attributes["transcript_id"] == "TCONS_00000026"

    # Test gtf parsers
    def test_parse_by_line(self):
        with open("test/short_jeq.gtf") as fd:
            gtf = [record for record in GTF.parse_by_line(fd)]
        assert len(gtf) == 15

    def test_parse(self):
        with open("test/short.CanFam3.gtf") as fd:
            gtf = GTF.parse(fd)
        # Nb of genes
        assert len(gtf) == 2
        # Nb of transcripts
        assert [len(gene.transcripts) for gene in gtf.values()] == [2, 1]
        # Exons by transcripts
        assert [
            len(transcript.exons)
            for gene in gtf.values()
            for transcript in gene.transcripts
        ] == [3, 2, 13]

    # Test gtf 3 level reconstruction
    def test_reconstruct(self):
        with open("test/short_jeq.gtf") as fd:
            genes = [gene for gene in GTF.reconstruct_full_gtf(fd)]
        # Nb of genes
        assert len(genes) == 4
        # Nb of transcripts
        assert [len(gene.transcripts) for gene in genes] == [1, 1, 3, 1]
        # Exons by transcripts
        assert [
            len(transcript.exons) for gene in genes for transcript in gene.transcripts
        ] == [2, 2, 3, 3, 3, 2]

    # Test gtf stats
    def test_stats(self):
        # Exon only
        with open("test/short_jeq.gtf") as fd:
            genes, transcripts, exons = GTF.stats(fd)
        assert (genes, transcripts, exons) == (4, 6, 15)

        # Full gtf
        with open("test/short.CanFam3.gtf") as fd:
            genes, transcripts, exons = GTF.stats(fd)
        assert (genes, transcripts, exons) == (2, 3, 18)
