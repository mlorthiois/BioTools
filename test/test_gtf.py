from ..GTF import Attributes, GtfRecord, GtfParent, GtfTranscript, GtfGene, GTF
import pytest


class TestAttributes:
    attr = 'gene_id "g1"; transcript_id "t1-3"; exon_number 16;'
    attr_error1 = 'gene_id "g1" ; transcript_id "t1"; exon_number "1";'
    attr_error2 = 'gene_id "g1"; transcript_id "t1"; exon_number "1"; '

    def test_init(self):
        attr = Attributes.from_str(self.attr)
        assert attr == {"gene_id": "g1", "transcript_id": "t1-3", "exon_number": "16"}

    def test_misformatted_line(self):
        with pytest.raises(Exception):
            Attributes.from_str(self.attr_error1)
            Attributes.from_str(self.attr_error2)

    def test_remove_attributes(self):
        attr = Attributes.from_str(self.attr)
        assert attr == {"gene_id": "g1", "transcript_id": "t1-3", "exon_number": "16"}
        attr.remove(["transcript_id"])
        assert attr == {"gene_id": "g1", "exon_number": "16"}

    def test_filter_attributes(self):
        attr = Attributes.from_str(self.attr)
        assert attr == {"gene_id": "g1", "transcript_id": "t1-3", "exon_number": "16"}
        attr.filter(["transcript_id"])
        assert attr == {"transcript_id": "t1-3"}


class TestGtfRecord:
    line = '1\tCufflinks\texon\t1\t76\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";'

    def test_init(self):
        record = GtfRecord.from_line(self.line)
        assert record.seqname == "1"
        assert record.source == "Cufflinks"
        assert record.feature == "exon"
        assert record.start == 1
        assert record.end == 76
        assert record.score == "."
        assert record.strand == "+"
        assert record.frame == "."
        assert record.attributes == {"gene_id": "g1", "transcript_id": "t1", "exon_number": "1"}
        assert record["gene_id"] == "g1"
        assert record["transcript_id"] == "t1"
        assert record["exon_number"] == "1"

    def test_all_specific_method(self):
        record = GtfRecord.from_line(self.line)
        assert (str(record)) == self.line
        assert len(record) == 75


class TestGtfParentRecord:
    ex1 = '1\tCufflinks\texon\t15\t36\t.\t+\t.\tgene_id "1"; transcript_id "1"; exon_number "1";'
    ex2 = '1\tCufflinks\texon\t70\t84\t.\t+\t.\tgene_id "1"; transcript_id "1"; exon_number "2";'
    par = '1\tCufflinks\texon\t15\t84\t.\t+\t.\tgene_id "1"; transcript_id "1"; exon_number "1";'

    def create_parent(self):
        rec_w_children = GtfParent()
        ex1 = GtfRecord.from_line(self.ex1)
        ex2 = GtfRecord.from_line(self.ex2)
        rec_w_children.add_child(ex1)
        rec_w_children.add_child(ex2)
        return rec_w_children, ex1, ex2

    def test_GtfParentRecord_init(self):
        rec_w_children = GtfParent()
        assert len(rec_w_children.children) == 0

    def test_GtfParentRecord_add_children(self):
        rec_w_children, ex1, ex2 = self.create_parent()
        assert len(rec_w_children.children) == 2
        assert rec_w_children.children[0] == ex1
        assert rec_w_children.children[1] == ex2

    def check_GtfParentRecord_features(self):
        rec_w_children, _, _ = self.create_parent()
        assert rec_w_children.start == 15
        assert rec_w_children.end == 84
        assert rec_w_children.seqname == "1"
        assert rec_w_children.source == "Cufflinks"
        assert rec_w_children.feature == "exon"
        assert rec_w_children.score == "."
        assert rec_w_children.strand == "+"
        assert rec_w_children.frame == "."
        assert len(rec_w_children) == 69
        assert "gene_id" in rec_w_children
        assert rec_w_children["gene_id"] == "1"

    def test_GtfParent_fromat_gtf(self):
        rec_w_children = GtfParent()
        rec_w_children.add_child(GtfRecord.from_line(self.ex1))
        rec_w_children.add_child(GtfRecord.from_line(self.ex2))
        assert rec_w_children.format_to_gtf(["exon"]) == f"{self.par}\n{self.ex1}\n{self.ex2}"


class TestTranscript:
    tx1 = '1\tCuff\ttranscript\t3\t80\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";'
    ex1 = '1\tCuff\texon\t3\t55\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; exon_number "1";'
    ex2 = '1\tCuff\texon\t70\t80\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; exon_number "2";'

    def create_transcript(self):
        tx = GtfTranscript()
        ex1 = GtfRecord.from_line(self.ex1)
        ex2 = GtfRecord.from_line(self.ex2)
        tx.add_child(ex1)
        tx.add_child(ex2)
        return tx, ex1, ex2

    def test_GtfTranscript_init(self):
        GtfTranscript()

    def test_GtfTranscript_add_child(self):
        tx, ex1, ex2 = self.create_transcript()
        assert len(tx.exons) == len(tx.children) == 2
        assert tx.exons[0] == ex1
        assert tx.exons[1] == ex2

    def test_GtfTranscript_specific_fields(self):
        tx, _, _ = self.create_transcript()
        assert tx.feature == "transcript"
        assert "exon_number" not in tx.attributes

    def test_GtfTranscript_format(self):
        tx, _, _ = self.create_transcript()
        assert str(tx.to_record()) == self.tx1


class TestGene:
    g = '1\tCufflinks\tgene\t3\t80\t.\t+\t.\tgene_id "g1";'

    t1 = '1\tCufflinks\ttranscript\t5\t80\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";'
    e1 = '1\tCufflinks\texon\t5\t80\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; exon_number "1";'

    t2 = '1\tCufflinks\ttranscript\t3\t76\t.\t+\t.\tgene_id "g1"; transcript_id "tx2";'
    e2 = '1\tCufflinks\texon\t3\t40\t.\t+\t.\tgene_id "g1"; transcript_id "tx2"; exon_number "1";'
    e3 = '1\tCufflinks\texon\t45\t76\t.\t+\t.\tgene_id "g1"; transcript_id "tx2"; exon_number "2";'

    def test_Gene_init(self):
        t1 = GtfTranscript()
        e1 = GtfRecord.from_line(self.e1)
        t1.add_child(e1)

        t2 = GtfTranscript()
        e2 = GtfRecord.from_line(self.e2)
        e3 = GtfRecord.from_line(self.e3)
        t2.add_child(e2, e2["transcript_id"])
        t2.add_child(e3)

        gene = GtfGene()
        gene.add_child(t1, t1["transcript_id"])
        gene.add_child(t2, t2["transcript_id"])

        assert list(gene.transcripts.values()) == [t1, t2]
        assert gene.exons == [e1, e2, e3]
        assert str(gene.to_record()) == self.g
        assert str(t1.to_record()) == self.t1
        assert str(t2.to_record()) == self.t2


class TestGTF:
    def test_parse_by_line(self):
        with open("test/short_jeq.gtf") as fd:
            gtf = [record for record in GTF.parse_by_line(fd)]
        assert len(gtf) == 15

    def test_parse(self):
        with open("test/short.CanFam3.gtf") as fd:
            gtf = GTF.parse(fd)

        assert len(gtf) == 2
        assert [len(gene.transcripts) for gene in gtf.values()] == [2, 1]
        assert [
            len(transcript.exons) for gene in gtf.values() for transcript in gene.transcripts
        ] == [3, 2, 13]

    def test_stats(self):
        with open("test/short_jeq.gtf") as fd:  # Exon only
            assert GTF.stats(fd) == (4, 6, 15)

        with open("test/short.CanFam3.gtf") as fd:  # Full gtf
            assert GTF.stats(fd) == (2, 3, 18)
