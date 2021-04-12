# BioTools

## GTF.py

Provide a class that allow you to parse [GTF files](https://www.ensembl.org/info/website/upload/gff.html), reconstruct GTF consisting of exons only, or calculate some basic statistics.

### Usage

#### From CLI

To reconstruct a full GTF (with 3 levels : Gene, Transcript, Exons) from exon only, use the command below and to reconstructed GTF will be print to STDOUT

```sh
GTF.py format {gtf_path}
```

To compute the number of genes, transcripts and exons in your GTF (GTF with exons only included), use :

```sh
GTF.py stats {gtf_path}
```

#### From Python script

First, import the GTF class. This class provide a static method to parse your file:

1. Your GTF is composed of 3 levels annotation. In that case, use:

```py
from GTF import GTF

for gene in GTF.parse({your GTF file}):
  print(gene)
```

2. Your GTF is composed of exons only. In that case, add the arg `by_line=True` which will parse your GTF line by line.

```py
from GTF import GTF

for record in GTF.parse({your GTF file}, by_line=True):
  print(record)
```

---
This parse method can also take as arg the feature of your choice, or the attribute of your choice, and return only the records that match your filters.

```py
for record in GTF.parse({your GTF file}, feature="transcript", attribute={"transcript_biotype" : "lncRNA"}):
  print(record)
```

##### GTFRecord
GTFRecord object provided by the iterator is an object with attributes (seqname, source, feature, start, end, score, strand, frame).
GTFRecord also provide a special `attribute` feature which is a dic with key, value from the last column of your GTF file.

1. `len(GTFRecord)` return the length of your record
2. `str(GTFRecord)` or `print(GTFRecord)` return the GTFRecord formatted as in a GTF file
3. `attribute in GTFRecord` return `True` if the attribute is in GTFRecord
4. `GTFRecord[attribute]` return the value of this attribute. Example: `GTFRecord["gene_biotype"] == "lncRNA"`

If your GTFRecord is a Gene object (in GTF.parse(by_line=False)):
1. `gene.transcripts` return a list of all the transcripts
2. `gene.exons` return a list that contains all the exons of all the transcripts

If your GTFRecord is a Transcript object:
1. `transcript.exons` return a list of all the exons


# To Do:
1. Merge GTF
