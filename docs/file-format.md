# File Structure Documentation

- [File Structure Documentation](#file-structure-documentation)
	- [1 Gene Variant Format](#1-gene-variant-format)
		- [1.1 File Metadata](#11-file-metadata)
		- [1.2 Point Mutation](#12-point-mutation)
		- [1.3 Fusion](#13-fusion)
		- [1.4 Alternative Splicing Site](#14-alternative-splicing-site)
		- [1.5 CircRNA](#15-circrna)
	- [2 Variant Peptide FASTA](#2-variant-peptide-fasta)


## 1 Gene Variant Format

In moPepGen we are interested in finding variant peptides caused by combinations of different types of variants, including single nucleotide substitution, INDEL, RNA editing site, gene fusion and alternative splicing. We are also interested in non-coding RNA and circRNA with unreported ORF or start codon gained from mutation.

The different mutation events are called by different algorithms with varying output formats. In moPepGen, data type and tool-specific parsers convert variant data from different sources to a standardized VCF-like format. They are used in the `moPepGen callVariant` command to create the transcript variant graph and call variant peptides.

In moPepGen, we define the GVF (Gene Variant Format) file format, extended and modified from the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file format, to represent the variant records. In a GVF file, each entry represents a variant associated with a transcript. The `CHROM` column is used to hold the gene ID, and the `POS` column indicates the variant position in reference to the transcript gene.

### 1.1 File Metadata

Each variant GVF file contains a metadata section, where each line starts with a double hashtag. The first line of the metadata must be the `fileformat` field in consistency with VCF standards. The line should read:

```
##fileformat=VCFv4.2
```

moPepGen-specific GVF metadata starts on the second line. Each line is a key-value pair separated by an equal sign (`=`). Keys follow the `snake_case`. See examples below.

```
##moPepGen_version=0.0.1
##parser=parseXXX
##reference_index=/path/to/reference-index
##genome_fasta=/path/to/genome.fasta
##annotation_gtf=/path/to/annotation.gtf
##source=SNP
```

+ `moPepGen_version`: the version number of moPepGen used to generate the variant file.
+ `parser`: the parser used to create the variant file.
+ `reference_index`: the path to the reference index used.
+ `genome_fasta`: The genome fasta file used.
+ `annotation_gtf`: The annotation GTF file used.
+ `source`: The source of variants (e.g., SNP, SNV, Fusion). This value is used in the splitDatabase subcommand.

If reference_index is not empty, `genome_fasta` and `annotation_gtf` refer to the files used to generate the index.

After the moPepGen metadata section, there should be a section for field information, that defines the fields used in either `ALT` or `INFO` column. For example:

```
##INFO=<ID=ID,Number=number,Type=type,Description="description">
```

### 1.2 Point Mutation

Below is an example of a GVF file for point mutation, including single nucleotide substitution and small INDEL.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseVEP
##reference_index=
##genome_fasta=
##annotation_gtf=
##source=SNP
##CHROM=<Description='Gene ID'>
##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
#CHROM    POS  ID           REF  ALT  QUAL  FILTER  INFO
ENSG0001  110  SNV-110-C-A  C    A    .     .       TRANSCRIPT_ID=ENST00011,ENST00012;GENE_SYMBOL=TP53;GENOMIC_POSITION="chr1:1000-1001"
ENSG0002  210  SNV-210-T-A  T    A    .     .       TRANSCRIPT_ID=ENST00021,ENST00022;GENE_SYMBOL=EGFR;GENOMIC_POSITION="chr1:1000-1001"
```

The `REF` and `ALT` must be explicit. The `INFO` column should contain the gene ID that the transcript belongs to. The `ID` column follows the pattern of '\<varant_type>-\<position>-\<ref>-\<alt>'.

### 1.3 Fusion

Below is an example of a GVF file for gene fusions.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseXXX
##reference_index=
##genome_fasta=
##annotation_gtf=
##source=Fusion
##CHROM=<Description='Gene ID'>
##ALT=<ID=FUSION,Description="Fusion">
##INFO=<ID=TRANSCRIP_ID,Number="+",Type=String,Description="5' Junction (Donor) Transcript ID">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=ACCEPTER_GENE_ID,Number=1,Type=String,Description="3' Junction (Accepter) Transcript's Gene ID">
##INFO=<ID=ACCEPTER_TRANSCRIPT_ID,Number=1,Type=String,Description="3' Junction (Accepter) Transcript's Transcript ID">
##INFO=<ID=ACCEPTER_POS,Number=1,Type=Integer,Description="Position of the break point of the accepter transcript">
##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
##INFO=<ID=DONOR_GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
##INFO=<ID=DONOR_GENE_SYMBOL,Number=1,Type=String,Description="5' Junction (Donor) Gene Symbol">
#CHROM    POS  ID                                 REF  ALT       QUAL  FILTER  INFO
ENSG0001  500  FUSION-ENSG0001:500-ENSG0011:1000  A    <FUSION>  .     .       TRANSCRIPT_ID=ENST0000;GENE_SYMBOL=SYMB1;ACCEPTER_GENE_ID=ENSG0010;ACCEPTER_TRANSCRIPT_ID=ENST0011;ACCEPTER_POSITION=1000;GENOMIC_POSITION=chr1:1000-1000;ACCEPTER_GENOMIC_POSITION=chr2:2000-2000;DONOR_GENE_SYMBOL=SYMB3
ENSG0002  500  FUSION-ENSG0002:500-ENSG0011:1000  A    <FUSION>  .     .       TRANSCRIPT_ID=ENST0000;GENE_SYMBOL=SYMB1;ACCEPTER_GENE_ID=ENSG0010;ACCEPTER_TRANSCRIPT_ID=ENST0011;ACCEPTER_POSITION=1000;GENOMIC_POSITION=chr1:1000-1000;ACCEPTER_GENOMIC_POSITION=chr2:2000-2000;DONOR_GENE_SYMBOL=SYMB3
ENSG0021  500  FUSION-ENSG0021:500-ENSG0031:1000  C    <FUSION>  .     .       TRANSCRIPT_ID=ENST0021;GENE_SYMBOL=SYMB2;ACCEPTER_GENE_ID=ENSG0030;ACCEPTER_TRANSCRIPT_ID=ENST0031;ACCEPTER_POSITION=1000;GENOMIC_POSITION=chr3:1000-1000;ACCEPTER_GENOMIC_POSITION=chr4:2000-2000;DONOR_GENE_SYMBOL=SYMB4
```

The `Info` column must contain the following fields:
+ `TRANSCRIPT_ID`: the transcript ID of the donor (upstream) transcript.
+ `ACCEPTER_GENE_ID`: the accepter (downstream) transcript gene ID.
+ `ACCEPTER_TRANSCRIPT_ID`: the accepter (downstream) transcript transcript ID.
+ `ACCEPTER_POSITION`: the position of the breakpoint of the ACCEPTER (downstream) transcript.
+ `GENOMIC_POSITION`: the genomic position of the donor (upstream) transcript, in the format of `<chrom name>:<breakpoint>-<breakpoint>`.
+ `ACCEPTER_GENOMIC_POSITION`: the genomic position of the accepter (downstream) transcript, in the format of `<chrom name>:<breakpoint>-<breakpoint>`.


### 1.4 Alternative Splicing Site

Alternative splicing sites called by [rMATS](http://rnaseq-mats.sourceforge.net/) have five types, *i.e.* skipped exon (SE), alternative 5' splice site (A5SS), alternative 3' splice site (A3SS), mutually exclusive exons (MXE), and retained intron (RI). Each alternative splicing event can be represented as a deletion, insertion or substitution.

SE is called when an exon is skipped given its upstream and downstream exons. It is represented as an **insertion** when the target transcript from the GTF file contains the exon, or is represented as a **deletion** when the target transcript is annotated without the exon.

A5SS and A3SS are called when an exon has two splicing sites that can generate either a long or a short version of the exon. When the longer version is annotated in the given transcript, the variant is represented as a **deletion**, and an **insertion** is used when the shorter version is annotated.

MXE is represented as the substitution of one exon with another exon.

RI is represented as an insertion of the intron sequence.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseRMATS
##reference_index=/path/to/reference/index
##genome_fasta=/path/to/genome.fasta
##annotation_gtf=/path/to/annotaion.gtf
##source=AlternativeSplicing
##CHROM=<Description='Gene ID'>
##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID">
##INFO=<ID=START,Number=1,Type=Integer,Description="Start Position">
##INFO=<ID=END,Number=1,Type=Integer,Description="End Position">
##INFO=<ID=DONOR_START,Number=1,Type=Integer,Description="Donor Start Position">
##INFO=<ID=DONOR_END,Number=1,Type=Integer,Description="Donor End Position">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
#CHROM    POS  ID            REF  ALT    QUAL  FILTER  INFO
ENSG0001  110  SE-300        C    <INS>  .     .       TRANSCRIPT_ID=ENST00011;DONOR_GENE_ID=ENSG0005;DONOR_START=300;DONOR_END=400;GENE_SYMBOL=TP53;GENOMIC_POSITION=chr1:1000-1001
ENSG0002  210  A5SS-210      T    <DEL>  .     .       TRANSCRIPT_ID=ENST00021;DONOR_GENE_ID=ENSG0006;DONOR_START=210;DONOR_END=400;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
ENSG0003  115  A3SS-320      T    <INS>  .     .       TRANSCRIPT_ID=ENST00031;DONOR_GENE_ID=ENSG0007;DONOR_START=320;DONOR_END=380;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
ENSG0004  277  MXE-477-1103  T    <SUB>  .     .       TRANSCRIPT_ID=ENST00041;DONOR_GENE_ID=ENSG0008;DONOR_START=477;DONOR_END=582;DONOR_START=1103;DONOR_END=1228;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

**Examples:**

```
ENSG0001	110	SE-300	C	<INS>	.	.	TRANSCRIPT_ID=ENST0001;DONOR_START=300;DONOR_END=400;GENE_SYMBOL=TP53;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents an SE (skipped exon), that the sequence of 300-400 of the gene ENSG0001 is inserted into the transcript of ENST0001 at position 110. In this case, all transcripts of the gene in the annotation GTF don't contain this exon.

```
ENSG0002	210	A5SS-210	T	<DEL>	.	.	TRANSCRIPT_ID=ENST0002;START=210;END=400;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents an A5SS (alternative 5' splicing site), where the nucleotides from positions 210 to 400 of the transcript ENST0002 are deleted. This A5SS is represented as a **deletion** because all transcripts of the gene in the annotation GTF have the longer version of the exon.

```
ENSG0003	115	MXE-320	T	<INS>	.	.	TRANSCRIPT_ID=ENST0003;START=320;END=380;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents an MXE (mutually exclusive exon), where the exon at position 320-380 of the gene ENSG0003 is retained in the transcript ENST0003 and resulted as an insertion at position 115 of the transcript. This MXE is represented as an **insertion** because none of the transcripts of this gene has the first exon retained and the second spliced in the GTF, but this transcript has both exons retained.

```
ENSG0004	277	MXE-477-1103	T	<SUB>	.	.	TRANSCRIPT_ID=ENST0004;START=477;END=582;DONOR_START=1103;DONOR_END=1228;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

This line above represents an MXE where the exon at position 447-582 in the gene (transcript ENST0004 position 277) is replaced with the exon at position 1103-1228 of the gene.

### 1.5 CircRNA

Circular RNAs are derived from back-spliced exons and introns. They exist as independent RNA molecules and have the potential to be translated into proteins. We are interested in finding the possible peptide sequences that could be the result of circRNA translation, with and without additional variants (SNP, INDEL, etc). In this case, circRNAs per se are rather new transcripts backbones than variants and are also recorded in GVF files in moPepGen. In such a GVF file, each row represents a circRNA, with the gene ID it is associated with, the start position at the gene coordinate, the offset and length of each segment, and their exon or intron indices. Normally each segment is an exon, but with intron-retained alternative splicing, they could be introns.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseCIRCexplorer
##reference_index=/path/to/reference/index
##genome_fasta=/path/to/genome.fasta
##annotation_gtf=/path/to/annotaion.gtf
##CHROM=<Description='Gene ID'>
##INFO=<ID=OFFSET,Number=+,Type=Integer,Description="Offsets of fragments (exons or introns)">
##INFO=<ID=LENGTH,Number=+,Type=Integer,Description="Length of fragments (exons or introns)">
##INFO=<ID=INTRON,Number=+,Type=Integer,Description="Indices of fragments that are introns">
##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description="Transcripts associated with this circRNA">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">
##POS=<Description="Gene coordinate of circRNA start">
#CHROM    POS  ID                            REF  ALT  QUAL  FILTER  INFO
ENSG00000128408.9	0	CIRC-ENST00000614167.2-0:464	.	.	.	.	OFFSET=0,323;LENGTH=323,82;INTRON=;TRANSCRIPT_ID=ENST00000614167.2;GENE_SYMBOL=RIBC2;GENOMIC_POSITION=chr22:0:464
ENSG00000099949.21	0	CIRC-ENST00000642151.1-0:197	.	.	.	.	OFFSET=0,98;LENGTH=78,42;INTRON=;TRANSCRIPT_ID=ENST00000642151.1;GENE_SYMBOL=LZTR1;GENOMIC_POSITION=chr22:4980:5177
ENSG00000099949.21	78	CIRC-ENST00000642151.1-78:98	.	.	.	.	OFFSET=0;LENGTH=20;INTRON=0;TRANSCRIPT_ID=ENST00000642151.1;GENE_SYMBOL=LZTR1;GENOMIC_POSITION=chr22:5058:5078
```

circRNAs are not variants that are added to the transcript variant graph, thus the `REF` and `ALT` columns should be kept empty as ".". The `INFO` column must contain the following fields.

+ **`OFFSET`**: The offset of each fragment after the `start` position of the gene. Each segment can be either an exon or an intron.
+ **`LENGTH`**: The length of each fragment.
+ **`INTRON`**: The indices of fragments that are introns.
+ **`TRANSCRIPT`** The transcript ID of a transcript that is able to generate this circRNA (e.g. contains all exons and introns of the circRNA).
+ **`GENE_SYMBOL`** The name of the gene.

circRNAs IDs follow the format CIRC-\<transcript_id>-<upstrea>:<downstream>, where \<transcript_id> is taken from the INFO column, and \<upstream> and \<downstream> represent the upstream and downstream gene coordinates of the backsplicing site. For example, `CIRC-ENST0001.1-78:135` refers to a circRNA derived from transcript ENST0001.1 with backsplicing occurring between gene positions 78 and 135.

## 2 Variant Peptide FASTA

In moPepGen, the headers of the final output variant peptide FASTA contain the transcript IDs and variants associated with this variant peptide. The header of a peptide record starts with the transcript backbone ID and is followed by the variant IDs that it is associated with, separated by '|'. The variant IDs are defined in the GVF files. An ORF ID is added if the peptide is found from a novel ORF of the transcript. In some cases, several non-canonical peptides from the same transcript may share the same variants. This is most common in cases of peptide miscleavages. In addition, a frameshifting variant may cause multiple non-canonical peptides. An integer index is thus always added to the end of each header entry to resolve redundancies.

If the same peptide is found in multiple transcripts, all are documented as separate entries in the fasta header, separated by space.

```
>ENST0001|SNV-110-C-A|1
XXXXXXXXXXXXXXXXXX
>ENST0002|SNV-210-T-A|SNV_220-G-C|1
XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
>ENST0003|INDEL-350-T-TACT|1 ENST0004|ENSG0004|SYMB4|SNV-55-C-G|1
XXXXXXXXXXXXXXXX
>ENST0005|INDEL-110-CAA-A|1
XXXXXXXXXXXXXXXXXX
>ENST0005|INDEL-110-CAA-A|2
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```

For circRNA, the FASTA headers follow this style: `<circRNA-ID>|<variant_id_1>|...|<variant_id_k>|<index>`

```
>CIRC-ENSG0001-15:68|1
XXXXXXXXXXXXX
>CIRC-ENSG0001-153:285|SNV-110-C-A|2
XXXXXXXXXXXXX
```
