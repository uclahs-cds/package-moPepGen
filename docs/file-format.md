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

In moPepGen we are interested in finding varianted peptide caused by combination of different types of variants, including single nucleotide substitution, INDEL, RNA editing site, gene fusion and alternative splicing. We are also interested in non-coding RNA and circRNA with unreported ORF, or start codon gained from mutation.

The different mutation events are called by different programs, and those files have different format. In moPepGen, for each type of data, we use data type and tool specific parsers to convert variant data from different sources to a standardized VCF-like format that moPepGen can use to create the transcript variant graph. They are then collected by `moPepGen callPeptide` command to call variant peptides.

In moPepGen, we define the GVF (Gene Variant Format) file format, that extended and modified from the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file format to represent the variant records. In a GVF file, each entry represents a variant associated with a transcript. The `CHROM` column is used to hold the gene ID, and the `POS` column indicates the position of the corresponding transcript.

### 1.1 File Metadata

Each variant file should contain a metadata section that each line starts with a double hashtag. The first line of the metadta must be the 'fileformat' field for VCF, to be consistant with VCF's standards. The link should read:

```
##fileformat=VCFv4.2
```

Starting from the second line should be moPepGen's metadata. Each line should be a key value pair separated by an equal sign ('='). Keys should follow the `snake_case`. See example below.

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
ENSG0001  500  FUSION-ENST0001:500-ENST0011:1000  A    <FUSION>  .     .       TRANSCRIPT_ID=ENST0000;GENE_SYMBOL=SYMB1;ACCEPTER_GENE_ID=ENSG0010;ACCEPTER_TRANSCRIPT_ID=ENST0011;ACCEPTER_POSITION=1000;GENOMIC_POSITION=chr1:1000-1000;ACCEPTER_GENOMIC_POSITION=chr2:2000-2000;DONOR_GENE_SYMBOL=SYMB3
ENSG0002  500  FUSION-ENST0002:500-ENST0011:1000  A    <FUSION>  .     .       TRANSCRIPT_ID=ENST0000;GENE_SYMBOL=SYMB1;ACCEPTER_GENE_ID=ENSG0010;ACCEPTER_TRANSCRIPT_ID=ENST0011;ACCEPTER_POSITION=1000;GENOMIC_POSITION=chr1:1000-1000;ACCEPTER_GENOMIC_POSITION=chr2:2000-2000;DONOR_GENE_SYMBOL=SYMB3
ENSG0021  500  FUSION-ENST0021:500-ENST0031:1000  C    <FUSION>  .     .       TRANSCRIPT_ID=ENST0021;GENE_SYMBOL=SYMB2;ACCEPTER_GENE_ID=ENSG0030;ACCEPTER_TRANSCRIPT_ID=ENST0031;ACCEPTER_POSITION=1000;GENOMIC_POSITION=chr3:1000-1000;ACCEPTER_GENOMIC_POSITION=chr4:2000-2000;DONOR_GENE_SYMBOL=SYMB4
```

The `Info` column must contain the following fields:
+ `TRANSCRIPT_ID`: the transcript ID of the donor (upstream) transcript.
+ `ACCEPTER_GENE_ID`: the accepter (downstream) transcript's gene ID.
+ `ACCEPTER_TRANSCRIPT_ID`: the accepter (downstream) transcript's transcript ID.
+ `ACCEPTER_POSITION`: the position of the break point of the ACCEPTER (downstream) transcript.
+ `GENOMIC_POSITION`: the genomic position of the donor (upstream) transcript, in the format of `<chrom name>:<breakpoint>:<breakpoint>`.
+ `ACCEPTER_GENOMIC_POSITION`: the genomic position of the accepter (downstream) transcript, in the format of `<chrom name>:<breakpoint>-<breakpoint>`.


### 1.4 Alternative Splicing Site

Alternative splicing site called by [rMATS](http://rnaseq-mats.sourceforge.net/) has five types, e.g. skipped exon (SE), alternative 5' splice site (A5SS), alternative 3' splice site (A3SS), mutually exclusive exons (MXE), and retained intron (RI). Each alternative splicing event can be represented as a deletion, insertion or a substitution.

SE is when a exon is skipped given its upstream and downstream exon. It is represented as a **insertion** when the target transcript from the GTF file contains the exon. And it is represented as a **deletion** when the target transcript is annotated without the exon.

A5SS and A3SS are when an exon has two splicing sites that can generate a longer and a short version. When the longer version is annotated in the given transcript, the variant is represented as a deletion, and a insertion when the shorter version is annotated.

MXE is represented as substitution of one exon with another exon.

RI is represented as an insertion or the intron sequence.

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
ENSG0004  277  MXE-477-1103  T    <SUB>  .     .       TRANSCRIPT_ID=ENST00041;DONOR_GENE_ID=ENSG0008;DONOR_START=477;DONOREND=582;DONOR_START=1103;DONOR_END=1228;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

**Examples:**

```
ENSG0001	110	SE-300	C	<INS>	.	.	TRANSCRIPT_ID=ENST0001;DONOR_START=300;DONOR_END=400;GENE_SYMBOL=TP53;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents an SE (skipped exon), that the sequence of 300-400 of the gene ENSG0001 is inserted to the t ranscript of ENST0001 at position 110. In this case, all transcripts of the gene in the annotation GTF don't contain this exon.

```
ENSG0002	210	A5SS-210	T	<DEL>	.	.	TRANSCRIPT_ID=ENST0002;START=210;END=400;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents a A5SS (alternative 5' splicing site), that the sequence from 210 to 400 of the transcript ENST0002 is deleted. In this case, all transcripts of the gene in the annotation GTF have the longer version of the exon.

```
ENSG0003	115	MXE-320	T	<INS>	.	.	TRANSCRIPT_ID=ENST0003;START=320;END=380;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents a MXE (mutually exclusive exon), that the exon of 320-380 of the gene ENSG0003 is retained in the transcript ENST0003 and resulted as an insertion at position 115 of the transcript. In this case, none of the transcripts of this gene has the first exon retained and second spliced at the same time. And this transcript has both exons retained.

```
ENSG0004	277	MXE-477-1103	T	<SUB>	.	.	TRANSCRIPT_ID=ENST0004;START=477;END=582;DONOR_START=1103;DONOR_END=1228;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

This line above represents a MXE that the exon 447-582 (transcript ENST0004 position 277) is replaced with exon 1103-1228 of the gene.

### 1.5 CircRNA

Circular RNAs are derived from back-spliced exons. They exist as individual RNA molecules and have the potential to be translated to proteins. We are then interested in finding the possible peptide sequences translated from circRNAs with and without variants (SNP, INDEL, etc). In this case, circRNAs per se are rather new transcripts than variants. Here we define a TSV file format to represent the circRNA molecules. In this TSV format, each row represent a circRNA, with the gene ID it is associated with, the start position at the gene, the offset and length of each segment, and IDs. Normally ach segment is an exon, but with intron retained alternative splicing, there could be introns.

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
ENSG0001  413  CIRC-ENST0001-E2-E3-E4        .    .    .     .       OFFSET=0,211,398;LENGTH=72,85,63;INTRON=;TRANSCRIPT=ENST0001;GENE_SYMBOL=SYMB1
ENSG0002  112  CIRC-ENST0001-E3-E4           .    .    .     .       OFFSET=0,175LENGTH=72,85;INTRON=;TRANSCRIPT=ENST0001;GENE_SYMBOL=SYMB2
ENSG0002  112  CIRC-ENST0001-E3-I3-E4        .    .    .     .       OFFSET=0,73,175;LENGTH=72,103,85;INTRON=;TRANSCRIPT=ENST0001,ENST0012;GENE_SYMBOL=SYMB2
ENSG0003  77   CIRC-ENST0003-E2-E3-E4        .    .    .     .       OFFSET=0,181,424;LENGTH=100,175,85;INTRON=;TRANSCRIPT=ENST0003;GENE_SYMBOL=SYMB3
ENSG0003  77   CIRC-ENST0003-E2-I2-E3-I3-E4  .    .    .     .       OFFSET=0,101,181,357,424;LENGTH=100,80,175,67,85;INTRON=;TRANSCRIPT=ENST0003;GENE_SYMBOL=SYMB3
ENSG0004  789  CI-ENST0004-I3                .    .    .     .       OFFSET=0;LENGTH=112;INTRON=1;TRNASCRIPT=ENST0004;GENE_SYMBOL=SYMB4
```

Technically, circRNAs are not variants that alters the gene/transcript sequence. We here still use the GVF file format to tr The `Info` column must contain the following fields:

+ **`OFFSET`**: The offset of each fragment after the `start` position of the gene. Each segment can be either an exon or intron.
+ **`LENGTH`**: The length of each fragmet.
+ **`INTRON`**: The indices of fragments that are introns.
+ **`TRANSCRIPT`** The transcript ID that are able to generate this circRNA (e.g. contains all exons and introns of the circRNA.)
+ **`GENE_SYMBOL`** The name of the gene.

The ID of circRNAs consist of two components. They all start with \<transcript_id>-circRNA or \<transcript_id>-ciRNA where `transcript_id` is the value from the `CHROM` column. Following that is the information for each fragment including E (exon) or I (intron) and the index of the fragment. For example,ENSG0001-circRNA-E2-I2-E3 is made up of the second exon, second intron, and the third exon of the gene ENSG0001.
## 2 Variant Peptide FASTA

In moPepGen, the headers of the final output variant peptide FASTA contains the transcript IDs and variants associated with this variant peptide. The header of a peptide record starts with the transcript ID, followed by the gene ID and gene symbol, and the variant IDs that it is associated with, separated by '|'. The Variant IDs are defined in the GVF files. In some cases, several non-canonical peptides from the same transcript may share the same variants. This is most common in cases of peptide miscleavages. In addition, a frameshifting variant may cause multiple non-canonical peptides. A integer index is thus always added to the end to resolve redundancies.

If the same peptide is found in multiple transcripts, the annotation is separated by space.

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
>ENSG0001-circRNA-E3-E4|1
XXXXXXXXXXXXX
>ENSG0001-circRNA-E3-E4|SNV-110-C-A|2
XXXXXXXXXXXXX
```
