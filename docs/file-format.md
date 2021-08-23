# File Structure Documentation
	
- [File Structure Documentation](#file-structure-documentation)
	- [1 Transcript Variant Format](#1-transcript-variant-format)
		- [1.1 File Metadata](#11-file-metadata)
		- [1.2 Point Mutation](#12-point-mutation)
		- [1.3 Fusion](#13-fusion)
		- [1.4 Alternative Splicing Site](#14-alternative-splicing-site)
	- [2 CircRNA TSV Format](#2-circrna-tsv-format)
	- [3 Variant Peptide FASTA](#3-variant-peptide-fasta)


## 1 Transcript Variant Format

In moPepGen we are interested in finding varianted peptide caused by combination of different types of variants, including single nucleotide substitution, INDEL, RNA editing site, gene fusion and alternative splicing. We are also interested in non-coding RNA and circRNA with unreported ORF, or start codon gained from mutation.

The different mutation events are called by different programs, and those files have different format. In moPepGen, for each type of data, we use data type and tool specific parsers to convert variant data from different sources to a standardized VCF-like format that moPepGen can use to create the transcript variant graph. They are then collected by `moPepGen callPeptide` command to call variant peptides.

In moPepGen, we define the TVF (Transcript Variant Format) file format, that extended and modified from the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file format. The a TVF file, each record represent a variant associated with a transcript. The `CHROM` column is used to hold the transcript_id, and the `POS` column is also representing the position of the corresponding transcript.

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
```

+ `moPepGen_version`: the version number of moPepGen used to generate the variant file.
+ `parser`: the parser used to create the variant file.
+ `reference_index`: the path to the reference index used.
+ `genome_fasta`: The genome fasta file used.
+ `annotation_gtf`: The annotation GTF file used.

If reference_index is not empty, `genome_fasta` and `annotation_gtf` refer to the files used to generate the index.

After the moPepGen metadata section, there should be a section for field information, that defines the fields used in either `ALT` or `INFO` column. For example:

```
##INFO=<ID=ID,Number=number,Type=type,Description="description">
```

### 1.2 Point Mutation

Below is an example of a TVF file for point mutation, including single nucleotide substitution and small INDEL.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseVEP
##reference_index=/path/to/reference-index
##genome_fasta=/path/to/genome.fasta
##annotation_gtf=/path/to/annotation.gtf
##CHROM=<Description='Transcript ID'>
##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Acceptor Transcript's Gene ID">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
ENST0001    110 SNV_110-C-A C   A   .   .   GENE_ID=ENSG0001;GENE_SYMBOL=TP53;GENOMIC_POSITION="chr1:1000-1001"
ENST0002    210 SNV_210-T-A T   A   .   .   GENE_ID=ENSG0002;GENE_SYMBOL=EGFR;GENOMIC_POSITION="chr1:1000-1001"
```

The `REF` and `ALT` must be explicit. The `INFO` column should contain the gene ID that the transcript belongs to. The `ID` column follows the pattern of '\<varant_type>_\<position>-\<ref>-\<alt>'. 

### 1.3 Fusion

Below is an example of a TVF file for gene fusions.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseXXX
##reference_index=/path/to/reference-index
##genome_fasta=/path/to/genome.fasta
##annotation_gtf=/path/to/annotation.gtf
##CHROM=<Description='Transcript ID'>
##ALT=<ID=FUSION,Description="Fusion">
##INFO=<ID=GENE_ID,Number=1,Type=String,Description="3' Junction (Acceptor) Transcript's Gene ID">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=DONOR_GENE_ID,Number=1,Type=String,Description="5' Junction (Donor) Transcript's Gene ID">
##INFO=<ID=DONOR_TRANSCRIPT_ID,Number=1,Type=String,Description="5' Junction (Donor) Transcript's Transcript ID">
##INFO=<ID=DONOR_POS,Number=1,Type=Integer,Description="Position of the break point of the donor transcript">
##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
##INFO=<ID=DONOR_GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
##INFO=<ID=DONOR_GENE_SYMBOL,Number=1,Type=String,Description="5' Junction (Donor) Gene Symbol">
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
ENST0001    500 FUSION_ENST0001:500-ENST0011:1000    A   <FUSION>   .   .   GENE_ID=ENSG0000;GENE_SYMBOL=SYMB1;DONOR_GENE_ID=ENSG0010;DONOR_TRANSCRIPT_ID=ENST0011;DONOR_POS=1000;GENOMIC_POSITION=chr1:1000-1000;DONOR_GENOMIC_POSITION=chr2:2000-2000;DONOR_GENE_SYMBOL=SYMB3
ENST0001    500 FUSION_ENST0001:500-ENST0012:1500    A   <FUSION>   .   .   GENE_ID=ENSG0000;GENE_SYMBOL=SYMB1;DONOR_GENE_ID=ENSG0010;DONOR_TRANSCRIPT_ID=ENST0012;DONOR_POS=1500;GENOMIC_POSITION=chr1:1000-1000;DONOR_GENOMIC_POSITION=chr2:2000-2000;DONOR_GENE_SYMBOL=SYMB3
ENST0002    500 FUSION_ENST0002:500-ENST0011:1000    A   <FUSION>   .   .   GENE_ID=ENSG0000;GENE_SYMBOL=SYMB1;DONOR_GENE_ID=ENSG0010;DONOR_TRANSCRIPT_ID=ENST0011;DONOR_POS=1000;GENOMIC_POSITION=chr1:1000-1000;DONOR_GENOMIC_POSITION=chr2:2000-2000;DONOR_GENE_SYMBOL=SYMB3
ENST0002    500 FUSION_ENST0002:500-ENST0012:1500    A   <FUSION>   .   .   GENE_ID=ENSG0000;GENE_SYMBOL=SYMB1;DONOR_GENE_ID=ENSG0010;DONOR_TRANSCRIPT_ID=ENST0012;DONOR_POS=1500;GENOMIC_POSITION=chr1:1000-1000;DONOR_GENOMIC_POSITION=chr2:2000-2000;DONOR_GENE_SYMBOL=SYMB3
ENST0021    500 FUSION_ENST0021:500-ENST0031:1000    C   <FUSION>   .   .   GENE_ID=ENSG0021;GENE_SYMBOL=SYMB2;DONOR_GENE_ID=ENSG0030;DONOR_TRANSCRIPT_ID=ENST0031;DONOR_POS=1000;GENOMIC_POSITION=chr3:1000-1000;DONOR_GENOMIC_POSITION=chr4:2000-2000;DONOR_GENE_SYMBOL=SYMB4
ENST0021    500 FUSION_ENST0021:500-ENST0032:1500    C   <FUSION>   .   .   GENE_ID=ENSG0021;GENE_SYMBOL=SYMB2;DONOR_GENE_ID=ENSG0030;DONOR_TRANSCRIPT_ID=ENST0032;DONOR_POS=1500;GENOMIC_POSITION=chr3:1000-1000;DONOR_GENOMIC_POSITION=chr4:2000-2000;DONOR_GENE_SYMBOL=SYMB4
ENST0022    500 FUSION_ENST0022:500-ENST0031:1000    C   <FUSION>   .   .   GENE_ID=ENSG0021;GENE_SYMBOL=SYMB2;DONOR_GENE_ID=ENSG0030;DONOR_TRANSCRIPT_ID=ENST0031;DONOR_POS=1000;GENOMIC_POSITION=chr3:1000-1000;DONOR_GENOMIC_POSITION=chr4:2000-2000;DONOR_GENE_SYMBOL=SYMB4
ENST0022    500 FUSION_ENST0022:500-ENST0032:1500    C   <FUSION>   .   .   GENE_ID=ENSG0021;GENE_SYMBOL=SYMB2;DONOR_GENE_ID=ENSG0030;DONOR_TRANSCRIPT_ID=ENST0032;DONOR_POS=1500;GENOMIC_POSITION=chr3:1000-1000;DONOR_GENOMIC_POSITION=chr4:2000-2000;DONOR_GENE_SYMBOL=SYMB4
```

The `Info` column must contain the following fields:
+ `GENE_ID`: the gene ID of the acceptor transcript.
+ `DONOR_GENE_ID`: the donor transcript's gene ID.
+ `DONOR_TRANSCRIPT_ID`: the donor transcript's transcript ID.
+ `DONOR_POS`: the position of the break poitn of the donor transcript. 
+ `GENOMIC_POSITION`: the genomic position, in the format of `<chrom name>:<breakpoint>:<breakpoint>`.
+ `DONOR_GENOMIC_POSITION`: the genomic position of the donor, in the format of `<chrom name>:<breakpoint>-<breakpoint>`.

In reality, gene fusion happens at the gene level. But in a TVF file, each line represents a transcript, so the same fusion events could appear multiple times, because both the acceptor and donor gene could have multiple transcript isoforms.

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
##CHROM=<Description='Transcript ID'>
##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Acceptor Transcript's Gene ID">
##INFO=<ID=START,Number=1,Type=Integer,Description="Start Position">
##INFO=<ID=END,Number=1,Type=Integer,Description="End Position">
##INFO=<ID=DONOR_START,Number=1,Type=Integer,Description="Donor Start Position">
##INFO=<ID=DONOR_END,Number=1,Type=Integer,Description="Donor End Position">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
ENST0001	110	SE-300	C	<INS>	.	.	GENE_ID=ENSG0001;START=300;END=400;GENE_SYMBOL=TP53;GENOMIC_POSITION=chr1:1000-1001
ENST0002	210	A5SS-210	T	<DEL>	.	.	GENE_ID=ENSG0002;START=210;END=400;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
ENST0003	115	A3SS-320	T	<INS>	.	.	GENE_ID=ENSG0003;START=320;END=380;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
ENST0003	115	MXE-320	T	<INS>	.	.	GENE_ID=ENSG0003;START=320;END=380;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
ENST0004	277	MXE-477-1103	T	<SUB>	.	.	GENE_ID=ENSG0004;START=477;END=582;DONOR_START=1103;DONOR_END=1228;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

**Examples:**

```
ENST0001	110	SE-300	C	<INS>	.	.	GENE_ID=ENSG0001;START=300;END=400;GENE_SYMBOL=TP53;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents an SE (skipped exon), that the sequence of 300-400 of the gene ENSG0001 is inserted to the t ranscript of ENST0001 at position 110. In this case, all transcripts of the gene in the annotation GTF don't contain this exon.

```
ENST0002	210	A5SS-210	T	<DEL>	.	.	GENE_ID=ENSG0002;START=210;END=400;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents a A5SS (alternative 5' splicing site), that the sequence from 210 to 400 of the transcript ENST0002 is deleted. In this case, all transcripts of the gene in the annotation GTF have the longer version of the exon.

```
ENST0003	115	MXE-320	T	<INS>	.	.	GENE_ID=ENSG0003;START=320;END=380;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

The line above represents a MXE (mutually exclusive exon), that the exon of 320-380 of the gene ENSG0003 is retained in the transcript ENST0003 and resulted as an insertion at position 115 of the transcript. In this case, none of the transcripts of this gene has the first exon retained and second spliced at the same time. And this transcript has both exons retained.

```
ENST0004	277	MXE-477-1103	T	<SUB>	.	.	GENE_ID=ENSG0004;START=477;END=582;DONOR_START=1103;DONOR_END=1228;GENE_SYMBOL=EGFR;GENOMIC_POSITION=chr1:1000-1001
```

This line above represents a MXE that the exon 447-582 (transcript ENST0004 position 277) is replaced with exon 1103-1228 of the gene.

## 2 CircRNA TSV Format

Circular RNAs are derived from back-spliced exons. They exist as individual RNA molecules and have the potential to be translated to proteins. We are then interested in finding the possible peptide sequences translated from circRNAs with and without variants (SNP, INDEL, etc). In this case, circRNAs per se are rather new transcripts than variants. Here we define a TSV file format to represent the circRNA molecules. In this TSV format, each row represent a circRNA, with the gene ID it is associated with, the start position at the gene, the offset and length of each segment, and IDs. Normally ach segment is an exon, but with intron retained alternative splicing, there could be introns.

```
##mopepgen_version=0.0.1
##parser=parseXXX
##reference_index=/path/to/reference-index
##genome_fasta=/path/to/genome.fasta
##annotation_gtf=/path/to/annotation.gtf
#gene_id	start	offset	length	intron	circ_id	transcript_id	gene_name
ENSG0001	413	0,211,398	72,85,63	.	CIRC-ENSG0001-E2-E3-E4	ENST0001,ENST0002	SYMB1
ENSG0002	112	0,175 	72,85	.	CIRC-ENSG0001-E3-E4	ENST0011,ENST0012	SYMB2
ENSG0002	112	0,73,175 	72,103,85	2	CIRC-ENSG0001-E3-I3-E4	ENST0011,ENST0012	SYMB2
ENSG0003	77	0,181,424	100,175,85	.	CIRC-ENSG0003-E2-E3-E4	ENST0021	SYMB3
ENSG0003	77	0,101,181,357,424	100,80,175,67,85	2,4	CIRC-ENSG0003-E2-I2-E3-I3-E4	ENST0021	SYMB3
ENSG0004	789	0	112	1	CI-ENSG0004-I3	ENST0041	SYMB4
```

The circRNA TSV file is defined here to represent all circRNAs to be passed to moPepGen to call for variant peptides. In the TSV file, each row represents a circRNA. The TSV file has the columns below:

+ **`gene_id`**: The gene ID of the gene where the circRNA is derived.
+ **`start`**: The start position of the circRNA at the gene.
+ **`offset`**: The offset of each fragment after the `start` position of the gene. Each segment can be either an exon or intron.
+ **`length`**: The length of each fragmet.
+ **`intron`**: The indices of fragments that are introns.
+ **`circ_id`**: The ID of circRNAs have two components. They all start with CIRC-\<gene_id> where `gene_id` is the value from the first column. Following that is the information for each fragment including E (exon) or I (intron) and the index of the fragment. For example, CIRC-ENSG0001-E2-I2-E3 is made up of the second exon, second intron, and the third exon of the gene ENSG0001.
+ **`transcript_id`** The transcript IDs that are able to generate this circRNA (e.g. contains all exons and introns of the circRNA.)
+ **`gene_name`** The name of the gene.

## 3 Variant Peptide FASTA

In moPepGen, the headers of the final output variant peptide FASTA contains the transcript IDs and variants associated with this variant peptide. The header of a peptide record starts with the transcript ID, followed by the gene ID and gene symbol, and the variant IDs that it is associated with, separated by '|'. The Variant IDs are defined in the TVF files. In some cases, several non-canonical peptides from the same transcript may share the same variants. This is most common in cases of peptide miscleavages. In addition, a frameshifting variant may cause multiple non-canonical peptides. A integer index is thus always added to the end to resolve redundancies.

If the same peptide is found in multiple transcripts, the annotation is separated by ';'.

```
>ENST0001|ENSG0001|SYMB1|SNV_110-C-A|1
XXXXXXXXXXXXXXXXXX
>ENST0002|ENSG0002|SYMB2|SNV_210-T-A|SNV_220-G-C|1
XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
>ENST0003|ENSG0003|SYMB3|INDEL_350-T-TACT|1||ENST0004|ENSG0004|SYMB4|SNV-55-C-G|1
XXXXXXXXXXXXXXXX
>ENST0005|ENSG0005|SYMB5|INDEL_110-CAA-A|1
XXXXXXXXXXXXXXXXXX
>ENST0005|ENSG0005|SYMB5|INDEL_110-CAA-A|2
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```
