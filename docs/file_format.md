# File Structure Documentation
  
- [File Structure Documentation](#file-structure-documentation)
  - [Transcript Variant Format](#transcript-variant-format)
    - [File Metadata](#file-metadata)
    - [Point Mutation](#point-mutation)
    - [Fusion](#fusion)
  - [CircRNA BED Format](#circrna-bed-format)
  - [Reference Index](#reference-index)
  - [Variant Peptide FASTA](#variant-peptide-fasta)


## Transcript Variant Format

In moPepGen we are interested in finding varianted peptide caused by combination of different types of variants, including single nucleotide substitution, INDEL, RNA editing site, gene fusion and alternative splicing. We are also interested in non-coding RNA and circRNA with unreported ORF, or start codon gained from mutation.

The different mutation events are called by different programs, and those files have different format. In moPepGen, for each type of data, we use data type and tool specific parsers to convert variant data from different sources to a standardized VCF-like format that moPepGen can use to create the transcript variant graph. They are then collected by `moPepGen callPeptide` command to call variant peptides.

In moPepGen, we define the TVF (Transcript Variant Format) file format, that extended and modified from the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file format. The a TVF file, each record represent a variant associated with a transcript. The `CHROM` column is used to hold the transcript_id, and the `POS` column is also representing the position of the corresponding transcript.

### File Metadata

Each variant file should contain a metadata section that each line starts with a double hashtag. The first line of the metadta must be the 'fileformat' field for VCF, to be consistant with VCF's standards. The link should read:

```
##fileformat=VCFv4.2
```

Starting from the second line should be moPepGen's metadata. Each line should be a key value pair separated by an equal sign ('='). Keys should follow the `snake_case`. See example below.

```
##moPepGen_version=0.0.1
##parser=parseXXX
##reference_index=
##genome_fasta=
##annotation_gtf=
```

+ `moPepGen_version`: the version number of moPepGen used to generate the variant file.
+ `parser`: the parser used to create the variant file.
+ `reference_index`: the path to the reference index used. Should be empty if not specified in the parser.
+ `genome_fasta`: The genome fasta file used. Should be empty if not used in the parser.
+ `annotation_gtf`: The annotation GTF file used. Should be empty if not specified in the parser.

After the moPepGen metadata section, there should be a section for field information, that defines the fields used in either `ALT` or `INFO` column. For example:

```
##INFO=<ID=ID,Number=number,Type=type,Description="description">
```

### Point Mutation

Below is an example of a TVF file for point mutation, including single nucleotide substitution and small INDEL.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseVEP
##reference_index=
##genome_fasta=
##annotation_gtf=
##CHROM=<Description='Transcript ID'>
##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Acceptor Transcript's Gene ID">
##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">
##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
ENST0001    110 SNV_110-C-A C   A   .   .   GENE_ID=ENSG0001;GENE_SYMBOL=TP53;GENOMIC_POSITION="chr1:1000-1001"
ENST0002    210 SNV_210-T-A T   A   .   .   GENE_ID=ENSG0002;GENE_SYMBOL=EGFR;GENOMIC_POSITION="chr1:1000-1001"
```

The `REF` and `ALT` must be explicit. The `INFO` column should contain the gene ID that the transcript belongs to. The `ID` column follows the pattern of '\<varant_type>_\<position>-\<ref>-\<alt>'. 

### Fusion

Below is an example of a TVF file for gene fusions.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseXXX
##reference_index=
##genome_fasta=
##annotation_gtf=
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

## CircRNA BED Format

Circular RNAs are derived from back-spliced exons. They exist as individual RNA molecules and have the potential to be translated to proteins. We are then interested in finding the possible peptide sequences translated from circRNAs with and without variants (SNP, INDEL, etc). In this case, circRNAs per se are rather new transcripts than variants. We then use a BED like format to represen the circRNA molecules. In the BED file, each row is an exon associated with a circRNA, represented as the transcript ID and the location (start and end) of the transcript. Each circRNA can have one or more exons, so the exons next to each other that share the same transcript ID and circRNA ID belong to the same circRNA molecule. The order of records represent the order of the exons apper on the circRNA. The same circRNA can be derived from more than one transcript of the same gene. In this case, each transcript will have its own serious of exons in the file.

```
##mopepgen_version=0.0.1
##parser=parseVEP
##reference_index=
##genome_fasta=
##annotation_gtf=
#transcript_id	start	end	id	gene_id	gene_name	
ENST0001	406	750	ENSG0001-3-4	ENSG0001	SYMB1
ENST0001	751	769	ENSG0001-3-4	ENSG0001	SYMB1
ENST0002	606	950	ENSG0001-3-4	ENSG0001	SYMB1	
ENST0002	951	969	ENSG0001-3-4	ENSG0001	SYMB1
```

## Reference Index


## Variant Peptide FASTA

In moPepGen, the headers of the final output variant peptide FASTA contains the transcript IDs and variants associated with this variant peptide. The header of a peptide record starts with the transcript ID, followed by the gene ID and gene symbol, and the variant IDs that it is associated with, separated by '|'. The Variant IDs are defined in the TVF files. In some cases, several non-carnonical from the same transcript may share the same variants. For example, a frameshifting variant may cause multiple non-carnonical peptides. A integer index is thus always added to the end to solve conflicts.

If a peptide is found in multiple transcripts, the information are separated by '||'.

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
