# File Structure Documentation
  
- [File Structure Documentation](#file-structure-documentation)
  - [Transcript Variant Format](#transcript-variant-format)
    - [File Metadata](#file-metadata)
    - [Point Mutation](#point-mutation)
    - [Fusion](#fusion)
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
##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Acceptor Transcript's Gene ID">
#CHROM      POS    ID             REF    ALT    QUAL    FILTER    INFO
ENST0001    110    SNV-110-C-A    C      A      .       .         GENE_ID=ENSG0001
ENST0002    210    SNV-210-T-A    T      A      .       .         GENE_ID=ENSG0002
```

The `REF` and `ALT` must be explicit. The `INFO` column should contain the gene ID that the transcript belongs to. The `ID` column follows the pattern of '\<varant_type>-\<position>-\<ref>-\<alt>'. 

### Fusion

Below is an example of a TVF file for gene fusions.

```
##fileformat=VCFv4.2
##mopepgen_version=0.0.1
##parser=parseXXX
##reference_index=
##genome_fasta=
##annotation_gtf=
##ALT=<ID=FUSION,Description="Fusion">
##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Acceptor Transcript's Gene ID">
##INFO=<ID=DONOR_GENE_ID,Number=1,Type=String,Description="Donor Transcript's Gene ID">
##INFO=<ID=DONOR_TRANSCRIPT_ID,Number=1,Type=String,Description="Donor Transcript's Transcript ID">
##INFO=<ID=DONOR_POS,Number=1,Type=Integer,Description="Position of the break point of the donor transcript">
#CHROM      POS    ID                                   REF    ALT         QUAL    FILTER    INFO
ENST0001    500    FUSION-ENST0001:500-ENST0011:1000    A      <FUSION>    .       .         GENE_ID=ENSG0000;DONOR_GENE_ID=ENSG0010;DONOR_TRANSCRIPT_ID=ENST0011;DONOR_POS=1000
ENST0001    500    FUSION-ENST0001:500-ENST0012:1500    A      <FUSION>    .       .         GENE_ID=ENSG0000;DONOR_GENE_ID=ENSG0010;DONOR_TRANSCRIPT_ID=ENST0012;DONOR_POS=1500
ENST0002    500    FUSION-ENST0002:500-ENST0011:1000    A      <FUSION>    .       .         GENE_ID=ENSG0000;DONOR_GENE_ID=ENSG0010;DONOR_TRANSCRIPT_ID=ENST0011;DONOR_POS=1000
ENST0002    500    FUSION-ENST0002:500-ENST0012:1500    A      <FUSION>    .       .         GENE_ID=ENSG0000;DONOR_GENE_ID=ENSG0010;DONOR_TRANSCRIPT_ID=ENST0012;DONOR_POS=1500
ENST0021    500    FUSION-ENST0021:500-ENST0031:1000    C      <FUSION>    .       .         GENE_ID=ENSG0021;DONOR_GENE_ID=ENSG0030;DONOR_TRANSCRIPT_ID=ENST0031;DONOR_POS=1000
ENST0021    500    FUSION-ENST0021:500-ENST0032:1500    C      <FUSION>    .       .         GENE_ID=ENSG0021;DONOR_GENE_ID=ENSG0030;DONOR_TRANSCRIPT_ID=ENST0032;DONOR_POS=1500
ENST0022    500    FUSION-ENST0022:500-ENST0031:1000    C      <FUSION>    .       .         GENE_ID=ENSG0021;DONOR_GENE_ID=ENSG0030;DONOR_TRANSCRIPT_ID=ENST0031;DONOR_POS=1000
ENST0022    500    FUSION-ENST0022:500-ENST0032:1500    C      <FUSION>    .       .         GENE_ID=ENSG0021;PDONOR_GENE_ID=ENSG0030;DONOR_TRANSCRIPT_ID=ENST0032;DONOR_POS=1500
```

The `Info` column must contain the following fields:
+ `GENE_ID`: the gene ID of the acceptor transcript.
+ `DONOR_GENE_ID`: the donor transcript's gene ID.
+ `DONOR_TRANSCRIPT_ID`: the donor transcript's transcript ID.
+ `DONOR_POS`: the position of the break poitn of the donor transcript. 

In reality, gene fusion happens at the gene level. But in a TVF file, each line represents a transcript, so the same fusion events could appear multiple times, because both the acceptor and donor gene could have multiple transcript isoforms.

## Reference Index


## Variant Peptide FASTA

In moPepGen, the headers of the final output variant peptide FASTA contains the transcript IDs and variants associated with this variant peptide. The header of a peptide record starts with the transcript ID, and followed by the variant IDs that it is associated separated by '|'. The Variant IDs are defined in the TVF files. If a peptide is found in multiple transcripts, the information are separated by '||'.

```
>ENST0001|SNV-110-C-A
XXXXXXXXXXXXXXXXXX
>ENST0002|SNV-210-T-A|SNV-220-G-C
XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
>ENST0003|INDEL-350-T-TACT||ENST0004|SNV-55-C-G
XXXXXXXXXXXXXXXX
```