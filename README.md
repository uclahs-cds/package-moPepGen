# <u>M</u>ulti-<u>O</u>mics <u>Pep</u>tide <u>Gen</u>erator

moPepGen (multi-omics peptide generator) uses data from multiple -omics experiments and calls variant peptides as custom database for proteogenomics library search.

moPepGen takes genomic variants such as single nucleotide variants (SNP or SNV), insertion/deletion (INDEL), gene fusion, and post transcriptional modifications such as RNA editing and alternative splicing, and detects variated peptides affected. 

## Installation

Install directly from github

```
pip install git+ssh://git@github.com:uclahs-cds/private-moPepGen.git
```

## Quick Start

### Build Reference Index

Creating index files for genome, proteome and GTF, and create a canonical peptide pool to be used in further steps. By default, the `generateIndex` subcommand uses trypsin, up to 2 miscleavages, and minimal molecular weight of 500. Da. However those settings can be adjusted using command line arguments. See `moPepGen generateIndex --help`.

```
moPepGen generateIndex \
    --genome-fasta path/to/genom.fasta \
    --annotation-gtf path/to/annotation.gtf \
    --proteome-fasta path/to/proteome.fasta \
    --output-dir path/to/index
```

### Parse VEP output

Convert the VEP output TXT files to a VCF-like format that moPepGen recognizes. It can take multiple VEP output files. For example, if SNP and INDEL are called separately, they can be parsed together.

```
moPepGen parseVEP \
    --vep-txt path/to/vep_1.txt path/to/vep_2.txt
    --index-dir path/to/index \
    --output-prefix path/to/parsed_vep
```

### Call Variant Peptides

Call variant peptides using the parsed variant table generated by `parseVEP`. `callPeptide` takes multiple variant tables. They must all be parsed by one of the parsers of `moPepGen`. Cleavage rule, miscleavage and minimal molecular weights are also set to be default as the values below.

```
moPepGen callPeptide \
    --input-variant path/to/vep.tvf \
    --index-dir path/to/index \
    --cleavage-rule trypsin \
    --miscleavage 2 \
    --min-mw 500  
    --output-fasta vep_moPepGen.fasta
```
