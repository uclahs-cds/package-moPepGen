## Quick Start

### Downloading Reference Files

moPepGen requires a coherent set of reference genome, proteome and Gene transfer format (GTF) files to run. The key is to ensure that the genome annotation version is consistent between the GTF and proteome FASTA and that the reference genome build version is the same between all three.

To ensure consistency, we recommend downloading the reference genome FASTA, reference proteome FASTA and genome annotation GTF from the same source. At this time reference files from GENCODE and Ensembl are supported. We do not foresee support for UniProt proteomes in the near future since there is no fail-safe method to map all IDs to formats found in commonly used genome annotation GTFs.

#### GENCODE

At the time of writing, the [GENCODE](https://www.gencodegenes.org/) [Human](https://www.gencodegenes.org/human/) release is at v41 (GRCh38.p13).

1. Under Fasta files, download the `Genome sequence (GRCh38.p13)` `ALL` `Fasta` file (`GRCh38.p13.genome.fa.gz`).
2. Under Fasta files, download the `Protein-coding transcript translation sequences` `CHR` `Fasta` file (`gencode.v41.pc_translations.fa.gz`) - Please ensure that this FASTA contains `amino acid` sequences, not nucleotide sequences.
3. Under GTF/GFF3 files, download the `Comprehensive gene annotation` `ALL` `GTF` file (`gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz`).

#### Ensembl

At the time of writing, the [Ensembl](https://www.ensembl.org/index.html) Release is at 107 (Jul 2022). We illustrate using the Human GRCh38.p13 genome from the FTP Download [page](https://www.ensembl.org/info/data/ftp/index.html).

1. Click on `DNA (FASTA)` and download `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
2. Click on `Protein sequence (FASTA)` and download `Homo_sapiens.GRCh38.pep.all.fa.gz`
3. Click on `Gene sets` `GTF` and download `Homo_sapiens.GRCh38.107.chr_patch_hapl_scaff.gtf.gz`

!!! warning
    **Do not mix and match** GENCODE and Ensembl reference files, the chromosome names and transcript IDs DO NOT MATCH

### Build Reference Index

Creating index files for the genome, proteome and GTF annotations, as well as a canonical peptide pool to be used in further steps. By default, the `generateIndex` subcommand uses trypsin as the protease, considering up to 2 miscleavages, and a minimal peptide molecular weight of 500 Da. However, those settings can be adjusted using command line arguments. See `moPepGen generateIndex --help`.

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
    -i path/to/vep_1.txt path/to/vep_2.txt
    --index-dir path/to/index \
    --output-prefix path/to/parsed_vep
```

### Call Variant Peptides

Call variant peptides using the parsed variant GVF file generated by `parseVEP`. `callVariant` takes multiple variant GVF files to combine variant types. They must all be parsed by `moPepGen` parsers. Cleavage rule, miscleavage and minimal molecular weight are also set to default values as below.

```
moPepGen callPeptide \
    -i path/to/vep.gvf \
    -o vep_moPepGen.fasta
    --index-dir path/to/index \
    --cleavage-rule trypsin \
    --miscleavage 2 \
    --min-mw 500
```
