{% set command = 'generateIndex' %}
# {{ command }}

::: moPepGen.cli.generate_index
	handler: python
    selection:
      members: false
    rendering:
      show_root_heading: false
      show_source: false

{% include 'partials/_caution_on_reference_version.md' %}

## Usage

```
{{ get_arg_usage(command) }}
```

## Arguments

{% with actions=get_arg_data(command) %}
{% include 'partials/_command_usage.md' %}
{% endwith %}

## Input

Three files are required for this command:

1. Reference genome FASTA file.
2. Genome annotation GTF file.
3. The protein sequence FASTA file.

All three files must be downloaded from the same release version of either [GENCODE](https://www.gencodegenes.org/) or [ENSEMBL](https://useast.ensembl.org/index.html). moPepGen does not support reference files in any other format (*e.g.* RefSeq) at this point.

For GENCODE, the primary assembly is recommended. The file name of the genome FASTA should look like 'GRCh38.primary_assembly.genome.fa'. For annotation, we recommend using the comprehensive gene annotation file for the primary assembly, which matches with the genome FASTA. The file name should look like 'gencode.v43.primary_assembly.annotation.gtf'. The protein sequence should be downloaded from the same source with the same version, and the file name should look like 'gencode.v43.pc_transcripts.fa'

Similarly, for ENSEMBL, we also recommend using the primary genome assembly and its annotation. The genome FASTA file should have a pattern of 'Homo_sapiens.GRCh38.dna.primary_assembly.fa', the GTF should look like 'Homo_sapiens.GRCh38.109.chr_patch_hapl_scaff.gtf', and the protein sequence file should look like 'Homo_sapiens.GRCh38.pep.all.fa'.

## Output

Users usually don't need to worry about the output files of this command. As long as the correct path is provided to the subsequent moPepGen commands, the correct index files will be recognized.

Files are created by this command:

| File Name | Description |
|:----------|:------------|
| `genome.pkl` | This file contains the entire reference genome. |
| `annotation.gtf` | A copy of the input annotation GTF file. If `--symlink-gtf` is used, this will be a symlink pointing to the input file. |
| `annotation_gene.idx` | A text file with the location of each gene in the GTF file. |
| `annotation_tx.idx` | A text file with the location of each transcript in the GTF file. |
| `proteome.pkl` | Contains all protein sequences. |
| `canonical_peptides.pkl` | Contains nonredundant canonical peptides from the input proteome. |
| `coding_transcripts.pkl` | Contains all coding transcripts. |
