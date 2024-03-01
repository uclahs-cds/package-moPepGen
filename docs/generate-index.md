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
3. Protein sequence FASTA file.

All three files must be downloaded from the same release version of either [GENCODE](https://www.gencodegenes.org/) or [ENSEMBL](https://useast.ensembl.org/index.html). moPepGen does not support reference files from other databases (*e.g.* RefSeq) at the moment.

For GENCODE, the primary assembly is recommended (*i.e.* 'GRCh38.primary_assembly.genome.fa'). For annotation, we recommend using the comprehensive gene annotation file for the primary assembly, which matches with the genome FASTA (*e.g.* 'gencode.vXX.primary_assembly.annotation.gtf'). The version-matched protein sequence should also be downloaded from GENCODE (*e.g.* 'gencode.vXX.pc_transcripts.fa').

Similarly, for ENSEMBL, we recommend using the primary genome assembly and its annotation. The genome FASTA file should resemble 'Homo_sapiens.GRCh38.dna.primary_assembly.fa', the GTF should look like 'Homo_sapiens.GRCh38.XX.chr_patch_hapl_scaff.gtf', and the protein sequence file should look like 'Homo_sapiens.GRCh38.pep.all.fa'.

## Output

Users usually don't need to worry about the output files of this command. As long as the correct path is provided to subsequent moPepGen commands, the correct index files will be recognized.

Files are created by this command:

| File Name | Description |
|:----------|:------------|
| `genome.pkl` | This file contains the entire reference genome. |
| `annotation.gtf` | A copy of the input annotation GTF file. If `--gtf-symlink` is used (default = `False`), this will be a symlink pointing to the input file. |
| `annotation_gene.idx` | A text file with the location of each gene in the GTF file. |
| `annotation_tx.idx` | A text file with the location of each transcript in the GTF file. |
| `proteome.pkl` | Contains all protein sequences. |
| `canonical_peptides_001.pkl` | Contains nonredundant canonical peptides from the input proteome. |
| `coding_transcripts.pkl` | Contains all coding transcripts. |
| `metadata.json` | Metadata including software versions as well as enzymes and cleavage parameters used to generate the canoincal peptide pool. |
