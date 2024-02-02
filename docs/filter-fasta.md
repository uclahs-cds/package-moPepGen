{% set command = 'filterFasta' %}
# {{ command }}

::: moPepGen.cli.filter_fasta
	handler: python
    selection:
      members: false
    rendering:
      show_root_heading: false
      show_source: false

## Usage

```
{{ get_arg_usage(command) }}
```

## Arguments

{% with actions=get_arg_data(command) %}
{% include 'partials/_command_usage.md' %}
{% endwith %}

## Examples

### Filter by Expression

The example below filters the variant peptide sequences based on their expression level. The expression table is given as TSV file, with the first column being the transcript ID, and the fourth column being the expression level. Peptides are removed if the transcript it is associated has the expression level smaller than 2. Any transcript quantitation value can be used, including read count, TPM, and FPKM.

```bash
moPepGen fitlerFasta \
  --input-path path/to/variant_peptides.fasta \
  --output-path path/to/variant_peptides_filter.fasta \
  --exprs-table path/to/expression.tsv \
  --delimiter '\t' \
  --tx-id-col 1 \
  --quant-col 4 \
  --quant-cutoff 2
```

### Filter by Expression and Miscleavages

This example is the same as above with the addition of a filter by miscleavages. Any peptides with more than 3 miscleavages will be dropped.

```bash
moPepGen fitlerFasta \
  --input-path path/to/variant_peptides.fasta \
  --output-path path/to/variant_peptides_filter.fasta \
  --exprs-table path/to/expression.tsv \
  --delimiter '\t' \
  --tx-id-col 1 \
  --quant-col 4 \
  --quant-cutoff 2 \
  --miscleavages 0:2 \
  --enzyme trypsin
```

### Filter Variant and Noncoding Peptides

This example takes both the variant peptide FASTA and the noncoding peptide FASTA and filters the peptides based on the expression level of the transcripts they are associated with.

```bash
moPepGen fitlerFasta \
  --input-path \
    path/to/variant_peptides.fasta \
    path/to/noncoding_peptides.fasta \
  --output-path path/to/variant_peptides_filter.fasta \
  --exprs-table path/to/expression.tsv \
  --delimiter '\t' \
  --tx-id-col 1 \
  --quant-col 4 \
  --quant-cutoff 2
```

### Filter by Denylist

This example here removes any peptide sequences that appear in the given denylist.

!!! warning:
When using noncoding peptides in a denylist, do not also pass the noncoding peptide FASTA as an input FASTA, because all peptides will be removed.

```bash
moPepGen fitlerFasta \
  --input-path path/to/variant_peptides.fasta \
  --output-path path/to/variant_peptides_filter.fasta \
  --denylist path/to/denylist.fasta
```

Use the `--keep-canonical` option to keep peptides that are called from canonical ORFs even if they are in the denylist. Canonical ORFs include coding transcripts with mutation(s) and fusion transcripts where the upstream transcript is coding. Peptides called from circRNAs are considered noncanonical ORFs.

```bash
moPepGen filterFasta \
  --input-path path/to/variant_peptides.fasta \
  --output-path path/to/variant_peptides_filter.fasta \
  --denylist path/to/denylist.fasta \
  --keep-canonical
```

### Complex Filtering

Sometimes we want a more complex filtering strategy. In the example below, we want to first remove any variant peptides that overlap with any noncoding peptides, and then filter again based on the expression level.

Remove variant peptides if they overlap with any noncoding peptide.

```path
moPepGen fitlerFasta \
  --input-path variant_peptides.fasta \
  --output-path variant_peptides_filter.fasta \
  --denylist noncoding.fasta
```

Filter again with both filtered variant peptides and noncoding peptides based on expression level.

```path
moPepGen fitlerFasta \
  --input-path \
    variant_peptides_filter.fasta \
    noncoding_peptides.fasta \
  --output-path combined_filter.fasta \
  --exprs-table expression.tsv \
  --delimiter '\t' \
  --tx-id-col 1 \
  --quant-col 4 \
  --quant-cutoff 2
```
