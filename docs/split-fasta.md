{% set command = 'splitFasta' %}
# {{ command }}

::: moPepGen.cli.split_fasta
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

## Examples

### Basic usage

A basic usage is below.

```bash
moPepGen splitFasta \
  --gvf \
    path/to/gSNP.gvf \
    path/to/gINDEL.gvf \
    path/to/reditools.gvf \
  --variant-peptides path/to/variant.fasta \
  --index-dir path/to/index \
  --max-source-groups 1 \
  --output-prefix path/to/split
```

The example above splits the variant peptide sequence database into four FASTA files, `split_gSNP.fasta`, `split_gINDEL.fasta`, `split_RNAEditing.fasta`, `split_Remaining.fasta`. Variant peptides are split into individual FASTA file of its variant source group based on the order of GVF files. Peptides with more than one variant sources are written into the `*_Remaining.fasta` because the `--max-source-groups` is set to 1.

### Group sources

Sometimes we want to group variant sources together. See example below.

```bash
moPepGen splitFasta \
  --gvf \
    path/to/gSNP.gvf \
    path/to/gINDEL.gvf \
    path/to/reditools.gvf \
  --variant-peptides path/to/variant.fasta \
  --index-dir path/to/index \
  --group-source Coding:gSNP,gINDEL \
  --order-source Coding,RNAEditing \
  --max-source-groups 1 \
  --output-prefix path/to/split
```

This example outputs three split FASTA files, `split_Coding`.fasta`, `split_RNAEditing.fasta`, and `split_Remaining.fasta`.

### Additional split

Additional split allows you to split the records with the source group specified that would otherwise be placed along with all other records exceeding `--max-source-groups` into `remaining.fasta`. See the example below.

```bash
moPepGen splitFasta \
  --gvf \
    path/to/gSNP.gvf \
    path/to/gINDEL.gvf \
    path/to/reditools.gvf \
  --variant-peptides path/to/variant.fasta \
  --index-dir path/to/index \
  --max-source-groups 1 \
  --additional-split gSNP-gINDL gSNP-RNAEditing \
  --output-prefix path/to/split
```

As a result, `split_gSNP.fasta` and `split_gINDEL.fasta` will be written. `split_gSNP-gINDEL.fasta` and `split_gSNP-RNAEditing.fasta` will also be written although the number of variant sources (2) is larger than the value specified through `--max-source-groups`.

## Arguments

{% with actions=get_arg_data(command) %}
{% include 'partials/_command_usage.md' %}
{% endwith %}
