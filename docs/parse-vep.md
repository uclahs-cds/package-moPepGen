{% set command = 'parseVEP' %}
# {{ command }}

::: moPepGen.cli.parse_vep
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

## Running VEP for moPepGen

moPepGen relies on Ensembl Variant Effect Predictor ([VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)) to annotate single nucleotide and small insertion and deletion (indel) variants, typically presented in VCF files. Note that VEP supports the annotation of multi-nucleotide variants and structural variants, but those are not currently supported by moPepGen and will be discarded by `parseVEP`.

Please refer to the official [VEP Tutorial](https://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html) for instructions on downloading, installing and running VEP. The key is to ensure that the annotation version used in `VEP` is the same as the one you would like to use with moPepGen.

`parseVEP` currently only supports the TSV output of `VEP`, please set the suffix of `--output_file` to `.tsv` to ensure the correct format is outputted by `VEP`.

### Using an Ensembl GTF

`VEP` automatically uses an Ensembl GTF for annotation. Please ensure that the release version and genome build version is consistent with the references files downloaded. This is available in the name of the cache in output headers (for example, `107_GRCh38`).

One can ensure that the correct version is used by downloading the cache following [instructions](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) and specifying the parameters:
`--species homo_sapiens --assembly GRCh38 --cache --cache_version 107`

### Using a Non-Ensembl GTF

If you have decided to use a set of reference files from GENCODE, it is important to supply VEP with the non-Ensembl GTF during annotation, so that chromosome names, transcript IDs and transcript coordinates match with your intended reference files. This can be specified in `VEP`.

As [instructed](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff) by `VEP`, your GTF file must be sorted in chromosomal order and indexed.
```
grep -v "#" PATH_TO_GTF | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > PATH_TO_GTF_GZ
tabix -p gff PATH_TO_GTF_GZ # gtf is not a tabix format option, gff works
```

To use the GTF for annotation, run VEP with the additional parameters
```
--custom PATH_TO_GTF,GENCODE,gtf --fasta PATH_TO_GENOME_FA
```

Followed by VEP, we run the `filter_vep` command to subset outputs to only those annotated using the GENCODE GTF, since VEP automatically outputs both native (Ensembl) and custom annotations.

```
filter_vep -i VEP_OUTPUT_TSV -o FILTERED_VEP_OUTPUT_TSV --filter Source = GENCODE
```

### Tips for Running VEP

We recommend the following parameters for running VEP for run time optimization, please select appropriate settings for your system.
- `--offline --cache`
- `--no_stats`
- `--fork`
- `--buffer_size`
- `--distance 0` - only annotation variants directly overlapping gene coordinates
- `--chr 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT` - focus on primary assembly
- `--no_intergenic` - ignore intergenic variants

We do not recommend adding any additional annotation flags since the information is not utilized by moPepGen and increases run time.
