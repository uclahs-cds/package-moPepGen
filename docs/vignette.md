# Vignette

Welcome to the vignette for moPepGen, a powerful Python package designed for generating custom proteomic databases for library search. This vignette aims to provide you with a comprehensive overview and step-by-step guide on utilizing moPepGen to create custom databases for non-canonical peptide detection in mass spectrometry-based proteomics experiments.

## Installation

moPepGen is a command line tool designated to execute in a unix-like environment. For MacOS and Linux users, moPepGen can be installed using the command below. For Windows users, we recommend installing and running moPepGen from WSL (Windows Subsystem for Linux).

```shell
pip install git+ssh://git@github.com:uclahs-cds/private-moPepGen.git
```

## Reference Data

moPepGen requires a set of reference files, including the reference genome, its annotation and the translated protein sequences. We currently support reference files downloaded from two sources, ENSEMBL and GENCODE. See [here](quick-start/#downloading-reference-files) for more details.

A simulated reference set is provided for demonstration. The demo reference set only contains 5 transcripts and is only about 40 KB in size so should run very easily on any computer. The demo reference set can be downloaded with the commands blew.

```shell
cd ~
mkdir -p moPepGen-demo
cd moPepGen-demo
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/genome.fasta
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/annotation.gtf
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/translate.fasta
```

Convert reference set into index files for quick access of moPepGen.

```shell
moPepGen generateIndex \
    -g genome.fasta \
    -a annotation.gtf \
    -p translate.fasta \
    -o ./index
```

## Parsing

moPepGen starts with parsing a variety of variant files into GVF, a TSV format derived from VCF, to be used by moPepGen to call variant peptides.

### SNP & INDEL

Single nucleotide variants (SNVs/SNPs) and small insertions/deletions (INDELs) called by variant callers (*e.g.* GATK and Mutect2) must be annotated by the Variant Effect Predictor (VEP) first to get the genes each variant is associated with. Here is a generic command we use. Note that, the VEP cache files must be downloaded prior to running VEP (see [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html)). The cache file of the correct version number should be used, however when running VEP, we recommend also providing a custom reference genome and annotation downloaded from eiither ENCODE or ENSEMBL. The exact genome FASTA and annotation GTF files should be used later when calling for variant peptides.

```shell
vep \
    --offline \
    --cache \
    --check_ref \
    --no_stats \
    --fork ${N_THREADS} \
    --buffer_size 10000 \
    --distance 0 \
    --assembly GRCh38 \
    --no_intergenic \
    --chr 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19.20,21,22,X,Y,M \
    -i ${INPUT_BED_FILE} \
    -o ${OUTPUT_FILE} \
    --fasta ${GENOME_FASTA} \
    --custom ${ANNOTATION_GTF},${REFERENCE_VERSION},gtf
```

If you provided a custom GTF file to VEP, please also use a filter command to only keep the annotation records from the custom reference.

```shell
filter_vep \
    --force_overwrite \
    -i ${INPUT_FILE} \
    -o ${OUTPUT_FILE} \
    --filter "Source = ${REFERENCE_VERSION}"
```

Download the demo VEP file TSV file.

```shell
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/vep/vep_snp.txt
```

The output VEP TSV file must be parsed by parseVEP into GVF format.

```shell
moPepGen parseVEP \
    -i vep_snp.txt \
    --index-dir index \
    -o vep_snp.gvf \
    --source SNV
```

The `--source` argument is used to specify the type of variants (*e.g.*, SNP, SNV, INDEL) parsed. The source names are used in later processing steps, and are required by all moPepGen parsers.

### Fusion

moPepGen provides parsers to three fusion callers, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion), [Arriba](https://github.com/suhrig/arriba) and [FusionCatcher](https://github.com/ndaniel/fusioncatcher).

As an example, we provide a STAR-Fusion TSV output file for demonstration.

```shell
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/fusion/star_fusion.txt
```

Parse it into GVF format.

```shell
moPepGen parseSTARFusion \
    -i star_fusion.txt \
    --index-dir index \
    --source Fusion \
    -o star_fusion.gvf
```

Be default, `parseSTARFusion` only keeps fusion events with minimal `est_J` value of 5. This can be altered by the `--min-est-j` argument.

### Alternative Splicing

moPepGen accepts alternative splicing (AS) events estimated by [rMATS](https://rnaseq-mats.sourceforge.net/). RMATS estimates five AS events: SE (skipped exon), A3SS (alternative 3' splicing), A5SS (alternative 5' splicing), MXE (mutually exclusive exons), and RI (retained introns), accepted by `parseRMATS` as separate input channels. Noted that only the *.JC.txt files are supported.

Example data:

```shell
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/alternative_splicing/rmats_se_case_1.txt
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/alternative_splicing/rmats_a3ss_case_1.txt
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/alternative_splicing/rmats_a5ss_case_1.txt
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/alternative_splicing/rmats_mxe_case_1.txt
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/alternative_splicing/rmats_ri_case_1.txt
```

Parse AS events output by rMATS into GVF format with `parseRMATS`.

```shell
moPepGen parseRMATS \
    --se rmats_se_case_1.txt \
    --a3ss rmats_a3ss_case_1 \
    --a5ss rmats_a5ss_case_1 \
    --mxe rmats_mxe_case_1 \
    --ri rmats_ri_case_1 \
    --index-dir index \
    -o alt_splice_rmats.gvf \
    --source AltSplice \
```

By default `parseRMATS` only accepts AS events with inclusion and exclusion junction counts at least 1. These cutoffs can be set by `--min-ijc` and `--min-sjc`. See [here](./parse-rmats) for a complete list of arguments.

### RNA Editing Sites

RNA editing sites are specific positions within mRNA molecules where nucleotides undergo post-transcriptional modifications. moPepGen supports RNA editing sites called by [REDItools](https://github.com/BioinfoUNIBA/REDItools). Noted that the REDItools output must be annotated by the `AnnotateTable.py` from the REDItools package prior to passing into `parseREDItools`.

Example data:

```shell
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/reditools/reditools_annotated.txt
```

Parse REDItools output to GVF:

```shell
moPepGen parseREDItools \
    -i reditools_annotated.txt \
    -o rna_editing_reditools.gvf \
    --index-dir index \
    --source RNAEditing
```

By default `parseREDItools` looks for the transcript ID at column in column 17. This can be changed with `--transcript-id-column`, which takes a 1-based column number. See [here](./parse-reditools) for a complete list of arguments.

### CircRNA

CircRNAs are commonly recognized as noncoding RNAs, but evidence has showed that they are potentially translatable. moPepGen accepts circRNA events called by [CIRCexplorer](https://circexplorer2.readthedocs.io/en/latest/).

Download demo data:

```shell
wget https://github.com/uclahs-cds/private-moPepGen/raw/main/test/files/circRNA/CIRCexplorer_circularRNA_known.txt
```

Parse it to GVF format.

```shell
moPepGen parseCIRCexplorer \
    -i CIRCexplorer_circularRNA_known.txt \
    -o circRNA_CIRCexplorer.gvf \
    --index-dir index \
    --source CircRNA
```

By default `parseCIRCexplorer` accepts the text file output by CIRCexplorer2, however CIRCexplorer3 is also supported with the `--circexplorer3` flag. We also provide a series of filtering parameters that can be found [here](./parse-circexplorer).

## Non-canonical Peptides Calling

moPepGen provides three commands for non-canonical peptides calling. `callVariant` for calling peptides from variant GVFs, `callNoncoding` for doing 3-frame translation from noncoding transcripts, and `callAltTranslation` for calling peptide that harbor alternative translation events such as selenocysteine termination and W > F substants.

### Variant Peptides

Variant peptides, peptides that harbor any variant, can be called using the `callVariant` command. `callVariant` must take one or more GVF files produced by moPepGen parsers.

```shell
moPepGen callVariant \
    -i vep_snp.gvf star_fusion.gvf alt_splice_rmats.gvf rna_editing_reditools.gvf circRNA_CIRCexplorer.gvf \
    --index-dir index \
    -o variant_peptides.fasta \
    --threads 4
```

`callVariant` supports multi-processing and the number of processors to use can be specified with `--threads` argument. The `--selenocysteine-termination` and `--w2f-reassignment` arguments can be used to calling variant peptides that also carry selenocysteine termination and W2F reassignment. By default, `callVariant` uses trypsin as the enzyme for *in silico* digestion and allows up to 2 miscleavages, and this can be specified with `--cleavage-rule` and `--miscleavage`. See [here](./call-variant) for a complete list of arguments supported by `callVariant`.

### Noncoding Peptides

Noncoding peptides, peptids that could potentially be translated from novel open reading frames in transcripts that are annotated as noncoding, can be called using `callNoncoding`. Note that `callNoncoding` does not take any variants as input but only works with the reference set of noncoding transcripts. There is no need to rerun `callNoncoding` unless you wish to use a different enzyme or reference set.

```shell
moPepGen callNoncoding \
    --index-dir index \
    -o noncoding_peptides.fasta
```

Similar to `callVariant`, trypsin is the default enzyme and the default maximal miscleavages to allow is 2. These can be specified with `--cleavage-rule` and `--miscleavage`. See [here](./call-noncoding) for a complete list of arguments supported by `callNoncoding`.

### Alternate Translation Peptides

Alternative translation peptides are those that harbor special events during translation, such as selenocysteine termination and W > F substants, the the genetic code are not altered (see [here](./call-alt-translation) for more details). Similiar to noncoding peptides, `callAltTranslation` only calls for peptides from the reference genome and annotation and do not take any GVF file as input.

```shell
moPepGen callAltTranslation \
    --index-dir index \
    -o alt_trans_peptides.fasta
```

And again, `callAltTranslation` also uses trypsin as the default enzyme, and up to 2 miscleavages by default. See [here](./call-alt-translation) for a complete list of arguments.

## Processing

moPepGen provides a series of processing commands that aims to deliver FASTA files ready for database searching. The processing tasks include summarization of a non-canonical database, splitting a detabase to separate databases, creating decoy databases, shortening fasta headers for search engines to handle, and merging multiple database files for multiplexing proteomic experiments.

### Summarizing

The first step after having a list of variant peptides is usually inspecting what you have. `summarizeFasta` takes a variant peptide fasta file output by `callVariant` and summarize the variant peptides by categories based on the value of `--source` you input when calling the parser commands. `summarizeFasta` must take all GVFs used when calling `callVariant`.

```shell
moPepGen summarizeFasta \
    --gvf vep_snp.gvf star_fusion.gvf alt_splice_rmats.gvf rna_editing_reditools.gvf circRNA_CIRCexplorer.gvf \
    --variant-peptides variant_peptides.fasta \
    --index-dir index
```

By default, `summarizeFasta` outputs to stdout, which can be saved to a text file using the `-o` argument. In case `--output-image` is given, a horizontal bar plot will be saved.

![summarize-fasta-bar-plot](img/summarize-fasta.png)

Because moPepGen calls for enzymatic cleaved peptides, there are chances that the same peptide can be called from different transcripts, or the same transcript with different or different combinations of variants. For example, the peptide below is called twice from two separate transcripts with different variants.

```
>ENST00000622235.5|SNV-100-G-T|4 ENST00000614167.2|RES-202-G-A|2
HETLFLLTFPR
```

To resolve the issue of collapsed peptides like this, we use the `--order-source` argument that takes the priority order of sources to be considered. It takes a comma separated format. For example `--order-source gSNP,RNAEditing` will prioritize gSNP over RNA editing events, thus the example peptide above will be assigned to the gSNP category. Noted that the values passed into `--order-source` must match with the values used for `--source` in the corresponding parser calls. If the `--order-source` is not provided, the source priority order will be inferred from the order of input GVF files.

Besides variant peptides called by `callVariant`, noncoding peptides and alternative translation peptides can also be passed to `summarizeFasta` with `--noncoding-peptides` and `--alt-translation-peptides`.

### Splitting

The `splitFasta` is provided to split a variant peptide database into several separate databases for tiered database searching.

```shell
mkdir -p split
moPepGen splitFasta \
    --variant-peptides variant_peptides.fasta \
    --gvf vep_snp.gvf star_fusion.gvf alt_splice_rmats.gvf rna_editing_reditools.gvf circRNA_CIRCexplorer.gvf \
    --output-prefix split/split \
    --index-dir index
```

Similar to `summarizeFasta`, `splitFasta` also takes a `--order-source` to specify the priority order of which category a peptide should be assigned to, and will be inferred from the input GVFs if not specified. `--group-source` is used to group sources as a super category. For example, `--group-source Germline:gSNP,gINDEL Somatic:sSNV,sINDEL` will group sources of `gSNP` and `gINDEL` together as `Germline`, and `sSNV` and `sINDEL` as `Somatic`.

Noted that, when assigning a peptide to a source category, it must carry exclusively the desired type(s) of variants. For example, a peptide of 'ENST00000622235.5|SNV-100-G-T|SNV-110-C-A|2' is assigned to `SNV`, while a peptide of 'ENST00000622235.5|SNV-100-G-T|RES-110-C-A|2' will be assigned to the category of `SNV-RNAEditing` but not `SNV`. `--max-source-groups` is used to specify the maximal number of source groups to split. The default value is 1, which means all peptides that contains two or more types of variants will not be written into their own FASTA file, but kept in the '\<prefix\>_Remaining.fasta' file.

Similar to `summarizeFasta`, noncoding and alternative translation peptides can be passed to `splitFasta` via `--noncoding-peptides` and `--alt-translation-peptides`.

See [here](./split-fasta) for a complete list of arguments.

### Target-Decoy Database

Most search engines expect a target-decoy database as input to estimate false discovery rate (FDR). We provide a `decoyFasta` command, that takes a variant peptide database and add decoy sequences with either `reverse` or `shuffle` algorithm.

```shell
moPepGen decoyFasta \
    -i split/split_gSNP_encode.fasta \
    -o split/split_gSNP_decoy.fasta \
    --method reverse \
    --non-shuffle-pattern K,R
```

`--non-shuffle-pattern` specifies the amino acid residues that should be fixed in the decoy sequence. By default, the N- and C-terminal residues are also fixed, which can be turned off by setting `--keep-peptide-nterm` or `--keep-peptide-cterm` to `false`. See [here](./decoy-fasta) for a complete list of arguments.

### Shortening FASTA headers

Some search engines have limits on how long the database FASTA headers can be. The headers of moPepGen's variant FASTA files could be very long, because the same peptide can be called from different transcripts. We provide a shortening approach by replacing all FASTA headers with UUIDs and storing the mapping information into a `.dict` file with the same prefix.

```shell
moPepGen encodeFasta \
    -i split/split_gSNP.fasta \
    -o split/split_gSNP_encode.fasta
```

Noted that for decoy peptides, the same UUID will be used as their corresponding target sequences, with the decoy prefix/suffix retained. The resulted `.dict` file can be used to map back to the original FASTA header.
