# File Structure Documentation
  
- [File Structure Documentation](#file-structure-documentation)
  - [Variant File](#variant-file)
    - [File header](#file-header)
    - [Point Mutation](#point-mutation)
    - [Structural Mutation](#structural-mutation)
  - [Reference Index](#reference-index)
  - [Variant Peptide FASTA](#variant-peptide-fasta)

This documentation describes the design of the intermediate files, including index files and the BED-like intermediate variant file, and the file FASTA output.


## [Variant File](#variant-file)

Variant files are BED-like files to hold the variant information.

### File header

All variant file should contain a file header section as below.

```
##MOPEPGEN_VERSION=0.0.1
##PARSER=parseXXX
##VARIANT_TYPE=SNP
##REFERENCE_INDEX=
##GENOME_FASTA=
##ANNOTATION_GTF=
##PROTEOME_FASTA=
##COMMAND=
```

Headers are key value paires that hold the metadata for this file.

+ `MOPEPGEN_VERSION`: the version number of moPepGen used to generate the variant file.
+ `PARSER`: the parser used, 

### Point Mutation

```
#transcript_id  start   end ref alt type    id
```

### Structural Mutation

```
#transcript_id  start end doner_transcript_id docker_start  doner_end
```

## [Reference Index](#reference-index)


## [Variant Peptide FASTA](#variant-peptide-fasta)


---

+ gene fusion: gene