## Codon Table & Start Codon

### Codon Table

The NCBI standard codon table is used by default, which is used for the majority of nulcear gene translation in eukaryote cells. See [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for a complete list of NCBI codon tables.

The default codon table can be override using `--codon-table`. For example:

```shell
--codon-table ’Ciliate Nuclear'
```

The `--chr-codon-table` can be used to specify the codon table used for a specific chomosome. The example below uses the 'Vertebrate Mitochondrial' (SGC1) codon table for genes from the mitochondria chomosome, and the standard codon table otherwise.

```shell
--codon-table Standard \
--chr-codon-table 'chrM:SGC1'
```

### Start Codons

Stard codons usually do not need to be specified. The standard start codon `ATG` is used by default, and it is translated as Methionine as start codon and in elongation. However, in some cases, for example, mitochondria, `ATA` and `ATT` may also be used as start codon. While `ATT` is translatted into Isoleucine during elongation, Methionine is still used as start codon.

Similar to codon table, the default codon table can be override using `--start-codons`.

```shell
--start-codons ATG
```

The `--chr-start-codon` can also be used to assign start codons to a specific chomosome. The example below assigns `ATG`, `ATA`, and `ATT` to the mitochondrial chromosome.

```shell
--chr-start-codons 'chrM:ATG,ATA,ATT'
```

### Default

The chromosome names must be specified correctly, same as what used in the genome fasta and annotation GTF file. By default, moPepGen infers the reference source of the annotation (*i.e.*, GENCODE or EMSEMBL), and uses the 'SGC1' codon table for mitochondirla chromosome. So the default is equivalent to:

```shell
--reference-source GENCODE \
--codon-table Standard \
--chr-codon-table 'chrM:SGC1' \
--start-codons 'ATG' \
--chr-start-codons 'chrM:ATG,ATA,ATT
```

or

```shell
--reference-source ENSEMBL \
--codon-table Standard \
--chr-codon-table 'MT:SGC1' \
--start-codons 'ATG' \
--chr-start-codons 'MT:ATG,ATA,ATT
```
