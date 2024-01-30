# <u>M</u>ulti-<u>O</u>mics <u>Pep</u>tide <u>Gen</u>erator

<!-- badges: start -->

[![Lifecycle:Maturing](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)
[![Tests](https://github.com/uclahs-cds/private-moPepGen/actions/workflows/tests.yaml/badge.svg)](https://github.com/uclahs-cds/private-moPepGen/actions/workflows/tests.yaml)
[![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=plastic&logo=docker&logoColor=white)](https://github.com/uclahs-cds/private-moPepGen/pkgs/container/mopepgen)
[![Documentation](https://img.shields.io/static/v1?style=plastic&message=ReadTheDocs&color=2C4AA8&logo=ReadTheDocs&logoColor=FFFFFF&label=Documentation)](https://uclahs-cds.github.io/private-moPepGen/)
[![License](https://img.shields.io/badge/License-GPL_v2-blue)](./LICENSE.txt)

<!-- badges: end -->

moPepGen (multi-omics peptide generator) uses data from one or more omics experiments and calls variant peptides as custom databases for proteogenomic library search.

moPepGen takes genomic and transcriptomic variants such as single nucleotide variants (SNPs or SNVs), small insertions/deletions (indels), gene fusion, alternative splicing, RNA circularization and RNA editing events, and generates noncanonical peptides affected by the variants.

## Installation

Install directly from GitHub

```shell
pip install git+ssh://git@github.com/uclahs-cds/private-moPepGen.git
```

## Documentation

See [here](https://uclahs-cds.github.io/private-moPepGen/index.html)
