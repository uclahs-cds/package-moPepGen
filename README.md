# <u>M</u>ulti-<u>O</u>mics <u>Pep</u>tide <u>Gen</u>erator

<!-- badges: start -->

[![Lifecycle:Stable](https://img.shields.io/badge/Lifecycle-Stable-97ca00)](https://img.shields.io/badge/Lifecycle-Stable-97ca00)
[![Tests](https://github.com/uclahs-cds/package-moPepGen/actions/workflows/tests.yaml/badge.svg)](https://github.com/uclahs-cds/package-moPepGen/actions/workflows/tests.yaml)
![PyPI - Version](https://img.shields.io/pypi/v/mopepgen?logo=Python)
[![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=plastic&logo=docker&logoColor=white)](https://github.com/uclahs-cds/package-moPepGen/pkgs/container/mopepgen)
[![Documentation](https://img.shields.io/static/v1?style=plastic&message=ReadTheDocs&color=2C4AA8&logo=ReadTheDocs&logoColor=FFFFFF&label=Documentation)](https://uclahs-cds.github.io/package-moPepGen/)
[![License](https://img.shields.io/badge/License-GPL_v2-blue)](./LICENSE.txt)
[![DOI:10.1101/2024.03.28.587261](https://zenodo.org/badge/DOI/10.1101/2024.03.28.587261.svg)](https://doi.org/10.1101/2024.03.28.587261)

<!-- badges: end -->

moPepGen (multi-omics peptide generator) uses data from one or more omics experiments and calls variant peptides as custom databases for proteogenomic library search.

moPepGen takes genomic and transcriptomic variants such as single nucleotide variants (SNPs or SNVs), small insertions/deletions (indels), gene fusion, alternative splicing, RNA circularization and RNA editing events, and generates noncanonical peptides affected by the variants.

## Installation

Install from PyPI

```shell
pip install mopepgen
```

Install directly from GitHub

```shell
pip install git+ssh://git@github.com/uclahs-cds/package-moPepGen.git
```

## Documentation

See [here](https://uclahs-cds.github.io/package-moPepGen/index.html)

## Citation

1. Zhu, C. et al. moPepGen: Rapid and Comprehensive Proteoform Identification. 2024.03.28.587261 Preprint at https://doi.org/10.1101/2024.03.28.587261 (2024).
