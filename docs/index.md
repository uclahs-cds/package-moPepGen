# moPepGen: <u>M</u>ulti-<u>O</u>mics <u>Pep</u>tide <u>Gen</u>erator

<!-- badges: start -->

[![Tests](https://github.com/uclahs-cds/private-moPepGen/actions/workflows/pull_request.yaml/badge.svg)](https://github.com/uclahs-cds/private-moPepGen/actions/workflows/pull_request.yaml)
[![Docs](https://github.com/uclahs-cds/private-moPepGen/actions/workflows/deploy.yaml/badge.svg)](https://github.com/uclahs-cds/private-moPepGen/actions/workflows/deploy.yaml)

<!-- badges: end -->

moPepGen (multi-omics peptide generator) uses data from multiple -omics experiments and calls variant peptides as custom database for proteogenomics library search.

moPepGen takes genomic variants such as single nucleotide variants (SNP or SNV), insertion/deletion (INDEL), gene fusion, and post transcriptional modifications such as RNA editing and alternative splicing, and detects variated peptides affected.

## Installation

Install directly from github

```
pip install git+ssh://git@github.com:uclahs-cds/private-moPepGen.git
```

## Overview

![graphic-abstract](img/graphic-abstract.png)