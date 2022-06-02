# Changelog

All notable changes to the tool_name Docker file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [0.6.1] - 2022-05-31

### Fixed

- Fixed the issue that large insertions, such as retained introns, that start right after the start codon were not converted to end inclusion successfully. #470

## [0.6.1] - 2022-05-31

### Fixed

- Deletion was mistreated as insertion when trying to convert to end-inclusion format. This only affects deletions that start on the third nucleotide of CDS. #468

---

## [0.6.0] - 2022-05-28

### Fixed

- `callVariant` very slow on certain cases when the transcript is large. #464

- `callVariant` failed to find start codon for transcript with fusion being the first variant and the donor breakpoint is intronic. #465

---

## [0.5.1] - 2022-05-22

### Added

- An option `--denylist` is added to `filterFasta` to accept a FASTA file to exclude peptide sequences from it.

---

## [0.5.0] - 2022-05-13

### Changed

- Fixed `callVariant` that it failed when the breakpoint of a fusion is right at the end of the start codon because it tried convert to end inclusion. #454

- Fixed `callVariant` that no variant peptides are called when the first variant is a fusion. #454

---

## [0.4.2] - 2022-05-11

### Changed

- A warning is raised in `parseVEP` when tryping to parse a MNV (multi-nucleotide variant) and skip the record instead of raising an error. #447

- Fixed `summarizeFasta` that source order of `Noncoding` was not recognized. #449

---

## [0.4.1] - 2022-04-27

### Changed

- Fixed the problem that in `summarizeFasta` output the order of variant sources in the same group is not consistent across runs. #428

- Argument `--ignore-missing-source` added to `summarizeFasta` so sources not present in any GVF can be ignored without raising any error. #436

- In `filterFasta`, when filter with expression table, changed to filter out peptides smaller than, instead of smaller or equal to, the value of `--quant-cutoff`.

- Fixed the issue that in `splitFasta`, variant sources are not grouped as they are specified by `--group-source` #439

### Added

- Resources usage including memory, CPU and time is now printed to stdout in the end of all command line programs.

### Fixed

- Fixed issue that `--additional-split` not recognized properly in `splitFasta`. #443

---

## [0.4.0] - 2022-03-17

### Added

- Added CLI command `summarizeFasta` to output a summary table of the variant peptide FASTA file output by `callVariant`.

### Changed

- Attribute key for transcript ID is fixed from 'TRANSCRIPT' to 'TRANSCRIPT_ID' in circRNA's GVF files output by `parseCIRCExplorer` to be the same as other GVF files.

- Genomic position for each record is added to the GVF file output by `parseCIRCExplorer`.

---

## [0.3.1] - 2022-03-01

### Changed

- Argument parameter `--decoy_string_position` is changed to `--decoy-string-position`. #417

---

## [0.3.0] - 2022-02-24

### Added

- Enable `filterFasta` to filter by number of miscleavages per peptide. #382

- Added CLI command `mergeFasta` to merge multiple variant peptide database Fasta files into one. This could be useful when working with multiplexed proteomic experiments such as TMT. [#380](https://github.com/uclahs-cds/private-moPepGen/issues/380)

- Added CLI command `decoyFasta` to generate decoy database by shuffling or reversing each sequence. [#386](https://github.com/uclahs-cds/private-moPepGen/issues/386)

- Added parameter `--min-coverage-rna` to `parseREDItools` to filter by total RNA reads at a given position. #392

- Added CLI command `encodeFasta` to replace the variant peptide headers with UUIDs. The original FASTA headers are stored in a text file together with the UUIDs. This is to make the FASTA header short enough for library search engines. #389

### Changed

- Donor and accepter transcript IDs are now explicitly included in the variant IDs of fusion in both GVFs and variaint peptide FASTA headers. Closed #376 via #377

- For fusion, `callVariant` now looks at the entire accepter sequence for potential variant peptides, rather than only the peptides that contains the breakpoint. #377

- `filterFasta` updated to support filter by number of miscleavages. #383

- In `parseVEP`, chromosome seqname for each record is now read directly from the gene annotation, to avoid the 'chr' prefix issue. #391

- The `--transcript-id-column` parameter of `parseREDItools` is changed to take 1-based index. #392

- Changed `splitDatabase` to `splitFasta` for consistency. #397

- Updated `generateIndex` to reduce the size of genomic annotation data and the memory usage when loaded. #395

---

## [0.2.0] - 2022-01-28

### Added

- Multi-threading is enabled for `callVariant` to run in parallel.

- CLI command `indexGVF` added to generate a index file for quickly access variant data from the corresponding GVF file. Noted that this command is not required to run.

### Changed

- To solve the complexity of subgraphs introduced by fusion and especially alternative splicing insertion and substitution, the `SubgraphTree` class is added to keep the graph-subgraph relationship between nodes.

- Variant records are now kept on disk rather than reading the entire GVF file(s) into memory, and only the file pointers to variant records are kept in memory. This significantly reduces the memory usage of `callVariant`.

- The command line arguments are standardized across all commands, for example '-i/--input-path' for inputs and '-o/--output-path' for outputs.

- `generateIndex` is changed to use compressed text format to store genomic annotation, because for some reason that we are not sure, when loading the pickled genomic annotation, the memory usage is almost doubled. [#394](https://github.com/uclahs-cds/private-moPepGen/issues/394)

---

## [0.1.0-beta.1] - 2021-12-23

### Added

- Initial beta release of moPepGen, with the three-frame graph based algorithm implemented to call noncanoinical peptides.
