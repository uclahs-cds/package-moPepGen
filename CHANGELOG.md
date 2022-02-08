# Changelog
All notable changes to the tool_name Docker file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added

- Enable `filterFasta` to filter by number of miscleavages per peptide. #382

- Added CLI command `mergeFasta` to merge multiple variant peptide database Fasta files into one. This could be useful when working on multiplexed proteomics experiments such as TMT. [#380](https://github.com/uclahs-cds/private-moPepGen/issues/380)

### Changed

- Donor and accepter transcript IDs are now explicitly included in the variant IDs of fusion in both GVFs and variaint peptide FASTA headers. Closed #376 via #377

- For fusion, `callVariant` now looks at the entire accepter sequence for potential variant peptides, rather than only the peptides that contains the breakpoint.

---

## [0.2.0] - 2021-01-28

### Added

- Multi-threading is enabled for `callVariant` to run in parallel.

- CLI command `indexGVF` added to generate a index file for quickly access variant data from the corresponding GVF file. Noted that this command is not required to run.

### Changed

- To solve the complexity of subgraphs introduced by fusion and especially alternative splicing insertion and substitution, the `SubgraphTree` class is added to keep the graph-subgraph relationship between nodes.

- Variant records are now kept on disk rather than reading the entire GVF file(s) into memory, and only the file pointers to variant records are kept in memory. This significantly reduces the memory usage of `callVariant`.

- The command line arguments are standardized across all commands, for example '-i/--input-path' for inputs and '-o/--output-path' for outputs.

---

## [0.1.0-beta.1] - 2021-12-23

### Added

- Initial beta release of moPepGen, with the three-frame graph based algorithm implemented to call noncanoinical peptides.
