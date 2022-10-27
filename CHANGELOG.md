# Changelog

All notable changes to the tool_name Docker file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [0.10.1] - 2022-10-27

### Fixed

- Transcriptional coordinate to genomic coordinate not converted successfully when it is the last nucleotide of the transcript. #592

- When creating the peptide cleavage graph, the end nodes of ate variant bubble with alt splice were collapsed with the reference node causing the graph cleavage process terminated too early resulting uncleaved nodes. #597

- in `callVariant` when filtering variants associated with the donor transcript, the left breakpoint coordinate not converted successfully if it is the end of the transcript. #598

## [0.10.0] - 2022-10-20

### Added

- Added support for fusion, alternative splicing and circRNA in `bruteForce`.

### Fixed

- Several issues of `bruteForce` were fixed for fusion, alternative splicing and circRNA to be consistent with `callVariant`.

- In `ThreeFrameTVG` and `PeptideVariantGraph`, large deletions (for alternative splicing) are no longer treated as subgraphs any more.

- Fixed the issue that the `subgraph_id` attributes of `TVGNode` and `PVGNode` are lost after nodes are merged. #566

- When expanding the aligned variant bubble, if the downstream node of the start node has multiple inbond nodes, nucleotides will be taken from the downstream node and added to each upstreams

- Fixed `callVariant` that when filtering variants that are compatible with fusion, the breakpoint site were not recognized correctly. #567

- For `ThreeFrameTVG', when aligning variant bubbles, if the end of the first variant is the start of the next (e.g. alternative splicing events that sharing the same splicing site), the merged bubble will then contain both variants

- Fixed `callVariant` that mutations are assigned as stop altering mutation when there is a start codon after it. #568

- Fixed `callVariant` that alternative splicing variants were not recognized as stop altering mutation correctly because their reference sequence from GVF is only the first nucleotide. #569

- Fixed `callVariant` that nodes being lost after an in-frame subgraph. #573

- Fixed `callVariant` that the actual fusion breakpoint was not found correctly when trying to tell whether a novel start site should be considered.

- Fixed `callVariant` that variant peptides were called with variants present in one loop but not in another. #576

- When collapsing the end nodes when creating the peptide cleavage graph, nodes that contains alt splice deletions are now separated from others. #580

- Fusion not inserted correctly when the breakpoint is intronic. #578

- When finding start altering variants from a node, wrong right position was used. #583

## [0.9.4] - 2022-09-07

### Fixed

- Fixed issue of alternative splicing deletion that starts at the third nucleotide of start codon. Those variants are now skipped. #560

## [0.9.3] - 2022-08-26

### Fixed

- Fixed issue that mouse reference genome and proteome were not parsed correctly. #546

- Fixed issue that stop-lost inframe deletions being missed during node collapsing. #549

- Fixed issue that variants missed by callVariant when multiple variants are causing the same sequence. #552

- Fixed issue that `cpop_collapsed` attribute was not retained after merging so peptides that don't end with cleavage sites were yield. #554

- Fixed problem caused by N in the reference DNA sequence. #556

---

## [0.9.2] - 2022-07-29

### Fixed

- For circRNA, each reading frame subgraph is now replicated for 3 times in order to catch all variant peptides that read through the junction site. #514

- For `decoyFasta`, overlapping decoy sequences are also counted and printed to the stdout when using reverse. #474

### Added

- Enzyme lysN is added. #523

---

## [0.9.1] - 2022-07-27

### Added

- UTIL command line tool `validateNoncodingCalling` added to validate the output of `callNoncoding` with `bruteForceNoncoding`. #524

## Fixed

- `parseVEP` updated to only parse the first 13 columns of the VEP TSV. #540

---

## [0.9.0] - 2022-07-18

### Fixed

- Fixed issue that peptides caused by start gain mutations are not called by `bruteForce`. #492

- Fixed issue that variant peptides caused by start gain mutations from noncoding transcripts are not called in some cases. Stop gain mutations from downstream are also missed in some cases. #495

- Fixed issue that some variant peptides were not called with in-frame deletions. #515

- Fixed issue that the shuffled decoy sequences produced by decoyFasta are not reproducible when input FASTA order is changed. #517

- Fixed issue that stop lost mutations are not recognized correctly. #519

- Fixed issue that start gain and stop lost mutations before the novel start site were not excluded. #520

- Fixed issue that frameshifting mutations on the same node of a novel ORF start site and is right after it was not carried over to downstream nodes. #526

- Fixed issue that stop lost mutations were not recognized for deletions. #527

- Fixed issue that in-frame deletion stop retaining mutations were not recognized #528

- Fixed issue that variant coordinates not handled correctly when the cleavage site is contained in the inserted sequence. #52

- Fixed issue that peptides are not called with overlapping deletions. #531

---

## [0.8.1] - 2022-06-24

### Fixed

- The attribute `global_variant` is added to the `PeptideVariantGraph` so for circRNA, if the entire variant peptide is called from an insertion sequence, the variant of circRNA is still going to be added. So the FASTA header will be written out correctly. #490

## [0.8.0] - 2022-06-22

### Fixed

- Fixed issue that variants are not filtered correctly on `ThreeFrameCVG` (some variants are discarded or retained incorrectly). #488

---

## [0.7.2] - 2022-06-20

### Changed

- Rolling back to not filtering candidate variant peptides that overlaps with canonical peptide pool at transcript level in `PeptideVariantGraph` because accessing the shared memory object is too slow (which could be optimized in the future).

---

## [0.7.1] - 2022-06-18

### Changed

- `callVariant` optimized so the memory usage and run time is reduced to better handle hyper-mutated regions. #483

- Switched from [pathos](https://pathos.readthedocs.io/en/latest/index.html) to [ray](https://www.ray.io/) because the latter has better support for shared memory.

---

## [0.7.0] - 2022-06-08

### Fixed

- Fixed the issue of long run time for transcripts with hypermutated regions by limiting number of additional variants per miscleavage (`--additional-variants-per-misc`). Closed #476 via #477

- Location was missing in variant labels output by `parseREDItools`. #478

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
