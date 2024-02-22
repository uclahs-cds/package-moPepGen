# Changelog

All notable changes to the tool_name Docker file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [1.2.2] - 2024-1-16

### Fixed:

- Adjacent variants were not merged as MNVs successfully. The function always exited with nothing.

- Because of the updating to on-disk GTF, the coding transcripts were not generated and saved successfully. `filterFasta` is the only command affected.

- Updated `splitFasta` and `summarizeFasta` to accept source combinations in `--order-source`.

- Fixed `parseCIRCexplorer` so the exon/intron indices in variant IDs are sorted correctly.

- Fixed `parseVEP` to handle insertions in start-inclusion. #840

- Fixed `callVariant` of 'no reference out node found'. #842

- Fixed `mergeFasta` to remove redundant FASTA header entries. #846

## [1.2.1] - 2023-10-05

### Add

- Added `--graph-output-dir` to save graph data in json.

- Added `--timeout-seconds` to callVariant.

## Fixed

- Fixed `summarizeFasta` that SEC and W2F on fusion peptides are ignored. #789

- Fixed `callVariant` that `variant_coordinates_to_gene` failed when the deletion is end inclusion and it overlaps with the last nucleotide of an exon. #793

- Fixed `splitFasta` that CodonReassign and SECT were not able to be grouped. #796

- Fixed `callVariant` that in-frame subgraphs not recognized when they are not in variant bubble.

- Fixed `callVariant` that peptides are falsely called if the last miscleaved node is missing a downstream cleavage altering variant. #800

- Fixed `TVGNode` that `get_max_subgraph_id` always returns the last subgraph ID. #802

- Fixed `callVariant` with altSplice insertion with intronic frameshift variant which is very closed to the end of the subgraph. #803

- Fixed `bruteForce` that accepter variants are skipped if the donor transcript has variants with the same coordinate. #810

- Fixed `splitFasta` that source and source group order gets overriden by GVF order. #805

- Fixed `summarizeFasta` and `splitFasta` being too slow. #795

- Fixed `splitFasta` to use top priority header for additional split

- Fixed `callVariant` that some accepter only ORFs maybe included for noncoding fusion transcripts.

- Fixed `callVariant` that when filtering variants for a given transcript/fusion/circRNA, coordinates of end inclusion insertions were not interpreted correctly.

- Fixed `bruteForce` that selenocysteine not fixed for deletion alt sequence.

- Fixed `callVariant` that nodes with fusion being treated as subgraph out or end node incorrectly.

- Fixed `callVariant` that upstream cleavage altering variants were affecting checks for whether a node is hybrid in circRNA.

- Fixed `callVariant` that two SNV at the same location in circRNA was affecting hybrid node identification.

- Fixed `callVariant`. When creating the cleavage graph, when a variant bubble is processed, the downstream node(s) needs to be identified for the next iteration, and only in-frame node should be used. However some nodes can span over two reading frames, so we should check the last reading frame index instead of the first.

- Fixed `callVariant` that accepter transcript variants very closed to the breakpoint were skipped.

## Added

- Added support for `--group-source` for `summarizeFasta`. #798

## [1.2.0] - 2023-08-03

### Fixed

- Fixed that reference source are not recognized. Switched to use upper case values. #758

- Fixed `generateIndex` that symlink was not created properly for GTF with `--gtf-symlink`.

- Fixed matplotlib warning message #742

- Fixed `parseVEP` that insertions not parsed successfully if the location is a single base. #766

- Fixed `callVariant` that peptides with upstream cleavage altering mutations were not called. #670

- Fixed fuzzTest to take multiple CPUs and use temporary directory.

- Fixed issue in callVariant of circRNA with a SNV being silent in the first loop but not in the second.

- Fixed issue that circRNA peptide nodes with and without a silent mutation were collapsed. #778

- Fixed that circRNA peptides that carry indels incompatible with the orf start node should not be called. #780

- Fixed that hybrid node was not identified if the variant is in the end of an exon. #782

- Fixed that in circRNA, cleavage gain from upstream node is added to the the wrong ORF. #783

- Fixed callVariant that `fit_into_codon` terminated early in fusion transcripts when there donor has a frameshift and accepter has a variant right after the breakpoint. #786

## [1.1.0] - 2023-06-28

### Fixed

- Reduced memory usage by mapping the GTF files into memory instead of reading it all at once. #371

## [1.0.0] - 2023-06-15

### Added

- Added the support for calling peptides of Selenocysteine terminated. #684

- Added the support for calling peptides of W > F codon reassignments. #484

- Added the support for calling peptides with adjacent SNP/INDEL. #691

- Updated fake.py and bruteForce to handle selenocysteine, W2F and MNV. #689

- `callAltTranslation` added to call peptides with alternative translation without any genomic or transcriptomic variations.

- Enabled `summarizeFasta` to create bar plot of the summary results.

### Fixed

- Fixed `fake` that simulated selenocysteine positions could be in introns.

- Fixed `fake` that the last exon was picked for A3SS or first exon for A5SS.

- Fixed fusion with very small intronic insertion. #707

- In ThreeFrameTVG when aligning variant bubbles and when nodes are merged, variants were not merged correctly.

- Fixed TVG that indel merged with downstream fusion treated as subgraph out. #708

- Fixed `parseRMATS` to handle more complex situations such as exons interjacent between splicing sites and exons spanning over the splicing site. #715, #716, #717, and PR #720

- Fixed `callVariant` that failed when there is a SNV very close to the end on a AltSplice insertion. #723

- Fixed `TranscriptAnnotationModel` for not recognizing transcripts with `mRNA_end_NF` correctly. #724

- Fixed `callVariant` issue of altSplice insertion carries an intronic indel that goes back to the original reading frame. #726

- Fixed `callVariant` to handle deletion that spans over an entire intron. #732

- Fixed `callVariant` to skip peptides earlier if they are either too long or too short to significantly improve efficiency. #736

- Fixed `callVariant` to handle hypermutated region with a dynamic cutoff. #738

- Fixed `decoyFasta` to make it as default to keep cleavage site amino acid residues unmodified. #750

- Fixed `SeqFeature` and `GTFSeqFeature` to remove the definition of strand and use `location.strand`. #616

- Refactored `util` so all functions are accessible. #749

## [0.11.5] - 2023-3-5

### Fixed

- Fixed `callVariant` that the command line argument `--max-variants-per-node` and `--additional-variants-per-misc` not passed to it and the default value was always used.

## [0.11.4] - 2023-2-23

### Fixed

- Fixed `callVariant` that variant peptides may get redundant labels with same information (transcript ID and variants). #679

## [0.11.3] - 2023-2-9

### Fixed

- `filterFasta` failed with the new noncoding peptide FASTA header. #675

## [0.11.2] - 2023-2-3

### Fixed

- Noncoding peptide headers not parsed successfully by summarizeFasta #672

## [0.11.1] - 2023-1-31

### Fixed

- Fixed `callVariant` that alt splice insertions were treated as stop altering when they are not.

- Argument `--orf-assignment` added to `callNoncoding` to allow choosing the min or max ORF. #667

## [0.11.0] - 2023-1-29

### Fixed

- When filtering variants for circRNA, those on fragments that are shorter than 3 nucleotides will not be included. #613

- When collapsing nodes with the same sequence, global variants (mostly circRNA) are no longer considered when comparing variants, so that nodes with no other variant won't be discarded mistakenly. #619

- Fixed issue that for hybrid nodes that span over a junction site on a circRNA, the location got lost and caused it to fail to identify whether the node is at least one loop downstream to an ORF start site. #621

- Fixed issue that circRNA with only 1 nucleotide was causing it to fail to filter variants. #623

- Fixed issue that peptide nodes on different subgraphs were collapsed after expanding the variant bubble, causing downstream nodes to be unprocessed and resulting in `*` in the final sequences. #625

- Fixed issue that the ORF start site position cannot be interpreted when checking whether it is at least one loop away, because it can be off by 1 when converting the location from the gene coordinate to amino acid. #630

- Fixed issue that the 'CHROM' attribute of GVF metadata not read in correctly. #629

- Fixed issue that when a frameshift insertion is on a alt splice frameshift substitution (or insertion), the node became disconnected after aligning the variant bubble. #635

- Fixed issue that when getting the stop altering mutations, location comparison was done incorrectly by 1. #636

- Node in circRNA missing downstream stop lost mutation called as variant peptide incorrectly. #637

- Silent mutation not excluded when it is very closed to anther mutation. #638

- Stop retaining mutation not excluded. #638

- Fusion with donor breakpoint smaller than 3 causing it fail to run. #633

- Alt splice insertion recognized as stop altering incorrectly. #640

- Fix that variants that overlap with last 3 nucleotides of the transcripts causing it to fail. #645

- Fixed `parseREDItools` that the ref and alt nucleotides were not set correctly for negative strands. #644

- Fixed `callVariant` that hybrid peptide sequences were called from circRNA. #653

- Fixed `callVariant` that peptides with variants of in-frame mutation causing deletion/insertion between two cleavage sites were missed. #655

- Fixed `callVariant` that when setting the max number of variants per peptide, the number of miscleavages was not used correctly. #657

### Changed

- The transcript trailing peptides (peptides at the end of the transcript sequence) are now excluded for transcripts with the `mRNA_end_NF` tag and circRNA regardless of it. Otherwise for transcripts (either coding or noncoding) that the mRNA end is confirmed (without the `mRNA_end_NF`) they are now included in the final FASTA. #649

- Change from ray back to pathos for parallelization. #643

- Gene ID of the transcript from which a noncoding peptide is called is added to the FASTA header. #662

## [0.10.1] - 2022-11-2

### Fixed

- Transcriptional coordinate to genomic coordinate not converted successfully when it is the last nucleotide of the transcript. #592

- When creating the peptide cleavage graph, the end nodes of ate variant bubble with alt splice were collapsed with the reference node causing the graph cleavage process terminated too early resulting uncleaved nodes. #597

- in `callVariant` when filtering variants associated with the donor transcript, the left breakpoint coordinate not converted successfully if it is the end of the transcript. #598

- Large deletion caused by alt splice raised the complexity drastically so had to go back to treat alt splice deletion as subgraph again. #600

- In circRNA, nodes that span over the backsplicing site with variant before the backsplicing site were not recognized correctly. #602

- Failed to find the downstream exon start when the breakpoint is the first nucleotide before the exon start of a transcript on the negative strand. #603

- In `bruteForce`, variants that start at the fusion breakpoint were not excluded. #604

- Variant peptides on circRNA are missed by `callVariant` when there are multiple ORF candidates. #606

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

- Added CLI command `mergeFasta` to merge multiple variant peptide database Fasta files into one. This could be useful when working with multiplexed proteomic experiments such as TMT. [#380](https://github.com/uclahs-cds/package-moPepGen/issues/380)

- Added CLI command `decoyFasta` to generate decoy database by shuffling or reversing each sequence. [#386](https://github.com/uclahs-cds/package-moPepGen/issues/386)

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

- `generateIndex` is changed to use compressed text format to store genomic annotation, because for some reason that we are not sure, when loading the pickled genomic annotation, the memory usage is almost doubled. [#394](https://github.com/uclahs-cds/package-moPepGen/issues/394)

---

## [0.1.0-beta.1] - 2021-12-23

### Added

- Initial beta release of moPepGen, with the three-frame graph based algorithm implemented to call noncanoinical peptides.
