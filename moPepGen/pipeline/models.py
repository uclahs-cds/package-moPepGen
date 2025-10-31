"""Data models for the variant peptide calling pipeline.

This module defines the core data structures used by the CallVariantOrchestrator
and its worker functions. These models facilitate communication between pipeline
components and support features like timeout tracking, step-down policies, and
result aggregation.

The models are designed to be lightweight with minimal imports to avoid
circular dependencies and improve import performance.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple, Any

# We keep types loosely annotated to avoid heavy imports at module import time.

class GraphPhase(str, Enum):
    """Enumeration of graph construction phases for timeout tracking.

    Used by TimeoutTracer to identify which phase of graph construction
    was active when a timeout occurred, helping diagnose performance issues.

    Attributes:
        TVG: Three-frame Transcript Variant Graph phase (DNA level)
        PVG: Peptide Variant Graph phase (protein level)
    """
    TVG = "TVG"
    PVG = "PVG"

@dataclass
class TimeoutTracer:
    """Tracks the execution context when timeouts occur.

    This tracer helps diagnose which part of the pipeline (main variants,
    fusions, or circRNA) and which graph construction phase (TVG or PVG)
    was running when a timeout happened. This information is logged to help
    users identify problematic transcripts or variants.

    Attributes:
        backbone: The processing backbone/pathway being executed.
            - 'main': Processing transcriptional variants (SNVs, indels, splicing)
            - 'fusion': Processing fusion gene variants
            - 'circRNA': Processing circular RNA variants
            - None: No specific backbone (initialization state)
        graph: The graph construction phase being executed.
            - GraphPhase.TVG: Building three-frame transcript variant graph (DNA)
            - GraphPhase.PVG: Building peptide variant graph (protein)
            - None: No specific phase (initialization state)
    """
    backbone: Optional[str] = None  # 'main' | 'fusion' | 'circRNA'
    graph: Optional[GraphPhase] = None  # GraphPhase

@dataclass
class Flags:
    """Boolean flags controlling pipeline behavior.

    These flags are extracted from command-line arguments and control
    various aspects of variant peptide calling behavior. They are passed
    to worker functions to maintain consistent behavior across all processing
    pathways.

    Attributes:
        noncanonical_transcripts: Whether to process non-protein-coding transcripts.
            When True, processes lncRNA, pseudogenes, etc. for potential novel ORFs.
        backsplicing_only: For circRNA, only report peptides spanning the backsplicing
            junction. When False, reports all variant peptides from circRNA.
        coding_novel_orf: Whether to search for alternative start sites in
            protein-coding transcripts. When False, only uses annotated start sites.
        skip_failed: Whether to skip transcripts that fail processing instead of
            raising exceptions. Useful for large-scale processing where some failures
            are acceptable.
        truncate_sec: Whether to consider selenocysteine (Sec/U) termination events.
            When True, treats UGA as both Sec and a stop codon, generating both variants.
        w2f_reassignment: Whether to consider tryptophan (W) to phenylalanine (F)
            reassignment. Relevant for some mitochondrial and alternative genetic codes.
    """
    noncanonical_transcripts: bool
    backsplicing_only: bool
    coding_novel_orf: bool
    skip_failed: bool
    truncate_sec: bool
    w2f_reassignment: bool

@dataclass
class Limits:
    """Numeric limits and thresholds for variant processing.

    These parameters control how the pipeline handles complex variant scenarios,
    particularly in hypermutated regions where combinatorial explosion of variant
    combinations could cause performance issues.

    Attributes:
        max_adjacent_as_mnv: Maximum distance (in bases) between adjacent variants
            to be merged as a multi-nucleotide variant (MNV). Larger values create
            more complex MNVs but may better capture variant phase information.
        max_variants_per_node: List of caps for maximum variants per graph node.
            Used with step-down policy - if a transcript times out, it's retried
            with progressively stricter caps from this list. Default: [7]
        additional_variants_per_misc: List of caps for additional variants allowed
            per miscleavage site when calling peptides. Used with step-down policy.
            Default: [2]
        in_bubble_cap_step_down: Number of steps to reduce the in-bubble variant
            cap when retrying after timeout. 0 means use the default cap without
            modification. Higher values are more aggressive in reducing complexity.
    """
    max_adjacent_as_mnv: int
    max_variants_per_node: List[int] = field(default_factory=lambda: [7])
    additional_variants_per_misc: List[int] = field(default_factory=lambda: [2])
    in_bubble_cap_step_down: int = 0

# Result graph container types: (main, {fusion_id: graph}, {circ_id: graph})
# These type aliases define the structure for storing DNA and peptide graphs
# from different variant sources.
DGraphs = Tuple[Any, Dict[str, Any], Dict[str, Any]]  # DNA-level graphs (TVG/CVG)
PGraphs = Tuple[Any, Dict[str, Any], Dict[str, Any]]  # Protein-level graphs (PVG)

@dataclass
class CallVariantDispatch:
    """Input parameters for variant peptide calling workers.

    This dataclass bundles all the parameters needed to call variant peptides
    for a single transcript. It's passed to worker functions which may be
    executed serially or in parallel.

    The design allows workers to be stateless functions that receive all needed
    context through this dispatch object, making them suitable for parallel
    execution via multiprocessing.

    Attributes:
        tx_id: Transcript identifier (e.g., 'ENST00000123456.2')
        variant_series: Collection of variants affecting this transcript, organized
            by source (main variants, fusions, circRNAs)
        tx_seqs: Dictionary mapping transcript IDs to their DNA sequences
            (DNASeqRecordWithCoordinates objects)
        gene_seqs: Dictionary mapping gene IDs to their DNA sequences, used for
            processing fusions
        reference_data: Container with reference data (genome, annotation, proteome,
            canonical peptides, codon tables)
        codon_table: Codon table information for this transcript, including
            translation table and start codons
        pool: Variant record pool for coordinate lookups during graph construction
        cleavage_params: Parameters controlling peptide cleavage (enzyme, miscleavage,
            length limits, flanking_size for archipel mode)
        flags: Boolean flags controlling processing behavior (see Flags class)
        limits: Numeric limits for variant handling (see Limits class)
        save_graph: Whether to save and return graph objects. When False, graphs
            are discarded after peptide calling to save memory.
        timeout_seconds: Maximum time allowed for processing this transcript.
            If exceeded, the transcript is retried with step-down parameters or
            skipped if skip_failed is True.
    """
    tx_id: str
    variant_series: Any  # VariantSeriesGroup (loosely typed to avoid imports)
    tx_seqs: Dict[str, Any]  # Dict[str, DNASeqRecordWithCoordinates]
    gene_seqs: Dict[str, Any]  # Dict[str, DNASeqRecordWithCoordinates]
    reference_data: Any  # ReferenceData
    codon_table: Any  # CodonTableInfo
    pool: Any  # VariantRecordPoolOnDisk
    cleavage_params: Any  # CleavageParams
    flags: Flags
    limits: Limits
    save_graph: bool
    timeout_seconds: int

@dataclass
class CallResult:
    """Results from variant peptide calling for a single transcript.

    This dataclass encapsulates all outputs from processing a transcript,
    including the called peptides, optionally the graph objects, and success
    flags for each processing pathway.

    The orchestrator collects these results from worker functions and aggregates
    them into the final output FASTA file and optional graph files.

    Attributes:
        tx_id: Transcript identifier that was processed
        peptide_anno: Dictionary mapping peptide sequences to their annotations.
            Keys are Bio.Seq.Seq objects (peptide sequences), values are lists
            of AnnotatedPeptideLabel objects containing variant information,
            genomic coordinates, and other metadata.
        dgraphs: Tuple of DNA-level graphs (TVG/CVG):
            - [0]: Main three-frame transcript variant graph (or None)
            - [1]: Dict mapping fusion IDs to their CVG graphs
            - [2]: Dict mapping circRNA IDs to their CVG graphs
            Set to (None, {}, {}) if save_graph=False
        pgraphs: Tuple of protein-level graphs (PVG):
            - [0]: Main peptide variant graph (or None)
            - [1]: Dict mapping fusion IDs to their PVGs
            - [2]: Dict mapping circRNA IDs to their PVGs
            Set to (None, {}, {}) if save_graph=False
        success: Dictionary tracking which processing pathways completed successfully:
            - 'variant': Main transcriptional variants (SNVs, indels, splicing)
            - 'fusion': Fusion gene variants
            - 'circRNA': Circular RNA variants
            Values are True if processing succeeded, False if it failed/timed out.
            Used for logging and determining whether to apply step-down retry.
    """
    tx_id: str
    peptide_anno: Dict[Any, List[Any]]  # Dict[Seq, List[AnnotatedPeptideLabel]]
    dgraphs: DGraphs
    pgraphs: PGraphs
    success: Dict[str, bool]  # keys: 'variant', 'fusion', 'circRNA'
