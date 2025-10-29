"""Orchestrator for parallel/serial variant peptide calling with step-down policy.

This module coordinates the processing of transcripts to call variant peptides
from genomic variants. It handles:
- Loading and organizing variant data by transcript
- Dispatching work to worker functions (serial or parallel)
- Step-down policy for handling timeouts (retry with reduced parameters)
- Aggregating results and writing output files

The orchestrator uses a step-down policy to handle transcripts that timeout
during processing. When a timeout occurs, it automatically retries with more
conservative parameters (fewer variants per node, stricter bubble caps) to
ensure the work can complete.

Key Components:
- CallVariantOrchestrator: Main class that coordinates the entire workflow
- StepDownPolicy: Manages parameter adjustments on timeout
- _process_with_stepdown: Wrapper for worker calls with retry logic
- _gather_dispatch: Prepares isolated data for each transcript
"""
from __future__ import annotations
from typing import TYPE_CHECKING
import argparse
from contextlib import ExitStack
import dataclasses
import os
import copy
from pathos.pools import ParallelPool

from moPepGen import svgraph, seqvar, gtf, params, get_logger
from moPepGen.svgraph import ThreeFrameTVG
from moPepGen.svgraph.VariantPeptideTable import (
    get_peptide_table_path,
    get_peptide_table_path_temp
)

from .models import CallVariantDispatch, CallResult, Flags, Limits, TimeoutTracer, GraphPhase
from .graph_writer import GraphWriter
from . import call_variant_workers as workers


if TYPE_CHECKING:
    from typing import List
    from seqvar import VariantRecordPoolOnDisk
    from svgraph import VariantPeptideTable

class StepDownPolicy:
    """Manages parameter step-down for timeout recovery.

    When a transcript times out during processing, this policy determines
    how to adjust parameters for the retry. It steps down parameters in
    this order:

    1. For TVG (DNA graph) timeouts: Increase in_bubble_cap_step_down
       to reduce variant combinations in the DNA graph.

    2. For PVG (peptide graph) timeouts: Reduce max_variants_per_node
       and additional_variants_per_misc to limit graph complexity.

    Each transcript gets its own StepDownPolicy instance, so timeouts
    in one transcript don't affect others.

    Attributes:
        max_variants_per_node: List of values to try for max variants per node.
        additional_variants_per_misc: List of values for additional misc variants.
        in_bubble_cap_step_down: Current bubble cap step-down level.
    """
    def __init__(self, limits: Limits):
        self.max_variants_per_node = list(limits.max_variants_per_node)
        self.additional_variants_per_misc = list(limits.additional_variants_per_misc)
        self.in_bubble_cap_step_down = limits.in_bubble_cap_step_down

    def next(self, tracer: TimeoutTracer, cleavage_params: params.CleavageParams, tx_id: str):
        """Generate next set of parameters after a timeout.

        Adjusts parameters based on which graph phase (TVG or PVG) timed out.
        Uses the tracer to determine where the timeout occurred.

        Args:
            tracer: TimeoutTracer indicating which graph phase timed out.
            cleavage_params: Current cleavage parameters to adjust.
            tx_id: Transcript ID (for error messages).

        Returns:
            Updated CleavageParams with stepped-down values.

        Raises:
            ValueError: If we've exhausted all step-down options.
        """
        p = dataclasses.replace(cleavage_params)
        if tracer.graph == GraphPhase.TVG:
            self.in_bubble_cap_step_down += 1
            if self.in_bubble_cap_step_down > len(ThreeFrameTVG.VARIANT_BUBBLE_CAPS):
                raise ValueError(
                    f"Failed to finish transcript: {tx_id} with in_bubble_cap_step_down: "
                    f"{self.in_bubble_cap_step_down}"
                )
            p.in_bubble_cap_step_down = self.in_bubble_cap_step_down
        else:
            if len(self.max_variants_per_node) > 1:
                self.max_variants_per_node.pop(0)
            else:
                # step down from current param until >0
                next_val = max(p.max_variants_per_node - 1, 0)
                if next_val <= 0:
                    raise ValueError(f"Failed to finish transcript: {tx_id}")
                self.max_variants_per_node = [next_val]

            if len(self.additional_variants_per_misc) > 1:
                self.additional_variants_per_misc.pop(0)
            else:
                self.additional_variants_per_misc = [0]

            p.max_variants_per_node = self.max_variants_per_node[0]
            p.additional_variants_per_misc = self.additional_variants_per_misc[0]
        return p


def _process_with_stepdown(dispatch: CallVariantDispatch):
    """Process a transcript with automatic retry on timeout.

    This is the main entry point for processing a single transcript. It
    calls the worker function and catches TimeoutError. On timeout, it
    uses StepDownPolicy to adjust parameters and retries.

    Each transcript gets its own StepDownPolicy instance, ensuring that
    timeout handling is isolated and doesn't affect other transcripts.

    Args:
        dispatch: CallVariantDispatch containing all data and parameters
                  needed to process one transcript.

    Returns:
        CallResult with peptide annotations and optional graph data.

    Raises:
        ValueError: If step-down policy exhausts all options and still
                    cannot complete the transcript.
    """
    logger = get_logger()
    policy = StepDownPolicy(dispatch.limits)
    # Prepare kwargs for worker
    base_kwargs = dict(
        tx_id=dispatch.tx_id,
        variant_series=dispatch.variant_series,
        tx_seqs=dispatch.tx_seqs,
        gene_seqs=dispatch.gene_seqs,
        reference_data=dispatch.reference_data,
        codon_table=dispatch.codon_table,
        pool=dispatch.pool,
        cleavage_params=dispatch.cleavage_params,
        noncanonical_transcripts=dispatch.flags.noncanonical_transcripts,
        max_adjacent_as_mnv=dispatch.limits.max_adjacent_as_mnv,
        truncate_sec=dispatch.flags.truncate_sec,
        w2f_reassignment=dispatch.flags.w2f_reassignment,
        backsplicing_only=dispatch.flags.backsplicing_only,
        save_graph=dispatch.save_graph,
        coding_novel_orf=dispatch.flags.coding_novel_orf,
        skip_failed=dispatch.flags.skip_failed,
        timeout=dispatch.timeout_seconds,
    )
    tracer = TimeoutTracer()
    base_kwargs['tracer'] = tracer

    while True:
        try:
            peptide_anno, tx_id, dgraphs, pgraphs, success_flags = \
                workers.call_variant_peptides_wrapper(**base_kwargs)
            return CallResult(
                tx_id=tx_id,
                peptide_anno=peptide_anno,
                dgraphs=dgraphs,
                pgraphs=pgraphs,
                success={
                    'variant': success_flags[0],
                    'fusion': success_flags[1],
                    'circRNA': success_flags[2]
                }
            )
        except TimeoutError as e:
            p = policy.next(tracer, base_kwargs['cleavage_params'], dispatch.tx_id)
            base_kwargs['cleavage_params'] = p
            tracer = TimeoutTracer()
            base_kwargs['tracer'] = tracer
            logger.warning(
                "Transcript %s timed out. Retry with --max-variants-per-node %i "
                "--additional-variants-per-misc %i and --in-bubble-cap-step-down %i",
                dispatch.tx_id, p.max_variants_per_node, p.additional_variants_per_misc,
                p.in_bubble_cap_step_down
            )


class CallVariantOrchestrator:
    """Orchestrates variant peptide calling across all transcripts.

    This is the main coordinator for the callVariant pipeline. It:
    1. Loads reference data (genome, annotation, canonical peptides)
    2. Creates variant pools from GVF files
    3. Sorts transcripts by rank for processing
    4. Dispatches work to workers (parallel or serial)
    5. Aggregates results and writes output files

    The orchestrator supports both parallel and serial execution modes.
    In parallel mode, it uses pathos.pools.ParallelPool to process multiple
    transcripts simultaneously. In serial mode, it processes one at a time.

    Graph saving is optional and controlled by graph_output_dir. When enabled,
    DNA graphs (TVG/CVG) and peptide graphs (PVG) are written to JSON files
    for debugging and visualization.

    Attributes:
        args: Command-line arguments from argparse.
        variant_files: List of GVF input files.
        output_path: Path for output FASTA file.
        graph_output_dir: Optional directory for saving graph JSON files.
        threads: Number of parallel threads (1 = serial).
        reference_data: Reference data (genome, annotation, canonical peptides).
        cleavage_params: Peptide cleavage parameters.
        flags: Boolean flags controlling behavior.
        limits: Numerical limits for graph construction.
    """
    def __init__(self, args: argparse.Namespace, reference_data: params.ReferenceData = None):
        self.args = args
        self.variant_files: List = args.input_path
        self.output_path = args.output_path
        self.graph_output_dir = args.graph_output_dir
        self.threads = args.threads
        self.logger = get_logger()
        self.reference_data = reference_data  # Can be provided or loaded later

        self.peptide_table_temp_path = get_peptide_table_path_temp(self.output_path)
        self.peptide_table_output_path = get_peptide_table_path(self.output_path)

        self.cleavage_params = params.CleavageParams(
            enzyme=args.cleavage_rule,
            exception=args.cleavage_exception,
            miscleavage=int(args.miscleavage),
            min_mw=float(args.min_mw),
            min_length=args.min_length,
            max_length=args.max_length,
            max_variants_per_node=args.max_variants_per_node[0],
            additional_variants_per_misc=args.additional_variants_per_misc[0],
            in_bubble_cap_step_down=args.in_bubble_cap_step_down,
            min_nodes_to_collapse=args.min_nodes_to_collapse,
            naa_to_collapse=args.naa_to_collapse,
        )

        self.flags = Flags(
            noncanonical_transcripts=args.noncanonical_transcripts,
            backsplicing_only=args.backsplicing_only,
            coding_novel_orf=args.coding_novel_orf,
            skip_failed=args.skip_failed,
            truncate_sec=args.selenocysteine_termination,
            w2f_reassignment=args.w2f_reassignment,
        )
        self.limits = Limits(
            max_adjacent_as_mnv=args.max_adjacent_as_mnv,
            max_variants_per_node=list(args.max_variants_per_node),
            additional_variants_per_misc=list(args.additional_variants_per_misc),
            in_bubble_cap_step_down=args.in_bubble_cap_step_down,
        )

        self.invalid_protein_as_noncoding = args.invalid_protein_as_noncoding

    def _load_reference(self):
        """Load reference data if not already provided via constructor.

        This method provides lazy loading of reference data, which includes:
        - Genome sequences
        - Gene/transcript annotations
        - Canonical peptides (for denylist)
        - Codon tables (per chromosome)

        If reference_data was provided to __init__, this is a no-op.
        Otherwise, it imports and calls common.load_references().
        """
        if self.reference_data is not None:
            return  # Already loaded

        # Lazy import only when needed (backward compatibility)
        from moPepGen.cli import common as _common
        self.reference_data = _common.load_references(
            args=self.args,
            invalid_protein_as_noncoding=self.invalid_protein_as_noncoding,
            cleavage_params=self.cleavage_params,
            load_codon_tables=True,
        )

    def _create_variant_pool(self):
        """Create variant record pool from GVF files.

        Builds a VariantRecordPoolOnDisk that indexes all variants by
        transcript ID. This allows efficient lookup of variants affecting
        each transcript during processing.
        """
        self.variant_record_pool = seqvar.VariantRecordPoolOnDisk(
            gvf_files=self.variant_files,
            anno=self.reference_data.anno,
            genome=self.reference_data.genome,
        )

    def _gather_dispatch(self, tx_id: str, pool: VariantRecordPoolOnDisk
            ) -> CallVariantDispatch | None:
        """Prepare dispatch object for a single transcript.

        Gathers all data needed to process one transcript in isolation:
        - Variant series for the transcript
        - Transcript and gene sequences
        - Related transcripts (for fusions)
        - Codon table for the chromosome
        - Isolated genomic annotation and variant pool

        This isolation is critical for parallel processing - each worker
        gets a self-contained dataset with no shared references to avoid
        race conditions.

        Args:
            tx_id: Transcript identifier to gather data for.
            pool: Variant record pool containing all variants.

        Returns:
            CallVariantDispatch object with all necessary data, or None
            if the transcript should be skipped (no variants, etc.).
        """
        ref = self.reference_data
        tx_ids = [tx_id]
        tx_model = ref.anno.transcripts[tx_id]
        codon_table = ref.codon_tables[tx_model.transcript.chrom]
        try:
            variant_series = pool[tx_id]
        except ValueError:
            if self.args.skip_failed:
                return None
            raise
        if variant_series.is_empty():
            return None
        if self.flags.noncanonical_transcripts and \
                not variant_series.has_any_noncanonical_transcripts():
            return None
        tx_ids += variant_series.get_additional_transcripts()

        gene_seqs = {}
        if variant_series.is_gene_sequence_needed():
            gene_id = tx_model.transcript.gene_id
            _chrom = tx_model.transcript.chrom
            gene_model = ref.anno.genes[gene_id]
            gene_seq = gene_model.get_gene_sequence(ref.genome[_chrom])
            gene_seqs[gene_id] = gene_seq

        tx_seqs = {}
        for _tx_id in tx_ids:
            _tx_model = ref.anno.transcripts[_tx_id]
            _chrom = _tx_model.transcript.chrom
            tx_seqs[_tx_id] = _tx_model.get_transcript_sequence(ref.genome[_chrom])
            _gene_id = _tx_model.transcript.gene_id
            if _gene_id not in gene_seqs:
                _gene_model = ref.anno.genes[_gene_id]
                _gene_seq = _gene_model.get_gene_sequence(ref.genome[_chrom])
                gene_seqs[_gene_id] = _gene_seq
        gene_ids = list({ref.anno.transcripts[x].transcript.gene_id for x in tx_ids})

        gene_models = {}
        for gene_id in gene_ids:
            dummy_gene_model = copy.deepcopy(ref.anno.genes[gene_id])
            dummy_gene_model.transcripts = [x for x in dummy_gene_model.transcripts if x in tx_ids]
            gene_models[gene_id] = dummy_gene_model

        dummy_anno = gtf.GenomicAnnotation(
            genes=gene_models,
            transcripts={tx_id: ref.anno.transcripts[tx_id] for tx_id in tx_ids},
            source=ref.anno.source,
        )

        reference_data = params.ReferenceData(
            genome=None,
            anno=dummy_anno,
            canonical_peptides=set(),
            codon_tables={},
        )
        dummy_pool = seqvar.VariantRecordPool()
        dummy_pool.anno = dummy_anno
        dummy_pool.data[tx_id] = variant_series
        for add_tx in tx_ids:
            if add_tx != tx_id:
                try:
                    dummy_pool[add_tx] = pool[add_tx]
                except KeyError:
                    continue

        return CallVariantDispatch(
            tx_id=tx_id,
            variant_series=variant_series,
            tx_seqs=tx_seqs,
            gene_seqs=gene_seqs,
            reference_data=reference_data,
            codon_table=codon_table,
            pool=dummy_pool,
            cleavage_params=self.cleavage_params,
            flags=self.flags,
            limits=self.limits,
            save_graph=self.graph_output_dir is not None,
            timeout_seconds=self.args.timeout_seconds,
        )

    def run(self):
        """Execute the complete variant peptide calling workflow.

        Main entry point that orchestrates the entire pipeline:

        1. Create variant pool from GVF files
        2. Open output files (FASTA and peptide table)
        3. Sort transcripts by rank for optimal processing order
        4. Process transcripts (parallel or serial based on self.threads)
        5. Write final FASTA and peptide table files
        6. Clean up temporary files

        Progress is logged every 1000 transcripts. The peptide table is
        written to a temporary file during processing, then sorted and
        moved to the final location.

        Raises:
            Various exceptions from worker functions if skip_failed=False.
        """
        self._create_variant_pool()
        ref = self.reference_data
        logger = self.logger

        with ExitStack() as stack:
            pool = stack.enter_context(
                seqvar.VariantRecordPoolOnDiskOpener(self.variant_record_pool)
            )
            seq_anno_handle = stack.enter_context(open(self.peptide_table_temp_path, 'w+'))

            peptide_table = svgraph.VariantPeptideTable(seq_anno_handle)
            peptide_table.write_header()

            for file in pool.gvf_files:
                logger.info("Using GVF file: %s", file)
            tx_rank = ref.anno.get_transcript_rank()
            tx_sorted = sorted(pool.pointers.keys(), key=lambda x: tx_rank[x])
            logger.info('Variants sorted')
            n_total = len(tx_sorted)
            logger.info('Start calling variant peptides..')

            graph_writer = GraphWriter(self.graph_output_dir) if self.graph_output_dir else None

            processed = 0
            if self.threads > 1:
                process_pool = ParallelPool(ncpus=self.threads)
                dispatches: List[CallVariantDispatch] = []
                for tx_id in tx_sorted:
                    dispatch = self._gather_dispatch(tx_id, pool)
                    if not dispatch:
                        continue
                    dispatches.append(dispatch)
                for result in process_pool.map(_process_with_stepdown, dispatches):
                    processed += 1
                    self._handle_result(result, peptide_table, ref, graph_writer)
                    if processed % 1000 == 0:
                        logger.info(
                            '%.2f %% ( %i / %i ) transcripts processed.',
                            processed / n_total * 100, processed, n_total
                        )
            else:
                # Serial execution
                for tx_id in tx_sorted:
                    dispatch = self._gather_dispatch(tx_id, pool)
                    if not dispatch:
                        continue
                    result = _process_with_stepdown(dispatch)
                    processed += 1
                    self._handle_result(result, peptide_table, ref, graph_writer)
                    if processed % 1000 == 0:
                        logger.info(
                            '%.2f %% ( %i / %i ) transcripts processed.',
                            processed / n_total * 100, processed, n_total
                        )

            logger.info('All transcripts processed. Preparing FASTA file..')
            peptide_table.write_fasta(self.output_path)
            logger.info('Variant peptide FASTA written.')
            peptide_table.sort_table(self.peptide_table_output_path)
            logger.info('Variant peptide table sorted.')

        os.remove(self.peptide_table_temp_path)

    def _handle_result(self, result: CallResult, peptide_table: VariantPeptideTable,
            ref: params.ReferenceData, graph_writer: GraphWriter | None):
        """Process results from a single transcript.

        Takes the peptide annotations and optional graphs from a completed
        transcript and:
        1. Filters peptides against canonical peptides (denylist)
        2. Validates peptides against cleavage parameters (length, MW)
        3. Writes valid peptides to the peptide table
        4. Optionally writes DNA and peptide graphs to JSON files

        Args:
            result: CallResult containing peptide annotations and graphs.
            peptide_table: VariantPeptideTable for writing peptides.
            ref: Reference data (for canonical peptides denylist).
            graph_writer: Optional GraphWriter for saving graph JSON files.
        """
        # Write peptides if valid relative to canonical pool
        for peptide in result.peptide_anno:
            is_valid = peptide_table.is_valid(
                seq=peptide,
                canonical_peptides=ref.canonical_peptides,
                cleavage_params=self.cleavage_params,
            )
            if is_valid:
                for seq_anno in result.peptide_anno[peptide]:
                    peptide_table.add_peptide(peptide, seq_anno)
        if graph_writer:
            graph_writer.write_dgraphs(result.tx_id, result.dgraphs)
            graph_writer.write_pgraphs(result.tx_id, result.pgraphs)
