""" `callVariant` is the core of moPepGen. It takes multiple GVF files, generated
by any moPepGen parser, and calls variant peptides caused by genomic variants
using a graph-based algorithm. For any transcript, it creates a three-frame
transcript variant graph by incorporating all variants from any sources (SNV,
INDEL, fusion, alternative splicing, RNA editing, and circRNA). The transcript
variant graph is then translated into a peptide variant graph, followed by
converting to a cleavage graph based on the enzymatic cleavage rule. The
variant peptide graph is than used to call for variant peptides that contains
at least one variant, and do not present in the canonical peptide pool.
 """
from __future__ import annotations
import argparse
from typing import List, Set, TYPE_CHECKING, Dict
from pathlib import Path
from pathos.pools import ParallelPool
from moPepGen import svgraph, aa, seqvar, logger, gtf
from moPepGen.cli.common import add_args_cleavage, add_args_quiet, \
    print_start_message, print_help_if_missing_args, add_args_reference, \
    load_references


if TYPE_CHECKING:
    from moPepGen import dna, circ


# pylint: disable=W0212
def add_subparser_call_variant(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callPeptides """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='callVariant',
        help='Call non-canonical peptides from genomic variants.',
        description="Genomic variant data must be generated by one of the"
        "moPepGen parser. See moPepGen --help",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        '-i', '--input-variant',
        type=Path,
        nargs='+',
        help='Path to input variant files. Must be generated by any of the'
        'moPepGen parser. This can be multiple.',
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '-o', '--output-fasta',
        type=Path,
        help='Filename for the output FASTA.',
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '--max-variants-per-node',
        type=int,
        help='Maximal number of variants per node. This argument can be useful'
        ' when there are local regions that are heavily mutated. When creating'
        ' the cleavage graph, nodes containing variants larger than this value'
        ' are skipped. Setting to -1 will avoid checking for this.',
        default=-1,
        metavar='<number>'
    )
    p.add_argument(
        '--noncanonical-transcripts',
        action='store_true',
        help='Process only noncanonical transcripts of fusion transcripts and'
        'circRNA. Canonical transcripts are skipped.'
    )
    p.add_argument(
        '--verbose-level',
        type=int,
        help='Level of verbose for logging.',
        default=1,
        metavar='<number>'
    )
    p.add_argument(
        '--threads',
        type=int,
        help='Set number of threads to be used.',
        default=1,
        metavar='<number>'
    )
    add_args_reference(p)
    add_args_cleavage(p)
    add_args_quiet(p)

    p.set_defaults(func=call_variant_peptide)
    print_help_if_missing_args(p)
    return p


class VariantPeptideCaller():
    """ Helper class to call variant peptides """
    def __init__(self, args:argparse.Namespace):
        """ constructor """
        self.args = args
        self.variant_files:List[Path] = args.input_variant

        for file in self.variant_files:
            if not file.exists():
                raise FileNotFoundError(f"GVF file not found: {file}")

        self.output_fasta:str = args.output_fasta
        self.quiet:bool = args.quiet
        self.rule:str = args.cleavage_rule
        self.miscleavage:int = int(args.miscleavage)
        self.min_mw:float = float(args.min_mw)
        self.exception = 'trypsin_exception' if self.rule == 'trypsin' else None
        self.min_length:int = args.min_length
        self.max_length:int = args.max_length
        self.max_variants_per_node:int = args.max_variants_per_node
        self.noncanonical_transcripts = args.noncanonical_transcripts
        self.verbose = args.verbose_level
        self.threads = args.threads
        if self.quiet is True:
            self.verbose = 0
        self.genome:dna.DNASeqDict = None
        self.anno:gtf.GenomicAnnotation = None
        self.canonical_peptides:Set[str] = None
        self.variant_peptides = aa.VariantPeptidePool()
        self.variant_record_pool:seqvar.VariantRecordPoolOnDisk = None

    def load_reference(self):
        """ load reference genome, annotation, and canonical peptides """
        self.genome, self.anno, _, self.canonical_peptides = load_references(self.args)

    def create_in_disk_vairant_pool(self):
        """ Create in disk variant pool """
        self.variant_record_pool = seqvar.VariantRecordPoolOnDisk(
            gvf_files=self.variant_files, anno=self.anno, genome=self.genome
        )

    def write_fasta(self):
        """ write variant peptides to fasta. """
        self.variant_peptides.write(self.output_fasta)
        if self.verbose >= 1:
            logger('Variant peptide FASTA file written to disk.')

def call_variant_peptides_wrapper(tx_id:str,
        variant_series:seqvar.TranscriptionalVariantSeries,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        anno:gtf.GenomicAnnotation, pool:seqvar.VariantRecordPool,
        rule:str, exception:str, miscleavage:int, max_variants_per_node:int,
        noncanonical_transcripts:bool) -> List[Set[aa.AminoAcidSeqRecord]]:
    """ wrapper function to call variant peptides """
    peptide_pool:List[Set[aa.AminoAcidSeqRecord]] = []
    if variant_series.transcriptional:
        try:
            if not noncanonical_transcripts or \
                    variant_series.has_any_alternative_splicing():
                peptides = call_peptide_main(
                    tx_id=tx_id, tx_variants=variant_series.transcriptional,
                    variant_pool=pool, anno=anno, genome=None, tx_seqs=tx_seqs,
                    gene_seqs=gene_seqs, rule=rule, exception=exception,
                    miscleavage=miscleavage,
                    max_variants_per_node=max_variants_per_node
                )
                peptide_pool.append(peptides)
        except:
            logger(f'Exception raised from {tx_id}')
            raise

    peptides = set()
    for variant in variant_series.fusion:
        try:
            _peptides = call_peptide_fusion(
                variant=variant, variant_pool=pool, anno=anno, genome=None,
                tx_seqs=tx_seqs, gene_seqs=gene_seqs, rule=rule,
                exception=exception, miscleavage=miscleavage,
                max_variants_per_node=max_variants_per_node
            )
        except:
            logger(f"Exception raised from fusion {variant.id}")
            raise

        peptides.update(_peptides)
    peptide_pool.append(peptides)
    peptides = set()
    for circ_model in variant_series.circ_rna:
        try:
            _peptides = call_peptide_circ_rna(
                record=circ_model, variant_pool=pool,
                gene_seqs=gene_seqs, rule=rule, exception=exception,
                miscleavage=miscleavage,
                max_variants_per_node=max_variants_per_node
            )
        except:
            logger(f"Exception raised from {circ_model.id}")
            raise
        peptides.update(_peptides)
    peptide_pool.append(peptides)

    return peptide_pool

def wrapper(dispatch):
    """ wrapper for ParallelPool """
    return call_variant_peptides_wrapper(*dispatch)

def call_variant_peptide(args:argparse.Namespace) -> None:
    """ Main entry point for calling variant peptide """
    caller = VariantPeptideCaller(args)
    print_start_message(args)
    caller.load_reference()
    caller.create_in_disk_vairant_pool()
    rule = caller.rule
    exception = caller.exception
    miscleavage = caller.miscleavage
    max_variants_per_node = caller.max_variants_per_node
    noncanonical_transcripts = caller.noncanonical_transcripts

    with seqvar.VariantRecordPoolOnDiskOpener(caller.variant_record_pool) as pool:
        tx_rank = caller.anno.get_transcirpt_rank()
        tx_sorted = sorted(pool.pointers.keys(), key=lambda x:tx_rank[x])
        # tx_sorted = pool.get_transcript_order()
        logger('Variants sorted')
        if caller.threads > 1:
            process_pool = ParallelPool(ncpus=caller.threads)
        dispatches = []
        i = 0
        for tx_id in tx_sorted:
            tx_ids = [tx_id]
            tx_model = caller.anno.transcripts[tx_id]
            variant_series = pool[tx_id]
            if variant_series.is_empty():
                continue
            if noncanonical_transcripts and \
                    not variant_series.has_any_noncanonical_transcripts():
                continue
            tx_ids += variant_series.get_additional_transcripts()
            tx_ids = set(tx_ids)

            gene_seqs = {}
            if variant_series.is_gene_sequence_needed():
                gene_id = tx_model.transcript.gene_id
                _chrom = tx_model.transcript.chrom
                gene_model = caller.anno.genes[gene_id]
                gene_seq = gene_model.get_gene_sequence(caller.genome[_chrom])
                gene_seqs[gene_id] = gene_seq
                tx_ids.update(gene_model.transcripts)

            tx_seqs = {}
            for _tx_id in tx_ids:
                _tx_model = caller.anno.transcripts[_tx_id]
                _chrom = _tx_model.transcript.chrom
                tx_seqs[_tx_id] = _tx_model.get_transcript_sequence(caller.genome[_chrom])
                _gene_id = _tx_model.transcript.gene_id
                if _gene_id not in gene_seqs:
                    _gene_model = caller.anno.genes[_gene_id]
                    _gene_seq = _gene_model.get_gene_sequence(caller.genome[_chrom])
                    gene_seqs[_gene_id] = _gene_seq

            gene_ids = list({caller.anno.transcripts[x].transcript.gene_id for x in tx_ids})
            for gene_id in gene_ids:
                tx_ids.update(caller.anno.genes[gene_id].transcripts)
            tx_ids = list(tx_ids)

            dummy_anno = gtf.GenomicAnnotation(
                genes={gene_id:caller.anno.genes[gene_id] for gene_id in gene_ids},
                transcripts={tx_id:caller.anno.transcripts[tx_id] for tx_id in tx_ids},
                source=caller.anno.source
            )

            dummy_pool = seqvar.VariantRecordPool()
            dummy_pool.anno = dummy_anno
            dummy_pool.data[tx_id] = variant_series
            for add_tx in tx_ids:
                if add_tx != tx_id:
                    try:
                        dummy_pool[add_tx] = caller.variant_record_pool[add_tx]
                    except KeyError:
                        continue

            dispatch = (
                tx_id, variant_series, tx_seqs, gene_seqs, dummy_anno,
                dummy_pool, rule, exception, miscleavage, max_variants_per_node,
                noncanonical_transcripts
            )
            dispatches.append(dispatch)

            reloaded = ((i + 1) % caller.threads == 0 or i + 1 == len(tx_sorted)) \
                and len(dispatches) > 0
            if reloaded:
                if caller.verbose >= 2:
                    logger([x[0] for x in dispatches])
                if caller.threads > 1:
                    results = process_pool.map(wrapper, dispatches)
                else:
                    results = [wrapper(dispatches[0])]

                for peptide_series in results:
                    for peptides in peptide_series:
                        for peptide in peptides:
                            caller.variant_peptides.add_peptide(
                                peptide=peptide,
                                canonical_peptides=caller.canonical_peptides,
                                min_mw=caller.min_mw, min_length=caller.min_length,
                                max_length=caller.max_length
                            )
                dispatches = []

            if caller.verbose >= 1:
                i += 1
                if i % 1000 == 0:
                    logger(f'{i} transcripts processed.')

    caller.write_fasta()


def call_peptide_main(tx_id:str, tx_variants:List[seqvar.VariantRecord],
        variant_pool:seqvar.VariantRecordPoolOnDisk, anno:gtf.GenomicAnnotation,
        genome:dna.DNASeqDict,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        rule:str, exception:str, miscleavage:int,
        max_variants_per_node:int) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for main variants (except cirRNA). """
    tx_model = anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        max_variants_per_node=max_variants_per_node
    )
    dgraph.init_three_frames()
    dgraph.create_variant_graph(
        variants=tx_variants, variant_pool=variant_pool, genome=genome,
        anno=anno, tx_seqs=tx_seqs, gene_seqs=gene_seqs,
    )
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()

    pgraph.create_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)

def call_peptide_fusion(variant:seqvar.VariantRecord,
        variant_pool:seqvar.VariantRecordPoolOnDisk, anno:gtf.GenomicAnnotation,
        genome:dna.DNASeqDict,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        rule:str, exception:str, miscleavage:int,
        max_variants_per_node:int) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for fusion """
    tx_id = variant.location.seqname
    tx_model = anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]
    tx_seq = tx_seq[:variant.location.start]

    if tx_id in variant_pool:
        tx_variants = [x for x in variant_pool[tx_id].transcriptional
            if x.location.end < variant.location.end]
    else:
        tx_variants = []

    tx_variants.append(variant)

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        max_variants_per_node=max_variants_per_node
    )
    dgraph.init_three_frames()
    dgraph.create_variant_graph(
        variants=tx_variants, variant_pool=variant_pool, genome=genome,
        anno=anno, tx_seqs=tx_seqs, gene_seqs=gene_seqs,
    )
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()
    pgraph.create_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)

def call_peptide_circ_rna(record:circ.CircRNAModel,
        variant_pool:seqvar.VariantRecordPool,
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        rule:str, exception:str, miscleavage:int, max_variants_per_node:int
        )-> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides from a given circRNA """
    gene_id = record.gene_id
    gene_seq = gene_seqs[gene_id]
    circ_seq = record.get_circ_rna_sequence(gene_seq)

    # Alternative splicing should not be included. Alternative splicing
    # represented as Insertion, Deletion or Substitution.
    exclusion_variant_types = ['Insertion', 'Deletion', 'Substitution']

    variant_records = variant_pool.filter_variants(
        tx_ids=[record.transcript_id], exclude_type=exclusion_variant_types,
        intron=False, segments=record.fragments
    )

    cgraph = svgraph.ThreeFrameCVG(
        circ_seq, _id=record.id, circ_record=record,
        max_variants_per_node=max_variants_per_node
    )
    cgraph.init_three_frames()
    cgraph.create_variant_circ_graph(variant_records)
    cgraph.extend_loop()
    cgraph.truncate_three_frames()
    cgraph.fit_into_codons()
    pgraph = cgraph.translate()
    pgraph.create_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)
