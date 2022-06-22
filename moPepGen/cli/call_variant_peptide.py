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
import copy
from typing import List, Set, TYPE_CHECKING, Dict
from pathlib import Path
import ray
from moPepGen import svgraph, aa, seqvar, logger, gtf, params
from moPepGen.cli import common


if TYPE_CHECKING:
    from moPepGen import dna, circ


INPUT_FILE_FORMATS = ['.gvf']
OUTPUT_FILE_FORMATS = ['.fasta', '.fa']

# pylint: disable=W0212
def add_subparser_call_variant(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callPeptides """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='callVariant',
        help='Call non-canonical peptides from genomic variants.',
        description='Genomic variant data must be generated by one of the'
        ' moPepGen parser. See moPepGen --help',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    common.add_args_input_path(
        parser=p,
        formats=INPUT_FILE_FORMATS,
        plural=True,
        message='File path to GVF files. Must be generated by any of the'
        ' moPepGen parsers.'
    )
    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    p.add_argument(
        '--max-variants-per-node',
        type=int,
        help='Maximal number of variants per node. This argument can be useful'
        ' when there are local regions that are heavily mutated. When creating'
        ' the cleavage graph, nodes containing variants larger than this value'
        ' are skipped. Setting to -1 will avoid checking for this.',
        default=7,
        metavar='<number>'
    )
    p.add_argument(
        '--additional-variants-per-misc',
        type=int,
        help='Additional variants allowed for every miscleavage. This argument'
        ' is used together with --max-variants-per-node to handle hypermutated'
        ' regions. Setting to -1 will avoid checking for this.',
        default=2,
        metavar='<number>'
    )
    p.add_argument(
        '--min-nodes-to-collapse',
        type=int,
        help='When making the cleavage graph, the minimal number of nodes'
        ' to trigger pop collapse.',
        default=30,
        metavar='<number>'
    )
    p.add_argument(
        '--naa-to-collapse',
        type=int,
        help='The number of bases used for pop collapse.',
        default=5,
        metavar='<number>'
    )
    p.add_argument(
        '--noncanonical-transcripts',
        action='store_true',
        help='Process only noncanonical transcripts of fusion transcripts and'
        ' circRNA. Canonical transcripts are skipped.'
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
    common.add_args_reference(p)
    common.add_args_cleavage(p)
    common.add_args_quiet(p)

    p.set_defaults(func=call_variant_peptide)
    common.print_help_if_missing_args(p)
    return p


class VariantPeptideCaller():
    """ Helper class to call variant peptides """
    def __init__(self, args:argparse.Namespace):
        """ constructor """
        self.args = args
        self.variant_files:List[Path] = args.input_path

        for file in self.variant_files:
            common.validate_file_format(file, INPUT_FILE_FORMATS, True)

        self.output_path:str = args.output_path
        common.validate_file_format(self.output_path, OUTPUT_FILE_FORMATS)

        self.quiet:bool = args.quiet
        self.cleavage_params = params.CleavageParams(
            enzyme=args.cleavage_rule,
            exception='trypsin_exception' if args.cleavage_rule == 'trypsin' else None,
            miscleavage=int(args.miscleavage),
            min_mw=float(args.min_mw),
            min_length=args.min_length,
            max_length=args.max_length
        )
        self.max_variants_per_node:int = args.max_variants_per_node
        self.additional_variants_per_misc:int = args.additional_variants_per_misc
        self.min_nodes_to_collapse:int = args.min_nodes_to_collapse
        self.naa_to_collapse:int = args.naa_to_collapse
        self.noncanonical_transcripts = args.noncanonical_transcripts
        self.invalid_protein_as_noncoding:bool = args.invalid_protein_as_noncoding
        self.verbose = args.verbose_level
        self.threads = args.threads
        if self.quiet is True:
            self.verbose = 0
        self.reference_data = None
        self.variant_peptides = aa.VariantPeptidePool()
        self.variant_record_pool:seqvar.VariantRecordPoolOnDisk = None

    def load_reference(self):
        """ load reference genome, annotation, and canonical peptides """
        genome, anno, _, canonical_peptides = common.load_references(
                args=self.args,
                invalid_protein_as_noncoding=self.invalid_protein_as_noncoding
            )
        self.reference_data = params.ReferenceData(
            genome=genome,
            anno=anno,
            canonical_peptides=canonical_peptides
        )

    def create_in_disk_variant_pool(self):
        """ Create in disk variant pool """
        self.variant_record_pool = seqvar.VariantRecordPoolOnDisk(
            gvf_files=self.variant_files,
            anno=self.reference_data.anno,
            genome=self.reference_data.genome
        )

    def write_fasta(self):
        """ write variant peptides to fasta. """
        self.variant_peptides.write(self.output_path)
        if self.verbose >= 1:
            logger('Variant peptide FASTA file written to disk.')

def call_variant_peptides_wrapper(tx_id:str,
        variant_series:seqvar.TranscriptionalVariantSeries,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        reference_data:params.ReferenceData,
        pool:seqvar.VariantRecordPool,
        cleavage_params:params.CleavageParams,
        noncanonical_transcripts:bool
        ) -> List[Set[aa.AminoAcidSeqRecord]]:
    """ wrapper function to call variant peptides """
    peptide_pool:List[Set[aa.AminoAcidSeqRecord]] = []
    if variant_series.transcriptional:
        try:
            if not noncanonical_transcripts or \
                    variant_series.has_any_alternative_splicing():
                peptides = call_peptide_main(
                    tx_id=tx_id, tx_variants=variant_series.transcriptional,
                    variant_pool=pool, ref=reference_data,
                    tx_seqs=tx_seqs, gene_seqs=gene_seqs,
                    cleavage_params=cleavage_params
                )
                peptide_pool.append(peptides)
        except:
            logger(f'Exception raised from {tx_id}')
            raise

    peptides = set()
    exclude_variant_types = ['Fusion', 'Insertion', 'Deletion', 'Substitution', 'circRNA']
    for variant in variant_series.fusion:
        filtered_variants = pool.filter_variants(
            tx_ids=[tx_id], start=0, end=variant.location.start,
            exclude_type=exclude_variant_types, intron=False,
            return_coord='transcript'
        )
        variant_pool = copy.copy(pool)
        variant_pool[tx_id] = copy.copy(pool[tx_id])
        variant_pool[tx_id].transcriptional = filtered_variants
        try:
            _peptides = call_peptide_fusion(
                variant=variant, variant_pool=variant_pool, ref=reference_data,
                tx_seqs=tx_seqs, gene_seqs=gene_seqs, cleavage_params=cleavage_params
            )
        except:
            logger(
                f"Exception raised from fusion {variant.id}, "
                f"donor:{tx_id}, accepter: {variant.attrs['ACCEPTER_TRANSCRIPT_ID']}"
            )
            raise

        peptides.update(_peptides)
    peptide_pool.append(peptides)
    peptides = set()
    for circ_model in variant_series.circ_rna:
        try:
            _peptides = call_peptide_circ_rna(
                record=circ_model, ref=reference_data,
                variant_pool=pool, gene_seqs=gene_seqs,
                cleavage_params=cleavage_params
            )
        except:
            logger(f"Exception raised from {circ_model.id}")
            raise
        peptides.update(_peptides)
    peptide_pool.append(peptides)

    return peptide_pool

@ray.remote
def wrapper_remote(dispatch):
    """ wrapper for ParallelPool """
    return call_variant_peptides_wrapper(*dispatch)

def wrapper_local(dispatch):
    """ wrapper for local """
    return call_variant_peptides_wrapper(*dispatch)

def call_variant_peptide(args:argparse.Namespace) -> None:
    """ Main entry point for calling variant peptide """
    caller = VariantPeptideCaller(args)
    common.print_start_message(args)
    caller.load_reference()
    caller.create_in_disk_variant_pool()
    noncanonical_transcripts = caller.noncanonical_transcripts
    ref = caller.reference_data

    with seqvar.VariantRecordPoolOnDiskOpener(caller.variant_record_pool) as pool:
        if caller.verbose >= 1:
            for file in pool.gvf_files:
                logger(f"Using GVF file: {file}")
        tx_rank = ref.anno.get_transcript_rank()
        tx_sorted = sorted(pool.pointers.keys(), key=lambda x:tx_rank[x])
        if caller.verbose >= 1:
            logger('Variants sorted')

        # Not providing canonical peptide pool to each task for now.
        canonical_peptides = set()
        if caller.threads > 1:
            ray.init(num_cpus=caller.threads)
            # ray.put(canonical_peptides)

        dispatches = []
        i = 0
        for tx_id in tx_sorted:
            tx_ids = [tx_id]
            tx_model = ref.anno.transcripts[tx_id]
            variant_series = pool[tx_id]
            if variant_series.is_empty():
                continue
            if noncanonical_transcripts and \
                    not variant_series.has_any_noncanonical_transcripts():
                continue
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
                dummy_gene_model.transcripts = [x for x in dummy_gene_model.transcripts
                    if x in tx_ids]
                gene_models[gene_id] = dummy_gene_model

            dummy_anno = gtf.GenomicAnnotation(
                genes=gene_models,
                transcripts={tx_id:ref.anno.transcripts[tx_id] for tx_id in tx_ids},
                source=ref.anno.source
            )
            reference_data = params.ReferenceData(
                genome=None,
                anno=dummy_anno,
                canonical_peptides=canonical_peptides
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
                tx_id, variant_series, tx_seqs, gene_seqs, reference_data,
                dummy_pool, caller.cleavage_params, noncanonical_transcripts
            )
            dispatches.append(dispatch)

            reloaded = ((i + 1) % caller.threads == 0 or i + 1 == len(tx_sorted)) \
                and len(dispatches) > 0
            if reloaded:
                if caller.verbose >= 2:
                    logger([x[0] for x in dispatches])
                if caller.threads > 1:
                    # results = process_pool.map(wrapper, dispatches)
                    results = ray.get([wrapper_remote.remote(d) for d in dispatches])
                else:
                    results = [wrapper_local(dispatches[0])]

                for peptide_series in results:
                    for peptides in peptide_series:
                        for peptide in peptides:
                            caller.variant_peptides.add_peptide(
                                peptide=peptide,
                                canonical_peptides=ref.canonical_peptides,
                                cleavage_params=caller.cleavage_params
                            )
                dispatches = []

            if caller.verbose >= 1:
                i += 1
                if i % 1000 == 0:
                    logger(f'{i} transcripts processed.')

    caller.write_fasta()


def call_peptide_main(tx_id:str, tx_variants:List[seqvar.VariantRecord],
        variant_pool:seqvar.VariantRecordPoolOnDisk,
        ref:params.ReferenceData,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        cleavage_params:params.CleavageParams
        ) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for main variants (except cirRNA). """
    tx_model = ref.anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        cleavage_params=cleavage_params
    )
    dgraph.init_three_frames()
    dgraph.create_variant_graph(
        variants=tx_variants, variant_pool=variant_pool, genome=ref.genome,
        anno=ref.anno, tx_seqs=tx_seqs, gene_seqs=gene_seqs,
    )
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()

    pgraph.create_cleavage_graph()
    return pgraph.call_variant_peptides(blacklist=ref.canonical_peptides)

def call_peptide_fusion(variant:seqvar.VariantRecord,
        variant_pool:seqvar.VariantRecordPoolOnDisk,
        ref:params.ReferenceData,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        cleavage_params:params.CleavageParams
        ) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for fusion """
    tx_id = variant.location.seqname
    tx_model = ref.anno.transcripts[tx_id]
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
        cleavage_params=cleavage_params
    )
    dgraph.init_three_frames()
    dgraph.create_variant_graph(
        variants=tx_variants, variant_pool=variant_pool, genome=ref.genome,
        anno=ref.anno, tx_seqs=tx_seqs, gene_seqs=gene_seqs,
    )
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()
    pgraph.create_cleavage_graph()
    return pgraph.call_variant_peptides(blacklist=ref.canonical_peptides)

def call_peptide_circ_rna(record:circ.CircRNAModel, ref:params.ReferenceData,
        variant_pool:seqvar.VariantRecordPool,
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        cleavage_params:params.CleavageParams
        )-> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides from a given circRNA """
    gene_id = record.gene_id
    gene_seq = gene_seqs[gene_id]
    record.fragments.sort()
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
        cleavage_params=cleavage_params
    )
    cgraph.init_three_frames()
    cgraph.create_variant_circ_graph(variant_records)
    cgraph.extend_loop()
    cgraph.truncate_three_frames()
    cgraph.fit_into_codons()
    pgraph = cgraph.translate()
    pgraph.create_cleavage_graph()
    return pgraph.call_variant_peptides(blacklist=ref.canonical_peptides)
