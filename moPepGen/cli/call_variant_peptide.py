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
import json
from typing import List, Set, TYPE_CHECKING, Dict, Tuple
from pathlib import Path
from pathos.pools import ParallelPool
from moPepGen import svgraph, aa, seqvar, logger, gtf, params
from moPepGen.cli import common
from moPepGen.SeqFeature import FeatureLocation, SeqFeature


if TYPE_CHECKING:
    from moPepGen import dna, circ
    from moPepGen.svgraph import ThreeFrameCVG, ThreeFrameTVG, PeptideVariantGraph


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
        required=True,
        message='File path to GVF files. Must be generated by any of the'
        ' moPepGen parsers.'
    )
    common.add_args_output_path(p, OUTPUT_FILE_FORMATS)
    p.add_argument(
        '--graph-output-dir',
        type=Path,
        required=False,
        default=None,
        help='Directory path that graph data are saved to. Graph data are not'
        ' saved if this is not given.',
        metavar='<file>'
    )
    p.add_argument(
        '--max-adjacent-as-mnv',
        type=int,
        help='Max number of adjacent variants that should be merged.',
        default=2
    )
    p.add_argument(
        '--selenocysteine-termination',
        action='store_true',
        help='Include peptides of selenoprotiens that the UGA is treated as '
        'termination instead of Sec.'
    )
    p.add_argument(
        '--w2f-reassignment',
        action='store_true',
        help='Include peptides with W > F (Tryptophan to Phenylalanine) '
        'reassignment.'
    )
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
            common.validate_file_format(
                file, INPUT_FILE_FORMATS, check_readable=True
            )

        self.output_path:str = args.output_path
        common.validate_file_format(
            self.output_path, OUTPUT_FILE_FORMATS, check_writable=True
        )
        self.graph_output_dir:Path = args.graph_output_dir
        if self.graph_output_dir is not None:
            common.validate_file_format(
                self.graph_output_dir, check_writable=True, is_directory=True
            )

        self.quiet:bool = args.quiet
        self.cleavage_params = params.CleavageParams(
            enzyme=args.cleavage_rule,
            exception=args.cleavage_exception,
            miscleavage=int(args.miscleavage),
            min_mw=float(args.min_mw),
            min_length=args.min_length,
            max_length=args.max_length,
            max_variants_per_node = args.max_variants_per_node,
            additional_variants_per_misc = args.additional_variants_per_misc,
            min_nodes_to_collapse = args.min_nodes_to_collapse,
            naa_to_collapse = args.naa_to_collapse
        )
        self.max_adjacent_as_mnv = args.max_adjacent_as_mnv
        self.truncate_sec:bool = args.selenocysteine_termination
        self.w2f_reassignment = args.w2f_reassignment
        self.noncanonical_transcripts = args.noncanonical_transcripts
        self.invalid_protein_as_noncoding:bool = args.invalid_protein_as_noncoding
        self.verbose = args.verbose_level
        self.threads = args.threads
        if self.quiet is True:
            self.verbose = 0
        self.reference_data = None
        self.variant_peptides = aa.VariantPeptidePool()
        self.variant_record_pool:seqvar.VariantRecordPoolOnDisk = None

        if not self.variant_files:
            if not (self.w2f_reassignment or self.truncate_sec):
                raise ValueError('Please provide at least a GVF file.')
            self.process_all_transcripts = True
            self.check_external_variants = False
        else:
            self.process_all_transcripts = False
            self.check_external_variants = True

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

    def write_dgraphs(self, tx_id:str, dgraphs:TypeDGraphs):
        """ Write ThreeFrameTVG data """

    def write_pgraphs(self, tx_id:str, pgraphs:TypePGraphs):
        """ Write PeptideVariantGraph data """
        if pgraphs[0]:
            data = pgraphs[0].jsonfy()
            with open(self.graph_output_dir/f"{tx_id}_main_PVG.json", 'wt') as handle:
                json.dump(data, handle)

        for var_id, pgraph in pgraphs[1].items():
            if pgraph is None:
                continue
            data = pgraph.jsonfy()
            with open(self.graph_output_dir/f"{tx_id}_Fusion_{var_id}_PVG.json", 'wt') as handle:
                json.dump(data, handle)

        for var_id, pgraph in pgraphs[2].items():
            if pgraph is None:
                continue
            data = pgraph.jsonfy()
            with open(self.graph_output_dir/f"{tx_id}_circRNA_{var_id}_PVG.json", 'wt') as handle:
                json.dump(data, handle)

TypeDGraphs = Tuple[
    svgraph.ThreeFrameTVG,
    Dict[str, svgraph.ThreeFrameTVG],
    Dict[str, svgraph.ThreeFrameCVG]
]
TypePGraphs = Tuple[
    svgraph.PeptideVariantGraph,
    Dict[str, svgraph.PeptideVariantGraph],
    Dict[str, svgraph.PeptideVariantGraph]
]
def call_variant_peptides_wrapper(tx_id:str,
        variant_series:seqvar.TranscriptionalVariantSeries,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        reference_data:params.ReferenceData,
        pool:seqvar.VariantRecordPool,
        cleavage_params:params.CleavageParams,
        noncanonical_transcripts:bool,
        max_adjacent_as_mnv:bool,
        truncate_sec:bool,
        w2f_reassignment:bool,
        save_graph:bool
        ) -> Tuple[Set[aa.AminoAcidSeqRecord], str, TypeDGraphs, TypePGraphs]:
    """ wrapper function to call variant peptides """
    peptide_pool:List[Set[aa.AminoAcidSeqRecord]] = []
    main_peptides = None
    denylist = call_canonical_peptides(
        tx_id=tx_id, ref=reference_data, tx_seq=tx_seqs[tx_id],
        cleavage_params=cleavage_params, truncate_sec=truncate_sec,
        w2f=w2f_reassignment
    )
    dgraphs:TypeDGraphs = (None, {}, {})
    pgraphs:TypePGraphs = (None, {}, {})
    if variant_series.transcriptional:
        try:
            if not noncanonical_transcripts or \
                    variant_series.has_any_alternative_splicing():
                main_peptides, dgraph, pgraph = call_peptide_main(
                    tx_id=tx_id, tx_variants=variant_series.transcriptional,
                    variant_pool=pool, ref=reference_data,
                    tx_seqs=tx_seqs, gene_seqs=gene_seqs,
                    cleavage_params=cleavage_params,
                    max_adjacent_as_mnv=max_adjacent_as_mnv,
                    truncate_sec=truncate_sec,
                    w2f=w2f_reassignment,
                    denylist=denylist,
                    save_graph=save_graph
                )
                peptide_pool.append(main_peptides)
                dgraphs = (dgraph, dgraphs[1], dgraphs[2])
                pgraphs = (pgraph, pgraphs[1], pgraphs[2])
        except:
            logger(f'Exception raised from {tx_id}')
            raise

    peptides = set()
    exclude_variant_types = ['Fusion', 'Insertion', 'Deletion', 'Substitution', 'circRNA']
    for variant in variant_series.fusion:
        try:
            donor_breakpoint_genomic = reference_data.anno.coordinate_transcript_to_genomic(
                index=variant.location.start - 1, transcript=tx_id
            )
            donor_breakpoint_gene = reference_data.anno.coordinate_genomic_to_gene(
                index=donor_breakpoint_genomic, gene=variant.attrs['GENE_ID']
            ) + 1
            filtered_variants = pool.filter_variants(
                tx_ids=[tx_id], start=0, end=donor_breakpoint_gene,
                exclude_type=exclude_variant_types, intron=False,
                return_coord='transcript'
            )
            variant_pool = copy.copy(pool)
            variant_pool[tx_id] = copy.copy(pool[tx_id])
            variant_pool[tx_id].transcriptional = filtered_variants
            _peptides, dgraph, pgraph = call_peptide_fusion(
                variant=variant, variant_pool=variant_pool, ref=reference_data,
                tx_seqs=tx_seqs, gene_seqs=gene_seqs, cleavage_params=cleavage_params,
                max_adjacent_as_mnv=max_adjacent_as_mnv, w2f_reassignment=w2f_reassignment,
                denylist=denylist, save_graph=save_graph
            )
            dgraphs[1][variant.id] = dgraph
            pgraphs[1][variant.id] = pgraph
        except:
            logger(
                f"Exception raised from fusion {variant.id}, "
                f"donor:{tx_id}, accepter: {variant.attrs['ACCEPTER_TRANSCRIPT_ID']}"
            )
            raise

        peptides.update(_peptides)
    peptide_pool.append(peptides)
    peptides = set()
    if main_peptides:
        denylist.update([str(x.seq) for x in main_peptides])
    for circ_model in variant_series.circ_rna:
        try:
            _peptides, cgraph, pgraph = call_peptide_circ_rna(
                record=circ_model, variant_pool=pool, gene_seqs=gene_seqs,
                cleavage_params=cleavage_params,
                max_adjacent_as_mnv=max_adjacent_as_mnv,
                w2f_reassignment=w2f_reassignment, denylist=denylist,
                save_graph=save_graph
            )

        except:
            logger(f"Exception raised from {circ_model.id}")
            raise
        peptides.update(_peptides)
        dgraphs[2][circ_model.id] = cgraph
        pgraphs[2][circ_model.id] = pgraph
    peptide_pool.append(peptides)

    return peptide_pool, tx_id, dgraphs, pgraphs

def wrapper(dispatch):
    """ wrapper for ParallelPool """
    return call_variant_peptides_wrapper(**dispatch)

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
            process_pool = ParallelPool(ncpus=caller.threads)

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

            dispatch = {
                'tx_id': tx_id,
                'variant_series': variant_series,
                'tx_seqs': tx_seqs,
                'gene_seqs': gene_seqs,
                'reference_data': reference_data,
                'pool': dummy_pool,
                'cleavage_params': caller.cleavage_params,
                'noncanonical_transcripts': noncanonical_transcripts,
                'max_adjacent_as_mnv': caller.max_adjacent_as_mnv,
                'truncate_sec': caller.truncate_sec,
                'w2f_reassignment': caller.w2f_reassignment,
                'save_graph': caller.graph_output_dir is not None
            }
            dispatches.append(dispatch)

            reloaded = ((i + 1) % caller.threads == 0 or i + 1 == len(tx_sorted)) \
                and len(dispatches) > 0
            if reloaded:
                if caller.verbose >= 2:
                    logger([x['tx_id'] for x in dispatches])
                if caller.threads > 1:
                    results = process_pool.map(wrapper, dispatches)
                else:
                    results = [wrapper(dispatches[0])]

                # pylint: disable=W0621
                for peptide_series, tx_id, _, pgraphs in results:
                    for peptides in peptide_series:
                        for peptide in peptides:
                            caller.variant_peptides.add_peptide(
                                peptide=peptide,
                                canonical_peptides=ref.canonical_peptides,
                                cleavage_params=caller.cleavage_params
                            )
                    if caller.graph_output_dir is not None:
                        caller.write_pgraphs(tx_id, pgraphs)
                dispatches = []

            if caller.verbose >= 1:
                i += 1
                if i % 1000 == 0:
                    logger(f'{i} transcripts processed.')

    caller.write_fasta()


def call_canonical_peptides(tx_id:str, ref:params.ReferenceData,
        tx_seq:dna.DNASeqRecordWithCoordinates,
        cleavage_params:params.CleavageParams,
        truncate_sec:bool, w2f:bool):
    """ Call canonical peptides """
    tx_model = ref.anno.transcripts[tx_id]
    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq, _id=tx_id,
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=False,
        cleavage_params=cleavage_params
    )
    dgraph.gather_sect_variants(ref.anno)
    dgraph.init_three_frames()
    pgraph = dgraph.translate()
    pgraph.create_cleavage_graph()
    peptides = pgraph.call_variant_peptides(
        check_variants=False, truncate_sec=truncate_sec, w2f=w2f,
        check_external_variants=False
    )
    return {str(x.seq) for x in peptides}

def call_peptide_main(tx_id:str, tx_variants:List[seqvar.VariantRecord],
        variant_pool:seqvar.VariantRecordPoolOnDisk,
        ref:params.ReferenceData,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        cleavage_params:params.CleavageParams,
        max_adjacent_as_mnv:bool, truncate_sec:bool, w2f:bool,
        denylist:Set[str], save_graph:bool
        ) -> Tuple[Set[aa.AminoAcidSeqRecord], ThreeFrameTVG, PeptideVariantGraph]:
    """ Call variant peptides for main variants (except cirRNA). """
    tx_model = ref.anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        cleavage_params=cleavage_params,
        max_adjacent_as_mnv=max_adjacent_as_mnv
    )
    dgraph.gather_sect_variants(ref.anno)
    dgraph.init_three_frames()
    dgraph.create_variant_graph(
        variants=tx_variants, variant_pool=variant_pool, genome=ref.genome,
        anno=ref.anno, tx_seqs=tx_seqs, gene_seqs=gene_seqs
    )
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()

    pgraph.create_cleavage_graph()
    peptides = pgraph.call_variant_peptides(
        denylist=denylist,
        truncate_sec=truncate_sec,
        w2f=w2f,
        check_external_variants=True
    )
    if not save_graph:
        dgraph, pgraph = None, None
    return peptides, dgraph, pgraph

def call_peptide_fusion(variant:seqvar.VariantRecord,
        variant_pool:seqvar.VariantRecordPoolOnDisk,
        ref:params.ReferenceData,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        cleavage_params:params.CleavageParams, max_adjacent_as_mnv:bool,
        w2f_reassignment:bool, denylist:Set[str], save_graph:bool
        ) -> Tuple[Set[aa.AminoAcidSeqRecord], ThreeFrameTVG, PeptideVariantGraph]:
    """ Call variant peptides for fusion """
    tx_id = variant.location.seqname
    tx_model = ref.anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]
    tx_seq = tx_seq[:variant.location.start]
    donor_tx_id = variant.accepter_transcript_id

    orf_start = tx_seq.orf.start + 3 if tx_seq.orf else 3
    if variant.location.start < orf_start:
        return [], None, None

    if tx_id in variant_pool:
        tx_variants = [x for x in variant_pool[tx_id].transcriptional
            if x.location.end < variant.location.start]
    else:
        tx_variants = []

    tx_variants.append(variant)

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=ref.anno.transcripts[donor_tx_id].is_mrna_end_nf(),
        cleavage_params=cleavage_params,
        max_adjacent_as_mnv=max_adjacent_as_mnv
    )
    dgraph.gather_sect_variants(ref.anno)
    dgraph.sect_variants = [v for v in dgraph.sect_variants
        if v.location.end < variant.location.start]
    dgraph.init_three_frames()
    dgraph.create_variant_graph(
        variants=tx_variants, variant_pool=variant_pool, genome=ref.genome,
        anno=ref.anno, tx_seqs=tx_seqs, gene_seqs=gene_seqs
    )
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()
    pgraph.create_cleavage_graph()
    peptides = pgraph.call_variant_peptides(
        denylist=denylist,
        w2f=w2f_reassignment,
        check_external_variants=True
    )
    if not save_graph:
        dgraph, pgraph = None, None
    return peptides, dgraph, pgraph

def call_peptide_circ_rna(record:circ.CircRNAModel,
        variant_pool:seqvar.VariantRecordPool,
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        cleavage_params:params.CleavageParams, max_adjacent_as_mnv:bool,
        w2f_reassignment:bool, denylist:Set[str], save_graph:bool
        )-> Tuple[Set[aa.AminoAcidSeqRecord], ThreeFrameCVG, PeptideVariantGraph]:
    """ Call variant peptides from a given circRNA """
    gene_id = record.gene_id
    gene_seq = gene_seqs[gene_id]
    record.fragments.sort()
    circ_seq = record.get_circ_rna_sequence(gene_seq)
    for loc in circ_seq.locations:
        loc.ref.seqname = record.id

    # Alternative splicing should not be included. Alternative splicing
    # represented as Insertion, Deletion or Substitution.
    exclusion_variant_types = ['Insertion', 'Deletion', 'Substitution']

    fragments = []
    for frag in record.fragments:
        if len(frag) <= 3:
            continue
        loc = FeatureLocation(
            start=frag.location.start + 3, end=frag.location.end
        )
        frag = SeqFeature(chrom=frag.chrom, location=loc, attributes=frag.attributes)
        fragments.append(frag)

    if not fragments:
        return set(), None, None

    variant_records = variant_pool.filter_variants(
        tx_ids=[record.transcript_id], exclude_type=exclusion_variant_types,
        intron=False, segments=fragments
    )

    cgraph = svgraph.ThreeFrameCVG(
        circ_seq, _id=record.id, circ_record=record,
        cleavage_params=cleavage_params,
        max_adjacent_as_mnv=max_adjacent_as_mnv
    )
    cgraph.init_three_frames()
    cgraph.create_variant_circ_graph(variant_records)
    cgraph.extend_loop()
    cgraph.truncate_three_frames()
    cgraph.fit_into_codons()
    pgraph = cgraph.translate()
    pgraph.create_cleavage_graph()
    peptides = pgraph.call_variant_peptides(
        denylist=denylist, circ_rna=record,
        w2f=w2f_reassignment, check_external_variants=True
    )
    if not save_graph:
        cgraph, pgraph = None, None
    return peptides, cgraph, pgraph
