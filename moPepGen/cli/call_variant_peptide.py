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
from typing import List, Set, TYPE_CHECKING, Dict, Tuple, Union
from pathlib import Path
import os
import dataclasses
from contextlib import ExitStack
from pathos.pools import ParallelPool
from moPepGen import svgraph, aa, seqvar, gtf, params, get_logger
from moPepGen.cli import common
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.svgraph import ThreeFrameTVG


if TYPE_CHECKING:
    from logging import Logger
    from Bio.Seq import Seq
    from moPepGen import dna, circ
    from moPepGen.svgraph import ThreeFrameCVG, PeptideVariantGraph
    from moPepGen.svgraph.VariantPeptideDict import AnnotatedPeptideLabel
    from moPepGen.gtf import GeneAnnotationModel
    from moPepGen.params import CodonTableInfo, CleavageParams


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
        '--backsplicing-only',
        action='store_true',
        help='For circRNA, only keep noncanonical peptides spaning the backsplicing site.'
    )
    p.add_argument(
        '--coding-novel-orf',
        action='store_true',
        help='Find alternative start site for coding transcripts.'
    )
    p.add_argument(
        '--max-variants-per-node',
        type=int,
        help='Maximal number of variants per node. This argument can be useful'
        ' when there are local regions that are heavily mutated. When creating'
        ' the cleavage graph, nodes containing variants larger than this value'
        ' are skipped. Setting to -1 will avoid checking for this.',
        default=(7,),
        nargs='+',
        metavar='<number>'
    )
    p.add_argument(
        '--additional-variants-per-misc',
        type=int,
        help='Additional variants allowed for every miscleavage. This argument'
        ' is used together with --max-variants-per-node to handle hypermutated'
        ' regions. Setting to -1 will avoid checking for this.',
        default=(2,),
        nargs='+',
        metavar='<number>'
    )
    p.add_argument(
        '--in-bubble-cap-step-down',
        type=int,
        default=0,
        metavar='<number>',
        help='In bubble variant caps default step down.'
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
        '--timeout-seconds',
        type=int,
        default=1800,
        help='Time out in seconds for each transcript.'
    )
    p.add_argument(
        '--threads',
        type=int,
        help='Set number of threads to be used.',
        default=1,
        metavar='<number>'
    )
    common.add_args_skip_failed(p)
    common.add_args_reference(p)
    common.add_args_cleavage(p)
    common.add_args_debug_level(p)

    p.set_defaults(func=call_variant_peptide)
    common.print_help_if_missing_args(p)
    return p


class TallyTable:
    """ Tally table for callVariant """
    def __init__(self, logger:Logger):
        """ constructor """
        self.n_transcripts_total = 0
        self.n_transcripts_processed = 0
        self.n_transcripts_invalid = 0
        self.n_transcripts_failed = {
            'variant': 0,
            'fusion': 0,
            'circRNA': 0
        }
        self.n_total_peptides = 0
        self.n_valid_peptides = 0
        self.logger = logger

    def log(self):
        """ log """
        self.logger.info("callVariant summary:")
        self.logger.info("Total transcripts with variants: %i", self.n_transcripts_total)
        self.logger.info("Total transcripts processed: %i", self.n_transcripts_processed)
        self.logger.info("Number of invalid transcripts: %i", self.n_transcripts_invalid)
        self.logger.info("Number of transcripts failed when calling for:")
        self.logger.info("    Variant peptides: %i", self.n_transcripts_failed['variant'])
        self.logger.info("    Fusion peptides: %i", self.n_transcripts_failed['fusion'])
        self.logger.info("    circRNA peptides: %i", self.n_transcripts_failed['circRNA'])
        self.logger.info(
            "Total variant peptides generated (including redundant): %i",
            self.n_total_peptides
        )
        self.logger.info("Total variant peptides saved: %i", self.n_valid_peptides)

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
        self.peptide_table_temp_path:Path = common.get_peptide_table_path_temp(
            path=Path(self.output_path)
        )
        self.peptide_table_output_path:Path = common.get_peptide_table_path(
            path=Path(self.output_path)
        )
        self.graph_output_dir:Path = args.graph_output_dir
        if self.graph_output_dir is not None:
            common.validate_file_format(
                self.graph_output_dir, check_writable=True, is_directory=True
            )

        self.cleavage_params = params.CleavageParams(
            enzyme=args.cleavage_rule,
            exception=args.cleavage_exception,
            miscleavage=int(args.miscleavage),
            min_mw=float(args.min_mw),
            min_length=args.min_length,
            max_length=args.max_length,
            max_variants_per_node = args.max_variants_per_node[0],
            additional_variants_per_misc = args.additional_variants_per_misc[0],
            in_bubble_cap_step_down=args.in_bubble_cap_step_down,
            min_nodes_to_collapse = args.min_nodes_to_collapse,
            naa_to_collapse = args.naa_to_collapse
        )
        self.max_adjacent_as_mnv = args.max_adjacent_as_mnv
        self.truncate_sec:bool = args.selenocysteine_termination
        self.backsplicing_only:bool = args.backsplicing_only
        self.w2f_reassignment = args.w2f_reassignment
        self.noncanonical_transcripts = args.noncanonical_transcripts
        self.invalid_protein_as_noncoding:bool = args.invalid_protein_as_noncoding
        self.threads = args.threads
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

        self.logger:Logger = None
        self.tally:TallyTable = None

    def load_reference(self):
        """ load reference genome, annotation, and canonical peptides """
        self.reference_data = common.load_references(
            args=self.args,
            invalid_protein_as_noncoding=self.invalid_protein_as_noncoding,
            cleavage_params=self.cleavage_params,
            load_codon_tables=True
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
        self.logger.info('Variant peptide FASTA file written to disk.')

    def write_dgraphs(self, tx_id:str, dgraphs:TypeDGraphs):
        """ Write ThreeFrameTVG data """
        if dgraphs[0]:
            data = dgraphs[0].jsonfy()
            with open(self.graph_output_dir/f"{tx_id}_main_TVG.json", 'wt') as handle:
                json.dump(data, handle)

        for var_id, pgraph in dgraphs[1].items():
            if pgraph is None:
                continue
            data = pgraph.jsonfy()
            with open(self.graph_output_dir/f"{tx_id}_Fusion_{var_id}_TVG.json", 'wt') as handle:
                json.dump(data, handle)

        for var_id, pgraph in dgraphs[2].items():
            if pgraph is None:
                continue
            data = pgraph.jsonfy()
            with open(self.graph_output_dir/f"{tx_id}_circRNA_{var_id}_TVG.json", 'wt') as handle:
                json.dump(data, handle)

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

    def gather_data_for_call_variant(self, tx_id:str, pool:seqvar.VariantRecordPoolOnDisk):
        """ Gather data for callVariant """
        ref = self.reference_data
        tx_ids = [tx_id]
        tx_model = ref.anno.transcripts[tx_id]
        codon_table = ref.codon_tables[tx_model.transcript.chrom]
        try:
            variant_series = pool[tx_id]
        except ValueError:
            if self.args.skip_failed:
                self.tally.n_transcripts_invalid += 1
                return None
            raise
        if variant_series.is_empty():
            return None
        if self.noncanonical_transcripts and \
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

        gene_models:dict[str,GeneAnnotationModel] = {}
        for gene_id in gene_ids:
            dummy_gene_model = copy.deepcopy(ref.anno.genes[gene_id])
            dummy_gene_model.transcripts = [x for x in dummy_gene_model.transcripts
                if x in tx_ids]
            gene_models[gene_id] = dummy_gene_model

        dummy_anno = gtf.GenomicAnnotation(
            genes=gene_models,
            transcripts={tx_id: ref.anno.transcripts[tx_id] for tx_id in tx_ids},
            source=ref.anno.source
        )

        reference_data = params.ReferenceData(
            genome=None,
            anno=dummy_anno,
            canonical_peptides=set(),    # Not providing canonical peptide at this point.
            codon_tables={}
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

        return {
            'tx_id': tx_id,
            'variant_series': variant_series,
            'tx_seqs': tx_seqs,
            'gene_seqs': gene_seqs,
            'reference_data': reference_data,
            'codon_table': codon_table,
            'pool': dummy_pool,
            'cleavage_params': self.cleavage_params,
            'noncanonical_transcripts': self.noncanonical_transcripts,
            'max_adjacent_as_mnv': self.max_adjacent_as_mnv,
            'truncate_sec': self.truncate_sec,
            'w2f_reassignment': self.w2f_reassignment,
            'save_graph': self.graph_output_dir is not None,
            'timeout': self.args.timeout_seconds,
            'max_variants_per_node': tuple(self.args.max_variants_per_node),
            'additional_variants_per_misc': tuple(self.args.additional_variants_per_misc),
            'in_bubble_cap_step_down': self.args.in_bubble_cap_step_down,
            'backsplicing_only': self.args.backsplicing_only,
            'coding_novel_orf': self.args.coding_novel_orf,
            'skip_failed': self.args.skip_failed
        }

if TYPE_CHECKING:
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
    CallVariantOutput = Tuple[
        Dict[Seq, List[AnnotatedPeptideLabel]],
        str,
        TypeDGraphs,
        TypePGraphs,
        Tuple[bool,bool,bool]
    ]
# pylint: disable=unused-argument
@common.timeout()
def call_variant_peptides_wrapper(tx_id:str,
        variant_series:seqvar.TranscriptionalVariantSeries,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        reference_data:params.ReferenceData,
        codon_table:CodonTableInfo,
        pool:seqvar.VariantRecordPool,
        cleavage_params:params.CleavageParams,
        noncanonical_transcripts:bool,
        max_adjacent_as_mnv:bool,
        truncate_sec:bool,
        w2f_reassignment:bool,
        backsplicing_only:bool,
        save_graph:bool,
        coding_novel_orf:bool,
        skip_failed:bool,
        tracer:TimeoutTracer,
        **kwargs
        ) -> CallVariantOutput:
    """ wrapper function to call variant peptides """
    # pylint: disable=W0702
    logger = get_logger()
    peptide_anno:Dict[Seq, Dict[str, AnnotatedPeptideLabel]] = {}

    def add_peptide_anno(x:Dict[Seq, List[AnnotatedPeptideLabel]]):
        for seq, seq_data in x.items():
            val = peptide_anno.setdefault(seq, {})
            for metadata in seq_data:
                if metadata.label not in val:
                    val[metadata.label] = metadata

    main_peptides = None
    denylist = call_canonical_peptides(
        tx_id=tx_id, ref=reference_data, tx_seq=tx_seqs[tx_id],
        cleavage_params=cleavage_params,
        codon_table=codon_table,
        truncate_sec=truncate_sec, w2f=w2f_reassignment
    )
    dgraphs:TypeDGraphs = (None, {}, {})
    pgraphs:TypePGraphs = (None, {}, {})
    success_flags = (True, True, True)
    tracer.backbone = 'main'
    if variant_series.transcriptional:
        try:
            if not noncanonical_transcripts or \
                    variant_series.has_any_alternative_splicing():
                peptide_map, dgraph, pgraph = call_peptide_main(
                    tx_id=tx_id, tx_variants=variant_series.transcriptional,
                    variant_pool=pool, ref=reference_data, codon_table=codon_table,
                    tx_seqs=tx_seqs, gene_seqs=gene_seqs,
                    cleavage_params=cleavage_params,
                    max_adjacent_as_mnv=max_adjacent_as_mnv,
                    truncate_sec=truncate_sec,
                    w2f=w2f_reassignment,
                    denylist=denylist,
                    save_graph=save_graph,
                    coding_novel_orf=coding_novel_orf,
                    tracer=tracer
                )
                main_peptides = set(peptide_map.keys())
                dgraphs = (dgraph, dgraphs[1], dgraphs[2])
                pgraphs = (pgraph, pgraphs[1], pgraphs[2])
                add_peptide_anno(peptide_map)
        except Exception as e:
            if  isinstance(e, TimeoutError):
                raise
            if skip_failed:
                logger.warning('Variant peptides calling failed from %s', tx_id)
                success_flags = (False, success_flags[1], success_flags[2])
            else:
                logger.error('Exception raised from %s', tx_id)
                raise

    tracer.backbone = 'fusion'
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
            peptide_map, dgraph, pgraph = call_peptide_fusion(
                variant=variant, variant_pool=variant_pool, ref=reference_data,
                codon_table=codon_table, tx_seqs=tx_seqs, gene_seqs=gene_seqs,
                cleavage_params=cleavage_params, max_adjacent_as_mnv=max_adjacent_as_mnv,
                w2f_reassignment=w2f_reassignment, denylist=denylist, save_graph=save_graph,
                coding_novel_orf=coding_novel_orf, tracer=tracer
            )
            dgraphs[1][variant.id] = dgraph
            pgraphs[1][variant.id] = pgraph
            add_peptide_anno(peptide_map)
        except Exception as e:
            if  isinstance(e, TimeoutError):
                raise
            if skip_failed:
                logger.warning(
                    'Variant peptides calling failed from %s with fusion: %s',
                    tx_id, variant.id
                )
                success_flags = (success_flags[0], False, success_flags[2])
            else:
                logger.error(
                    "Exception raised from fusion %s, donor: %s, accepter: %s",
                    variant.id, tx_id, variant.attrs['ACCEPTER_TRANSCRIPT_ID']
                )
                raise

    tracer.backbone = 'circRNA'
    if main_peptides:
        denylist.update([str(x) for x in main_peptides])
    for circ_model in variant_series.circ_rna:
        try:
            peptide_map, cgraph, pgraph = call_peptide_circ_rna(
                record=circ_model, variant_pool=pool, gene_seqs=gene_seqs,
                cleavage_params=cleavage_params, codon_table=codon_table,
                max_adjacent_as_mnv=max_adjacent_as_mnv,
                w2f_reassignment=w2f_reassignment, denylist=denylist,
                save_graph=save_graph, backsplicing_only=backsplicing_only,
                tracer=tracer
            )

        except Exception as e:
            if  isinstance(e, TimeoutError):
                raise
            if skip_failed:
                logger.warning(
                    'Variant peptides calling failed from %s with circRNA: %s',
                    tx_id, circ_model.id
                )
                success_flags = (success_flags[0], success_flags[1], False)
            else:
                logger.error("Exception raised from %s", circ_model.id)
                raise
        dgraphs[2][circ_model.id] = cgraph
        pgraphs[2][circ_model.id] = pgraph
        add_peptide_anno(peptide_map)

    peptide_anno = {k:list(v.values()) for k,v in peptide_anno.items()}

    return peptide_anno, tx_id, dgraphs, pgraphs, success_flags

@dataclasses.dataclass
class TimeoutTracer:
    """ Tracking time out """
    backbone:str = None
    graph:str = None

def caller_reducer(dispatch):
    """ wrapper for ParallelPool. Also reduces the complexity if the run is timed out. """
    max_variants_per_node = dispatch['max_variants_per_node']
    additional_variants_per_misc = dispatch['additional_variants_per_misc']
    in_bubble_cap_step_down = dispatch['in_bubble_cap_step_down']
    tx_id = dispatch['tx_id']
    tracer = TimeoutTracer()
    dispatch['tracer'] = tracer
    while True:
        try:
            return call_variant_peptides_wrapper(**dispatch)
        except TimeoutError as e:
            new_dispatch = copy.copy(dispatch)
            p:CleavageParams = copy.copy(new_dispatch['cleavage_params'])
            if tracer.graph == 'TVG':
                in_bubble_cap_step_down += 1
                if in_bubble_cap_step_down > len(ThreeFrameTVG.VARIANT_BUBBLE_CAPS):
                    raise ValueError(
                        f"Failed to finish transcript: {tx_id} "
                        f"with in_bubble_cap_step_down: {in_bubble_cap_step_down}"
                    ) from e
            else:
                max_variants_per_node = max_variants_per_node[1:]
                if len(max_variants_per_node) == 0:
                    max_variants_per_node = (p.max_variants_per_node - 1, )
                    if max_variants_per_node[0] <= 0:
                        raise ValueError(f"Failed to finish transcript: {tx_id}") from e
                additional_variants_per_misc = additional_variants_per_misc[1:]
                if len(additional_variants_per_misc) == 0:
                    additional_variants_per_misc = (0,)
            p.max_variants_per_node = max_variants_per_node[0]
            p.additional_variants_per_misc = additional_variants_per_misc[0]
            p.in_bubble_cap_step_down = in_bubble_cap_step_down
            new_dispatch['cleavage_params'] = p
            new_dispatch['tracer'] = TimeoutTracer()
            dispatch = new_dispatch
            get_logger().warning(
                "Transcript %s timed out. Retry with "
                "--max-variants-per-node %i --additional-variants-per-misc %i"
                "and `--in-bubble-cap-step-down` %i",
                tx_id, p.max_variants_per_node, p.additional_variants_per_misc,
                in_bubble_cap_step_down
            )

def call_variant_peptide(args:argparse.Namespace) -> None:
    """ Main entry point for calling variant peptide """
    caller = VariantPeptideCaller(args)
    common.print_start_message(args)
    caller.load_reference()
    caller.create_in_disk_variant_pool()
    ref = caller.reference_data
    logger = get_logger()
    caller.logger = logger
    caller.tally = TallyTable(logger)

    with ExitStack() as stack:
        pool = stack.enter_context(seqvar.VariantRecordPoolOnDiskOpener(caller.variant_record_pool))
        seq_anno_handle = stack.enter_context(open(caller.peptide_table_temp_path, 'w+'))

        peptide_table = svgraph.VariantPeptideTable(seq_anno_handle)
        peptide_table.write_header()

        for file in pool.gvf_files:
            logger.info("Using GVF file: %s", file)
        tx_rank = ref.anno.get_transcript_rank()
        tx_sorted = sorted(pool.pointers.keys(), key=lambda x:tx_rank[x])
        logger.info('Variants sorted')
        caller.tally.n_transcripts_total = len(tx_sorted)

        if caller.threads > 1:
            process_pool = ParallelPool(ncpus=caller.threads)

        logger.info('Start calling variant peptides..')

        dispatches = []
        i = 0
        for tx_id in tx_sorted:
            dispatch = caller.gather_data_for_call_variant(tx_id, pool)
            if not dispatch:
                continue
            dispatches.append(dispatch)

            caller.tally.n_transcripts_processed += 1

            reloaded = ((i + 1) % caller.threads == 0 or i + 1 == len(tx_sorted)) \
                and len(dispatches) > 0
            if reloaded:
                logger.debug([x['tx_id'] for x in dispatches])
                if caller.threads > 1:
                    results = process_pool.map(
                        caller_reducer,
                        dispatches
                    )
                else:
                    results = [ caller_reducer(dispatches[0]) ]

                # pylint: disable=W0621
                for peptide_anno, tx_id, dgraphs, pgraphs, success_flags in results:
                    caller.tally.n_total_peptides += len(peptide_anno)
                    for peptide in peptide_anno:
                        is_valid = peptide_table.is_valid(
                            seq=peptide,
                            canonical_peptides=ref.canonical_peptides,
                            cleavage_params=caller.cleavage_params
                        )
                        if is_valid:
                            for seq_anno in peptide_anno[peptide]:
                                peptide_table.add_peptide(peptide, seq_anno)
                    if caller.graph_output_dir is not None:
                        caller.write_dgraphs(tx_id, dgraphs)
                        caller.write_pgraphs(tx_id, pgraphs)
                    if not success_flags[0]:
                        caller.tally.n_transcripts_failed['variant'] += 1
                    if not success_flags[1]:
                        caller.tally.n_transcripts_failed['fusion'] += 1
                    if not success_flags[2]:
                        caller.tally.n_transcripts_failed['circRNA'] += 1
                dispatches = []

            i += 1
            if i % 1000 == 0:
                logger.info(
                    '%.2f %% ( %i / %i ) transcripts processed.',
                    i/len(tx_sorted) * 100, i, len(tx_sorted)
                )
        caller.tally.n_valid_peptides = len(peptide_table.index)
        logger.info('All transcripts processed. Preparing FASTA file..')
        peptide_table.write_fasta(caller.output_path)
        logger.info('Variant peptide FASTA written.')
        peptide_table.sort_table(caller.peptide_table_output_path)
        logger.info('Variant peptide table soted.')
    os.remove(caller.peptide_table_temp_path)
    caller.tally.log()

def call_canonical_peptides(tx_id:str, ref:params.ReferenceData,
        tx_seq:dna.DNASeqRecordWithCoordinates,
        cleavage_params:params.CleavageParams,
        codon_table:CodonTableInfo, truncate_sec:bool, w2f:bool):
    """ Call canonical peptides """
    tx_model = ref.anno.transcripts[tx_id]
    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq, _id=tx_id,
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=False,
        cleavage_params=cleavage_params,
        coordinate_feature_type='transcript',
        coordinate_feature_id=tx_id
    )
    dgraph.gather_sect_variants(ref.anno)
    dgraph.init_three_frames()
    pgraph = dgraph.translate(
        table=codon_table.codon_table,
        start_codons=codon_table.start_codons
    )
    pgraph.create_cleavage_graph()
    peptide_map = pgraph.call_variant_peptides(
        check_variants=False, truncate_sec=truncate_sec, w2f=w2f,
        check_external_variants=False
    )
    return {str(x) for x in peptide_map}

if TYPE_CHECKING:
    TypeCallPeptideReturnData = Tuple[
        Dict[Seq, List[AnnotatedPeptideLabel]],
        Union[ThreeFrameTVG, ThreeFrameCVG],
        PeptideVariantGraph
    ]

def call_peptide_main(tx_id:str, tx_variants:List[seqvar.VariantRecord],
        variant_pool:seqvar.VariantRecordPoolOnDisk,
        ref:params.ReferenceData, codon_table:CodonTableInfo,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        cleavage_params:params.CleavageParams,
        max_adjacent_as_mnv:bool, truncate_sec:bool, w2f:bool,
        denylist:Set[str], save_graph:bool, coding_novel_orf:bool,
        tracer:TimeoutTracer
        ) -> TypeCallPeptideReturnData:
    """ Call variant peptides for main variants (except cirRNA). """
    tx_model = ref.anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]

    tracer.graph = 'TVG'

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        cleavage_params=cleavage_params,
        max_adjacent_as_mnv=max_adjacent_as_mnv,
        coordinate_feature_type='transcript',
        coordinate_feature_id=tx_id
    )
    dgraph.gather_sect_variants(ref.anno)
    dgraph.init_three_frames()
    dgraph.create_variant_graph(
        variants=tx_variants, variant_pool=variant_pool, genome=ref.genome,
        anno=ref.anno, tx_seqs=tx_seqs, gene_seqs=gene_seqs
    )
    dgraph.fit_into_codons()
    pgraph = dgraph.translate(
        table=codon_table.codon_table,
        start_codons=codon_table.start_codons
    )

    tracer.graph = 'PVG'

    pgraph.create_cleavage_graph()

    if tx_model.is_protein_coding:
        peptide_map = pgraph.call_variant_peptides(
            denylist=denylist,
            truncate_sec=truncate_sec,
            w2f=w2f,
            check_external_variants=True,
            check_orf=False
        )
    else:
        peptide_map = {}

    if not tx_model.is_protein_coding or coding_novel_orf:
        peptide_novel_orf = pgraph.call_variant_peptides(
            denylist=denylist,
            truncate_sec=truncate_sec,
            w2f=w2f,
            check_external_variants=True,
            check_orf=True
        )
        for peptide, labels in peptide_novel_orf.items():
            if peptide not in peptide_map:
                peptide_map[peptide] = labels

    if not save_graph:
        dgraph, pgraph = None, None
    return peptide_map, dgraph, pgraph

def call_peptide_fusion(variant:seqvar.VariantRecord,
        variant_pool:seqvar.VariantRecordPoolOnDisk,
        ref:params.ReferenceData, codon_table:CodonTableInfo,
        tx_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        cleavage_params:params.CleavageParams, max_adjacent_as_mnv:bool,
        w2f_reassignment:bool, denylist:Set[str], save_graph:bool,
        coding_novel_orf:bool, tracer:TimeoutTracer
        ) -> TypeCallPeptideReturnData:
    """ Call variant peptides for fusion """
    tx_id = variant.location.seqname
    tx_model = ref.anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]
    tx_seq = tx_seq[:variant.location.start]
    donor_tx_id = variant.accepter_transcript_id

    tracer.graph = 'TVG'

    orf_start = tx_seq.orf.start + 3 if tx_seq.orf else 3
    if variant.location.start < orf_start:
        return {}, None, None

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
        max_adjacent_as_mnv=max_adjacent_as_mnv,
        coordinate_feature_type='transcript',
        coordinate_feature_id=tx_id
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
    pgraph = dgraph.translate(
        table=codon_table.codon_table,
        start_codons=codon_table.start_codons
    )

    tracer.graph = 'PVG'

    pgraph.create_cleavage_graph()

    if tx_model.is_protein_coding:
        peptide_map = pgraph.call_variant_peptides(
            denylist=denylist,
            w2f=w2f_reassignment,
            check_external_variants=True,
            check_orf=False
        )
    else:
        peptide_map = {}

    if not tx_model.is_protein_coding or coding_novel_orf:
        peptide_novel_orf = pgraph.call_variant_peptides(
            denylist=denylist,
            w2f=w2f_reassignment,
            check_external_variants=True,
            check_orf=True
        )
        for peptide, labels in peptide_novel_orf.items():
            if peptide not in peptide_map:
                peptide_map[peptide] = labels

    if not save_graph:
        dgraph, pgraph = None, None
    return peptide_map, dgraph, pgraph

def call_peptide_circ_rna(record:circ.CircRNAModel,
        variant_pool:seqvar.VariantRecordPool,
        gene_seqs:Dict[str, dna.DNASeqRecordWithCoordinates],
        codon_table:CodonTableInfo,
        cleavage_params:params.CleavageParams, max_adjacent_as_mnv:bool,
        backsplicing_only:bool, w2f_reassignment:bool, denylist:Set[str],
        save_graph:bool, tracer:TimeoutTracer
        )-> TypeCallPeptideReturnData:
    """ Call variant peptides from a given circRNA """
    gene_id = record.gene_id
    gene_seq = gene_seqs[gene_id]
    record.fragments.sort()
    circ_seq = record.get_circ_rna_sequence(gene_seq)
    for loc in circ_seq.locations:
        loc.ref.seqname = record.id

    tracer.graph = 'TVG'

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
        return {}, None, None

    variant_records = variant_pool.filter_variants(
        tx_ids=[record.transcript_id], exclude_type=exclusion_variant_types,
        intron=False, segments=fragments
    )

    cgraph = svgraph.ThreeFrameCVG(
        circ_seq, _id=record.id, circ_record=record,
        cleavage_params=cleavage_params,
        max_adjacent_as_mnv=max_adjacent_as_mnv,
        coordinate_feature_type='gene',
        coordinate_feature_id=gene_id
    )
    cgraph.init_three_frames()
    cgraph.create_variant_circ_graph(variant_records)
    cgraph.extend_loop()
    cgraph.truncate_three_frames()
    cgraph.fit_into_codons()
    pgraph = cgraph.translate(
        table=codon_table.codon_table,
        start_codons=codon_table.start_codons
    )

    tracer.graph = 'PVG'

    pgraph.create_cleavage_graph()
    peptide_map = pgraph.call_variant_peptides(
        denylist=denylist, circ_rna=record, backsplicing_only=backsplicing_only,
        w2f=w2f_reassignment, check_external_variants=True
    )
    if not save_graph:
        cgraph, pgraph = None, None
    return peptide_map, cgraph, pgraph
