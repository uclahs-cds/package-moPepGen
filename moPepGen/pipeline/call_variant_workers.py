from __future__ import annotations
from typing import Dict, List, Set, Tuple, TYPE_CHECKING
from moPepGen import svgraph, get_logger, params
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.util.timeout import timeout
from .models import TimeoutTracer, GraphPhase

if TYPE_CHECKING:
    from Bio.Seq import Seq
    from moPepGen import seqvar, dna, circ
    from moPepGen.params import CodonTableInfo
    from moPepGen.svgraph import ThreeFrameTVG, ThreeFrameCVG, PeptideVariantGraph
    from moPepGen.svgraph.VariantPeptideDict import AnnotatedPeptideLabel

# Type aliases
TypeDGraphs = Tuple[
    'ThreeFrameTVG',
    Dict[str, 'ThreeFrameTVG'],
    Dict[str, 'ThreeFrameCVG']
]
TypePGraphs = Tuple[
    'PeptideVariantGraph',
    Dict[str, 'PeptideVariantGraph'],
    Dict[str, 'PeptideVariantGraph']
]
TypeCallPeptideReturnData = Tuple[
    Dict['Seq', List['AnnotatedPeptideLabel']],
    'ThreeFrameTVG | ThreeFrameCVG',
    'PeptideVariantGraph'
]


def call_canonical_peptides(tx_id: str, ref: params.ReferenceData,
        tx_seq: 'dna.DNASeqRecordWithCoordinates',
        cleavage_params: params.CleavageParams,
        codon_table: 'CodonTableInfo', truncate_sec: bool, w2f: bool):
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


def call_peptide_main(tx_id: str, tx_variants: List['seqvar.VariantRecord'],
        variant_pool: 'seqvar.VariantRecordPoolOnDisk',
        ref: params.ReferenceData, codon_table: 'CodonTableInfo',
        tx_seqs: Dict[str, 'dna.DNASeqRecordWithCoordinates'],
        gene_seqs: Dict[str, 'dna.DNASeqRecordWithCoordinates'],
        cleavage_params: params.CleavageParams,
        max_adjacent_as_mnv: int, truncate_sec: bool, w2f: bool,
        denylist: Set[str], save_graph: bool, coding_novel_orf: bool,
        tracer: TimeoutTracer
        ) -> TypeCallPeptideReturnData:
    tx_model = ref.anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]

    tracer.graph = GraphPhase.TVG

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

    tracer.graph = GraphPhase.PVG

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


def call_peptide_fusion(variant: 'seqvar.VariantRecord',
        variant_pool: 'seqvar.VariantRecordPoolOnDisk',
        ref: params.ReferenceData, codon_table: 'CodonTableInfo',
        tx_seqs: Dict[str, 'dna.DNASeqRecordWithCoordinates'],
        gene_seqs: Dict[str, 'dna.DNASeqRecordWithCoordinates'],
        cleavage_params: params.CleavageParams, max_adjacent_as_mnv: int,
        w2f_reassignment: bool, denylist: Set[str], save_graph: bool,
        coding_novel_orf: bool, tracer: TimeoutTracer
        ) -> TypeCallPeptideReturnData:
    tx_id = variant.location.seqname
    tx_model = ref.anno.transcripts[tx_id]
    tx_seq = tx_seqs[tx_id]
    tx_seq = tx_seq[:variant.location.start]
    donor_tx_id = variant.accepter_transcript_id

    tracer.graph = GraphPhase.TVG

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

    tracer.graph = GraphPhase.PVG

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


def call_peptide_circ_rna(record: 'circ.CircRNAModel',
        variant_pool: 'seqvar.VariantRecordPool',
        gene_seqs: Dict[str, 'dna.DNASeqRecordWithCoordinates'],
        codon_table: 'CodonTableInfo',
        cleavage_params: params.CleavageParams, max_adjacent_as_mnv: int,
        backsplicing_only: bool, w2f_reassignment: bool, denylist: Set[str],
        save_graph: bool, tracer: TimeoutTracer
        ) -> TypeCallPeptideReturnData:
    gene_id = record.gene_id
    gene_seq = gene_seqs[gene_id]
    record.fragments.sort()
    circ_seq = record.get_circ_rna_sequence(gene_seq)
    for loc in circ_seq.locations:
        loc.ref.seqname = record.id

    tracer.graph = GraphPhase.TVG

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

    tracer.graph = GraphPhase.PVG

    pgraph.create_cleavage_graph()
    peptide_map = pgraph.call_variant_peptides(
        denylist=denylist, circ_rna=record, backsplicing_only=backsplicing_only,
        w2f=w2f_reassignment, check_external_variants=True
    )
    if not save_graph:
        cgraph, pgraph = None, None
    return peptide_map, cgraph, pgraph


# pylint: disable=unused-argument
@timeout()
def call_variant_peptides_wrapper(tx_id: str,
        variant_series: 'seqvar.TranscriptionalVariantSeries',
        tx_seqs: Dict[str, 'dna.DNASeqRecordWithCoordinates'],
        gene_seqs: Dict[str, 'dna.DNASeqRecordWithCoordinates'],
        reference_data: params.ReferenceData,
        codon_table: 'CodonTableInfo',
        pool: 'seqvar.VariantRecordPool',
        cleavage_params: params.CleavageParams,
        noncanonical_transcripts: bool,
        max_adjacent_as_mnv: int,
        truncate_sec: bool,
        w2f_reassignment: bool,
        backsplicing_only: bool,
        save_graph: bool,
        coding_novel_orf: bool,
        skip_failed: bool,
        tracer: TimeoutTracer,
        **kwargs
        ) -> Tuple[
            Dict['Seq', List['AnnotatedPeptideLabel']],
            str,
            TypeDGraphs,
            TypePGraphs,
            Tuple[bool, bool, bool]
        ]:
    logger = get_logger()
    peptide_anno: Dict['Seq', Dict[str, 'AnnotatedPeptideLabel']] = {}

    def add_peptide_anno(x: Dict['Seq', List['AnnotatedPeptideLabel']]):
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
    dgraphs: TypeDGraphs = (None, {}, {})
    pgraphs: TypePGraphs = (None, {}, {})
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
            if isinstance(e, TimeoutError):
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
            variant_pool = pool.__class__()
            variant_pool.anno = pool.anno
            variant_pool.data = dict(pool.data)
            variant_pool[tx_id] = pool[tx_id].__class__(**pool[tx_id].__dict__)
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
            if isinstance(e, TimeoutError):
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
            if isinstance(e, TimeoutError):
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

    peptide_anno = {k: list(v.values()) for k, v in peptide_anno.items()}

    return peptide_anno, tx_id, dgraphs, pgraphs, success_flags
