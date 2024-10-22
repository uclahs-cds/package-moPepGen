""" Test module for moPepGen """
from __future__ import annotations
from typing import Dict, List, Tuple, Union, Set, TYPE_CHECKING
import copy
from pathlib import Path
import pickle
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen import params, svgraph, dna, seqvar, gtf, aa
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from moPepGen.seqvar import VariantRecordPool
from moPepGen.svgraph.SubgraphTree import SubgraphTree


GENOME_DATA = {
    'chr1':
    'ATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCTATGTACTGGTCCTTCTGCCT'
    'CCTCCCAATAAAGTCGAATTTTGGAACCGAATTCCCTTTTTTCGGGAAAAGCTACTAGGG'
}
ANNOTATION_ATTRS = [
    {
        'gene': {
            'gene_id'  : 'ENSG0001',
            'gene_name': 'SYMBO1'
        },
        'transcripts': [{
            'transcript_id': 'ENST0001.1',
            'gene_id'      : 'ENSG0001',
            'protein_id'   : 'ENSP0001',
            'gene_name'    : 'SYMBO1'
        }]
    }, {
        'gene': {
            'gene_id'  : 'ENSG0002',
            'gene_name': 'SYMBO2'
        },
        'transcripts': [{
            'transcript_id': 'ENST0002.1',
            'gene_id'      : 'ENSG0002',
            'protein_id'   : 'ENSP0002',
            'gene_name'    : 'SYMBO2'
        }]
    }
]
ANNOTATION_DATA = {
    'genes': [{
        'gene_id': ANNOTATION_ATTRS[0]['gene']['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (0, 40, ANNOTATION_ATTRS[0]['gene']),
        'transcripts': ['ENST0001.1']
    }, {
        'gene_id': ANNOTATION_ATTRS[1]['gene']['gene_id'],
        'chrom': 'chr1',
        'strand': 1,
        'gene': (60, 100, ANNOTATION_ATTRS[1]['gene']),
        'transcripts': ['ENST0002.1']
    }],
    'transcripts': [{
        'transcript_id': ANNOTATION_ATTRS[0]['transcripts'][0]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (5, 35, ANNOTATION_ATTRS[0]['transcripts'][0]),
        'exon': [
            (5,  12, ANNOTATION_ATTRS[0]['transcripts'][0]),
            (17, 23, ANNOTATION_ATTRS[0]['transcripts'][0]),
            (27, 35, ANNOTATION_ATTRS[0]['transcripts'][0])
        ]
    }, {
        'transcript_id': ANNOTATION_ATTRS[1]['transcripts'][0]['transcript_id'],
        'chrom': 'chr1',
        'strand': 1,
        'transcript': (65, 95, ANNOTATION_ATTRS[1]['transcripts'][0]),
        'exon': [
            (65, 72, ANNOTATION_ATTRS[1]['transcripts'][0]),
            (77, 83, ANNOTATION_ATTRS[1]['transcripts'][0]),
            (87, 95, ANNOTATION_ATTRS[1]['transcripts'][0])
        ]
    }]
}

def load_references(base_dir:Path=None, index:bool=False
        ) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotation]:
    """ Load reference files """
    genome, anno = None, None
    if index:
        with open(base_dir/'genome.pkl', 'rb') as handle:
            genome = pickle.load(handle)
        with open(base_dir/'annotation.pkl', 'rb') as handle:
            anno = pickle.load(handle)
    else:
        genome = dna.DNASeqDict()
        genome.dump_fasta(base_dir/'genome.fasta')
        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(base_dir/'annotation.gtf')
    return genome, anno


Type = Tuple[Union[svgraph.ThreeFrameTVG, svgraph.ThreeFrameCVG],
        Dict[int, svgraph.TVGNode]]
def create_three_frame_tvg(nodes:Dict[int,list], seq:str, graph_id:str='') -> Type:
    """ Create a ThreeFrameTVG

    Args:
        nodes (Dict[int,list]): Data for nodes must be a dict with keys are int
            and values are list of attributes of building the node. The value
            can hold up to x attributes.

            0 (str): the sequence of the node.
            1 (List[int|str]): List of int or str values to indicate the inbond
                node. If it is a string, it must follow the pattern of
                'RF[0-2]' to indicate that the node is outbond to the root of
                x reading frame.
            2 (List[tuple]): List of attributes to build variants.
            3 (int): Optional for the start coordinate of the sequence.
    """
    node_list:Dict[int,svgraph.TVGNode] = {}
    raw_seq = Seq(seq)
    seq = dna.DNASeqRecordWithCoordinates(raw_seq, locations=[])
    graph = svgraph.ThreeFrameTVG(seq, _id=graph_id)
    graph.cleavage_params = params.CleavageParams()
    for edge in copy.copy(graph.root.out_edges):
        graph.remove_edge(edge)
    graph.reading_frames[0] = svgraph.TVGNode(None, reading_frame_index=0,
        subgraph_id=graph.id)
    graph.reading_frames[1] = svgraph.TVGNode(None, reading_frame_index=1,
        subgraph_id=graph.id)
    graph.reading_frames[2] = svgraph.TVGNode(None, reading_frame_index=2,
        subgraph_id=graph.id)

    for key,val in nodes.items():
        _seq = Seq(val[0])

        in_nodes = []
        is_head_node = False
        for i in val[1]:
            if isinstance(i, str) and i.startswith('RF'):
                in_nodes.append(graph.reading_frames[int(i[-1])])
                is_head_node = True
            else:
                in_nodes.append(node_list[i])

        for node in in_nodes:
            if not node.variants:
                # upstream is the in_node that has no variants
                upstream = node
                break

        orf_idx = upstream.reading_frame_index

        if len(val) == 4:
            node_start = val[3]
        elif is_head_node:
            node_start = 0
        else:
            for in_node in in_nodes:
                if in_node.seq.locations:
                    loc = in_node.seq.locations[-1]
                    if loc.query.end != len(in_node.seq.seq):
                        continue
                    node_start = loc.ref.end
                    break
            else:
                raise ValueError('cannot infer the node start')

        variants:List[seqvar.VariantRecordWithCoordinate] = []
        for var_data in val[2]:
            var_start = node_start + var_data[0]
            var_end = var_start + len(var_data[1])
            var_record = seqvar.VariantRecord(
                location=FeatureLocation(
                    start=var_start, end=var_end,
                    reading_frame_index=orf_idx
                ),
                ref=var_data[1],
                alt=var_data[2],
                _type=var_data[3],
                _id=var_data[4]
            )
            var_start = var_data[5] if len(var_data) >= 6 else var_data[0]
            var_end = var_data[6] if len(var_data) >= 7 else \
                    var_start + len(var_data[1])
            var_location = FeatureLocation(
                start=var_start, end=var_end, reading_frame_index=orf_idx,
                seqname=graph_id
            )
            variant = seqvar.VariantRecordWithCoordinate(
                variant=var_record,
                location=var_location
            )
            variants.append(variant)

        left = 0
        seq_locations = []
        for variant in variants:
            right = variant.location.start
            if right > left:
                ref_start = left + node_start
                ref_end = right  + node_start
                seq_location = MatchedLocation(
                    query=FeatureLocation(
                        start=left, end=right, reading_frame_index=orf_idx
                    ),
                    ref=FeatureLocation(
                        start=ref_start, end=ref_end, reading_frame_index=orf_idx,
                        seqname=graph_id
                    )
                )
                seq_locations.append(seq_location)
            left = variant.location.end

        if left < len(_seq):
            right = len(_seq)
            ref_start = left + node_start
            ref_end = right + node_start
            seq_location = MatchedLocation(
                query=FeatureLocation(start=left, end=right,
                    reading_frame_index=orf_idx),
                ref=FeatureLocation(start=ref_start, end=ref_end,
                    reading_frame_index=orf_idx, seqname=graph_id)
            )
            seq_locations.append(seq_location)

        seq = dna.DNASeqRecordWithCoordinates(_seq, locations=seq_locations)
        node = svgraph.TVGNode(seq, variants,
            reading_frame_index=orf_idx, subgraph_id=graph.id)
        node_list[key] = node

        # add edges to in_nodes
        for in_node in in_nodes:
            if not in_node.variants and not node.variants:
                edge_type = 'reference'
            elif in_node.variants and not node.variants:
                edge_type = 'variant_end'
            else:
                edge_type = 'variant_start'
            graph.add_edge(in_node, node, edge_type)

    return graph, node_list

if TYPE_CHECKING:
    VariantData = Tuple[int, int, str, str, str, str, int, int, bool]
    PGraphData = Dict[int, Tuple[
        str,
        List[int],
        List[VariantData],
        List[Tuple[Tuple[int,int], Tuple[int,int]]],
        int
    ]]

def create_pgraph(data:PGraphData, _id:str, known_orf:List[int]=None,
        ) -> Tuple[svgraph.PeptideVariantGraph,Dict[int, svgraph.PVGNode]]:
    """ Create a peptide variant graph from data """
    root = svgraph.PVGNode(None, None, subgraph_id=_id)
    if not known_orf:
        known_orf = [None, None]
    cleavage_params = params.CleavageParams(
        enzyme='trypsin', exception = 'trypsin',
        miscleavage=0, min_mw=0, min_length=0
    )
    graph = svgraph.PeptideVariantGraph(root, _id, known_orf, cleavage_params)
    graph.subgraphs = SubgraphTree()
    graph.subgraphs.add_root(
        _id, feature_type='transcript', feature_id=_id, variant=None
    )
    node_list:Dict[int,svgraph.PVGNode] = {0: root}
    for key, val in data.items():
        locs = []
        for (query_start, query_end), (ref_start, ref_end) in val[3]:
            loc = MatchedLocation(
                query=FeatureLocation(
                    start=query_start, end=query_end, reading_frame_index=val[4]
                ),
                ref=FeatureLocation(start=ref_start, end=ref_end, seqname=_id)
            )
            locs.append(loc)

        seq = aa.AminoAcidSeqRecordWithCoordinates(
            Seq(val[0]),
            _id='ENST00001',
            transcript_id='ENST00001',
            locations=locs
        )
        seq.__class__ = aa.AminoAcidSeqRecordWithCoordinates
        variants:List[seqvar.VariantRecordWithCoordinate] = []

        for it in val[2]:
            if it is None:
                continue
            location_transcript = FeatureLocation(start=it[0], end=it[1])
            location_peptide = FeatureLocation(
                start=it[6], end=it[7], reading_frame_index=val[4]
            )
            var_record = seqvar.VariantRecord(
                location=location_transcript,
                ref=it[2],
                alt=it[3],
                _type=it[4],
                _id=it[5]
            )
            variant = seqvar.VariantRecordWithCoordinate(
                variant=var_record,
                location=location_peptide
            )

            variants.append(variant)

        node = svgraph.PVGNode(
            seq, val[4], variants=variants, subgraph_id=_id,
            left_cleavage_pattern_end=1,
            right_cleavage_pattern_start=len(seq) - 1
        )

        node_list[key] = node
        for i in val[1]:
            node_list[i].add_out_edge(node)
        if 0 in val[1]:
            graph.reading_frames[val[4]] = node

    return graph, node_list


def create_variant(start:int, end:int, ref:str, alt:str, _type:str, _id:str,
        attrs:dict=None, seqname:str=None) -> seqvar.VariantRecord:
    """ Helper function to create a VariantRecord """
    location = FeatureLocation(start=start, end=end, seqname=seqname)
    return seqvar.VariantRecord(
        location=location, ref=ref, alt=alt,
        _type=_type, _id=_id, attrs=attrs
    )

def create_variants(data) -> List[seqvar.VariantRecord]:
    """ Helper function to create a list of VariantRecord """
    return [create_variant(*x) for x in data]

def create_variant_with_coordinate(query_start:int, query_end:int, start:int, end:int,
        ref:str, alt:str, _type:str, _id:str, attrs:dict=None, seqname:str=None):
    """ Create VariantWithCoordinate """
    return seqvar.VariantRecordWithCoordinate(
        variant = create_variant(
            start=start, end=end, ref=ref, alt=alt, _type=_type, _id=_id,
            attrs=attrs, seqname=seqname
        ),
        location=FeatureLocation(
            start=query_start, end=query_end
        )
    )

def create_dgraph1(seq, variants, has_known_orf:bool=True,
        variant_pool:VariantRecordPool=None,
        genome:dna.DNASeqDict=None, anno:gtf.GenomicAnnotation=None
        ) -> svgraph.ThreeFrameTVG:
    """ Create a dna transcript graph from sequence and variants """
    location = MatchedLocation(
        query=FeatureLocation(start=0, end=len(seq)),
        ref=FeatureLocation(start=0, end=len(seq))
    )
    seq = dna.DNASeqRecordWithCoordinates(
        seq=Seq(seq),
        locations=[location],
        orf=FeatureLocation(start=0, end=len(seq))
    )
    graph = svgraph.ThreeFrameTVG(seq, '', has_known_orf=has_known_orf)
    records = create_variants(variants)
    graph.create_variant_graph(records, variant_pool, genome, anno)
    return graph

def create_transcript_model(data:dict) -> gtf.TranscriptAnnotationModel:
    """ Create a transcript model from data.
    """
    chrom = data['chrom']
    strand = data['strand']
    entry = data['transcript']
    location = FeatureLocation(
        start=entry[0], end=entry[1], seqname=chrom, strand=strand
    )
    transcript = GTFSeqFeature(chrom=chrom, location=location,
        attributes=entry[2], source='GENCODE')
    exons = []
    for entry in data['exon']:
        location = FeatureLocation(start=entry[0], end=entry[1], strand=strand)
        feature = GTFSeqFeature(
            chrom=chrom, location=location, attributes=entry[2], source='GENCODE'
        )
        exons.append(feature)
    cds = []
    if 'cds' in data:
        for entry in data['cds']:
            location = FeatureLocation(start=entry[0], end=entry[1], strand=strand)
            feature = GTFSeqFeature(
                chrom=chrom, location=location,
                attributes=entry[2], source='GENCODE', frame=0
            )
            cds.append(feature)
    model = gtf.TranscriptAnnotationModel(transcript, cds, exons)
    return model

def create_gene_model(data:dict) -> gtf.GeneAnnotationModel:
    """ Create a gene model from data from testing """
    chrom = data['chrom']
    strand = data['strand']
    entry = data['gene']
    location = FeatureLocation(start=entry[0], end=entry[1], seqname=chrom, strand=strand)
    return gtf.GeneAnnotationModel(
        chrom=chrom, location=location, attributes=entry[2],
        transcripts=data['transcripts']
    )

def create_genomic_annotation(data:dict) -> gtf.GenomicAnnotation:
    """ Create a GenomicAnnotation obejct from data for testing """
    anno = gtf.GenomicAnnotation()
    for entry in data['genes']:
        anno.genes[entry['gene_id']] = create_gene_model(entry)
    for entry in data['transcripts']:
        anno.transcripts[entry['transcript_id']] = create_transcript_model(entry)
    for tx_model in anno.transcripts.values():
        tx_model.is_protein_coding = len(tx_model.cds) > 0
    return anno

def create_dna_record_dict(data:dict) -> dna.DNASeqDict:
    """ Create a DNASeqDict as genome for testing """
    genome = dna.DNASeqDict()
    for key, val in data.items():
        genome[key] = dna.DNASeqRecord(Seq(val))
    return genome

T = List[Tuple[Tuple[int,int], Tuple[int,int]]]
def create_dna_seq_with_coordinates(seq:str, locations:T=None, orf=Tuple[int,int]
        ) -> dna.DNASeqRecordWithCoordinates:
    """ Create a dna.DNASeqRecordWithCoordinates instance """
    if not locations:
        locations = [((0,len(seq)), (0, len(seq)))]

    locs = []
    for (a,b),(c,d) in locations:
        loc = MatchedLocation(
            query=FeatureLocation(start=a, end=b),
            ref=FeatureLocation(start=c, end=d)
        )
        locs.append(loc)
    return dna.DNASeqRecordWithCoordinates(
        seq=Seq(seq),
        locations=locs,
        orf=orf
    )

def create_aa_seq_with_coordinates(seq:str, locations=T, orf=Tuple[int,int]
        ) -> aa.AminoAcidSeqRecordWithCoordinates:
    """ creat amino acid sequence with coordinates """
    locs = []
    for (a,b),(c,d) in locations:
        loc = MatchedLocation(
            query=FeatureLocation(start=a, end=b),
            ref=FeatureLocation(start=c, end=d)
        )
        locs.append(loc)
    return aa.AminoAcidSeqRecordWithCoordinates(
        seq=Seq(seq),
        locations=locs,
        orf=orf
    )

def create_aa_record(seq:str, description:str):
    """ Create a AminoAcidSeqRecord """
    seq = Seq(seq)
    return aa.AminoAcidSeqRecord(seq, _id=description, name=description,
        description=description)

def get_tx2gene_and_coding_tx(anno:gtf.GenomicAnnotation) -> Tuple[Dict[str,str], Set[str]]:
    """ extract tx to gene map and all coding transcripts """
    tx2gene = {}
    coding_tx = set()
    for tx_id in anno.transcripts:
        tx_model = anno.transcripts[tx_id]
        tx2gene[tx_id] = tx_model.transcript.gene_id
        if tx_model.is_protein_coding:
            coding_tx.add(tx_id)
    return tx2gene, coding_tx

def  create_pvg_node(seq:str, reading_frame_index=0, subgraph_id='TEST001', variants=None,
        ) -> svgraph.PVGNode:
    """ Create a PVGNode """
    if variants is None:
        variants = []
    return svgraph.PVGNode(
        seq = aa.AminoAcidSeqRecordWithCoordinates(
            seq=Seq(seq),
            locations=[MatchedLocation(
                query=FeatureLocation(start=0, end=len(seq)),
                ref=FeatureLocation(start=0, end=len(seq))
            )],
        ),
        reading_frame_index=reading_frame_index,
        subgraph_id=subgraph_id,
        variants = [create_variant_with_coordinate(*x) for x in variants]
    )
