""" Test module for moPepGen """
from typing import Dict, List, Tuple, Union
import copy
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen import svgraph, dna, seqvar, gtf

Type = Tuple[Union[svgraph.TranscriptVariantGraph, svgraph.CircularVariantGraph],
        Dict[int, svgraph.TVGNode]]
def create_dgraph2(data:dict, circular:bool=False, cds_start_nf:bool=False) -> Type:
    """ Create DNA transcript graph from node individuals.
    """
    node_list:Dict[int, svgraph.TVGNode] = {}
    graph = None
    for key, val in data.items():
        _seq = Seq(val[0]) if val[0] else None

        if not graph:
            if _seq:
                seq_location = dna.MatchedLocation(
                    query=FeatureLocation(start=0, end=len(_seq)),
                    ref=FeatureLocation(start=0, end=len(_seq))
                )
                seq = dna.DNASeqRecordWithCoordinates(_seq, [seq_location])
            else:
                seq = None
            graph = svgraph.TranscriptVariantGraph(seq, 'ENST0001', cds_start_nf)
            node_list[key] = graph.root
            continue

        in_nodes = [node_list[i] for i in val[1]]
        for node in in_nodes:
            if not node.variants:
                # upstream is the in_node that has no variants
                upstream = node
        node_start = upstream.seq.locations[-1].ref.end if upstream.seq else 0
        variants = []
        frameshifts = copy.copy(upstream.frameshifts)
        for var_data in val[2]:
            var_start = node_start + var_data[0]
            var_end = var_start + len(var_data[1])
            var_record = seqvar.VariantRecord(
                location=FeatureLocation(start=var_start, end=var_end),
                ref=var_data[1],
                alt=var_data[2],
                _type=var_data[3],
                _id=var_data[4]
            )
            var_location = FeatureLocation(
                start=var_data[0],
                end=var_data[0] + len(var_data[1])
            )
            variant = seqvar.VariantRecordWithCoordinate(
                variant=var_record,
                location=var_location
            )
            variants.append(variant)
            if variant.variant.is_frameshifting():
                frameshifts.add(variant.variant)

        left = 0
        seq_locations = []
        variant:seqvar.VariantRecordWithCoordinate
        for variant in variants:
            right = variant.location.start
            if right > left:
                ref_start = left + node_start
                ref_end = right  + node_start
                seq_location = dna.MatchedLocation(
                    query=FeatureLocation(start=left, end=right),
                    ref=FeatureLocation(start=ref_start, end=ref_end)
                )
                seq_locations.append(seq_location)
            left = variant.location.end

        if left < len(_seq):
            right = len(_seq)
            ref_start = left + node_start
            ref_end = right + node_start
            seq_location = dna.MatchedLocation(
                query=FeatureLocation(start=left, end=right),
                ref=FeatureLocation(start=ref_start, end=ref_end)
            )
            seq_locations.append(seq_location)

        # create a node
        seq = dna.DNASeqRecordWithCoordinates(_seq, seq_locations)
        node = svgraph.TVGNode(seq, variants, frameshifts)
        if len(val) == 4:
            node.branch = val[3]
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

        # update graph.seq if the current node has no variants
        if not node.variants:
            root_seq = graph.seq.seq + _seq if graph.seq else _seq
            root_seq_location = dna.MatchedLocation(
                query=FeatureLocation(start=0, end=len(root_seq)),
                ref=FeatureLocation(start=0, end=len(root_seq))
            )
            graph.seq = dna.DNASeqRecordWithCoordinates(root_seq,
                [root_seq_location])
    if circular:
        for i in data[1][1]:
            graph.add_edge(node_list[i], graph.root, 'reference')
            graph.__class__ = svgraph.CircularVariantGraph

    return graph, node_list


def create_variant(start:int, end:int, ref:str, alt:str, _type:str, _id:str,
        attrs:dict=None) -> seqvar.VariantRecord:
    """ Helper function to create a VariantRecord """
    location = FeatureLocation(start=start, end=end)
    return seqvar.VariantRecord(
        location=location, ref=ref, alt=alt,
        _type=_type, _id=_id, attrs=attrs
    )

def create_variants(data) -> List[seqvar.VariantRecord]:
    """ Helper function to create a list of VariantRecord """
    return [create_variant(*x) for x in data]

def create_dgraph1(seq, variants) -> svgraph.TranscriptVariantGraph:
    """ Create a dna transcript graph from sequence and variants """
    location = dna.MatchedLocation(
        query=FeatureLocation(start=0, end=len(seq)),
        ref=FeatureLocation(start=0, end=len(seq))
    )
    seq = dna.DNASeqRecordWithCoordinates(
        seq=Seq(seq),
        locations=[location],
        orf=FeatureLocation(start=0, end=len(seq))
    )
    graph = svgraph.TranscriptVariantGraph(seq, '')
    records = create_variants(variants)
    graph.create_variant_graph(records)
    return graph


def create_transcript_model(data:dict) -> gtf.TranscriptAnnotationModel:
    """ Create a transcript model from data.
    """
    chrom = data['chrom']
    strand = data['strand']
    entry = data['transcript']
    location = FeatureLocation(start=entry[0], end=entry[1])
    transcript = SeqFeature(chrom=chrom, location=location,
        attributes=entry[2], strand=strand)
    exons = []
    for entry in data['exon']:
        location = FeatureLocation(start=entry[0], end=entry[1])
        exons.append(SeqFeature(chrom=chrom, location=location,
            attributes=entry[2], strand=strand))
    cds = []
    if 'cds' in data:
        for entry in data['cds']:
            location = FeatureLocation(start=entry[0], end=entry[1])
            cds.append(SeqFeature(chrom=chrom, location=location,
                attributes=entry[2], strand=strand))
    model = gtf.TranscriptAnnotationModel(transcript, cds, exons)
    return model

def create_gene_model(data:dict) -> gtf.GeneAnnotationModel:
    """ Create a gene model from data from testing """
    chrom = data['chrom']
    strand = data['strand']
    entry = data['gene']
    location = FeatureLocation(start=entry[0], end=entry[1])
    return gtf.GeneAnnotationModel(chrom=chrom, location=location,
        attributes=entry[2], strand=strand, transcripts=data['transcripts'])


def create_genomic_annotation(data:dict) -> gtf.GenomicAnnotation:
    """ Create a GenomicAnnotation obejct from data for testing """
    anno = gtf.GenomicAnnotation()
    for entry in data['genes']:
        anno.genes[entry['gene_id']] = create_gene_model(entry)
    for entry in data['transcripts']:
        anno.transcripts[entry['transcript_id']] = create_transcript_model(entry)
    return anno

def create_dna_record_dict(data:dict) -> dna.DNASeqDict:
    """ Create a DNASeqDict as genome for testing """
    genome = dna.DNASeqDict()
    for key, val in data.items():
        genome[key] = dna.DNASeqRecord(Seq(val))
    return genome

def create_dna_seq_with_coordinates(seq, start=None, end=None):
    """ Create a dna.DNASeqRecordWithCoordinates instance """
    if not start:
        start = 0
    if not end:
        end = len(seq)
    location = dna.MatchedLocation(
        query=FeatureLocation(start=0, end=len(seq)),
        ref=FeatureLocation(start=start, end=end)
    )
    return dna.DNASeqRecordWithCoordinates(seq=Seq(seq), locations=[location])
