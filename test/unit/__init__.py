""" Test module for moPepGen """
from typing import Dict, List, Tuple, Union
import copy
from pathlib import Path
import pickle
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation, SeqFeature, MatchedLocation
from moPepGen import svgraph, dna, seqvar, gtf, aa
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from moPepGen.seqvar import VariantRecordPool


def load_references(base_dir:Path=None, index:bool=False
        ) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotation]:
    """ Load reference files """
    genome, anno = None, None
    if index:
        with open(base_dir/'genome.pickle', 'rb') as handle:
            genome = pickle.load(handle)
        with open(base_dir/'annotation.pickle', 'rb') as handle:
            anno = pickle.load(handle)
    else:
        genome = dna.DNASeqDict()
        genome.dump_fasta(base_dir/'genome.fasta')
        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(base_dir/'annotation.gtf')
    return genome, anno


Type = Tuple[Union[svgraph.ThreeFrameTVG, svgraph.ThreeFrameCVG],
        Dict[int, svgraph.TVGNode]]
def create_three_frame_tvg(nodes:Dict[int,list], seq:str) -> Type:
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
    seq = dna.DNASeqRecordWithCoordinates(raw_seq, [])
    graph = svgraph.ThreeFrameTVG(seq, _id='')
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
                location=FeatureLocation(start=var_start, end=var_end),
                ref=var_data[1],
                alt=var_data[2],
                _type=var_data[3],
                _id=var_data[4]
            )
            var_start = var_data[5] if len(var_data) >= 6 else var_data[0]
            var_end = var_data[6] if len(var_data) >= 7 else \
                    var_start + len(var_data[1])
            var_location = FeatureLocation(start=var_start, end=var_end)
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
                    query=FeatureLocation(start=left, end=right),
                    ref=FeatureLocation(start=ref_start, end=ref_end)
                )
                seq_locations.append(seq_location)
            left = variant.location.end

        if left < len(_seq):
            right = len(_seq)
            ref_start = left + node_start
            ref_end = right + node_start
            seq_location = MatchedLocation(
                query=FeatureLocation(start=left, end=right),
                ref=FeatureLocation(start=ref_start, end=ref_end)
            )
            seq_locations.append(seq_location)

        seq = dna.DNASeqRecordWithCoordinates(_seq, seq_locations)
        orf_idx = upstream.reading_frame_index
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
    location = FeatureLocation(start=entry[0], end=entry[1], seqname=chrom)
    transcript = GTFSeqFeature(chrom=chrom, location=location,
        attributes=entry[2], strand=strand, source='GENCODE')
    exons = []
    for entry in data['exon']:
        location = FeatureLocation(start=entry[0], end=entry[1])
        exons.append(GTFSeqFeature(chrom=chrom, location=location,
            attributes=entry[2], strand=strand, source='GENCODE'))
    cds = []
    if 'cds' in data:
        for entry in data['cds']:
            location = FeatureLocation(start=entry[0], end=entry[1])
            cds.append(GTFSeqFeature(chrom=chrom, location=location,
                attributes=entry[2], strand=strand, source='GENCODE',
                frame=0))
    model = gtf.TranscriptAnnotationModel(transcript, cds, exons)
    return model

def create_gene_model(data:dict) -> gtf.GeneAnnotationModel:
    """ Create a gene model from data from testing """
    chrom = data['chrom']
    strand = data['strand']
    entry = data['gene']
    location = FeatureLocation(start=entry[0], end=entry[1], seqname=chrom)
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

T = List[Tuple[Tuple[int,int], Tuple[int,int]]]
def create_dna_seq_with_coordinates(seq:str, locations=T, orf=Tuple[int,int]
        ) -> dna.DNASeqRecordWithCoordinates:
    """ Create a dna.DNASeqRecordWithCoordinates instance """
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
