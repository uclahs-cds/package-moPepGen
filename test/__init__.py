""""""
from typing import Dict, List, Tuple
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import svgraph, dna, seqvar


def create_dgraph(data:dict
        ) -> Tuple[svgraph.TranscriptVariantGraph, Dict[int, svgraph.DNANode]]:
    """
    example data: {
        1: ('ATGTGGC', [0], []),
        2: ('A', [1], [(0, 'C', 'A', 'SNV', '')])
    }
    """
    node_list:Dict[int, svgraph.DNANode] = {}
    graph = None
    for key, val in data.items():
        _seq = Seq(val[0])
        
        if not graph:
            seq_location = dna.MatchedLocation(
                query=FeatureLocation(start=0, end=len(_seq)),
                ref=FeatureLocation(start=0, end=len(_seq))
            )
            seq = dna.DNASeqRecordWithCoordinates(_seq, [seq_location])
            graph = svgraph.TranscriptVariantGraph(seq, 'ENST0001')
            node_list[key] = graph.root
            continue

        in_nodes = [node_list[i] for i in val[1]]
        for node in in_nodes:
            if not node.variants:
                # upstream is the in_node that has no variants
                upstream = node
        variants = []
        frameshifts = upstream.frameshifts
        for var_data in val[2]:
            node_start = upstream.seq.locations[-1].ref.end
            var_start = node_start + var_data[0]
            var_end = var_start + len(var_data[1])
            var_record = seqvar.VariantRecord(
                location=FeatureLocation(start=var_start, end=var_end),
                ref=var_data[1],
                alt=var_data[2],
                type=var_data[3],
                id=var_data[4]
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
                ref_start = left + upstream.seq.locations[-1].ref.end
                ref_end = right  + upstream.seq.locations[-1].ref.end
                seq_location = dna.MatchedLocation(
                    query=FeatureLocation(start=left, end=right),
                    ref=FeatureLocation(start=ref_start, end=ref_end)
                )
                seq_locations.append(seq_location)
            left = variant.location.end
        
        if left < len(_seq):
            right = len(_seq)
            ref_start = left + upstream.seq.locations[-1].ref.end
            ref_end = right + upstream.seq.locations[-1].ref.end
            seq_location = dna.MatchedLocation(
                query=FeatureLocation(start=left, end=right),
                ref=FeatureLocation(start=ref_start, end=ref_end)
            )
            seq_locations.append(seq_location)
        
        # create a node
        seq = dna.DNASeqRecordWithCoordinates(_seq, seq_locations)
        node = svgraph.DNANode(seq, variants, frameshifts)
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
            root_seq = graph.seq.seq + _seq
            root_seq_location = dna.MatchedLocation(
                query=FeatureLocation(start=0, end=len(root_seq)),
                ref=FeatureLocation(start=0, end=len(root_seq))
            )
            graph.seq = dna.DNASeqRecordWithCoordinates(root_seq,
                [root_seq_location])

    return graph, node_list