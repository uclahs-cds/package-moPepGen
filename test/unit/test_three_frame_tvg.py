""" Test module for ThreeFrameTVG """
import copy
from typing import Dict, List, Tuple, Union
import unittest
from test.unit import create_variant, create_variants, \
    create_genomic_annotation, create_dna_record_dict
from test.unit.test_vep_parser import GENOME_DATA, ANNOTATION_DATA
from Bio.Seq import Seq
from moPepGen.SeqFeature import MatchedLocation, FeatureLocation
from moPepGen.svgraph.ThreeFrameTVG import ThreeFrameTVG
from moPepGen.svgraph import TVGNode
from moPepGen import dna, seqvar


def create_three_frame_tvg(data:Dict[Union[int,str],List]
        ) -> Tuple[ThreeFrameTVG,Dict[int,TVGNode]]:
    node_list:Dict[int,TVGNode] = {}
    raw_seq = Seq(data['seq'])
    seq = dna.DNASeqRecordWithCoordinates(raw_seq, [])
    graph = ThreeFrameTVG(seq, _id='')
    for edge in copy.copy(graph.root.out_edges):
        graph.remove_edge(edge)
    node_list[0] = graph.root
    node_data:Dict[int,list] = data['nodes']
    for key,val in node_data.items():
        _seq = Seq(val[0])

        in_nodes = [node_list[i] for i in val[1]]
        for node in in_nodes:
            if not node.variants:
                # upstream is the in_node that has no variants
                upstream = node
                break

        if 0 in val[1]:
            node_start = data['node_start'][key]
        else:
            upstream.seq.locations[-1].ref.end

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
            var_start = var_data[5] if len(var_data) >= 6 else var_data[0]
            var_end = var_data[6] if len(var_data) >= 7 else \
                    var_start + len(var_data[1])
            var_location = FeatureLocation(start=var_start, end=var_end)
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
        node = TVGNode(seq, variants, frameshifts)
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

    for i, key in enumerate(data['reading_frames']):
        graph.reading_frames[i] = node_list[key]

    return graph, node_list


class TestCaseThreeFrameTVG(unittest.TestCase):
    """ Test case for ThreeFrameTVG """
    def test_apply_variant_inframe(self):
        """ apply_variant """
        data = {
            'seq': 'AAAA',
            'reading_frames': [1,2,3],
            'node_start': {1: 0, 2:1, 3:2},
            'nodes': {
                1: ['AAAAAT', [0], []],
                2: ['AAAAT', [0], []],
                3: ['AAAT', [0], []]
            }
        }

        graph, nodes = create_three_frame_tvg(data)
        variant = create_variant(3, 4, 'A', 'T', 'SNV', '')
        source, target = graph.apply_variant(nodes[1], nodes[1], variant)
        self.assertIs(source, target)
        out_nodes = {str(e.out_node.seq.seq) for e in source.out_edges}
        expected = {'A', 'T'}
        self.assertEqual(out_nodes, expected)

    def test_apply_variant_deletion(self):
        """ apply_variant """
        data = {
            'seq': 'AAAA',
            'reading_frames': [1,2,3],
            'node_start': {1: 0, 2:1, 3:2},
            'nodes': {
                1: ['AAAAAT', [0], []],
                2: ['AAAAT', [0], []],
                3: ['AAAT', [0], []]
            }
        }

        graph, nodes = create_three_frame_tvg(data)
        variant = create_variant(3, 5, 'AA', 'A', 'INDEL', '')
        source, target = graph.apply_variant(nodes[1], nodes[2], variant)
        self.assertIn(source, graph.reading_frames)
        self.assertIn(target, graph.reading_frames)
        self.assertEqual(str(source.seq.seq), 'AAA')
        self.assertEqual(str(target.seq.seq), 'AAAA')
        out_nodes = {str(e.out_node.seq.seq) for e in source.out_edges}
        self.assertEqual(out_nodes, {'A', 'AAT'})

    def test_create_variant_graph_in_frame(self):
        """ apply_variant """
        data = {
            'seq': 'AAATAAATAAAT',
            'reading_frames': [1,2,3],
            'node_start': {1: 0, 2:1, 3:2},
            'nodes': {
                1: ['AAATAAATAAAT', [0], []],
                2: ['AATAAATAAAT', [0], []],
                3: ['ATAAATAAAT', [0], []]
            }
        }

        var_data = [
            (3, 4, 'T', 'A', 'SNV', ''),
            (6, 7, 'A', 'T', 'SNV', '')
        ]

        graph, _ = create_three_frame_tvg(data)
        variants = create_variants(var_data)
        graph.create_variant_graph(variants)
        received = {str(n.seq.seq) for n in graph.reading_frames}
        expected = {'AAA', 'AA', 'A'}
        self.assertEqual(received, expected)
        received = {str(e.out_node.seq.seq) for e in graph.reading_frames[0].out_edges}
        expected = {'A', 'T'}
        self.assertEqual(received, expected)

    def test_create_variant_graph_frameshifting(self):
        """ apply_variant """
        data = {
            'seq': 'AAATAAATAAAT',
            'reading_frames': [1,2,3],
            'node_start': {1: 0, 2:1, 3:2},
            'nodes': {
                1: ['AAATAAATAAAT', [0], []],
                2: ['AATAAATAAAT', [0], []],
                3: ['ATAAATAAAT', [0], []]
            }
        }

        var_data = [
            (3, 5, 'TA', 'T', 'INDEL', ''),
            (7, 8, 'T', 'TC', 'INDEL', '')
        ]

        graph, _ = create_three_frame_tvg(data)
        variants = create_variants(var_data)
        graph.create_variant_graph(variants)
        received = {str(n.seq.seq) for n in graph.reading_frames}
        expected = {'AAA', 'AA', 'A'}
        self.assertEqual(received, expected)

    def test_apply_insertion(self):
        """ apply_insertion """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        genome = create_dna_record_dict(GENOME_DATA)
        data = {
            'seq': 'AAATAAATAAAT',
            'reading_frames': [1,2,3],
            'node_start': {1: 0, 2:0, 3:0},
            'nodes': {
                1: ['AAATAAATAAAT', [0], []],
                2: ['AAATAAATAAAT', [0], []],
                3: ['AAATAAATAAAT', [0], []]
            }
        }
        gene_id = 'ENSG0001'
        ins_attrs = {
            'DONOR_GENE_ID': 'ENSG0001',
            'DONOR_START': 17,
            'DONOR_END': 23
        }
        ins_var = create_variant(3,4,'T','<INS>','Insertion','',ins_attrs)

        var_data = [(10, 11, 'A', 'T', 'SNV', '')]
        vars = create_variants(var_data)
        var_pool = seqvar.VariantRecordPool(genetic={gene_id: vars})

        graph, nodes = create_three_frame_tvg(data)
        cursors = [nodes[1], nodes[2], nodes[3]]
        graph.apply_insertion(cursors, ins_var, var_pool, genome, anno)

        received = {str(n.seq.seq) for n in graph.reading_frames}
        self.assertEqual(received, {'AAA'})

