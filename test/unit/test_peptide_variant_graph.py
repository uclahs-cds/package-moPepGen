""" Module to test PeptideVariantGraph """
from typing import Tuple, Dict, List
import unittest
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation
from moPepGen import aa, params, svgraph, seqvar
from moPepGen.svgraph.PVGOrf import PVGOrf
from moPepGen.svgraph.PeptideVariantGraph import PVGTraversal, PVGCursor
from moPepGen.svgraph.VariantPeptideDict import VariantPeptideDict


def create_pgraph(data:dict, _id:str, known_orf:List[int]=None,
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

        node = svgraph.PVGNode(seq, val[4], variants=variants, subgraph_id=_id)

        node_list[key] = node
        for i in val[1]:
            node_list[i].add_out_edge(node)
        if 0 in val[1]:
            graph.reading_frames[val[4]] = node

    return graph, node_list


class TestPeptideVariantGraph(unittest.TestCase):
    """ Test case for peptide variant graph """
    def test_find_reference_next(self):
        """ Test that the reference next and prev can be found """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        data = {
            1: ('SSSR', [0], [None], [((0,4), (0,4))], 0),
            2: ('SSSV', [1], [variant_1], [((0,3), (4,7))], 0),
            3: ('SSSG', [1], [None], [((0,4), (4,8))], 0),
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        self.assertIs(nodes[1].find_reference_next(), nodes[3])
        graph.add_stop(nodes[2])
        self.assertIs(nodes[2].find_reference_next(), graph.stop)

    def test_find_reference_next_chimeric(self):
        """ When nodes are chimeric, such as combined because of frameshifting
        mutations. """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        variant_2 = (5, 8, 'ACC', 'A', 'INDEL', '0:ACC-A', 0, 1, True)
        data = {
            1: ('SSSR', [0], [variant_1], [((1,4),(1,4))], 0),
            2: ('SSSV', [1],[variant_1, variant_2], [((1,4),(1,4))], 0),
            3: ('SSSG', [1], [variant_1], [((1,4),(1,4))], 0),
        }
        _, nodes = create_pgraph(data, 'ENST0001')
        self.assertIs(nodes[1].find_reference_next(), nodes[3])

    def test_find_reference_prev_chimeric(self):
        """ When nodes are chimeric, such as combined because of frameshifting
        mutations. """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        variant_2 = (5, 8, 'ACC', 'A', 'INDEL', '0:ACC-A', 0, 1, True)
        data = {
            1: ('SSST', [0], [], [((0,4),(0,4))], 0),
            2: ('SSSR', [1], [variant_1], [((0,4),(0,4))], 0),
            3: ('SSSV', [1],[variant_1, variant_2], [((0,4),(0,4))], 0),
            4: ('SSSG', [2,3], [variant_1], [((0,4),(0,4))], 0),
        }
        _, nodes = create_pgraph(data, 'ENST0001')
        self.assertIs(nodes[4].find_reference_prev(), nodes[2])

    def test_find_orf_stop_gain(self):
        """ test finding orf start position of a stop gain """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('MSSSR', [0], [None], locations, 0),
            2: ('GSS', [1],[variant_1], [], 0),
            3: ('GSSSK', [1], [None], locations, 0)
        }
        graph,nodes = create_pgraph(data, '')
        graph.add_stop(nodes[2])
        self.assertEqual(nodes[2].get_orf_start(), -1)

    def test_find_orf_varinats_outnode(self):
        """ test finding orf start two variant nodes in a row """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        variant_2 = (0, 1, 'T', 'TCT', 'INDEL', '0:T-TCT', 1, 2, True)
        locations1 = [((0,4),(0,4))]
        locations2 = [((0,1), (5,6)),((2,4),(7,9))]
        data = {
            1: ('MSSSR', [0], [None], locations1, 0),
            2: ('GSSR', [1],[variant_1], [], 0),
            3: ('GSSSK', [1], [None], locations1, 0),
            4: ('SSSG', [2], [variant_2], locations2, 0)
        }
        _,nodes = create_pgraph(data, '')
        self.assertEqual(nodes[2].get_orf_start(), 15)

    def test_node_is_subgraph_end_case1(self):
        """ Test case that the node should not be identified as subgraph end
        if it's next node is the stop node """
        data = {
            1:  ('VLR',          [0],   [None], [((0,3),(0,3))],    1),
            2:  ('NG',           [1],   [None], [((0,2),(3,5))],    1),
            3:  ('NGALT',        [2],   [None], [((0,5),(3,8))],    1)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.add_stop(nodes[3])
        nodes[2].subgraph_id = 'ENSG0002'
        nodes[3].subgraph_id = 'ENSG0002'
        self.assertFalse(graph.node_is_subgraph_end(nodes[3]))

    def test_expand_backward_no_inbound(self):
        r""" > Test expand backward when the node has no inbound nodes.

                  V-P
                 /
            MNAAC-VC-PRN
                 \  /
                  DC
        """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('MNAAC', [0], [None], locations, 0),
            2: ('V', [1],[variant_1], locations, 0),
            3: ('P', [2], [None], locations, 0),
            4: ('VC', [1], [None], locations, 0),
            5: ('DC', [1], [variant_1], locations, 0),
            6: ('PRN', [4, 5], [None], locations, 0)
        }

        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        graph.expand_backward(nodes[1])

        seqs = {str(node.seq.seq) for node in graph.root.out_nodes}
        self.assertEqual(seqs, {'MNAACV', 'MNAACVC', 'MNAACDC'})
        self.assertEqual(len(nodes[1].in_nodes), 0)
        self.assertEqual(len(nodes[1].out_nodes), 0)
        self.assertNotIn(nodes[1], graph.root.out_nodes)

    def test_expand_backward_one_inbound(self):
        r""" > Test expand backward when the node has no inbound nodes.

                       V-P
                      /
            MNAAC-GCVV-VC-PRN
                      \  /
                       DC
        """
        variant_1 = (0,3,'TCT','T','INDEL', '0:TCT-T', 0, 1, True)
        variant_2 = (0, 1, 'T', 'A', 'SNV', '0:T-A', 0, 1, False)
        locations = [((0,4),(0,4))]
        data = {
            1: ('MNAAC', [0], [None], locations, 0),
            2: ('GCVV', [1], [None], locations, 0),
            3: ('V', [2],[variant_1], locations, 0),
            4: ('P', [3], [None], locations, 0),
            5: ('VC', [2], [None], locations, 0),
            6: ('DC', [2], [variant_2], locations, 0),
            7: ('PRN', [5, 6], [None], locations, 0)
        }

        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        branches,_ = graph.expand_backward(nodes[2])

        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        self.assertEqual(seqs, {'GCVVV', 'GCVVVC', 'GCVVDC'})
        self.assertEqual(len(nodes[2].in_nodes), 0)
        self.assertEqual(len(nodes[2].out_nodes), 0)
        self.assertNotIn(nodes[2], nodes[1].out_nodes)

        self.assertEqual(len(branches), 2)
        seqs = {str(branch.seq.seq) for branch in branches}
        self.assertEqual(seqs, {'GCVVV', 'PRN'})

    def test_expand_forward(self):
        r""" > expand forward
                  GCVVV-P
                 /
            MNAAR-GCVVVC-PDNK-DCAGP
                 \      /
                  GCVVDC
        """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '', 0, 1, True)
        variant_2 = (0, 1, 'T', 'A', 'SNV', '', 0, 1, False)
        locations = locations = [((0,4),(0,4))]
        data = {
            1: ('MNAAR', [0], [None], locations, 0),
            2: ('GCVVV', [1], [variant_1], locations, 0),
            3: ('P', [3], [None], locations, 0),
            4: ('GCVVVC', [1], [None], locations, 0),
            5: ('GCVVDC', [1], [variant_2], locations, 0),
            6: ('PDNK', [4, 5], [None], locations, 0),
            7: ('DCAGP', [6], [None], locations, 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        graph.expand_forward(nodes[6])
        seqs = {str(node.seq.seq) for node in nodes[7].in_nodes}
        self.assertEqual(seqs, {'GCVVVCPDNK', 'GCVVDCPDNK'})
        self.assertEqual(len(nodes[6].in_nodes), 0)
        self.assertEqual(len(nodes[6].out_nodes), 0)
        self.assertNotIn(nodes[6], nodes[4].out_nodes)
        self.assertNotIn(nodes[6], nodes[5].out_nodes)
        self.assertNotIn(nodes[6], nodes[7].in_nodes)

    def test_expand_forward_with_new_cleave_site(self):
        r""" > expand forward when there is a new cleave site
                  GCVVV-P
                 /
            MNAAR-GCVVVC-PDNK-DCAGP
                 \      /
                  GCVVRC
        """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '', 0, 1, True)
        variant_2 = (0, 1, 'T', 'A', 'SNV', '', 0, 1, False)
        locations = locations = [((0,4),(0,4))]
        data = {
            1: ('MNAAR', [0], [None], locations, 0),
            2: ('GCVVV', [1], [variant_1], locations, 0),
            3: ('P', [3], [None], locations, 0),
            4: ('GCVVVC', [1], [None], locations, 0),
            5: ('GCVVRC', [1], [variant_2], locations, 0),
            6: ('PDNK', [4, 5], [None], locations, 0),
            7: ('DCAGP', [6], [None], locations, 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        graph.expand_forward(nodes[6])
        seqs = {str(node.seq.seq) for node in nodes[7].in_nodes}
        self.assertEqual(seqs, {'GCVVVCPDNK', 'CPDNK'})
        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        self.assertIn('GCVVR', seqs)
        self.assertEqual(len(nodes[6].in_nodes), 0)
        self.assertEqual(len(nodes[6].out_nodes), 0)
        self.assertNotIn(nodes[6], nodes[4].out_nodes)
        self.assertNotIn(nodes[6], nodes[5].out_nodes)
        self.assertNotIn(nodes[6], nodes[7].in_nodes)

    def test_merge_join(self):
        r"""
                                       NCWHSTQV
                                      /        \
            NCWV     V               | NCWVSTQV |
           /    \   / \              |/        \|
        AER-NCWH-STQ-Q-PK    ->    AER-NCWHSTQQ-PK
                                      \        /
                                       NCWVSTQQ
        """
        variant_1 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        variant_2 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        locations = [((0,4),(0,4))]
        data = {
            1: ('AER', [0], [None], locations, 0),
            2: ('NCWH', [1], [None], locations, 0),
            3: ('NCWV', [1], [variant_1], locations, 0),
            4: ('STQ', [2,3], [], locations, 0),
            5: ('Q', [4], [], locations, 0),
            6: ('V', [4], [variant_2], locations, 0),
            7: ('PK', [5,6], [], locations, 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        branches,_ = graph.merge_join(nodes[4])
        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        seqs_expect = {'NCWHSTQV', 'NCWVSTQV', 'NCWHSTQQ', 'NCWVSTQQ'}
        self.assertEqual(seqs, seqs_expect)
        self.assertEqual(len(nodes[4].in_nodes), 0)
        self.assertEqual(len(nodes[4].out_nodes), 0)
        seqs = {str(node.seq.seq) for node in branches}
        expected = {'PK'}
        self.assertEqual(seqs, expected)

    def test_merge_join_with_cleave(self):
        r"""
                                       NCWHSTQV-
                                      /         \
            NCWR     V               | NCWR-STQV |
           /    \   / \              |/         \|
        AER-NCWH-STQ-Q-PK    ->    AER-NCWHSTQQ--PK
                                      \         /
                                       NCWR-STQQ
        """
        variant_1 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        variant_2 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        locations = [((0,4),(0,4))]
        data = {
            1: ('AER', [0], [None], locations, 0),
            2: ('NCWH', [1], [None], locations, 0),
            3: ('NCWR', [1], [variant_1], locations, 0),
            4: ('STQ', [2,3], [], locations, 0),
            5: ('Q', [4], [], locations, 0),
            6: ('V', [4], [variant_2], locations, 0),
            7: ('PK', [5,6], [], locations, 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        branches,_ = graph.merge_join(nodes[4])
        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        seqs_expect = {'NCWHSTQV', 'NCWR', 'NCWHSTQQ', 'NCWR'}
        self.assertEqual(seqs, seqs_expect)
        self.assertEqual(len(nodes[4].in_nodes), 0)
        self.assertEqual(len(nodes[4].out_nodes), 0)
        seqs = {str(node.seq.seq) for node in branches}
        expected = {'PK'}
        self.assertEqual(seqs, expected)
        seqs = {str(node.seq.seq) for node in nodes[7].in_nodes}
        seqs_expect = {'NCWHSTQV', 'STQV', 'NCWHSTQQ', 'STQQ'}
        self.assertEqual(seqs, seqs_expect)

    def test_cross_join(self):
        r"""
            NCWV           V                 NCWVSTEEK-LPAQV
           /    \         / \               /         X     \
        AER-NCWH-STEEKLPAQ-Q-PK    ->    AER-NCWHSTEEK-LPAQQ-PK
        """
        variant_1 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        variant_2 = (0, 1, 'A', 'T', 'SNV', '', 0, 1, False)
        locations = [((0,4),(0,4))]
        data = {
            1: ('AER', [0], [None], locations, 0),
            2: ('NCWH', [1], [None], locations, 0),
            3: ('NCWV', [1], [variant_1], locations, 0),
            4: ('STEEKLPAQ', [2,3], [None], locations, 0),
            5: ('Q', [4], [None], locations, 0),
            6: ('V', [4], [variant_2], locations, 0),
            7: ('PK', [5,6], [None], locations, 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        branches,_ = graph.cross_join(nodes[4], 5)

        seqs = {str(node.seq.seq) for node in branches}
        expected = {'PK'}
        self.assertEqual(seqs, expected)

        seqs = {str(node.seq.seq) for node in nodes[1].out_nodes}
        seqs_expected = {'NCWVSTEEK', 'NCWHSTEEK'}
        self.assertEqual(seqs, seqs_expected)

        seqs = {str(node.seq.seq) for node in nodes[7].in_nodes}
        seqs_expected = {'LPAQV', 'LPAQQ'}
        self.assertEqual(seqs, seqs_expected)

        self.assertTrue(len(nodes[2].out_nodes) == 0)
        self.assertTrue(len(nodes[2].in_nodes) == 0)
        self.assertTrue(len(nodes[3].out_nodes) == 0)
        self.assertTrue(len(nodes[3].in_nodes) == 0)

    def test_fit_into_cleavage_multiple_upstream_case1(self):
        """ Test case when a node has multiple splicing site and multiple
        upstream

                            WHRTWHR-CYKLTLT
                                   /
                  AGGL-DKQLRAARQYQQ
                 /    /
            NGALR-ALGL
           /
        VLR-NG

        """
        data = {
            1:  ('VLR',          [0],   [None], [((0,3),(0,3))],    1),
            2:  ('NG',           [1],   [None], [((0,2),(3,5))],    1),
            3:  ('NGALT',        [1],   [None], [((0,5),(3,8))],    1),
            4:  ('ALGL',         [3],   [None], [((0,4),(8,12))],   1),
            5:  ('AGGL',         [3],   [None], [((0,4),(8,12))],   1),
            6:  ('DKQLRAARQYQQ', [4,5], [None], [((0,12),(12,24))], 1),
            10: ('WHRWHR',       [0],   [None], [], 0),
            11: ('CYKLTLT',      [6,10],[None], [], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.rule = 'trypsin'
        downstreams,_ = graph.fit_into_cleavage_multiple_upstream(nodes[6])
        self.assertTrue(any(x.seq.seq == 'QYQQ' for x in nodes[11].in_nodes))
        self.assertEqual(len(downstreams), 0)

    def test_call_variant_peptides_ambiguous_amino_acid(self):
        """ test call peptides with ambiguous amino acid """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations, 0),
            2: ('SSSSR', [1],[variant_1], locations, 0),
            3: ('SSSIR', [1], [None], locations, 0),
            4: ('SSXPK', [2,3], [None], locations, 0)
        }
        graph, _ = create_pgraph(data, 'ENST0001', known_orf=[0,45])
        graph.cleavage_params.miscleavage = 1
        peptides = graph.call_variant_peptides()
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_miscleavages(self):
        """ test micleavages is handled correctly """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations, 0),
            2: ('SSSSR', [1],[variant_1], locations, 0),
            3: ('SSSIR', [1], [None], locations, 0),
            4: ('SSSPK', [2,3], [None], locations, 0),
            5: ('SSSVK', [4], [None], locations, 0)
        }
        graph, _ = create_pgraph(data, 'ENST0001', known_orf=[0,60])
        graph.cleavage_params.miscleavage = 1
        peptides = graph.call_variant_peptides()
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR', 'SSSSRSSSPK'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_miscleavages_not_in_cds(self):
        """ test micleavages is handled correctly """
        variant_1 = (10, 11, 'C', 'T', 'INDEL', '0:TCT-T', 4, 5, True)
        data = {
            1: ('SSSSSSSSIR', [0], [None], [((0,5),(0,5))], 0),
            2: ('SSSSR', [0], [variant_1], [((0,4),(0,4))], 0),
            3: ('SSSIR', [2], [None], [((0,5),(5,10))], 0),
            4: ('SSSSS', [1,3], [None], [((0,5),(10,15))], 0),
        }
        graph, _ = create_pgraph(data, 'ENST0001', known_orf=[0,60])
        graph.cleavage_params.miscleavage = 1
        peptides = graph.call_variant_peptides()
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSR', 'SSSSRSSSIR', 'SSSIR', 'SSSIRSSSSS'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_stop_gain(self):
        """ test stop gain mutation is included """
        variant_1 = (25, 26, 'T', 'A', 'SNV', '0:T-A', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], [((0,5),(0,5))], 0),
            2: ('SS', [1], [None], [((0,2), (5,7))], 0),
            3: ('*', [2], [variant_1], [], 0),
            4: ('SR', [3], [None], [((3,5),(8,10))], 0),
            # 2: ('SS*SR', [1],[variant_1], [((0,2), (5,7)),((3,5),(8,10))], 0),
            5: ('SSSIR', [1], [None], locations, 0),
            6: ('SSSPK', [4,5], [None], locations, 0),
            7: ('SSSVK', [6], [None], locations, 0)
        }
        graph, _ = create_pgraph(data, 'ENST0001', known_orf=[0,60])
        graph.cleavage_params.miscleavage = 1
        peptides = graph.call_variant_peptides()
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSS', 'SS'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_truncated(self):
        """ test tuncated peptide is handled properly """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations, 0),
            2: ('SSSSR', [1],[variant_1], locations, 0),
            3: ('SSSIR', [1], [None], locations, 0),
            4: ('SSSPK', [2,3], [None], locations, 0),
            5: ('SSSVK', [4], [None], locations, 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001', known_orf=[0,60])
        graph.cleavage_params.miscleavage = 2
        nodes[5].truncated = True
        peptides = graph.call_variant_peptides()
        received = {str(x.seq) for x in peptides}
        expected = {'SSSSKSSSSR', 'SSSSR', 'SSSSKSSSSRSSSPK', 'SSSSRSSSPK'}
        self.assertEqual(received, expected)

    def test_call_variant_peptides_small_orf(self):
        """ test very small ORF is handled """
        variant_1 = (10, 11, 'C', 'T', 'SNV', '1:C-T', 0, 1, True)
        data = {
            1: ('SSSSK', [0], [None], [((0,5),(0,5))], 0),
            2: ('SMSS', [1], [None], [((0,4),(5,9))], 0),
            3: ('*', [2], [variant_1], [], 0),
            4: ('R', [3], [None], [((5,6),(10,11))], 0),
            # 2: ('SMSS*R', [1],[variant_1], [((0,3),(5,8)), ((4,6),(9,11))], 0),
            5: ('SMSI', [1], [None], [((0,4),(5,9))], 0),
            6: ('*', [5], [None], [((0,1),(9,10))], 0),
            7: ('R', [6], [None], [((0,1),(10,11))], 0),
            #3: ('SMSI*R', [1], [None], [((0,6),(5,11))], 0),
            8: ('SSSPK', [4,7], [None], [((0,5),(11,16))], 0),
            9: ('SSSVK', [8], [None], [((0,5),(16,21))], 0)
        }
        graph, _ = create_pgraph(data, 'ENST0001', known_orf=[18, 27])
        peptides = graph.call_variant_peptides()
        received = {str(x.seq) for x in peptides}
        expected = {'MSS', 'SS'}
        self.assertEqual(received, expected)

    def test_call_peptides_check_orf(self):
        """ test calling peptides when checking for ORF """
        locations = [((0,4),(0,4))]
        data = {
            1: ('MSSSK', [0], [], locations, 0),
            2: ('MSSYK', [0], [], locations, 1),
            3: ('SSSSR', [1], [], locations, 0),
            4: ('SSSSR', [2], [], locations, 1)
        }
        graph_id = 'ENST0001'
        graph, _ = create_pgraph(data, graph_id, )
        peptides = graph.call_variant_peptides(check_variants=False, check_orf=True)

        received = {str(x.seq) for x in peptides}
        expected = {'MSSSK', 'MSSYK', 'SSSSR', 'SSYK', 'SSSK'}
        self.assertEqual(received, expected)

        expected = {
            f'{graph_id}|ORF1|1', f'{graph_id}|ORF2|1',
            f'{graph_id}|ORF1|2', f'{graph_id}|ORF2|2',
            f'{graph_id}|ORF1|3', f'{graph_id}|ORF2|3'
        }
        received = set()
        for peptide in peptides:
            received.update(peptide.description.split(' '))
        self.assertEqual(expected, received)

    def test_call_peptides_id(self):
        """ test peptide ids are named correctly. """
        variant_1 = (0, 3, 'TCT', 'T', 'INDEL', '0:TCT-T', 0, 1, True)
        locations = [((0,4),(0,4))]
        data = {
            1: ('SSSSK', [0], [None], locations, 0),
            2: ('SSSSR', [1],[variant_1], locations, 0),
            3: ('SSSIR', [1], [None], locations, 0),
            4: ('SSXPK', [2,3], [None], locations, 0)
        }
        graph, _ = create_pgraph(data, 'ENST0001')
        peptides = graph.call_variant_peptides(1)
        for peptide in peptides:
            self.assertTrue(peptide.description.startswith(graph.id))

    def test_call_peptides_empty_graph(self):
        """ When the graph is empty """
        data = {}
        graph, _ = create_pgraph(data, 'ENST0001')
        graph.add_stop(graph.root)
        peptides = graph.call_variant_peptides(2)
        self.assertEqual(len(peptides), 0)

    def test_call_and_stage_known_orf_incds(self):
        """ Test call and stage known orf """
        variant_1 = (8, 9, 'T', 'A', 'INDEL', '8:TCT-T', 0, 1, True)
        data = {
            1: ('SSSSK', [0], [None], [((0,5),(0,5))], 0),
            2: ('SSSSR', [1],[variant_1], [((1,5),(6,10))], 0),
            3: ('SSXPK', [2], [None], [((0,5),(10,15))], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.known_orf = [0,30]
        pool = VariantPeptideDict(graph.id)
        traversal = PVGTraversal(True, False, pool, (0,30), (0,10))
        orf = PVGOrf([0, None])
        cursor = PVGCursor(nodes[1], nodes[2], True, [orf])
        graph.call_and_stage_known_orf(cursor,  traversal)

        received = {str(x.seq) for x in traversal.pool.peptides.keys()}
        expected = {str(nodes[2].seq.seq)}
        self.assertEqual(received, expected)

        self.assertTrue(traversal.queue[-1].in_cds)

    def test_call_and_stage_known_orf_start_and_frameshift(self):
        """ Test when a frameshift mutation is in the same node as start codon
        """
        variant_1 = (8, 11, 'TCT', 'T', 'INDEL', '8:TCT-T', 2, 3, True)
        data = {
            1: ('SSSSK', [0], [None], [((0,5),(0,5))], 0),
            2: ('*MSSR', [1],[variant_1], [((0,3),(5,8)), ((4,5),(9,10))], 0),
            3: ('*MSSK', [1], [None], [((0,5),(5,10))], 0),
            4: ('SSSPK', [2], [None], [((0,5),(10,15))], 0),
            5: ('SSSPK', [4], [None], [((0,5),(15,20))], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        nodes[2].reading_frame_index = 0
        nodes[4].reading_frame_index = 2
        graph.known_orf = [18,60]
        pool = VariantPeptideDict(graph.id)
        traversal = PVGTraversal(True, False, pool, (18,60), (6,20))
        orf = PVGOrf([0, None], set())
        cursor = PVGCursor(nodes[1], nodes[2], False, [orf])
        graph.call_and_stage_known_orf(cursor,  traversal)
        self.assertEqual(len(traversal.queue), 1)
        self.assertTrue(traversal.queue[-1].in_cds)
        self.assertEqual(len(traversal.queue[-1].orfs[0].start_gain), 1)

    def test_call_and_stage_known_orf_multiple_methionine(self):
        """ Test when a mutation is in the same cleavage node as start codon,
        and there is another M before the actual start codon.
        """
        variant_1 = (27, 28, 'T', 'TCCC', 'INDEL', '8:T-TCCC', 9, 10, True)
        data = {
            1: ('SSSSK', [0], [None], [((0,5),(0,5))], 0),
            2: ('SMSMRIK', [1],[variant_1], [((0,5),(5,10)), ((6,7),(10,11))], 0),
            3: ('SMSMRK', [1], [None], [((0,6),(5,11))], 0),
            4: ('SSSPK', [2], [None], [((0,5),(11,16))], 0),
            5: ('SSSPK', [4], [None], [((0,5),(16,21))], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.known_orf = [24,90]
        pool = VariantPeptideDict(graph.id)
        traversal = PVGTraversal(True, False, pool, (24,90), (8,30))
        cursor = PVGCursor(nodes[1], nodes[2], False, [0, None], [])
        graph.call_and_stage_known_orf(cursor,  traversal)
        self.assertEqual(len(pool.peptides), 2)
        seqs = {str(x.seq) for x in pool.peptides.keys()}
        expected = {'MRIK', 'RIK'}
        self.assertEqual(seqs, expected)

    def test_call_and_stage_known_orf_start_altering(self):
        """ Test when the transcript is not cds_start_NF, and the mutation is
        start altering and.
        """
        variant_1 = (7, 8, 'T', 'TCCC', 'INDEL', '8:T-TCCC', 2, 3, True)
        data = {
            1: ('SSMSK', [0], [None], [((0,5), (0,5))], 0),
            2: ('SSSSK', [0], [variant_1], [], 0),
            3: ('SMSMRK', [1, 2], [None], [((0,6),(5,11))], 0),
            4: ('SSSPK', [3], [None], [((0,5),(11,16))], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.cds_start_nf = True
        graph.known_orf = [6,90]
        pool = VariantPeptideDict(graph.id)
        traversal = PVGTraversal(True, False, pool, (6,90), (2,30))
        cursor = PVGCursor(graph.root, nodes[2], False, [0, None], [])
        graph.call_and_stage_known_orf(cursor,  traversal)
        self.assertEqual(len(pool.peptides), 0)

    def test_call_and_stage_known_orf_cds_stop_gain(self):
        """ Test when the transcript is cds_start_NF, and the mutation is
        start altering.
        """
        variant_1 = (4, 5, 'T', 'A', 'SNV', '8:T-A', 0, 1, True)
        data = {
            1: ('SSSSK', [0], [None], [((0,5), (0,5))], 0),
            2: ('SMSMRK', [1], [None], [((0,6),(5,11))], 0),
            3: ('SSSPK', [2], [None], [((0,5),(11,16))], 0),
            4: ('SS', [2], [None], [((0,2), (11,13))], 0),
            5: ('*', [4], [variant_1], [], 0),
            6: ('PK', [5], [None], [((0,2), (14,16))], 0)
            # 4: ('SS*PK', [2], [variant_1], [((0,2),(11,13)), ((3,5),(14,16))], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        nodes[5].variants[0].is_stop_altering = True
        graph.known_orf = [0,90]
        pool = VariantPeptideDict(graph.id)
        traversal = PVGTraversal(True, False, pool, (0,90), (0,30))
        orf = PVGOrf([0, None])
        cursor = PVGCursor(nodes[2], nodes[4], True, [orf])
        graph.call_and_stage_known_orf(cursor,  traversal)
        self.assertEqual(len(pool.peptides), 1)
        seqs = {str(x.seq) for x in pool.peptides.keys()}
        expected = {'SS'}
        self.assertEqual(seqs, expected)

    def test_call_and_stage_known_orf_cds_stop_lost(self):
        """ Test when the transcript is cds_start_NF, and the mutation is
        start altering.
        """
        variant_1 = (41, 42, 'G', 'C', 'SNV', '41:G-C', 2, 3, True)
        data = {
            1: ('SSSSK', [0], [None], [((0,5), (0,5))], 0),
            2: ('SMSMRK', [1], [None], [((0,6),(5,11))], 0),
            3: ('SS', [2], [None], [((0,2), (11,13))], 0),
            4: ('*', [3], [None], [((0,1),(13,14))], 0),
            5: ('PK', [4], [None], [((0,2), (14,16))], 0),
            6: ('SSSPK', [2], [variant_1], [((0,2),(11,13)), ((3,5),(14,16))], 0),
            7: ('SPPPK', [5,6], [None], [((0,5),(17,22))], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.known_orf = [0,39]
        pool = VariantPeptideDict(graph.id)
        traversal = PVGTraversal(True, False, pool, (0,42), (0,14))
        orf = PVGOrf([0,None])
        cursor = PVGCursor(nodes[4], nodes[5], False, [orf], [])
        graph.call_and_stage_known_orf(cursor,  traversal)
        orf = PVGOrf([0, None])
        cursor = PVGCursor(nodes[2], nodes[6], True, [orf])
        graph.call_and_stage_known_orf(cursor,  traversal)
        self.assertEqual(len(traversal.queue[0].orfs[0].start_gain), 1)

    def test_call_and_stage_known_orf_not_in_cds_no_start_altering(self):
        """ Test when the transcript doesn't have any start altering variants,
        M not going to be treated as new start codon.
        """
        variant_1 = (66, 67, 'G', 'C', 'SNV', '41:G-C', 2, 3, True)
        data = {
            1: ('SSSSK', [0], [None], [((0,5), (0,5))], 0),
            2: ('SMSMRK', [1], [None], [((0,6),(5,11))], 0),
            3: ('SS*PGGK', [2], [None], [((0,5),(11,16))], 0),
            4: ('JIWIK', [3], [None], [((0,5),(17,22))], 0),
            5: ('JIMIT', [3], [variant_1], [((0,2),(17,19)),((3,5),(20,22))], 0),
            6: ('VSSSP', [4,5], [None], [((0,5),(22,27))], 0),
            7: ('JIWQT', [6], [None], [((0,5),(27,32))], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.known_orf = [18,39]
        pool = VariantPeptideDict(graph.id)
        traversal = PVGTraversal(True, False, pool, (6,13), (18,39))
        cursor = PVGCursor(nodes[5], nodes[6], False, [0, None], [])
        graph.call_and_stage_known_orf(cursor,  traversal)
        self.assertFalse(traversal.queue[0].in_cds)

    def test_call_and_stage_known_orf_cleavage_gain_downstream(self):
        """ Test when a cleavage gain was caused by downstream node.
        """
        variant_1 = (0, 1, 'C', 'T', 'SNV', '0:C-T', 0, 1, True)
        data = {
            1: ('SSSSK', [0], [None], [((0,5),(0,5))], 0),
            2: ('SSSSSSSSSR', [1],[None], [((0,10),(5,15))], 0),
            3: ('SSSSK', [1],[None], [((0,5),(5,10))], 0),
            4: ('SSSSR', [3], [variant_1], [((1,5),(11,15))], 0),
            5: ('SSSPK', [2,4], [None], [((0,5), (15,20))], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        graph.known_orf = [0,39]
        pool = VariantPeptideDict(graph.id)
        traversal = PVGTraversal(True, False, pool, (6,13), (18,39))
        orf = PVGOrf([0, None])
        cursor = PVGCursor(nodes[1], nodes[3], True, [orf])
        graph.call_and_stage_known_orf(cursor,  traversal)
        label = list(list(traversal.pool.peptides.values())[0])[0].label
        self.assertEqual(label.count('|'), 1)

    def test_fit_into_cleavage_bridge_node_needs_merge(self):
        r""" Test the fit into cleavage for bridge node that needs to be merged
        forward
                     Q                       WSPYQ-
                    /                       /      \
        0      W-SPY-QT               0     -WSPYQT-
                /                                  /|
              NG              ->            NGSPYQ- |
             /                             /       /
        1 VLR-NGALT                       | NGSPYQT
                                          |/
                                      1 VLR-NGALT
        """
        data = {
            1: ('VLR',  [0],  [None], [((0,3),( 0, 3))], 1),
            2: ('NG',   [1],  [None], [((0,1),( 3, 4))], 1),
            3: ('NGALT',[1],  [None], [((0,5),( 3, 8))], 1),
            4: ('W',    [0],  [None], [((0,1),(0,1))], 0),
            5: ('SPY',  [4,2],[None], [((0,3),(1,4))], 0),
            6: ('QT',   [5],  [None], [((0,2),(4,6))], 0),
            7: ('Q',    [5],  [None], [], 0)
        }
        graph, nodes = create_pgraph(data, 'ENST0001')
        branches,_ = graph.fit_into_cleavages_single_upstream(nodes[4])
        self.assertEqual(len(branches), 2)
        actual = {str(x.seq.seq) for x in nodes[1].out_nodes}
        expect = {'NGSPYQT', 'NGSPYQ'}
        self.assertTrue(expect.issubset(actual))
