""" Test Module for VariantPeptideDict """
import unittest
from test.unit import create_variants
from Bio.Seq import Seq
from moPepGen import params
from moPepGen.svgraph.PVGPeptideFinder import PVGCandidateNodePaths, PVGPeptideFinder, \
    PVGPeptideMetadata
from moPepGen.svgraph.PVGOrf import PVGOrf
import moPepGen.aa.VariantPeptideIdentifier as vpi
from test.unit import create_pgraph


def create_pvg_peptide_finder(tx_id, data) -> PVGPeptideFinder:
    """ create a VariantPeptideDict """
    peptides = {}
    for x,y in data:
        seq = Seq(x)
        metadatas = {}
        for it in y:
            variants = set(create_variants(it[0]))
            label = vpi.create_variant_peptide_id(tx_id, variants, None)
            is_pure_circ_ran = len(variants) == 1 and list(variants)[0].is_circ_rna()
            metadata = PVGPeptideMetadata(label, it[1], is_pure_circ_ran)
            metadata.has_variants = bool(variants)
            metadatas[metadata.label] = metadata
        peptides[seq] = metadatas
    return PVGPeptideFinder(tx_id=tx_id, peptides=peptides)

class TestCasePVGPeptideFinder(unittest.TestCase):
    """ Test cases for VariantPeptideDict """
    def test_get_peptide_sequences(self):
        """ Get peptide sequence """
        tx_id = 'ENST0001'
        data = [
            (
                'SSSSSSSSSR', [
                    (
                        [(100, 101, 'A', 'T', 'SNV', 'SNV-100-A-T', {}, 'ENST0001')],
                        [0, None]
                    ), (
                        [(200, 201, 'G', 'C', 'SNV', 'SNV-200-G-C', {}, 'ENST0001')],
                        [0, None]
                    )
                ]
            )
        ]
        finder = create_pvg_peptide_finder(tx_id, data)
        res = finder.get_peptide_sequences()
        self.assertEqual({str(x) for x in res}, {'SSSSSSSSSR'})
        seq_data = list(res.values())[0]
        peptide_id = {f'{tx_id}|SNV-100-A-T|1', f'{tx_id}|SNV-200-G-C|1'}
        self.assertEqual({x.label for x in seq_data}, peptide_id)

    def test_get_peptide_sequences_circ_rna(self):
        """ Get peptide sequence with circRNA """
        tx_id = 'ENST0001'
        variants = [
            (100, 101, 'A', '<circRNA>', 'circRNA', 'CIRC-ENST0001-E1-E2-E3'),
            (100, 101, 'A', '<circRNA>', 'circRNA', 'CIRC-ENST0001-E1-E2-E3')
        ]
        data = [
            (
                'SSSSSSSSSR', [
                    ([variants[0]], [0, None]),
                    ([variants[1]], [1, None])
                ]
            )
        ]
        finder = create_pvg_peptide_finder(tx_id, data)
        seqs = finder.get_peptide_sequences()
        self.assertEqual({str(x) for x in seqs}, {'SSSSSSSSSR'})
        self.assertEqual(list(seqs.values())[0][0].label, 'CIRC-ENST0001-E1-E2-E3|1')

    def test_find_candidate_node_paths_archipel_1(self):
        """ Find candidate node path for archipel peptides """
        gene_id = 'ENSG0001'
        tx_id = 'ENST0001'
        v1 = (41, 42, 'G', 'CTTTT', 'INDEL', '41:G-CTTTT', 0, 3, True)
        data = {
            1: ('S',  [0],   [None], [((0,1),(0,1))], 0),
            2: ('K',  [1],   [None], [((0,1),(1,2))], 0),
            3: ('L',  [2],   [None], [((0,1),(2,3))], 0),
            4: ('HY', [3],   [v1  ], [],              0),
            5: ('V',  [3],   [None], [((0,1),(3,4))], 0),
            6: ('C',  [4,5], [None], [((0,1),(4,5))], 0),
            7: ('W',  [6],   [None], [((0,1),(5,6))], 0),
            8: ('I',  [7],   [None], [((0,1),(6,7))], 0),
            9: ('*',  [8],   [None], [((0,1),(7,8))], 0)
        }
        orf = PVGOrf(orf = [0, 7])
        _, nodes = create_pgraph(data, 'ENST0001')
        finder = PVGPeptideFinder(tx_id)
        paths = finder.find_candidate_node_paths_archipel(
            node=nodes[1], orfs=[orf], cleavage_params=None,
            tx_id=tx_id, gene_id=gene_id, leading_node=None, subgraphs=None,
            is_circ_rna=False, backsplicing_only=False, is_start_codon=False,
            flanking_size=3
        )
        self.assertEqual(len(paths.data), 1)
        self.assertEqual(
            ''.join([str(x.seq.seq) for x in paths.data[0].nodes]),
            'SKLHYCWI'
        )

    def test_find_candidate_node_paths_archipel_2(self):
        """ Find candidate node path for archipel peptides """
        gene_id = 'ENSG0001'
        tx_id = 'ENST0001'
        v1 = (41, 42, 'G', 'CTTTT', 'INDEL', '41:G-CTTTT', 0, 3, True)
        data = {
            1: ('S',  [0],   [None], [((0,1),(0,1))], 0),
            2: ('K',  [1],   [None], [((0,1),(1,2))], 0),
            3: ('L',  [2],   [None], [((0,1),(2,3))], 0),
            4: ('HY', [3],   [v1  ], [],              0),
            5: ('V',  [3],   [None], [((0,1),(3,4))], 0),
            6: ('C',  [4,5], [None], [((0,1),(4,5))], 0),
            7: ('W',  [6],   [None], [((0,1),(5,6))], 0),
            8: ('*',  [7],   [None], [((0,1),(6,7))], 0)
        }
        orf = PVGOrf(orf = [0, 7])
        _, nodes = create_pgraph(data, 'ENST0001')
        finder = PVGPeptideFinder(tx_id)
        paths = finder.find_candidate_node_paths_archipel(
            node=nodes[1], orfs=[orf], cleavage_params=None,
            tx_id=tx_id, gene_id=gene_id, leading_node=None, subgraphs=None,
            is_circ_rna=False, backsplicing_only=False, is_start_codon=False,
            flanking_size=3
        )
        self.assertEqual(len(paths.data), 1)
        self.assertEqual(
            ''.join([str(x.seq.seq) for x in paths.data[0].nodes]),
            'SKLHYCW'
        )

    def test_find_candidate_node_paths_archipel_3(self):
        """ Find candidate node path for archipel peptides """
        gene_id = 'ENSG0001'
        tx_id = 'ENST0001'
        v1 = (41, 42, 'G', 'CTTTT', 'INDEL', '41:G-CTTTT', 0, 3, True)
        data = {
            1: ('M',  [0],   [None], [((0,1),(1,2))], 0),
            2: ('L',  [1],   [None], [((0,1),(2,3))], 0),
            3: ('HY', [2],   [v1  ], [],              0),
            4: ('V',  [2],   [None], [((0,1),(3,4))], 0),
            5: ('C',  [3,4], [None], [((0,1),(4,5))], 0),
            6: ('W',  [5],   [None], [((0,1),(5,6))], 0),
            7: ('I',  [6],   [None], [((0,1),(6,7))], 0),
            8: ('*',  [7],   [None], [((0,1),(7,8))], 0)
        }
        orf = PVGOrf(orf = [0, 7])
        _, nodes = create_pgraph(data, 'ENST0001')
        finder = PVGPeptideFinder(tx_id)
        paths = finder.find_candidate_node_paths_archipel(
            node=nodes[1], orfs=[orf], cleavage_params=None,
            tx_id=tx_id, gene_id=gene_id, leading_node=None, subgraphs=None,
            is_circ_rna=False, backsplicing_only=False, is_start_codon=True,
            flanking_size=3
        )
        self.assertEqual(len(paths.data), 1)
        self.assertEqual(
            ''.join([str(x.seq.seq) for x in paths.data[0].nodes]),
            'MLHYCWI'
        )

    def test_find_candidate_node_paths_archipel_4(self):
        """ Find candidate node path for archipel peptides """
        gene_id = 'ENSG0001'
        tx_id = 'ENST0001'
        v1 = (41, 42, 'G', 'CTTTT', 'INDEL', '41:G-CTTTT', 0, 3, True)
        data = {
            1: ('MHY', [0],  [v1  ], [],              0),
            2: ('V',  [0],   [None], [((0,1),(3,4))], 0),
            3: ('C',  [1,2], [None], [((0,1),(4,5))], 0),
            4: ('W',  [3],   [None], [((0,1),(5,6))], 0),
            5: ('I',  [4],   [None], [((0,1),(6,7))], 0),
            6: ('*',  [5],   [None], [((0,1),(7,8))], 0)
        }
        orf = PVGOrf(orf = [0, 7])
        _, nodes = create_pgraph(data, 'ENST0001')
        finder = PVGPeptideFinder(tx_id)
        paths = finder.find_candidate_node_paths_archipel(
            node=nodes[1], orfs=[orf], cleavage_params=None,
            tx_id=tx_id, gene_id=gene_id, leading_node=None, subgraphs=None,
            is_circ_rna=False, backsplicing_only=False, is_start_codon=True,
            flanking_size=3
        )
        self.assertEqual(len(paths.data), 1)
        self.assertEqual(
            ''.join([str(x.seq.seq) for x in paths.data[0].nodes]),
            'MHYCWI'
        )

    def test_find_candidate_node_paths_archipel_5(self):
        """ Find candidate node path for archipel peptides """
        gene_id = 'ENSG0001'
        tx_id = 'ENST0001'
        v1 = (41, 42, 'G', 'CTTTT', 'INDEL', '41:G-CTTTT', 0, 3, True)
        v2 = (51, 52, 'G', 'CTTTT', 'INDEL', '51:G-CTTTT', 0, 3, True)
        data = {
             1: ('S',  [0],   [None], [((0,1),(0,1))], 0),
             2: ('K',  [1],   [None], [((0,1),(1,2))], 0),
             3: ('L',  [2],   [None], [((0,1),(2,3))], 0),
             4: ('HY', [3],   [v1  ], [],              0),
             5: ('V',  [3],   [None], [((0,1),(3,4))], 0),
             6: ('C',  [4,5], [None], [((0,1),(4,5))], 0),
             7: ('W',  [6],   [None], [((0,1),(5,6))], 0),
             8: ('WE', [6],   [v2],   [],              0),
             9: ('I',  [7,8],   [None], [((0,1),(6,7))], 0),
            10: ('Q',  [9],   [None], [((0,1),(6,7))], 0),
            11: ('N',  [10],   [None], [((0,1),(6,7))], 0),
            12: ('R',  [11],   [None], [((0,1),(6,7))], 0),
            13: ('*',  [12],   [None], [((0,1),(7,8))], 0)
        }
        orf = PVGOrf(orf = [0, 7])
        _, nodes = create_pgraph(data, 'ENST0001')
        finder = PVGPeptideFinder(tx_id)
        paths = finder.find_candidate_node_paths_archipel(
            node=nodes[1], orfs=[orf], cleavage_params=None,
            tx_id=tx_id, gene_id=gene_id, leading_node=None, subgraphs=None,
            is_circ_rna=False, backsplicing_only=False, is_start_codon=False,
            flanking_size=3
        )
        self.assertEqual(len(paths.data), 2)
        self.assertEqual(
            {''.join([str(x.seq.seq) for x in path.nodes]) for path in paths.data},
            {'SKLHYCWI', 'SKLHYCWEIQN'}
        )

class TestCasePVGCandidateNodePaths(unittest.TestCase):
    """ Test cases for MiscleavedNodes """
    def test_is_valid_x(self):
        """ Test that when X is found in the sequence, it is recognized as an
        invalid sequence. """
        cleavage_params = params.CleavageParams(enzyme='trypsin')
        misc_nodes = PVGCandidateNodePaths([], cleavage_params)
        self.assertFalse(misc_nodes.is_valid_seq('AAAAXAAA', set(), set()))
