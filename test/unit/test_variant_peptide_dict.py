""" Test Module for VariantPeptideDict """
import unittest
from test.unit import create_aa_record, create_variants
from Bio.Seq import Seq
from moPepGen import params
from moPepGen.svgraph.VariantPeptideDict import MiscleavedNodes, VariantPeptideDict, \
    VariantPeptideMetadata
import moPepGen.aa.VariantPeptideIdentifier as vpi


def create_variant_peptide_dict(tx_id, data) -> VariantPeptideDict:
    """ create a VariantPeptideDict """
    peptides = {}
    for x,y in data:
        seq = Seq(x)
        metadatas = set()
        for it in y:
            variants = set(create_variants(it[0]))
            label = vpi.create_variant_peptide_id(tx_id, variants, None)
            is_pure_circ_ran = len(variants) == 1 and list(variants)[0].is_circ_rna()
            metadata = VariantPeptideMetadata(label, it[1], is_pure_circ_ran)
            metadata.has_variants = bool(variants)
            metadatas.add(metadata)
        peptides[seq] = metadatas
    return VariantPeptideDict(tx_id=tx_id, peptides=peptides)

class TestCaseVariantPeptideDict(unittest.TestCase):
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
        pool = create_variant_peptide_dict(tx_id, data)
        res = pool.get_peptide_sequences()
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
        pool = create_variant_peptide_dict(tx_id, data)
        seqs = pool.get_peptide_sequences()
        self.assertEqual({str(x) for x in seqs}, {'SSSSSSSSSR'})
        self.assertEqual(list(seqs.values())[0][0].label, 'CIRC-ENST0001-E1-E2-E3|1')


class TestCaseMiscleavedNodes(unittest.TestCase):
    """ Test cases for MiscleavedNodes """
    def test_is_valid_x(self):
        """ Test that when X is found in the sequence, it is recognized as an
        invalid sequence. """
        cleavage_params = params.CleavageParams(enzyme='trypsin')
        misc_nodes = MiscleavedNodes([], cleavage_params)
        self.assertFalse(misc_nodes.is_valid_seq('AAAAXAAA', set(), set()))
