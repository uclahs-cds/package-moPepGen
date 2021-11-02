""" Test Module for VariantPeptideDict """
import unittest
from test.unit import create_aa_record, create_variants
from moPepGen.svgraph.VariantPeptideDict import VariantPeptideDict, \
    VariantPeptideMetadata


def create_variant_peptide_dict(tx_id, data) -> VariantPeptideDict:
    """ create a VariantPeptideDict """
    peptides = {}
    for x,y in data:
        seq = create_aa_record(*x)
        metadatas = set()
        for it in y:
            variants = set(create_variants(it[0]))
            metadatas.add(VariantPeptideMetadata(variants, it[1]))
        peptides[seq] = metadatas
    return VariantPeptideDict(tx_id=tx_id, peptides=peptides)

class TestCaseVariantPeptideDict(unittest.TestCase):
    """ Test cases for VariantPeptideDict """
    def test_get_peptide_sequences(self):
        """ Get peptide sequence """
        tx_id = 'ENST0001'
        data = [
            (
                ('SSSSSSSSSR', tx_id), [
                    ([(100, 101, 'A', 'T', 'SNV', 'SNV-100-A-T')], [0, None]),
                    ([(200, 201, 'G', 'C', 'SNV', 'SNV-200-G-C')], [0, None])
                ]
            )
        ]
        pool = create_variant_peptide_dict(tx_id, data)
        seqs = pool.get_peptide_sequences()
        self.assertEqual({str(x.seq) for x in seqs}, {'SSSSSSSSSR'})
        seq = list(seqs)[0]
        peptide_id = {f'{tx_id}|SNV-100-A-T|1', f'{tx_id}|SNV-200-G-C|1'}
        self.assertEqual(set(seq.description.split(' ')), peptide_id)

    def test_get_peptide_sequences_circ_rna(self):
        """ Get peptide sequence with circRNA """
        tx_id = 'ENST0001'
        variants = [
            (100, 101, 'A', '<circRNA>', 'circRNA', 'CIRC-ENST0001-E1-E2-E3'),
            (100, 101, 'A', '<circRNA>', 'circRNA', 'CIRC-ENST0001-E1-E2-E3')
        ]
        data = [
            (
                ('SSSSSSSSSR', tx_id), [
                    ([variants[0]], [0, None]),
                    ([variants[1]], [1, None])
                ]
            )
        ]
        pool = create_variant_peptide_dict(tx_id, data)
        seqs = pool.get_peptide_sequences()
        self.assertEqual({str(x.seq) for x in seqs}, {'SSSSSSSSSR'})
        self.assertEqual(list(seqs)[0].description, 'CIRC-ENST0001-E1-E2-E3|1')