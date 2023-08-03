""" Test module for VariantPeptideLabel """
import copy
import unittest
from test.unit import create_genomic_annotation
from Bio.Seq import Seq
from moPepGen import aa
from moPepGen.aa.VariantPeptideLabel import (
    VariantPeptideInfo, LabelSourceMapping, VariantSourceSet
)
from moPepGen.aa.VariantPeptideLabel import (
    SOURCE_SEC_TERMINATION, SOURCE_CODON_REASSIGNMENT
)
from test.unit.test_peptide_pool_splitter import (
    ANNOTATION_DATA, LABEL_MAP1, SOURCE_ORDER
)


class TestVariantPeptideInfo(unittest.TestCase):
    """ Test case for VariantPeptideInfo """
    def test_from_variant_peptide_fusion_only(self):
        """ Fusion only """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno = create_genomic_annotation(anno_data)
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        fusion_id = 'FUSION-ENST0001:1050-ENST0003:3090'
        peptide = aa.AminoAcidSeqRecord(
            seq=Seq('TTTTTR'),
            description=f"{fusion_id}|1"
        )
        levels = copy.copy(SOURCE_ORDER)
        VariantSourceSet.set_levels(levels)
        labels = VariantPeptideInfo.from_variant_peptide(
            peptide=peptide, anno=anno, label_map=label_map,
            check_source=True
        )
        self.assertEqual(labels[0].sources, {'Fusion'})

    def test_from_variant_peptide_fusion_w2f(self):
        """ Fusion + W2F """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno = create_genomic_annotation(anno_data)
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        fusion_id = 'FUSION-ENST0001:1050-ENST0003:3090'
        peptide = aa.AminoAcidSeqRecord(
            seq=Seq('TTTTTR'),
            description=f"{fusion_id}|W2F-2|1"
        )
        levels = copy.copy(SOURCE_ORDER)
        levels[SOURCE_SEC_TERMINATION] = max(levels.values()) + 1
        levels[SOURCE_CODON_REASSIGNMENT] = max(levels.values()) + 1
        VariantSourceSet.set_levels(levels)
        labels = VariantPeptideInfo.from_variant_peptide(
            peptide=peptide, anno=anno, label_map=label_map,
            check_source=True
        )
        self.assertEqual(labels[0].sources, {'Fusion', SOURCE_CODON_REASSIGNMENT})

    def test_from_variant_peptide_fusion_sect(self):
        """ Fusion + SECT """
        anno_data = copy.deepcopy(ANNOTATION_DATA)
        anno = create_genomic_annotation(anno_data)
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        fusion_id = 'FUSION-ENST0001:1050-ENST0003:3090'
        peptide = aa.AminoAcidSeqRecord(
            seq=Seq('TTTTTR'),
            description=f"{fusion_id}|SECT|1"
        )
        levels = copy.copy(SOURCE_ORDER)
        levels[SOURCE_SEC_TERMINATION] = max(levels.values()) + 1
        levels[SOURCE_CODON_REASSIGNMENT] = max(levels.values()) + 1
        VariantSourceSet.set_levels(levels)
        labels = VariantPeptideInfo.from_variant_peptide(
            peptide=peptide, anno=anno, label_map=label_map,
            check_source=True
        )
        self.assertEqual(labels[0].sources, {'Fusion', SOURCE_SEC_TERMINATION})
