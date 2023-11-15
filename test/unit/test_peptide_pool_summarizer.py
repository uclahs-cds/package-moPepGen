""" Module for testing PeptidePoolSummarizer """
import copy
import io
from contextlib import redirect_stdout
import unittest
from test.unit import create_aa_record, create_genomic_annotation, get_tx2gene_and_coding_tx
from test.unit.test_peptide_pool_splitter import (
    LABEL_MAP1, SOURCE_ORDER, ANNOTATION_DATA,
)
from moPepGen.aa.PeptidePoolSummarizer import PeptidePoolSummarizer
from moPepGen.aa.PeptidePoolSplitter import LabelSourceMapping
from moPepGen.aa import VariantPeptidePool


SOURCE_PARSER_MAP = {
    'gSNP': 'parseVEP',
    'gINDEL': 'parseVEP',
    'sSNV': 'parseVEP',
    'sINDEL': 'parseVEP',
    'altSplice': 'parseRMATS',
    'Fusion': 'parseSTARFusion',
    'circRNA': 'parseCIRCExplorer'
}

class TestPeptidePoolSummarizer(unittest.TestCase):
    """ Test cases for PeptidePoolSummarizer """
    def test_summarize_fasta_case1(self):
        """ basic test """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [[ 'SSSSSSSR', 'ENST0001|SNV-1001-T-A|1' ]]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        summarizer = PeptidePoolSummarizer(
            peptides, order=copy.copy(SOURCE_ORDER), label_map=label_map,
        )
        summarizer.count_peptide_source(
            tx2gene=tx2gene,
            coding_tx=coding_tx,
            enzyme='trypsin'
        )
        self.assertEqual(set(summarizer.summary_table.data.keys()), {frozenset(['gSNP'])})

    def test_summarize_fasta_source_comb_order(self):
        """ When source combination is present in --order-source """
        anno = create_genomic_annotation(ANNOTATION_DATA)
        anno.transcripts['ENST0005'] = copy.deepcopy(anno.transcripts['ENST0002'])
        anno.transcripts['ENST0005'].is_protein_coding = False
        tx2gene, coding_tx = get_tx2gene_and_coding_tx(anno)
        peptides_data = [
            [
                'SSSSSSSR',
                'CIRC-ENST0002-E1-E2|1 ENST0005|SE-2100|1'
            ]
        ]
        peptides = VariantPeptidePool({create_aa_record(*x) for x in peptides_data})
        label_map = LabelSourceMapping(copy.copy(LABEL_MAP1))
        # order = copy.copy(SOURCE_ORDER)
        order = {
            'altSplice': 1,
            frozenset(['altSplice', 'Noncoding']): 2,
            'Noncoding': 3,
            'circRNA': 4
        }
        source_parser_map = copy.deepcopy(SOURCE_PARSER_MAP)
        summarizer = PeptidePoolSummarizer(
            peptides, order=order, label_map=label_map, source_parser_map=source_parser_map
        )
        summarizer.count_peptide_source(
            tx2gene=tx2gene,
            coding_tx=coding_tx,
            enzyme='trypsin'
        )
        self.assertEqual(
            set(summarizer.summary_table.data.keys()),
            {frozenset(['altSplice', 'Noncoding'])}
        )

        handle = io.StringIO()
        with redirect_stdout(handle):
            summarizer.write_summary_table(handle)
