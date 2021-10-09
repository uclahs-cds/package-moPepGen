""" Module to test CircRNAModel """
import io
import unittest
from moPepGen import circ


GVF_DATA = [
    'ENSG00000128408.9\t0\tENSG00000128408.9-circRNA-E3-E4\t.\t.\t.\t.\t'
        'OFFSET=0,323;LENGTH=323,82;INTRON=;TRANSCRIPT=ENST00000614167.2;'
        'GENE_SYMBOL=RIBC2',
    'ENSG00000099949.21\t0\tENSG00000099949.21-circRNA-E1-E2\t.\t.\t.\t.\t'
        'OFFSET=0,98;LENGTH=78,42;INTRON=;TRANSCRIPT=ENST00000642151.1;'
        'GENE_SYMBOL=LZTR1',
    'ENSG00000099949.21\t0\tENSG00000099949.21-circRNA-E1-I1-E2\t.\t.\t.\t.\t'
        'OFFSET=0,78,98;LENGTH=78,20,42;INTRON=2;TRANSCRIPT=ENST00000642151.1;'
        'GENE_SYMBOL=LZTR1',
    'ENSG00000244486.9\t1376\tENSG00000244486.9-circRNA-E7-E8\t.\t.\t.\t.\t'
        'OFFSET=0,118;LENGTH=118,116;INTRON=;TRANSCRIPT=ENST00000622235.5;'
        'GENE_SYMBOL=SCARF2'
]

class TestCircRNA(unittest.TestCase):
    """ Test case for circ RNA """
    def test_parse_circ_rna_gvf(self):
        """ Test to parse circRNA gvf file """
        with io.BytesIO('\n'.join(GVF_DATA).encode('utf8')) as binary_file:
            with io.TextIOWrapper(binary_file, encoding='utf8') as handle:
                records = list(circ.io.parse(handle))
        for record in records:
            self.assertIsInstance(record, circ.io.CircRNAModel)
        self.assertEqual(len(records[0].fragments), 2)
