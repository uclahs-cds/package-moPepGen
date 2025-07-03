""" Module for testing seqvar.io """
import unittest
import shutil
from pathlib import Path
from test.unit import create_variant
from moPepGen import seqvar


BASE_DIR = Path(__file__).parent.absolute()
WORK_DIR = BASE_DIR / 'work_dir'

class TestSeqvarIO(unittest.TestCase):
    """ Test case for parsing seqvar files. """
    def setUp(self):
        """ set up working directory """
        super().setUp()
        shutil.rmtree(WORK_DIR, ignore_errors=True)
        WORK_DIR.mkdir(parents=False, exist_ok=True)

    def tearDown(self):
        """ remove working files """
        super().tearDown()
        shutil.rmtree(WORK_DIR, ignore_errors=True)

    def test_seqvar_parse(self):
        """ Test parsing seqvar files. """
        gvf_path = 'test/files/vep/vep_gSNP.gvf'
        i = 0
        for record in seqvar.io.parse(gvf_path):
            i += 1
            self.assertIsInstance(record, seqvar.VariantRecord)
            if i > 5:
                break

    def test_coordinate_base_io(self):
        """ Test the writer is using 1-base position """
        attrs = {'GENE_ID': 'ENSG0001', 'START': '10', 'ACCEPTER_START': '10',
            'DONOR_START': '10'}
        variant = create_variant(10, 11, 'A', 'T', 'SNV', 'SNV-1', attrs)
        variant.location.seqname = 'ENST0001'
        output_file = WORK_DIR/'test.tvf'
        metadata = seqvar.GVFMetadata('parseXXX', source='XXX', chrom='Gene ID')
        seqvar.io.write([variant], output_file, metadata)
        with open(output_file, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                self.assertEqual(fields[1], '11')
                for item in fields[7].split(';'):
                    key,val =  item.split('=')
                    if 'START' in key:
                        self.assertEqual(int(val), 11)

        record = list(seqvar.io.parse(output_file))[0]
        self.assertEqual(record.location.start, 10)
        for key, val in record.attrs.items():
            if 'START' in key:
                self.assertEqual(int(val), 10)

    def test_parse_metadata(self):
        """ Test the metadata header is parsed correctly """
        gvf_path = 'test/files/vep/vep_gSNP_UCLA0001.gvf'
        with open(gvf_path, 'rt') as handle:
            metadata = seqvar.GVFMetadata.parse(handle)
            records = list(seqvar.io.parse(handle))
            for record in records:
                if 'PHASE_SETS' in record.attrs:
                    self.assertIsInstance(record.attrs['PHASE_SETS'], list)
        self.assertEqual(metadata.source, 'gSNP')
        self.assertEqual(metadata.parser, 'parseVEP')
        self.assertIn('TRANSCRIPT_ID', metadata.info)
        self.assertIn('PHASE_SETS', metadata.info)
        self.assertEqual(len(records), 16)
