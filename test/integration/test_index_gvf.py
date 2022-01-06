""" Test indexGVF """
import argparse
import shutil
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseIndexGVF(TestCaseIntegration):
    """ Test cases for moPepGen indexGVF """
    def test_index_gvf(self):
        """ Test parseFusionCatcher """
        filename = 'vep_gSNP.gvf'
        gvf_file = self.data_dir/'vep'/filename
        input_file = self.work_dir/filename
        shutil.copy2(gvf_file, input_file)
        args = argparse.Namespace()
        args.command = 'indexGVF'
        args.input_gvf = input_file
        args.quiet = True
        cli.index_gvf(args)
        files = {str(file.name) for file in self.work_dir.glob('*')
            if str(file.name) != filename}
        expected = {f"{filename}.idx"}
        self.assertEqual(files, expected)
