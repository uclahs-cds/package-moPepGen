""" Test indexGVF """
import argparse
import shutil
import subprocess as sp
import sys
from test.integration import TestCaseIntegration
from moPepGen import cli


class TestParseIndexGVF(TestCaseIntegration):
    """ Test cases for moPepGen indexGVF """
    def test_index_gvf_cli(self):
        """ Test indexGVF cli """
        filename = 'vep_gSNP.gvf'
        gvf_file = self.data_dir/'vep'/filename
        input_file = self.work_dir/filename
        shutil.copy2(gvf_file, input_file)
        cmd = f"""
        {sys.executable} -m moPepGen.cli indexGVF -i {input_file}
        """
        res = sp.run(cmd, shell=True, check=False, capture_output=True)
        try:
            self.assertEqual(res.returncode, 0)
        except:
            print(cmd)
            print(res.stderr.decode('utf-8'))
            raise

    def test_index_gvf(self):
        """ Test parseFusionCatcher """
        filename = 'vep_gSNP.gvf'
        gvf_file = self.data_dir/'vep'/filename
        input_file = self.work_dir/filename

        shutil.copy2(gvf_file, input_file)
        args = argparse.Namespace()
        args.command = 'indexGVF'
        args.input_path = input_file
        args.quiet = True
        cli.index_gvf(args)
        files = {str(file.name) for file in self.work_dir.glob('*')
            if str(file.name) != filename}
        expected = {f"{filename}.idx"}
        self.assertEqual(files, expected)
