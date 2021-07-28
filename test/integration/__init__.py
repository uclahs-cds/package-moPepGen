""" Integration Test """
from pathlib import Path
import shutil
import unittest


class TestCaseIntegration(unittest.TestCase):
    """ Test cases for the command line interface """

    def __init__(self, *args, **kwargs):
        """"""
        super().__init__(*args, **kwargs)
        self.base_dir = Path(__file__).parent.parent.absolute()
        self.work_dir = self.base_dir/'work_dir'
        self.data_dir = self.base_dir/'files'

    def setUp(self):
        """ set up working directory """
        super().setUp()
        shutil.rmtree(self.work_dir, ignore_errors=True)
        self.work_dir.mkdir(parents=False, exist_ok=True)

    def tearDown(self):
        """ remove working files """
        super().tearDown()
        shutil.rmtree(self.work_dir, ignore_errors=True)
