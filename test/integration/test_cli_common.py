""" Test the cli.common module """
from test.integration import TestCaseIntegration
from moPepGen.cli import common


class TestCliCommon(TestCaseIntegration):
    """ Test cases for cli.common """
    def test_validate_file_format_readable_successful(self):
        """ check readable successful """
        file = self.work_dir/'test.txt'
        file.touch(0o666)
        common.validate_file_format(file, ['.txt'], check_readable=True)

    def test_validate_file_format_readable_permission_error(self):
        """ check readable failed because of permission error """
        file = self.work_dir/'test.txt'
        file.touch(mode=0o222)
        with self.assertRaises(PermissionError):
            common.validate_file_format(file, ['.txt'], check_readable=True)

    def test_validate_file_format_writable_successful(self):
        """ check writable successful """
        file = self.work_dir/'test.txt'
        file.touch(mode=0o666)
        common.validate_file_format(file, ['.txt'], check_writable=True)

    def test_validate_file_format_writable_file_permission_error(self):
        """ check writable failed because file is not writable """
        file = self.work_dir/'test.txt'
        file.touch(mode=0o555)
        with self.assertRaises(PermissionError):
            common.validate_file_format(file, ['.txt'], check_writable=True)

    def test_validate_file_format_writable_dir_permission_error(self):
        """ check writable failed because dir is not writable """
        dir_path = self.work_dir/'test'
        dir_path.mkdir(mode=0o555)
        file = dir_path/'test.txt'
        with self.assertRaises(PermissionError):
            common.validate_file_format(file, ['.txt'], check_writable=True)
