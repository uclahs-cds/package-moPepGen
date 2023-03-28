""" Version """
from __future__ import annotations
import sys
from typing import Tuple
import Bio
from moPepGen import __version__


MINIMAL_VERSION = '0.12.0'

class MetaVersion():
    """ Versions """
    def __init__(self):
        """ constructor """
        self.python = sys.version_info[:3]
        self.biopython = Bio.__version__
        self.mopepgen = __version__

    def __eq__(self, other:MetaVersion):
        """ equal to """
        return self.python == other.python and \
            self.biopython == other.biopython and \
            self.mopepgen == other.mopepgen

    @staticmethod
    def get_semver(version:str) -> Tuple[int, int, int]:
        """ Get major, minor and patch number from semver """
        return tuple(int(x) for x in version.split('-')[0].split('.'))

    def is_valid_mpg_version(self, version:str) -> bool:
        """ Check if the given moPepGen version is valid """
        that = self.get_semver(version)
        minimal = self.get_semver(MINIMAL_VERSION)
        return that >= minimal

    def is_valid(self, version:MetaVersion) -> bool:
        """ Check if the given version is valid """
        return self.python == version.python and \
            self.biopython == version.biopython and \
            self.is_valid_mpg_version(version.mopepgen)

    def __ne__(self, other:MetaVersion):
        """ not equal to """
        return not self == other

    def __repr__(self) -> str:
        """ str representation """
        python = '.'.join([str(x) for x in self.python])
        return f"python={python}, biopython={self.biopython}, moPepGen={self.mopepgen}"
