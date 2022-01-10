""" Version """
from __future__ import annotations
import sys
import Bio
from moPepGen import __version__


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

    def __ne__(self, other:MetaVersion):
        """ not equal to """
        return not self == other

    def __repr__(self) -> str:
        """ str representation """
        python = '.'.join([str(x) for x in self.python])
        return f"python={python}, biopython={self.biopython}, moPepGen={self.mopepgen}"
