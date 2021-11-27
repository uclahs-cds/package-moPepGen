""" Version """
from __future__ import annotations
import sys
import Bio


__version__ = '0.0.1'


class MetaVersion():
    """ Versions """
    def __init__(self):
        """ constructor """
        self.python = sys.version_info[:3]
        self.biopython = Bio.__version__
        self.mopepgen = __version__

    def __eq__(self, other:MetaVersion):
        """ equal to """
        self.python == other.python
        self.biopython == other.biopython
        self.mopepgen == other.mopepgen

    def __ne__(self, other:MetaVersion):
        """ not equal to """
        return not self == other

    def __repr__(self) -> str:
        """ str representation """
        python = '.'.join(self.python)
        return f"python={python}, biopython={self.biopython}, moPepGen={self.mopepgen}"
