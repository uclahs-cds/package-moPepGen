""" """
from typing import IO
from moPepGen.version import MetaVersion

MOPEPGEN_GTF_INDEX_METADATA_PREFIX = '##'

class GTFIndexMetadata():
    """ """
    def __init__(self, source:str=None, version:MetaVersion=None):
        """ """
        self.source = source
        self.version = version or MetaVersion()

    def write(self, handle:IO):
        """ """
        prefix = MOPEPGEN_GTF_INDEX_METADATA_PREFIX
        handle.write(f"{prefix}source={self.source}\n")
        handle.write(f"{prefix}python={self.version.python}\n")
        handle.write(f"{prefix}moPepGen={self.version.mopepgen}\n")
        handle.write(f"{prefix}biopython={self.version.biopython}\n")

    @classmethod
    def parse(cls, handle:IO):
        """ """
        pos = handle.tell()
        it = handle.readline()
        metadata = {}
        version = {}
        prefix = MOPEPGEN_GTF_INDEX_METADATA_PREFIX
        while it and it.startswith(prefix):
            key, val = it.lstrip(prefix).rstrip().split('=', 1)
            if val == '':
                val = None
            if key in ('python', 'moPepGen', 'biopython'):
                version[key.lower()] = val
            else:
                metadata[key] = val

            pos = handle.tell()
            it = handle.readline()

        handle.seek(pos)

        version = MetaVersion(**version)

        return cls(source=metadata['source'], version=version)
