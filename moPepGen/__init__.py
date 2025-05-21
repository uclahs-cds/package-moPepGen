""" moPepGen """
from __future__ import annotations
from datetime import datetime
import hashlib
import itertools
from typing import Iterable, IO
import logging
from . import constant


__version__ = '1.4.6'

## Error messages
ERROR_INDEX_IN_INTRON = 'The genomic index seems to be in an intron'
ERROR_NO_TX_AVAILABLE = 'No transcripts available'
ERROR_VARIANT_NOT_IN_GENE_COORDINATE = 'Variant is not in gene coordinates.'
ERROR_REF_LENGTH_NOT_MATCH_WITH_LOCATION = \
    'Length of ref does not match with the location.'
# global variables

# Used in the moPepGen.aa.VariantPeptidePool module to serve as a delimiter of
# variant peptide seqname if it is found in multiple transcripts. E.g.:
# ENSG0001|SNV-10-T-A|1 ENSG0002|INDEL-20-CCCC-C|2
VARIANT_PEPTIDE_SOURCE_DELIMITER = ' '

# Used in moPepGen.aa.PeptidePoolSplitter to server as the separater in split
# database keys. E.g.:
# gSNP-gINDEL
SPLIT_DATABASE_KEY_SEPARATER = '-'

# GVF header
GVF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

class _CaptureEq:
    """Object wrapper that remembers "other" for successful equality tests.
    Adopted from https://code.activestate.com/recipes/499299/
    """
    def __init__(self, obj):
        self.obj = obj
        self.match = obj

    def __eq__(self, other):
        result = (self.obj == other)
        if result:
            self.match = other
        return result

    def __getattr__(self, name):  # support hash() or anything else needed by __contains__
        return getattr(self.obj, name)

    def __hash__(self):
        return hash(self.obj)

def get_equivalent(container, item, default=None):
    '''Gets the specific container element matched by: "item in container".

    Adopted from https://code.activestate.com/recipes/499299/

    Useful for retreiving a canonical value equivalent to "item".  For example, a
    caching or interning application may require fetching a single representative
    instance from many possible equivalent instances).

    >>> get_equivalent(set([1, 2, 3]), 2.0)             # 2.0 is equivalent to 2
    2
    >>> get_equivalent([1, 2, 3], 4, default=0)
    0

    NOTE: Use this function with cautious because the direction of == can go
    in different directions, and sometimes the _CaptureEq.__eq__ is not called.
    moPepGen.aa.AminoAcidSeqRecord.__eq__ is an example of a work around.
    '''
    t = _CaptureEq(item)

    if t in container:
        return t.match
    return default


def logger(message:str) -> None:
    """ Print message to the stdout. """
    print(
        f'[ {datetime.now().strftime(format="%Y-%m-%d %H:%M:%S")} ] {message}',
        flush=True
    )

def get_logger() -> logging.Logger:
    """ get logger """
    return logging.getLogger(constant.PROG_NAME)

def all_equal(iterable:Iterable) -> bool:
    """ Check if all elements are equal """
    it = itertools.groupby(iterable)
    return next(it, True) and not next(it, False)

def check_sha512(handle:IO):
    """ check sha512sum. Handle must be readable """
    sum_val = hashlib.sha512()
    # pylint: disable=W0640
    for byte_block in iter(lambda: handle.read(4096), b""):
        sum_val.update(byte_block)
    return sum_val.hexdigest()
