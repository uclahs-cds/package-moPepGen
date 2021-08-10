""" moPepGen """
from datetime import datetime


__version__ = '0.0.1'

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
