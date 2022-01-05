""" In disk class for VariantRecordPool """
from __future__ import annotations
from typing import Dict, IO, Iterable, TYPE_CHECKING, Union
from moPepGen import circ
from . import io


# To avoid circular import
if TYPE_CHECKING:
    from moPepGen.circ.CircRNA import CircRNAModel
    from moPepGen.seqvar import VariantRecord


class GVFPointer():
    """ File pointer """
    def __init__(self, key:str, start:int=None, end:int=None):
        """ constructor """
        self.key = key
        self.start = start
        self.end = end

    def __len__(self) -> int:
        """ length """
        return self.end - self.start

    def to_line(self) -> str:
        """ to line """
        return f"{self.key}\t{int(self.start)}\t{str(len(self))}"

def iterate_pointer(handle, is_circ_rna:bool) -> Iterable[GVFPointer]:
    """ """
    cur_key = None
    line_end = 0
    pointer = None
    for line in handle:

        line_start = line_end
        line_end += len(line)

        if line.startswith('#'):
            continue

        if is_circ_rna:
            record = circ.io.line_to_circ_model(line)
            key = record.id
        else:
            record = io.line_to_variant_record(line)
            if record.type == 'Fusion':
                key = record.id
            else:
                key = record.transcript_id

        if cur_key != key:
            if not cur_key is None:
                yield pointer
            cur_key = key
            pointer = GVFPointer(key=cur_key, start=line_start, end=line_end)
        else:
            pointer.end = line_end

    if pointer is not None:
        yield pointer

class GVFIndex():
    """ This class holds the index of a GVF file. """
    def __init__(self, handle:IO, pointers:Dict[str, GVFPointer]=None,
            is_circ_rna:bool=False):
        """ Constructor """
        self.handle = handle
        self.pointers = pointers or {}
        self.is_circ_rna = is_circ_rna

    def generate_indices(self):
        """ Generate indices from the GVF file that the file handle is pointing
        to. """
        for pointer in iterate_pointer(self.handle, self.is_circ_rna):
            self.pointers[pointer.key] = pointer

    def write(self, handle:IO):
        """ Write indices to file """
        for key, pointer in self.pointers.items():
            line = f"{key}\t{int(pointer.start)}\t{len(pointer)}\n"
            handle.write(line)

    def iterate_records(self, key:str) -> Iterable[Union[VariantRecord, CircRNAModel]]:
        """ Load and iterate variant records """
        pointer = self.pointers[key]
        self.handle.seek(pointer.start)
        for line in self.handle.read(len(pointer)):
            if self.is_circ_rna:
                yield circ.io.line_to_circ_model(line)
            else:
                yield io.line_to_variant_record(line)
