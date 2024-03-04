""" Varaint Peptide Table """
from __future__ import annotations
from typing import TYPE_CHECKING, Dict, Set, IO, List, Tuple
from Bio import SeqUtils
from Bio.SeqIO import FastaIO
from moPepGen import VARIANT_PEPTIDE_SOURCE_DELIMITER, aa


if TYPE_CHECKING:
    from pathlib import Path
    from Bio.Seq import Seq
    from moPepGen.svgraph.VariantPeptideDict import AnnotatedPeptideLabel
    from moPepGen.params import CleavageParams

VARIANT_PEPTIDE_TABLE_HEADERS = [
    '#sequence', 'header', 'start', 'end', 'feature_type', 'feature_id',
    'ref_start', 'ref_end', 'start_offset', 'end_offset', 'variant'
]

class VariantPeptideTable:
    """ Variant Peptide Segment Annotation """
    def __init__(self, handle:IO, index:Dict[Seq, List[Tuple[int,int]]]=None):
        self.handle = handle
        self.index = index or {}
        self.header_delimeter = VARIANT_PEPTIDE_SOURCE_DELIMITER

    def write_header(self):
        """ Write header """
        self.handle.write('\t'.join(VARIANT_PEPTIDE_TABLE_HEADERS) + '\n')

    def is_valid(self, seq:Seq, canonical_peptides:Set[str],
            cleavage_params:CleavageParams=None):
        """ Check if the peptide is valid """
        min_mw = cleavage_params.min_mw
        min_length = cleavage_params.min_length
        max_length = cleavage_params.max_length
        if SeqUtils.molecular_weight(seq, 'protein') < min_mw:
            return False
        if len(seq) < min_length or len(seq) > max_length:
            return False
        if str(seq) in canonical_peptides:
            return False
        return True

    def add_peptide(self, seq:Seq, peptide_anno:AnnotatedPeptideLabel):
        """ Write """
        start = self.handle.tell()
        for line in peptide_anno.to_lines():
            self.handle.write(f"{str(seq)}\t{line}\n")
        end = self.handle.tell()
        cur = (start, end)
        if seq in self.index:
            self.index[seq].append(cur)
        else:
            self.index[seq] = [cur]

    def load_pepitde(self, seq:Seq):
        """ Load peptide from table """
        labels = set()
        for start, end in self.index[seq]:
            self.handle.seek(start, 0)
            buffer:str = self.handle.read(end - start)
            for line in buffer.rstrip().split('\n'):
                fields = line.split('\t')
                if seq != fields[0]:
                    raise ValueError(
                        f"Peptide ({seq}) do not match with table record ({fields[0]})."
                    )
                labels.add(fields[1])
        label = self.header_delimeter.join(labels)

        return aa.AminoAcidSeqRecord(
            seq=seq,
            description=label,
            name=label
        )

    def write_fasta(self, path:Path):
        """ write fasta """
        with open(path, 'wt') as handle:
            record2title = lambda x: x.description
            writer = FastaIO.FastaWriter(handle, record2title=record2title)
            for seq in self.index:
                peptide = self.load_pepitde(seq)
                writer.write_record(peptide)
