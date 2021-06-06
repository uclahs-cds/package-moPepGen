""" This module defines the class for a single GTF record
"""
from __future__ import annotations
from typing import Dict


class GTFRecord():
    """ A GTFRecord object holds an entry from the a GTF file.

    The schema of GTFRecord follows the definition from ENSMEBL (https://uswest
    .ensembl.org/info/website/upload/gff.html).

    Attributes:
        seqname (str): The reference sequence name of the chromosome or
            scaffold.
        source (str): Name of the program that generated this feature, or the
            data source (database or project name).
        feature (str): Feature type name, e.g. gene, transcript, exone.
        start (int): Start position* of the feature, with sequence numbering
            starting at 1.
        end (int): End position* of the feature, with sequence numbering
            starting at 1.
        score (float): A floating point value.
        strand (str): Defined as + (forward) or - (reverse).
        frame (int): One of '0', '1' or '2'. '0' indicates that the first base
            of the feature is the first base of a codon, '1' that the second
            base is the first base of a codon, and so on.
        attributes (dict[str, str]): additional information about each feature,
            e.g., gene_id, transcript_id, gene_type.
    """
    def __init__(
            self, seqname: str, source: str, feature: str, start: int,
            end: int, score: float, strand: str, frame: int,
            attributes: Dict[str, str]):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = attributes

    @property
    def biotype(self):
        """ biotype """
        return self.attributes['gene_type']

    def __eq__(self, other:GTFRecord)->bool:
        """ equal to """
        return self.start == other.start and self.end == other.end

    def __ne__(self, other:GTFRecord)->bool:
        """ not equal """
        return not self == other

    def __gt__(self, other:GTFRecord)->bool:
        """ greater than """
        if self.start > other.start:
            return True
        if self.start < other.start:
            return False
        if self.end > other.end:
            return True
        if self.end < other.end:
            return False
        return False

    def __ge__(self, other:GTFRecord)->bool:
        """ greater or equal """
        return self > other or self == other

    def __lt__(self, other:GTFRecord)->bool:
        """ less than """
        return not (self > other or self == other)

    def __le__(self, other:GTFRecord)->bool:
        """" less or equal """
        return self < other or self == other
