""" Module for errors """
from __future__ import annotations
from typing import TYPE_CHECKING
from moPepGen import get_logger

if TYPE_CHECKING:
    from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
    from moPepGen.version import MetaVersion
    from moPepGen.params import CleavageParams

class VariantSourceNotFoundError(Exception):
    """ Error to be raised when the variant source of a peptide is not found """
    def __init__(self, gene_id:str=None, variant:str=None):
        """ constructor """
        message = 'Variant source not found '
        if gene_id:
            message += f"transcript [{gene_id}] "
        if variant:
            message += f"variant [{variant}]"
        message += '. Please verify all GVF files are imported'
        super().__init__(message)

class TranscriptionStopSiteMutationError(Exception):
    """ Error to to used when there is a variant altering the transcriptional
    stop site. """
    def __init__(self, transcript_id:str=None, variant:str=None):
        """ constructor """
        message = 'The variant alters the transcription stop site. '
        if variant:
            message += f"variant [{variant}]"
        if transcript_id:
            message += f"transcript [{transcript_id}] "
        super().__init__(message)

class TranscriptionStartSiteMutationError(Exception):
    """ Error to to used when there is a variant altering the transcriptional
    start site. """
    def __init__(self, transcript_id:str=None, variant:str=None):
        """ constructor """
        message = 'The variant alters the transcription start site. '
        if variant:
            message += f"variant [{variant}]"
        if transcript_id:
            message += f"transcript [{transcript_id}] "
        super().__init__(message)

class FailedToFindVariantBubbleError(Exception):
    """ Error to handle when variant bubble was not found correctly. """

class MuteableErrorMixin():
    """ Mixin for muteable errors """
    @classmethod
    def mute(cls):
        """ Mute it """
        cls.raised = True

class ReferenceSeqnameNotFoundError(MuteableErrorMixin, Exception):
    """ Error to be raised when the seqname is not found from the reference
    genome """
    raised = False
    def __init__(self, seqname:str):
        """ constructor """
        msg = f"This seqname is not found from the reference genome: {seqname}"
        super().__init__(msg)

class ExonNotFoundError(Exception):
    """ Error to be raised when an exon is not found """
    def __init__(self, gene_id:str, feature:GTFSeqFeature):
        """ constructor """
        msg = f"The feature {str(feature.location)} is not a valid exon from"+\
        f" gene {gene_id}"
        super().__init__(msg)

class IntronNotFoundError(Exception):
    """ Error to be raised when an exon is not found """
    def __init__(self, gene_id:str, feature:GTFSeqFeature):
        """ constructor """
        msg = f"The feature {str(feature.location)} is not a valid intron" +\
            f"from gene {gene_id}"
        super().__init__(msg)

class InvalidIndexError(Exception):
    """ Error to be raised when the index version does not match with the
    current environment """
    def __init__(self, this_index:MetaVersion, other_index:MetaVersion):
        """ constructor """
        msg = "Current runtime environment or cleavage params do not match with the index." +\
            f"Version: current: {this_index}; index: {other_index}"
        super().__init__(msg)

class GeneNotFoundError(Exception):
    """ Error to be raised when a gene is not found from GTF """
    def __init__(self, gene_id:str):
        """ constructor """
        msg = f"Gene {gene_id} not found."
        super().__init__(msg)

class MNVParsingError(Exception):
    """ Error to be raised when trying to parse MNVs (multi-nucleotide variant). """
    def __init__(self):
        """ constructor """
        msg = "Trying to parse a MNV, which is currently unsupported."
        super().__init__(msg)

def warning(msg:str) -> None:
    """ print a warning message """
    get_logger().warning("[ !!! moPepGen WARNING !!! ] %s", msg)

class MoPepGenWarning():
    """ Base warning class """
    def __init__(self, msg):
        """ constructor """
        self.msg = msg
        self.warn()

    def warn(self):
        """ Print warning message """
        warning(self.msg)

class HypermutatedRegionWarning(MoPepGenWarning):
    """ Warning to be printed when hyper mutated region is detected. """
    def __init__(self, graph_id:str, max_variants_per_node:int,
            additional_variants_per_misc:int):
        """ constructor """
        self.max_variants_per_node = max_variants_per_node
        msg = f"Hypermutated region detected from graph: '{graph_id}'. The" +\
            f" argument max_variants_per_node = {max_variants_per_node} and" +\
            f" additional_variants_per_misc = {additional_variants_per_misc}" +\
            ' was used to reduce complexity.'
        super().__init__(msg)
