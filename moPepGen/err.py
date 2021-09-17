""" Module for errors """
from __future__ import annotations
from typing import TYPE_CHECKING
from moPepGen import logger


if TYPE_CHECKING:
    from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature

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

def warning(msg:str) -> None:
    """ print a warning message """
    logger(f"[ !!! moPepGen WARNING !!! ] {msg}")
