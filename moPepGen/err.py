""" Module for errors """
from moPepGen import logger

class VariantSourceNotFoundError(Exception):
    """ Error to be raised when the variant source of a peptide is not found """
    def __init__(self, transcript_id:str=None, variant:str=None):
        """ constructor """
        message = 'Variant source not found '
        if transcript_id:
            message += f"transcript [{transcript_id}] "
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

class ReferenceSeqnameNotFoundError(Exception):
    """ Error to be raised when the seqname is not found from the reference
    genome """
    raised = False
    def __init__(self, seqname:str):
        """ constructor """
        msg = f"This seqname is not found from the reference genome: {seqname}"
        super().__init__(msg)

    @classmethod
    def mute(cls):
        """ Mute it """
        cls.raised = True

def warning(msg:str) -> None:
    """ print a warning message """
    logger(f"[ !!! moPepGen WARNING !!! ] {msg}")
