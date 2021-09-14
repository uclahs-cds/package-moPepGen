""" Module for errors """

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
