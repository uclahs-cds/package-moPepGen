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
