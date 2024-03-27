""" moPepGen constants """
from enum import Enum


class VariantPrefix(Enum):
    """ Varinat Prefix"""
    SNV    = 'SNV'
    INDEL  = 'INDEL'
    MNV    = 'MNV'
    RES    = 'RES'
    SE     = 'SE'
    RI     = 'RI'
    A3SS   = 'A3SS'
    A5SS   = 'A5SS'
    MXE    = 'MXE'
    W2F    = 'W2F'
    SECT   = 'SECT'
    CIRC   = 'CIRC'
    CI     = 'CI'
    FUSION = 'FUSION'

    def __str__(self):
        """ str """
        return self.value

    @classmethod
    def ctbv(cls):
        """ Conserved Transcript Backbone Variants """
        return [
            cls.SNV,
            cls.INDEL,
            cls.MNV,
            cls.RES,
            cls.SE,
            cls.RI,
            cls.A3SS,
            cls.A5SS,
            cls.MXE,
            cls.W2F,
            cls.SECT
        ]

    @classmethod
    def ntbv(cls):
        """ Novel Transcript Backbone Variants """
        return [
            cls.CIRC,
            cls.CI,
            cls.FUSION
        ]

    @classmethod
    def alt_translation(cls):
        """ Alternative translation variants """
        return [
            cls.W2F,
            cls.SECT
        ]

# Variant related constants
SINGLE_NUCLEOTIDE_SUBSTITUTION = ['SNV', 'SNP', 'INDEL', 'MNV', 'RNAEditingSite']
ATTRS_POSITION = ['START', 'DONOR_START', 'ACCEPTER_START', 'ACCEPTER_POSITION']
ALTERNATIVE_SPLICING_TYPES = ['Insertion', 'Deletion', 'Substitution']
RMATS_TYPES = ['SE', 'RI', 'A3SS', 'A5SS', 'MXE']
CODON_REASSIGNMENTS_TYPES = ['W2F']
SEC_TERMINATION_TYPE = 'SECT'

# Variant sources
SOURCE_NOVEL_ORF = 'NovelORF'
SOURCE_CODON_REASSIGNMENT = 'CodonReassign'
SOURCE_SEC_TERMINATION = 'SECT'

# CLI
PROG_NAME = 'moPepGen'
