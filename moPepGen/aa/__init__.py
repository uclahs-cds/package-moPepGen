""" Module for Amino Acid """
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord,\
    AminoAcidSeqRecordWithCoordinates
from moPepGen.aa.AminoAcidSeqDict import AminoAcidSeqDict
from moPepGen.aa.VariantPeptidePool import VariantPeptidePool
from moPepGen.aa.PeptidePoolSplitter import PeptidePoolSplitter, \
    VariantSourceSet, VariantPeptideInfo
from moPepGen.aa.VariantPeptideIdentifier import create_variant_peptide_id, \
    parse_variant_peptide_id
