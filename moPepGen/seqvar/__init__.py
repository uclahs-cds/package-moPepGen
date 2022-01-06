""" Module for sequence variant """
from moPepGen.seqvar.VariantRecord import VariantRecord, \
    SINGLE_NUCLEOTIDE_SUBSTITUTION
from moPepGen.seqvar.VariantRecordWithCoordinate import \
    VariantRecordWithCoordinate
from moPepGen.seqvar import io
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from moPepGen.seqvar.VariantRecordPool import VariantRecordPool
from moPepGen.seqvar.VariantRecordPoolInDisk import VariantRecordPoolInDisk, \
    TranscriptionalVariantSeries, CircRNAModelSeries
