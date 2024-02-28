""" Module for sequence variant """
from moPepGen.seqvar.VariantRecord import VariantRecord
from moPepGen.seqvar.VariantRecordWithCoordinate import \
    VariantRecordWithCoordinate
from moPepGen.seqvar import io
from moPepGen.seqvar.GVFMetadata import GVFMetadata
from moPepGen.seqvar.VariantRecordPool import VariantRecordPool
from moPepGen.seqvar.VariantRecordPoolOnDisk import VariantRecordPoolOnDisk, \
    TranscriptionalVariantSeries, VariantRecordPoolOnDiskOpener
from moPepGen.seqvar.VariantRecord import create_variant_sect, create_variant_w2f, \
    create_mnv_from_adjacent, find_mnvs_from_adjacent_variants
from moPepGen.seqvar.SplicingJunction import SpliceJunction, \
    SpliceJunctionTranscriptAlignment
