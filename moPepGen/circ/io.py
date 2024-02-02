""" IO """
from __future__ import annotations
from typing import Iterable, List, IO, TYPE_CHECKING
from moPepGen import GVF_HEADER
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from .CircRNA import CircRNAModel


if TYPE_CHECKING:
    from moPepGen.seqvar import GVFMetadata

def parse(handle:IO) -> Iterable[CircRNAModel]:
    """ Parse a circRNA TSV file and returns an iterable of CircRNAModel
    object """
    for line in handle:
        if line.startswith('#'):
            continue
        yield line_to_circ_model(line)

def line_to_circ_model(line) -> CircRNAModel:
    """ Convert GVF line to CircRNAModel """
    fields = line.rstrip().split('\t')

    gene_id = fields[0]
    start = int(fields[1])
    circ_id = fields[2]

    attrs = {}
    for attr in fields[7].split(';'):
        key,val = attr.split('=')
        if key in ['OFFSET', 'LENGTH']:
            val = [int(x) for x in val.split(',')]
        elif key == 'INTRON':
            val = [] if val == '' else [int(x) for x in val.split(',')]
        attrs[key] = val

    offsets = attrs['OFFSET']
    lengths = attrs['LENGTH']
    introns = attrs['INTRON']
    tx_id = attrs['TRANSCRIPT_ID']
    gene_name = attrs['GENE_SYMBOL']
    genomic_location = attrs.get('GENOMIC_LOCATION', '')

    fragments:List[SeqFeature] = []
    for j, (position, length) in enumerate(zip(offsets, lengths)):
        start_j = start + position
        end_j = start_j + length
        location = FeatureLocation(seqname=gene_id, start=start_j, end=end_j)
        frag = SeqFeature(
            chrom=gene_id, location=location, attributes={},
            type='intron' if j+1 in introns else 'exon'
        )
        fragments.append(frag)

    return CircRNAModel(tx_id, fragments, introns, circ_id, gene_id, gene_name,
            genomic_location)

def write(records:Iterable[CircRNAModel], metadata:GVFMetadata, handle:IO):
    """ Write circRNA records to file. """
    metadata.add_info('circRNA')
    for line in metadata.to_strings():
        handle.write(line + '\n')
    handle.write('#' + '\t'.join(GVF_HEADER) + '\n')
    for record in records:
        line = record.to_string()
        handle.write(line + '\n')
