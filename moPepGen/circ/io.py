""" IO """
from typing import Iterable, List, IO
from moPepGen import GVF_HEADER
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.seqvar import GVFMetadata
from .CircRNA import CircRNAModel


def parse(handle:IO) -> Iterable[CircRNAModel]:
    """ Parse a circRNA TSV file and returns an iterable of CircRNAModel
    object """
    for i, line in enumerate(handle):
        if line.startswith('#'):
            continue

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
        tx_id = attrs['TRANSCRIPT']
        gene_name = attrs['GENE_SYMBOL']

        if len(offsets) != len(lengths):
            raise ValueError(
                f'Number of postions and lengths not match in line {i}'
            )

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

        yield CircRNAModel(tx_id, fragments, introns, circ_id,
            gene_id, gene_name)


def write(records:Iterable[CircRNAModel], metadata:GVFMetadata, handle:IO):
    """ Write circRNA records to file. """
    metadata.add_info('circRNA')
    for line in metadata.to_strings():
        handle.write(line + '\n')
    handle.write('#' + '\t'.join(GVF_HEADER) + '\n')
    for record in records:
        line = record.to_string()
        handle.write(line + '\n')
