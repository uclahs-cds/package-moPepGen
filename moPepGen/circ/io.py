""" IO """
from typing import Iterable, List
from pathlib import Path
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.seqvar import GVFMetadata
from .CircRNA import CircRNAModel


def parse(path:Path) -> Iterable[CircRNAModel]:
    """ Parse a circRNA TSV file and returns an iterable of CircRNAModel
    object """
    with open(path, 'r') as handle:
        i = 0
        line = next(handle, None)

        while line:
            i += 1
            if line:
                if line.startswith('#'):
                    line = next(handle, None)
                    continue
                fields = line.rstrip().split('\t')

            gene_id = fields[0]
            start = int(fields[1])

            positions = fields[2].split(',')
            lengths = fields[3].split(',')

            if fields[4] == '.':
                introns = []
            else:
                introns = [int(x) for x in fields[4].split(',')]

            circ_id = fields[5]
            transcript_ids = fields[6].split(',')
            gene_name = fields[7]

            if len(positions) != len(lengths):
                raise ValueError(
                    f'Number of postions and lengths not match in line {i}'
                )

            fragments:List[SeqFeature] = []
            for j, (position, length) in enumerate(zip(positions, lengths)):
                position = int(position)
                length = int(length)
                start_j = start + position
                end_j = start_j + length
                location = FeatureLocation(seqname=gene_id, start=start_j, end=end_j)
                frag = SeqFeature(
                    chrom=gene_id, location=location, attributes={},
                    type='intron' if j+1 in introns else 'exon'
                )
                fragments.append(frag)

            yield CircRNAModel(gene_id, fragments, introns, circ_id,
                transcript_ids, gene_name)
            line = next(handle, None)


def write(records:Iterable[CircRNAModel], metadata:GVFMetadata, path:Path):
    """ Write circRNA records to file. """
    headers = ['gene_id','start','offset','length','circ_id','gene_name']
    with open(path, 'wt') as handle:
        for line in metadata.to_strings():
            handle.write(line + '\n')

        handle.write('#' + ','.join(headers) + '\n')

        for record in records:
            line = record.to_string() + '\n'
            handle.write(line)
