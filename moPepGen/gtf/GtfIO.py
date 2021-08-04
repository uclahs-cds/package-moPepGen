""" Module for GTF IO """
from typing import IO, Union, Iterable
from pathlib import Path
from Bio.SeqIO.Interfaces import SequenceIterator
from moPepGen.SeqFeature import SeqFeature, FeatureLocation
from moPepGen.gtf import GenomicAnnotation


class GtfIterator(SequenceIterator):
    """ GTF Iterator """
    def __init__(self, source:Union[IO, str], mode='t'):
        """ Constructor """
        super().__init__(source=source, mode=mode, fmt='GTF')

    def parse(self, handle:IO[str]) -> Iterable[SeqFeature]:
        """ parse
        """
        records = self.iterate(handle)
        return records

    @staticmethod
    def iterate(handle:IO[str]) -> Iterable[SeqFeature]:
        """ Iterate through a GTF file and yield a record each time. """
        for line in handle:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            fields = line.split('\t')

            try:
                strand={'+':1, '-':-1, '?':0}[fields[6]]
            except KeyError:
                strand=None

            attributes = {}
            attribute_list = [field.strip().split(' ') for field in \
                fields[8].rstrip(';').split(';')]
            for key,val in attribute_list:
                val = val.strip('"')
                if key == 'tag':
                    if key not in attributes:
                        attributes[key] = []
                    attributes[key].append(val)
                else:
                    attributes[key] = val

            location=FeatureLocation(
                seqname=fields[0],
                start=int(fields[3])-1,
                end=int(fields[4])
            )

            record = SeqFeature(
                chrom=fields[0],
                attributes=attributes,
                location=location,
                type=fields[2],
                strand=strand
            )
            # record.id = record.attributes['transcript_id']
            yield record

def parse(handle:Union[IO[str], str]) -> GtfIterator:
    """ Parser for GTF files.
    """
    return GtfIterator(handle)

def to_gtf_record(record:SeqFeature) -> str:
    """ Convert a SeqFeature object to a GTF record """
    if record.strand == 1:
        strand = '+'
    elif record.strand == -1:
        strand = '-'
    else:
        strand = '.'

    attrs = ""
    for key, val in record.attributes.items():
        if isinstance(val, list):
            for vali in val:
                attrs += f" {key} {vali};"
        else:
            attrs += f" {key} {val};"

    record_data = [
        record.chrom, '.', record.type, str(int(record.location.start)),
        str(int(record.location.end)), '.', strand, '0', attrs
    ]
    return '\t'.join(record_data)

def write(path:Path, anno:GenomicAnnotation) -> None:
    """ Write an GenomicAnnotation as a GTF file. """
    with open(path, 'wt') as handle:
        for gene_model in anno.genes.values():
            handle.write(to_gtf_record(gene_model) + '\n')
            for tx_id in gene_model.transcripts:
                tx_model = anno.transcripts[tx_id]
                handle.write(to_gtf_record(tx_model.transcript) + '\n')
                records = tx_model.cds + tx_model.exon
                records.sort()
                for record in records:
                    handle.write(to_gtf_record(record) + '\n')
