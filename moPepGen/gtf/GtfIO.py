""" Module for GTF IO """
from __future__ import annotations
from typing import IO, Union, Iterable, TYPE_CHECKING
from Bio.SeqIO.Interfaces import SequenceIterator
from moPepGen.SeqFeature import FeatureLocation
from .GTFSeqFeature import GTFSeqFeature


if TYPE_CHECKING:
    from moPepGen.gtf import GenomicAnnotation

class GtfIterator(SequenceIterator):
    """ GTF Iterator """
    def __init__(self, source:Union[IO, str], mode='t'):
        """ Constructor """
        super().__init__(source=source, mode=mode, fmt='GTF')

    def parse(self, handle:IO[str]) -> Iterable[GTFSeqFeature]:
        """ parse
        """
        records = self.iterate(handle)
        return records

    @staticmethod
    def iterate(handle:IO[str]) -> Iterable[GTFSeqFeature]:
        """ Iterate through a GTF file and yield a record each time. """
        for line in handle:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            record = line_to_seq_feature(line)

            yield record

def line_to_seq_feature(line:str) -> GTFSeqFeature:
    """ Conver line to SeqFeature """
    line = line.rstrip()
    fields = line.split('\t')

    try:
        strand={'+':1, '-':-1, '?':0}[fields[6]]
    except KeyError:
        strand=None

    attributes = {}
    attributes_to_keep = ['gene_id', 'transcript_id', 'protein_id',
        'gene_name', 'gene_type', 'gene_biotype', 'tag', 'is_protein_coding']
    attribute_list = [field.strip().split(' ', 1) for field in \
        fields[8].rstrip(';').split(';')]

    for key,val in attribute_list:
        if key not in attributes_to_keep:
            continue
        val = val.strip('"')
        if key == 'tag':
            if key not in attributes:
                attributes[key] = []
            attributes[key].append(val)
        else:
            attributes[key] = val

    location = FeatureLocation(
        seqname=fields[0],
        start=int(fields[3])-1,
        end=int(fields[4]),
        strand=strand,
    )

    frame = None if fields[7] == '.' else int(fields[7])

    return GTFSeqFeature(
        chrom=fields[0],
        attributes=attributes,
        location=location,
        type=fields[2],
        frame=frame
    )

def parse(handle:Union[IO[str], str]) -> GtfIterator:
    """ Parser for GTF files.
    """
    return GtfIterator(handle)

def to_gtf_record(record:GTFSeqFeature, is_protein_coding:bool=None) -> str:
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

    if is_protein_coding is not None:
        is_protein_coding = 'true' if is_protein_coding is True else 'false'
        attrs += f" is_protein_coding {is_protein_coding};"

    frame = '.' if record.frame is None else str(record.frame)
    record_data = [
        record.chrom, '.', record.type, str(int(record.location.start)+1),
        str(int(record.location.end)), '.', strand, frame, attrs
    ]
    return '\t'.join(record_data)

def write(handle:IO, anno:GenomicAnnotation) -> None:
    """ Write an GenomicAnnotation as a GTF file. """
    for gene_model in anno.genes.values():
        handle.write(to_gtf_record(gene_model) + '\n')
        for tx_id in gene_model.transcripts:
            tx_model = anno.transcripts[tx_id]
            record = to_gtf_record(tx_model.transcript, tx_model.is_protein_coding)
            handle.write(record + '\n')
            records = tx_model.cds + tx_model.exon
            records.sort()
            records.extend(tx_model.utr)
            records = tx_model.selenocysteine + records
            for record in records:
                handle.write(to_gtf_record(record) + '\n')
