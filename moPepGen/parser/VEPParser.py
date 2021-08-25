""" This module defines the class for a single VEP record
"""
from __future__ import annotations
from typing import List, Tuple, Iterable
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import seqvar, dna, gtf


def parse(path:str) -> Iterable[VEPRecord]:
    """ Parse a VEP output text file and return as an iterator.

    Args:
        path (str): Path to the REDItools output table.

    Return:
        A iterable of VEPRecord.
    """
    with open(path, 'r') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            fields = line.split('\t')

            consequences = fields[6].split(',')

            amino_acids = tuple(aa for aa in fields[10].split('/'))
            if len(amino_acids) == 1:
                amino_acids = (amino_acids[0], '')

            codons = tuple(codon for codon in fields[11].split('/'))
            if len(codons) == 1:
                codons = (codons[0], '')

            extra = {}
            for field in fields[13].split(';'):
                key, val = field.split('=')
                extra[key] = val

            yield VEPRecord(
                uploaded_variation=fields[0],
                location=fields[1],
                allele=fields[2],
                gene=fields[3],
                feature=fields[4],
                feature_type=fields[5],
                consequences=consequences,
                cdna_position=fields[7],
                cds_position=fields[8],
                protein_position=fields[9],
                amino_acids=amino_acids,
                codons=codons,
                existing_variation='' if fields[12] == '-' else fields[12],
                extra=extra
            )

class VEPRecord():
    """ A VEPRecord object holds the an entry from the VEP output. The VEP
    output is defined at https://uswest.ensembl.org/info/docs/tools/vep/
    vep_formats.html#output

    Attributes:
        uploaded_variation (str): as chromosome_start_alleles
        location (str): in standard coordinate format (chr:start or
            chr:start-end)
        allele (str): the variant allele used to calculate the consequence
        gene (str): Ensembl stable ID of affected gene
        feature (str): Ensembl stable ID of feature
        feature_type (str): type of feature. Currently one of Transcript,
            RegulatoryFeature, MotifFeature.
        consequence (List[str]): consequence type of this variant.
            See: https://uswest.ensembl.org/info/genome/variation/prediction/
            predicted_data.html#consequences
        cdna_position (str): relative position of base pair in cDNA sequence
        cds_position (str): relative position of base pair in coding sequence
        protein_position (str): relative position of amino acid in protein
        amino_acids (Tuple[str]): only given if the variant affects the
            protein-coding sequence
        codons (Tuple[str]) the alternative codons with the variant base in
            upper case
        existing_variation (str) known identifier of existing variant
        extra (dict): this column contains extra information.
    """
    def __init__(
            self, uploaded_variation: str, location: str, allele: str,
            gene: str, feature: str, feature_type:str,
            consequences: List[str], cdna_position: str, cds_position: str,
            protein_position: str, amino_acids: Tuple[str, str],
            codons: Tuple[str, str], existing_variation: str, extra: dict):
        """ Construct a VEPRecord object. """
        self.uploaded_variation = uploaded_variation
        self.location = location
        self.allele = allele
        self.gene = gene
        self.feature = feature
        self.feature_type=feature_type
        self.consequences = consequences
        self.cdna_position = cdna_position
        self.cds_position = cds_position
        self.protein_position = protein_position
        self.amino_acids = amino_acids
        self.codons = codons
        self.existing_variation = existing_variation
        self.extra = extra

    def __repr__(self)->str:
        """Return representation of the VEP record."""
        consequences = '|'.join(self.consequences)
        return f"< {self.feature}, {consequences}, {self.location} >"

    def convert_to_variant_record(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> seqvar.VariantRecord:
        """ Convert a VepRecord to a generic VariantRecord object.

        Args:
            seq (dna.DNASeqRecord): The DNA sequence of the transcript.
        """
        chrom_seqname, alt_position = self.location.split(':')
        if alt_position.find('-') > -1:
            alt_start_genomic, alt_end_genomic = \
                [int(x) for x in alt_position.split('-')]
        else:
            alt_start_genomic = int(alt_position)
            alt_end_genomic = alt_start_genomic
        alt_start_genomic -= 1

        gene_model = anno.genes[self.gene]
        strand = gene_model.strand
        seq = gene_model.get_gene_sequence(genome[chrom_seqname])

        alt_start = anno.coordinate_genomic_to_gene(alt_start_genomic, self.gene)
        alt_end = anno.coordinate_genomic_to_gene(alt_end_genomic - 1, self.gene)
        if strand == -1:
            alt_start, alt_end = alt_end, alt_start
        alt_end += 1

        if self.allele == '-':
            if alt_start == 0:
                ref = str(seq.seq[alt_start:alt_end])
                alt = str(genome[chrom_seqname].seq[alt_start_genomic - 1])
            else:
                alt_start -= 1
                ref = str(seq.seq[alt_start:alt_end])
                alt = str(seq.seq[alt_start])
        else:
            allele = self.allele
            if strand == -1:
                allele = str(Seq(allele).reverse_complement())
            if alt_end - alt_start == 1:
                if len(allele) > 1:
                    raise ValueError('Could not recognize the VEP record')
                ref = str(seq.seq[alt_start])
                alt = allele
            elif alt_end - alt_start == 2:
                alt_end -= 1
                ref = str(seq.seq[alt_start])
                alt = ref + allele
            else:
                raise ValueError('Could not recognize the VEP record')

        _type = 'SNV' if len(ref) == 1 and len(alt) == 1 else 'INDEL'
        _id = f'{_type}-{alt_start + 1}-{ref}-{alt}'

        attrs = {
            'TRANSCRIPT_ID': self.feature,
            'GENOMIC_POSITION': self.location,
            'GENE_SYMBOL': gene_model.attributes['gene_name']
        }

        try:
            return seqvar.VariantRecord(
                location=FeatureLocation(
                    seqname=self.gene,
                    start=alt_start,
                    end=alt_end
                ),
                ref=ref,
                alt=alt,
                _type=_type,
                _id=_id,
                attrs=attrs
            )
        except ValueError as e:
            raise ValueError(e.args[0] + f' [{self.feature}]') from e
