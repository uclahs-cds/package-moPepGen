""""""
from __future__ import annotations
from typing import List
from moPepGen import __version__, seqvar


ALT_DESCRIPTION = {
    'DEL': 'Deletion relative to the reference',
    'INS': 'Insertion of novel sequence relative to the reference',
    'DUP': 'Region of elevated copy number relative to the reference',
    'INV': 'Inversion of reference sequence',
    'DUP:TANDEM': 'Tandem duplication',
    'DEL:ME': 'Deletion of mobile element relative to the reference',
    'INS:ME': 'Insertion of a mobile element relative to the reference',
    'FUSION': 'Gene fusion'
}

class TVFMetadata():
    """ A TVFMetadata object contains the metadata section to be written into
    a TVF file.

    Attributes:
        parser (str): The moPepGen parser
        reference_index (str): Path to the reference index directory.
        genome_fasta (str): Path to the genome fasta file.
        annotation_gtf (str): Path to the annotation GTF file.
        alt (dict): Alternative allele fields.
        info (dict): Information fields.
    """
    def __init__(self, parser:str, reference_index:str=None,
            genome_fasta:str=None, annotation_gtf:str=None):
        """ Construct a TVFMetadata object. """
        self.parser = parser
        self.reference_index = reference_index
        self.genome_fasta = genome_fasta
        self.annotation_gtf = annotation_gtf
        self.alt = {}
        self.info = {}

    def add_info(self, variant_type:str) -> None:
        """ Add a INFO field to the metadata. """
        self.add_alt(variant_type)
        if variant_type in seqvar.SINGLE_NUCLEOTIDE_SUBSTITUTION:
            self.add_info_single_nucleotide()
        elif variant_type == 'Fusion':
            self.add_info_fusion()
        else:
            raise ValueError(f'Unknown variant type: {variant_type}')

    def add_alt(self, variant_type:str) -> None:
        """ Add a ALT field to the metadata. """
        if variant_type in ALT_DESCRIPTION and variant_type not in self.alt:
            self.alt[variant_type] = ALT_DESCRIPTION[variant_type]

    def add_info_single_nucleotide(self):
        """ Add a INFO field of a single nucleotide substitution. """
        if 'GENE_ID' not in self.info:
            self.info['GENE_ID'] = {
                'Number': 1,
                'Type': 'String',
                'Description': 'Transcript\'s Gene ID'
            }

    def add_info_fusion(self):
        """ Add a INFO fild of a gene fusion event. """
        if 'GENE_ID' not in self.info:
            self.info['GENE_ID'] = {
                'Number': 1,
                'Type': 'String',
                'Description': 'Transcript\'s Gene ID'
            }

        if 'DONOR_GENE_ID' not in self.info:
            self.info['DONOR_GENE_ID'] = {
                'Number': 1,
                'Type': 'String',
                'Description': 'Donor Transcript\'s Gene ID'
            }

        if 'DONOR_TRANSCRIPT_ID' not in self.info:
            self.info['DONOR_TRANSCRIPT_ID'] = {
                'Number': 1,
                'Type': 'String',
                'Description': 'Donor Transcript\'s Transcript ID'
            }

        if 'DONOR_POS' not in self.info:
            self.info['DONOR_POS'] = {
                'Number': 1,
                'Type': 'Integer',
                'Description': 'Position of the break point of the donor transcript'
            }

    def to_strings(self) -> List[str]:
        """ Convert metadata to a list of strings to be written to a TVF file. """
        info_lines = []
        for key,val in self.info.items():
            number = val['Number']
            _type = val['Type']
            desc = val['Description']
            line = f'##INFO=<ID={key},Number={number},Type={_type},Description="{desc}">'
            info_lines.append(line)
        ref_index = self.reference_index if self.reference_index else ''
        genome_fasta = self.genome_fasta if self.genome_fasta else ''
        annotation_gtf = self.annotation_gtf if self.annotation_gtf else ''
        return [
            '##fileformat=VCFv4.2',
            f'##mopepgen_version={__version__}',
            f'##parser={self.parser}',
            f'##reference_index={ref_index}',
            f'##genome_fasta={genome_fasta}',
            f'##annotation_gtf={annotation_gtf}',
            *[f'##ALT=<ID={key},Description="{val}">' for key, val in self.alt],
            *info_lines,
        ]
