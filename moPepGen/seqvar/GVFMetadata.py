""" Module for GVF metadata """
from __future__ import annotations
from typing import TYPE_CHECKING
import copy
from moPepGen import __version__, constant
from .GVFMetadataInfo import GVF_METADATA_INFO, GVF_METADATA_ADDITIONAL


if TYPE_CHECKING:
    from typing import List, IO, Tuple, Set

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

class GVFMetadata():
    """ A GVFMetadata object contains the metadata section to be written into
    a GVF file.

    Attributes:
        parser (str): The moPepGen parser
        reference_index (str): Path to the reference index directory.
        genome_fasta (str): Path to the genome fasta file.
        annotation_gtf (str): Path to the annotation GTF file.
        alt (dict): Alternative allele fields.
        info (dict): Information fields.
        additional (dict): Additional metadata fields.
        version (str): Version of moPepGen.
        phase_groups (List[Set[str,str]]): Phase groups for phased variants.
        source (str): Source of the variants, e.g., 'VEP', 'CIRCexplorer'.
        chrom (str): Chromosome or gene ID.
    """
    def __init__(self, parser:str, source:str, chrom:str, reference_index:str=None,
            genome_fasta:str=None, annotation_gtf:str=None, version:str=None, info=None,
            additional=None, phase_groups:List[Set[str,str]]=None):
        """ Construct a TVFMetadata object. """
        self.parser = parser
        self.source = source
        self.chrom = chrom
        self.reference_index = reference_index
        self.genome_fasta = genome_fasta
        self.annotation_gtf = annotation_gtf
        self.alt = {}
        self.info = info or GVF_METADATA_INFO['Base']
        self.added_types = []
        self.additional = additional
        self.version = version or __version__
        self.phase_groups = phase_groups or list()

    def add_info(self, variant_type:str, is_phased:bool=False) -> None:
        """ Add a INFO field to the metadata. """
        if is_phased:
            self.info.update({
                'PHASE_SETS': GVF_METADATA_ADDITIONAL['Base']['PHASE_SETS']
            })

        self.add_alt(variant_type)
        if variant_type in self.added_types:
            return
        if variant_type in constant.SINGLE_NUCLEOTIDE_SUBSTITUTION:
            return
        if variant_type == 'Fusion':
            self.info.update(GVF_METADATA_INFO['Fusion'])
            self.added_types.append('Fusion')
        elif variant_type == 'Insertion':
            self.info.update(GVF_METADATA_INFO['Insertion'])
            self.added_types.append('Insertion')
        elif variant_type == 'Deletion':
            self.info.update(GVF_METADATA_INFO['Deletion'])
            self.added_types.append('Deletion')
        elif variant_type == 'Substitution':
            self.info.update(GVF_METADATA_INFO['Substitution'])
            self.added_types.append('Substitution')
        elif variant_type == 'circRNA':
            self.info.update(GVF_METADATA_INFO['circRNA'])
            self.added_types.append('circRNA')
            self.additional = GVF_METADATA_ADDITIONAL['circRNA']
        else:
            raise ValueError(f'Unknown variant type: {variant_type}')

    def add_alt(self, variant_type:str) -> None:
        """ Add a ALT field to the metadata. """
        if variant_type in ALT_DESCRIPTION and variant_type not in self.alt:
            self.alt[variant_type] = ALT_DESCRIPTION[variant_type]

    def to_strings(self) -> List[str]:
        """ Convert metadata to a list of strings to be written to a GVF file. """
        info_lines = []
        for key,val in self.info.items():
            number = val['Number']
            _type = val['Type']
            desc = val['Description']
            line = f'##INFO=<ID={key},Number={number},Type={_type},Description="{desc}">'
            info_lines.append(line)

        additional = []
        if self.additional:
            for key,val in self.additional.items():
                values = []
                for k,v in val.items():
                    values.append(f"{k}={v}")
                additional.append(f"##{key}=<{','.join(values)}>")

        ref_index = self.reference_index if self.reference_index else ''
        genome_fasta = self.genome_fasta if self.genome_fasta else ''
        annotation_gtf = self.annotation_gtf if self.annotation_gtf else ''
        res = [
            '##fileformat=VCFv4.2',
            f'##mopepgen_version={self.version}',
            f'##parser={self.parser}',
            f'##reference_index={ref_index}',
            f'##genome_fasta={genome_fasta}',
            f'##annotation_gtf={annotation_gtf}',
            f'##source={self.source}'
        ]
        if self.phase_groups:
            phase_groups = sorted(self.phase_groups, key=lambda x: x[0])
            phase_groups = ','.join([f'{ps1}|{ps2}' for ps1,ps2 in phase_groups])
            res.append(f'##phase_groups={phase_groups}')
        res += [
            f'##CHROM=<Description="{self.chrom}">',
            *[f'##ALT=<ID={key},Description="{val}">' for key, val in self.alt],
            *info_lines,
            *additional
        ]
        return res

    def is_circ_rna(self) -> bool:
        """ checks if this is a circRNA """
        return self.parser in ['parseCIRCexplorer']

    @classmethod
    def parse(cls, handle:IO) -> GVFMetadata:
        """ Parse the metadata from a GVF file """
        pos = handle.tell()
        it = handle.readline()
        metadata = {}
        while it and it.startswith('##'):
            key, val = it.lstrip('##').rstrip().split('=', 1)
            if val == '':
                val = None
            if key == 'phase_groups':
                metadata[key] = [set(ps.split('|')) for ps in val.split(',')]
            elif key == 'CHROM':
                key = key.lower()
                for it in val.strip('<>').split(','):
                    k,v = it.split('=')
                    v = v.strip('"').strip("'")
                    if k == 'Description':
                        metadata[key] = v
            elif key == 'INFO':
                key = key.lower()
                if key not in metadata:
                    metadata[key] = {}
                info_key = None
                info_val = {}
                for it in val.strip('<>').split(','):
                    k,v = it.split('=')
                    v = v.strip('"')
                    if k == 'ID':
                        info_key = v
                    else:
                        info_val[k] = v
                if info_key is None:
                    raise ValueError('Could not find the key.')
                metadata[key][info_key] = info_val
            elif key in ['ALT', 'POS']:
                pass
            elif key == 'mopepgen_version':
                metadata['version'] = val
            elif key != 'fileformat':
                metadata[key] = val
            pos = handle.tell()
            it = handle.readline()
        handle.seek(pos)
        return cls(**metadata)

    def combine_phase_groups(self, phase_groups:List[Set[str]]) -> List[Set[str]]:
        """ Merge the phase groups from this metadata to the given phase groups. """
        phase_groups = copy.deepcopy(phase_groups)
        existing_phases = set().union(*phase_groups)
        incoming_phases = set()
        for phase_set in self.phase_groups:
            if any(x in existing_phases for x in phase_set):
                raise ValueError(
                    "Phase sets in the GVF metadata are not disjoint with existing"
                    " phase sets. Please check the GVF file."
                )
            if any(x in incoming_phases for x in phase_set):
                raise ValueError(
                    "Phase sets in the GVF metadata are not disjoint. "
                    "Please check the GVF file."
                )
            incoming_phases.update(phase_set)
            phase_groups.append(phase_set)
        return phase_groups
