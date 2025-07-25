""" This module defines the class for a single VEP record
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Literal
import copy
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.err import TranscriptionStopSiteMutationError, \
    TranscriptionStartSiteMutationError
from moPepGen import seqvar, dna, gtf


if TYPE_CHECKING:
    from typing import List, Tuple, Iterable, IO, Dict, Set

def parse(handle:IO, format:str=Literal['tsv', 'vcf'], samples:List[str]=None,
        current_phase_sets:Set[Tuple[str, str]]=None) -> Iterable[VEPRecord]:
    """ Parse a VEP output file and return as an iterator.

    Args:
        handle (IO): A file-like object containing the VEP output.
        format (str): The format of the VEP output, either 'tsv' or 'vcf'.
        current_phase_sets (Set[Tuple[str, str]]): A set of phase sets that have
            been used.

    Return:
        An iterable of VEPRecord.
    """
    if format == 'tsv':
        return parse_tsv(handle)
    return parse_vcf(handle, samples=samples, current_phase_sets=current_phase_sets)

def parse_tsv(handle:IO) -> Iterable[VEPRecord]:
    """ Parse a VEP output text file and return as an iterator.

    Args:
        path (str): Path to the REDItools output table.

    Return:
        A iterable of VEPRecord.
    """
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

def parse_vcf(handle:IO, samples:List[str]=None, current_phase_sets:Set[Tuple[str, str]]=None
        ) -> Iterable[VEPRecord]:
    """ Parse a VEP output VCF file and return as an iterator.

    Args:
        handle (IO): A file-like object containing the VEP output in VCF format.
        samples (List[str]): A list of sample names from the VCF file to be parsed.
            If None, all samples will be parsed.
        current_phase_sets (Set[Tuple[str, str]]): A set of phase sets that have
            been used.

    Return:
        A iterable of VEPRecord.
    """
    if current_phase_sets is None:
        current_phase_sets = set()
    max_phase_set = max(int(ps[1][2:]) for ps in current_phase_sets) \
        if current_phase_sets else 0
    phase_sets = (f"PS{max_phase_set + 1}", f"PS{max_phase_set + 2}")
    for line in handle:
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            headers = line.rstrip().split('\t')
            sample_index:Dict[str, int] = {}
            for i, header in enumerate(headers):
                if i <= 8:
                    continue
                if not samples or header in samples:
                    sample_index[header] = i
            # Checks if all samples are present in the VEP VCF headers.
            if samples:
                for sample in samples:
                    if sample not in sample_index:
                        raise ValueError(
                            f"Sample '{sample}' not found in VEP VCF headers."
                        )
                continue
        fields = line.rstrip().split('\t')
        info = {}
        item:str
        for item in fields[7].split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info[key] = value
            else:
                info[item] = True
        if 'CSQ' not in info:
            continue
        vep_records = []
        for csq in info['CSQ'].split(','):
            attrs = csq.split('|')
            record = VEPRecord(
                uploaded_variation=fields[2],
                location=fields[0] + ':' + fields[1],
                allele=attrs[0],
                gene=attrs[4],
                feature=attrs[6],
                feature_type=attrs[5],
                consequences=attrs[1].split('&'),
                cdna_position=attrs[12],
                cds_position=attrs[13],
                protein_position=attrs[14],
                amino_acids=tuple(attrs[15].split('/')),
                codons=tuple(attrs[16].split('/')),
                existing_variation=attrs[17] if attrs[17] != '-' else '',
                extra={
                    'IMPACT': attrs[2]
                }
            )
            vep_records.append(record)

        sample_genotypes = {}
        keys = fields[8].split(':')
        for sample, i in sample_index.items():
            values = fields[i].split(':')
            genotype = dict(zip(keys, values))
            gt = genotype['GT']
            is_phased = '|' in gt
            is_detected = gt not in ['.', './.', '.|.', '0|0', '0/0']
            cur_phase_set = []
            if is_phased:
                current_phase_sets.add(phase_sets)
                if gt.startswith('1|'):
                    cur_phase_set.append(phase_sets[0])
                if gt.endswith('|1'):
                    cur_phase_set.append(phase_sets[1])
            sample_genotypes[sample] = {
                'GT': gt,
                'is_phased': is_phased,
                'is_detected': is_detected,
                'phase_sets': cur_phase_set
            }

        for record in vep_records:
            for sample, genotype in sample_genotypes.items():
                if not genotype['is_detected']:
                    # No variant for this sample, so we can skip it.
                    continue
                cur_record = copy.deepcopy(record)
                cur_record.extra.update({
                    'SAMPLE': sample,
                    'PHASE_SETS': genotype['phase_sets']
                })
                yield cur_record

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
        return f"< {self.feature}, {consequences}, {self.location}, {self.allele} >"

    def convert_to_variant_record(self, anno:gtf.GenomicAnnotation,
            genome:dna.DNASeqDict) -> seqvar.VariantRecord:
        """ Convert a VepRecord to a generic VariantRecord object.

        Args:
            seq (dna.DNASeqRecord): The DNA sequence of the transcript.
        """
        _, alt_position = self.location.split(':')
        if alt_position.find('-') > -1:
            alt_start_genomic, alt_end_genomic = \
                [int(x) for x in alt_position.split('-')]
        else:
            alt_start_genomic = int(alt_position)
            alt_end_genomic = alt_start_genomic
        alt_start_genomic -= 1

        gene_model = anno.genes[self.gene]
        chrom_seqname = gene_model.chrom
        strand = gene_model.strand
        seq = gene_model.get_gene_sequence(genome[chrom_seqname])

        alt_start = anno.coordinate_genomic_to_gene(alt_start_genomic, self.gene)
        alt_end = anno.coordinate_genomic_to_gene(alt_end_genomic - 1, self.gene)

        tx_id = self.feature
        tx_model = anno.transcripts[tx_id]

        if strand == 1:
            tx_first = tx_model.transcript.location.start
            tx_last = tx_model.transcript.location.end - 1
            tx_start_genetic = anno.coordinate_genomic_to_gene(tx_first, self.gene)
            tx_end_genetic = anno.coordinate_genomic_to_gene(tx_last, self.gene) + 1
        else:
            tx_first = tx_model.transcript.location.end - 1
            tx_last = tx_model.transcript.location.start
            tx_start_genetic = anno.coordinate_genomic_to_gene(tx_first, self.gene)
            tx_end_genetic = anno.coordinate_genomic_to_gene(tx_last, self.gene) + 1

        if strand == -1:
            alt_start, alt_end = alt_end, alt_start
        alt_end += 1

        if alt_start < tx_start_genetic or (alt_start == tx_start_genetic \
                and not tx_model.is_cds_start_nf()):
            raise TranscriptionStartSiteMutationError(tx_id)

        if alt_end > tx_end_genetic:
            raise TranscriptionStopSiteMutationError(tx_id)

        if self.allele == '-':
            if alt_start == tx_start_genetic:
                alt_end += 1
                ref = str(seq.seq[alt_start:alt_end])
                alt = str(seq.seq[alt_end-1:alt_end])
            else:
                alt_start -= 1
                ref = str(seq.seq[alt_start:alt_end])
                alt = str(seq.seq[alt_start])
        else:
            allele = self.allele
            if strand == -1:
                allele = str(Seq(allele).reverse_complement())
            if alt_end - alt_start == 1:
                if len(allele) > 1: # insertion
                    # Sometimes insertions are reported by VEP in the end-inclusion
                    # way (e.g., C -> TACC), which needs to be converted into
                    # start-inclusion (A -> ATAC) for variants on + strand genes.
                    ref = str(seq.seq[alt_start])
                    if ref == allele[-1]:
                        alt_start -= 1
                        alt_end = alt_start + 1
                        ref = str(seq.seq[alt_start])
                        alt = ref + allele[:-1]
                    elif ref == allele[0]:
                        ref = str(ref)
                        alt = allele
                    else:
                        raise ValueError(f"Don't know how to process this variant: {self}")
                else: # SNV
                    ref = str(seq.seq[alt_start])
                    alt = allele
            elif alt_end - alt_start == 2:  # insertion
                alt_end -= 1
                ref = str(seq.seq[alt_start])
                alt = ref + allele
            else:
                ref = str(seq.seq[alt_start:alt_end])
                alt = allele

        if len(ref) == len(alt) == 1:
            _type = 'SNV'
        elif len(ref) == 1 or len(ref) == 1:
            _type = 'INDEL'
        else:
            _type = 'MNV'
        _id = f'{_type}-{alt_start + 1}-{ref}-{alt}'

        attrs = {
            'TRANSCRIPT_ID': self.feature,
            'GENOMIC_POSITION': self.location,
            'GENE_SYMBOL': gene_model.gene_name
        }

        if 'PHASE_SETS' in self.extra and self.extra['PHASE_SETS']:
            attrs['PHASE_SETS'] = ','.join(self.extra['PHASE_SETS'])

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
