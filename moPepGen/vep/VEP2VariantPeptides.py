""" The module defines the main logic that calls variant peptides from the VEP
output.
"""
from __future__ import annotations
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from typing import List, Set
from moPepGen import get_equivalent
from moPepGen.dna import DNASeqDict
from moPepGen.aa import AminoAcidSeqDict
from moPepGen.gtf import TranscriptGTFDict
from moPepGen.vep.TranscriptVariantDict import TranscriptVariantDict, \
    VEPVariantRecords
from moPepGen.vep.VEPVariantRecord import VEPVariantRecord
from moPepGen.vep.TranscriptVariantGraph import TranscriptVariantGraph
from moPepGen.vep import VepIO
from moPepGen.gtf.TranscriptGTFDict import TranscriptAnnotationModel


class VEP2VariantPeptides():
    """ Defines the variant peptide caller from VEP output.

    Attributes:
        vep (TranscriptVariantDict): A dict-like object that contains the
            variation effects predicted by VEP. Keys are transcript IDs and
            values are list of VEPRecord.
        annotation (TranscriptGTFDict): A dict-like object that contains the
            annotation parsed from a GTF file. Keys are transcript IDs and
            values are liks of GTFRecord.
        genome (Dict[str, [SeqRecord]])
        proteome (Dict[str, [SeqRecord]])
        canonical_peptides (set(SeqRecord))
    """
    def __init__(self,
            vep: TranscriptVariantDict,
            annotation: TranscriptGTFDict,
            genome: DNASeqDict,
            proteome: AminoAcidSeqDict,
            canonical_peptides: Set[Seq]=None):
        """ Create a VEP2VariantPeptides object.

        Args:
            vep (TranscriptVariantDict): A dict-like object that contains the
                variation effects predicted by VEP. Keys are transcript IDs
                and values are list of VEPRecord.
            annotation (TranscriptGTFDict): A dict-like object that contains
                the annotation parsed from a GTF file. Keys are transcript IDs
                and values are liks of GTFRecord.
            genome (Dict[str, [SeqRecord]])
            proteome (Dict[str, [SeqRecord]])
            canonical_peptides (set(SeqRecord))
        """
        self.vep = vep if vep is not None else TranscriptVariantDict()
        self.annotation = annotation if annotation is not None \
            else TranscriptGTFDict()
        self.genome = genome if genome is not None else DNASeqDict()
        self.proteome = proteome if proteome is not None else \
            AminoAcidSeqDict()
        self.canonical_peptides = canonical_peptides
    
    def dump_data_files(self, vep_path:List[str], gtf_path:str,
            genome_path:str, proteome_path:str)->None:
        """ Load VEP, GTF, and FASTA files from disk.
        """
        self.annotation.dump_gtf(gtf_path)
        self.genome.dump_fasta(genome_path)
        self.proteome.dump_fasta(proteome_path)
        for path in vep_path:
            self.dump_vep(path)
        for transcript in self.vep.values():
            transcript.sort()

    def dump_vep(self, path:str)->None:
        """ Load a VEP file and only keep the variant information, including
        location, ref and alt of the transcript sequence.
        
        Args:
            path (str): Path to the VEP output file.
        """
        for record in VepIO.parse(path):
            transcript_id = record.feature
            if transcript_id not in self.vep.keys():
                self.vep[transcript_id] = VEPVariantRecords()

            chrom_seqname = record.location.split(':')[0]   
            transcript_seq = self.annotation[record.feature]\
                .get_transcript_sequence(self.genome[chrom_seqname])
            
            variant_record = VEPVariantRecord.from_vep_record(
                vep=record,
                seq=transcript_seq,
                transcript_id=transcript_id
            )

            self.vep[transcript_id].append(variant_record)

    @staticmethod
    def set_up(vep_path:List[str], gtf_path:str, genome_path:str,
            proteome_path:str)->VEP2VariantPeptides:
        """ This should be the main starting point of running the
        VEP2VariantPeptides module.

        Args:
            vep_path (str): Path to the VEP file.
            gtf_path (str): Path to the annotation GTF.
            genome_path (str): Path to the genome assembly FASTA.
            proteome_path (str): path to the translated FASTA of the coding
                transcripts of the genome assembly.
        
        Example:
            adapter = VEP2VariantPaptides.set_up(
                vep_path='path/to/vep.txt',
                gtf_path='path/to/annotation.gtf',
                genome_path='path/to/genome_assembly.fasta',
                proteome_path='path/to/proteins.fasta'
            )
        """
        adapter = VEP2VariantPeptides(
            vep=None,
            annotation=None,
            genome=None,
            proteome=None
        )
        adapter.dump_data_files(
            vep_path=vep_path,
            gtf_path=gtf_path,
            genome_path=genome_path,
            proteome_path=proteome_path
        )
        return adapter

    def call_variant_peptides(self, rule:str, exception:str=None,
            miscleavage:int=2, min_mw:float=500.):
        """
        Args:
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
            start (int): Index to start searching.
            min_mw (float): Minimal molecular weight of the peptides to report.
                Defaults to 500.
        """
        carnonical_peptides = self.proteome.create_unique_peptide_pool(
            rule=rule, exception=exception,
            miscleavage=miscleavage, min_mw=min_mw
        )
        variant_peptides = set()
        variants: List[VEPVariantRecord]
        for transcript_id, variants in self.vep.items():
            anno:TranscriptAnnotationModel = self.annotation[transcript_id]
            chrom = anno.transcript.location.seqname
            transcript_seq = anno.get_transcript_sequence(self.genome[chrom])
            graph = TranscriptVariantGraph(
                seq=transcript_seq, transcript_id=transcript_id)
            graph.create_variant_graph(variants=variants)
            peptides = graph.walk_and_splice('trypsin')
            peptide:AminoAcidSeqRecord
            for peptide in peptides:
                if str(peptide.seq) in carnonical_peptides:
                    continue
                same_peptide = get_equivalent(variant_peptides, peptide)
                if same_peptide:
                    same_peptide.id += '||' + peptide.id
                    same_peptide.name = same_peptide.id
                    same_peptide.description = same_peptide.id
                    continue
                variant_peptides.add(peptide)
        return variant_peptides


if __name__ == '__main__':
    # from moPepGen.vep.VEP2VariantPeptides import VEP2VariantPeptides, TranscriptVariantGraph
    file_dir = 'test/files/downsampled_set'
    adapter = VEP2VariantPeptides.set_up(
        vep_path=[f'{file_dir}/CPCG0100_gencode_v34_snp_chr22.tsv',
            f'{file_dir}/CPCG0100_gencode_v34_indel_chr22.tsv'],
        gtf_path=f'{file_dir}/gencode_v34_chr22.gtf',
        genome_path=f'{file_dir}/gencode_v34_genome_chr22.fasta',
        proteome_path=f'{file_dir}/gencode_v34_translations_chr22.fasta'
    )
    vep = adapter.vep
    gtf = adapter.annotation
    genome = adapter.genome
    # ENST00000400588.5: stop codon gained mutation
    # ENST00000643316.1: edge issue
    # ENST00000651146.1: *
    transcript_id = 'ENST00000651146.1'
    graph = TranscriptVariantGraph(gtf[transcript_id]\
        .get_transcript_sequence(genome['chr22']), transcript_id)
    graph.create_variant_graph(vep[transcript_id])
    graph.walk_and_splice('trypsin')
    # peptides = adapter.call_variant_peptides(rule='trypsin')
