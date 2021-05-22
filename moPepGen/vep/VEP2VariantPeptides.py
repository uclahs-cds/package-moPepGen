""" The module defines the main logic that calls variant peptides from the VEP
output.
"""
from __future__ import annotations
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from typing import List, Set
from pathos.pools import ProcessPool
from Bio import SeqIO
from moPepGen import get_equivalent
from moPepGen import dna
from moPepGen import aa
from moPepGen import gtf
from moPepGen.vep.TranscriptVariantDict import TranscriptVariantDict, \
    VEPVariantRecords
from moPepGen.vep.VEPVariantRecord import VEPVariantRecord
from moPepGen.vep.TranscriptVariantGraph import TranscriptVariantGraph
from moPepGen.vep import VepIO


_VALID_BIOTYPES = ['protein_coding']

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
            annotation: gtf.TranscriptGTFDict,
            genome: dna.DNASeqDict,
            proteome: aa.AminoAcidSeqDict,
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
            else gtf.TranscriptGTFDict()
        self.genome = genome if genome is not None else dna.DNASeqDict()
        self.proteome = proteome if proteome is not None else \
            aa.AminoAcidSeqDict()
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

            biotype = self.annotation[transcript_id].transcript.biotype
            
            if biotype not in _VALID_BIOTYPES:
                continue

            if transcript_id not in self.vep.keys():
                self.vep[transcript_id] = VEPVariantRecords()

            chrom_seqname = record.location.split(':')[0]
            try:
                transcript_seq = self.annotation[transcript_id]\
                    .get_transcript_sequence(self.genome[chrom_seqname])
            except ValueError:
                print(transcript_id)
                raise
            
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

    def create_carnonical_peptide_pool(self, rule:str, exception:str=None,
            miscleavage:int=2, min_mw:float=500.):
        """ Creates the carnonical peptide pool with the protein sequences.
        
        Args:
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
            start (int): Index to start searching.
            min_mw (float): Minimal molecular weight of the peptides to report.
                Defaults to 500.
        """
        self.canonical_peptides = self.proteome.create_unique_peptide_pool(
            rule=rule, exception=exception,
            miscleavage=miscleavage, min_mw=min_mw
        )

    def call_variant_peptides(self, rule:str, exception:str=None,
            miscleavage:int=2, min_mw:float=500., ncpus:int=1):
        """
        Args:
            rule (str): The rule for enzymatic cleavage, e.g., trypsin.
            exception (str): The exception for cleavage rule.
            start (int): Index to start searching.
            min_mw (float): Minimal molecular weight of the peptides to report.
                Defaults to 500.
        """
        pool = ProcessPool(nodes=ncpus)

        if self.canonical_peptides is None:
            self.create_carnonical_peptide_pool(
                rule=rule,
                exception=exception,
                miscleavage=miscleavage,
                min_mw=min_mw
            )
        self.variant_peptides = set()
        variants: List[VEPVariantRecord]

        def wrapper(transcript_seq, transcript_id, protein_seq, variants):
            
            if protein_seq.seq.startswith('N'):
                transcript_seq, protein_seq = transcript_seq\
                    .validate_orf(protein_seq)

            graph = TranscriptVariantGraph(
                seq=transcript_seq,
                transcript_id=transcript_id,
                protein=protein_seq
            )
            graph.create_variant_graph(variants=variants)
            
            return graph.walk_and_splice(
                rule=rule,
                exception=exception,
                miscleavage=miscleavage,
                min_mw=min_mw
            )

        transcript_seqs = []
        transcript_ids = []
        protein_seqs = []
        variants_list = []

        for transcript_id, variants in self.vep.items():
            anno:gtf.TranscriptAnnotationModel = self.annotation[transcript_id]
            chrom = anno.transcript.location.seqname
            
            transcript_seq = anno.get_transcript_sequence(self.genome[chrom])
            protein_seq:AminoAcidSeqRecord = self.proteome[transcript_id]

            transcript_seqs.append(transcript_seq)
            transcript_ids.append(transcript_id)
            protein_seqs.append(protein_seq)
            variants_list.append(variants)

            if len(transcript_ids) == ncpus:
                results = pool.map(wrapper, transcript_seqs, transcript_ids,
                    protein_seqs, variants_list)
                transcript_seqs = []
                transcript_ids = []
                protein_seqs = []
                variants_list = []

                for peptides in results:
                    for peptide in peptides:
                        if str(peptide.seq) in self.canonical_peptides:
                            continue
                        same_peptide = get_equivalent(
                            self.variant_peptides, peptide)
                        if same_peptide:
                            same_peptide.id += '||' + peptide.id
                            same_peptide.name = same_peptide.id
                            same_peptide.description = same_peptide.id
                            continue
                        self.variant_peptides.add(peptide)
        
    def write_peptides(self, output_path:str):
        """ Write the variant peptides to FASTA file. """
        with open(output_path, 'w') as handle:
            SeqIO.write(self.variant_peptides, handle, 'fasta')


if __name__ == '__main__':
    from Bio import SeqIO
    # from moPepGen.vep.VEP2VariantPeptides import VEP2VariantPeptides, TranscriptVariantGraph
    # file_dir = 'test/files/downsampled_set'
    # adapter = VEP2VariantPeptides.set_up(
    #     vep_path=[f'{file_dir}/CPCG0100_gencode_v34_snp_chr22.tsv',
    #         f'{file_dir}/CPCG0100_gencode_v34_indel_chr22.tsv'],
    #     gtf_path=f'{file_dir}/gencode_v34_chr22.gtf',
    #     genome_path=f'{file_dir}/gencode_v34_genome_chr22.fasta',
    #     proteome_path=f'{file_dir}/gencode_v34_translations_chr22.fasta'
    # )
    gtf_path='/hot/ref/reference/GRCh38-EBI-GENCODE34/gencode.v34.chr_patch_hapl_scaff.annotation.gtf'
    genome_path='/hot/ref/reference/GRCh38-EBI-GENCODE34/GRCh38.p13.genome.fa'
    proteome_path='/hot/ref/reference/GRCh38-EBI-GENCODE34/gencode.v34.pc_translations.fa'
    adapter = VEP2VariantPeptides.set_up(               
        vep_path=[
            '/hot/projects/cpcgene/noncanonical_peptides/Mutation/gencodev34_grch38/VEP/germline/filtered_snv/CPCG0100.gencode.aa.tsv',
            '/hot/projects/cpcgene/noncanonical_peptides/Mutation/gencodev34_grch38/VEP/germline/filtered_indel/CPCG0100.gencode.aa.tsv'
        ],
        gtf_path=gtf_path,
        genome_path=genome_path,
        proteome_path=proteome_path,
    )
    # vep = adapter.vep
    # gtf = adapter.annotation
    # genome = adapter.genome
    # transcript_id = 'ENST00000519684.5'
    # transcript = gtf[transcript_id].get_transcript_sequence(genome['chr22'])
    # protein = adapter.proteome[transcript_id]
    # graph = TranscriptVariantGraph(transcript, transcript_id, protein)
    # graph.create_variant_graph(vep[transcript_id])
    # graph.walk_and_splice('trypsin')
    peptides = adapter.call_variant_peptides(rule='trypsin', ncpus=8)
