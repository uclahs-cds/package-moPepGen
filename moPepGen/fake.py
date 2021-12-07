""" Module to generate fast data """
import random
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.dna.DNASeqDict import DNASeqDict
from moPepGen.gtf.GenomicAnnotation import GenomicAnnotation
from moPepGen.seqvar import VariantRecord


DNA_ALPHABET = ['A', 'T', 'G', 'C']

def create_random_dna_sequence(size:int) -> str:
    """ Create a randeom DNA sequence in str """
    return ''.join(random.choices(DNA_ALPHABET, k=size))

def fake_variant_record(anno:GenomicAnnotation, genome:DNASeqDict,
        tx_id:str, max_size:int, exonic_only:bool)->VariantRecord:
    """ Create a fake VariantRecord object """
    tx_model = anno.transcripts[tx_id]
    gene_id = tx_model.transcript.gene_id
    gene_model = anno.genes[gene_id]
    chrom = gene_model.chrom
    gene_seq = gene_model.get_gene_sequence(genome[chrom])

    var_size = random.randint(1, max_size)
    var_type = 'SNV' if var_size == 1 else 'INDEL'
    tx_start = tx_model.transcript.location.start
    tx_end = tx_model.transcript.location.end
    if tx_model.transcript.strand == -1:
        tx_start, tx_end = tx_end, tx_start
    tx_start = anno.coordinate_genomic_to_gene(tx_start, gene_id)
    tx_end = anno.coordinate_genomic_to_gene(tx_end - 1, gene_id) + 1

    start = random.randint(tx_start, tx_end - var_size)
    end = start + var_size
    start_genomic = anno.coordinate_gene_to_genomic(start, gene_id)
    end_genomic = anno.coordinate_gene_to_genomic(end - 1, gene_id) + 1

    if exonic_only:
        while not tx_model.is_exonic(start_genomic) or not tx_model.is_exonic(end_genomic):
            start = random.randint(tx_start, tx_end - var_size)
            end = start + var_size
            start_genomic = anno.coordinate_gene_to_genomic(start, gene_id)
            end_genomic = anno.coordinate_gene_to_genomic(end - 1, gene_id) + 1

    location = FeatureLocation(
        start=start,
        end=start + 1,
        seqname=gene_id
    )

    ref_seq = str(gene_seq.seq[start])
    alt_seq = create_random_dna_sequence(var_size)
    while alt_seq == ref_seq:
        alt_seq = create_random_dna_sequence(var_size)

    var_id = f"{gene_id}-{start}:{end}"

    genomic_start = anno.coordinate_gene_to_genomic(start, gene_id)
    genomic_end = anno.coordinate_gene_to_genomic(end - 1, gene_id) + 1
    genomic_position = f"{chrom}-{genomic_start}:{genomic_end}"

    attrs = {
        'TRANSCRIPT_ID': tx_id,
        'GENOMIC_POSITION': genomic_position,
        'GENE_SYMBOL': gene_model.gene_name
    }

    return VariantRecord(
        location=location,
        ref=ref_seq,
        alt=alt_seq,
        _type=var_type,
        _id=var_id,
        attrs=attrs
    )
