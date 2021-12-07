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

    frames_shifted = random.randint(-max_size + 1, max_size - 1)
    var_type = 'SNV' if frames_shifted == 0 else 'INDEL'
    tx_start = tx_model.transcript.location.start
    tx_end = tx_model.transcript.location.end
    if tx_model.transcript.strand == -1:
        tx_start, tx_end = tx_end, tx_start
    tx_start = anno.coordinate_genomic_to_gene(tx_start, gene_id)
    tx_end = anno.coordinate_genomic_to_gene(tx_end - 1, gene_id) + 1

    while True:
        if frames_shifted >= 0:
            start = random.randint(tx_start, tx_end - 1)
            end = start + 1
        else:
            start = random.randint(tx_start, tx_end - frames_shifted - 1)
            end = start - frames_shifted + 1
        start_genomic = anno.coordinate_gene_to_genomic(start, gene_id)
        end_genomic = anno.coordinate_gene_to_genomic(end - 1, gene_id) + 1
        if not exonic_only:
            break
        if all(tx_model.is_exonic(x) for x in [start_genomic, end_genomic]):
            break

    location = FeatureLocation(
        start=start,
        end=end,
        seqname=gene_id
    )

    ref_seq = str(gene_seq.seq[start:end])

    while True:
        if frames_shifted > 0:
            alt_seq = create_random_dna_sequence(frames_shifted)
            alt_seq = ref_seq + alt_seq
        elif frames_shifted == 0:
            alt_seq = create_random_dna_sequence(1)
        else:
            alt_seq = ref_seq[0]

        if alt_seq != ref_seq:
            break

    var_id = f"{gene_id}-{start}-{ref_seq}-{alt_seq}"

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
