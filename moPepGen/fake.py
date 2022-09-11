""" Module to generate fake data """
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

def fake_fusion(anno:GenomicAnnotation, genome:DNASeqDict, tx_id:str) -> VariantRecord:
    """ Create a fake fusion """
    if len(anno.transcripts) <=1 :
        raise ValueError("Not enough transcripts to create fake fusion.")
    donor_tx_model = anno.transcripts[tx_id]
    donor_chrom = donor_tx_model.transcript.chrom
    donor_tx_seq = donor_tx_model.get_transcript_sequence(genome[donor_chrom])
    donor_gene_id = donor_tx_model.gene_id
    donor_gene_model = anno.genes[donor_gene_id]
    donor_gene_seq = donor_gene_model.get_gene_sequence(genome[donor_chrom])

    if donor_tx_model.is_protein_coding:
        cds_start_tx = donor_tx_model.get_cds_start_index()
        cds_start_genomic = anno.coordinate_transcript_to_genomic(cds_start_tx, tx_id)
        cds_start_gene = anno.coordinate_genomic_to_gene(cds_start_genomic, donor_gene_id)
        cds_end_tx = donor_tx_model.get_cds_end_index(donor_tx_seq, donor_tx_seq.orf.start)
        cds_end_genomic = anno.coordinate_transcript_to_genomic(cds_end_tx, tx_id)
        cds_end_gene = anno.coordinate_genomic_to_gene(cds_end_genomic, donor_gene_id)
    else:
        cds_start_genomic = anno.coordinate_gene_to_genomic(0, tx_id)
        cds_start_gene = anno.coordinate_genomic_to_gene(cds_start_genomic, donor_gene_id)
        cds_end_tx = len(donor_tx_seq)
        cds_end_genomic = anno.coordinate_transcript_to_genomic(cds_end_tx, tx_id)
        cds_end_gene = anno.coordinate_genomic_to_gene(cds_end_genomic, donor_gene_id)
    donor_breakpoint = random.randint(cds_start_gene, cds_end_gene - 1)
    donor_breakpoint_genomic = anno.coordinate_gene_to_genomic(donor_breakpoint, donor_gene_id)
    ref_seq = donor_gene_seq.seq[donor_breakpoint]

    while True:
        accepter_tx_id = random.choice(list(anno.transcripts.keys()))
        if accepter_tx_id != tx_id:
            break
    accepter_tx_model = anno.transcripts[accepter_tx_id]
    accepter_chrom = accepter_tx_model.transcript.chrom
    accepter_tx_seq = accepter_tx_model.get_transcript_sequence(genome[accepter_chrom])
    accepter_gene_id = accepter_tx_model.gene_id
    cds_start_tx = accepter_tx_model.get_cds_start_index()
    cds_start_genomic = anno.coordinate_transcript_to_genomic(cds_start_tx, accepter_tx_id)
    cds_start_gene = anno.coordinate_genomic_to_gene(cds_start_genomic, accepter_gene_id)
    cds_end_tx = accepter_tx_model.get_cds_end_index(accepter_tx_seq, accepter_tx_seq.orf.start)
    cds_end_genomic = anno.coordinate_transcript_to_genomic(cds_end_tx, accepter_tx_id)
    cds_end_gene = anno.coordinate_genomic_to_gene(cds_end_genomic, accepter_gene_id)

    accepter_breakpoint = random.randint(cds_start_gene, cds_end_gene - 1)

    accepter_breakpoint_genomic = anno.coordinate_gene_to_genomic(
        accepter_breakpoint, accepter_gene_id)

    location = FeatureLocation(
        seqname=donor_gene_id,
        start=donor_breakpoint,
        end=donor_breakpoint + 1
    )

    attrs = {
        'TRANSCRIPT_ID': tx_id,
        'GENE_SYMBOL': donor_tx_model.gene_name,
        'GENOMIC_POSITION': donor_breakpoint_genomic,
        'ACCEPTER_GENE_ID': accepter_gene_id,
        'ACCEPTER_TRANSCRIPT_ID': accepter_tx_id,
        'ACCEPTER_SYMBOL': accepter_tx_model.gene_name,
        'ACCEPTER_POSITION': accepter_breakpoint,
        'ACCEPTER_GENOMIC_POSITION': accepter_breakpoint_genomic
    }

    fusion_id = f"FUSION-{tx_id}:{donor_breakpoint}-{accepter_tx_id}:{accepter_breakpoint}"

    return VariantRecord(
        location=location,
        ref=ref_seq,
        alt='<FUSION>',
        _type='Fusion',
        _id=fusion_id,
        attrs=attrs
    )
