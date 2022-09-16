""" Module to generate fake data """
import random
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.dna.DNASeqDict import DNASeqDict
from moPepGen.gtf.GenomicAnnotation import GenomicAnnotation
from moPepGen.seqvar import VariantRecord
from moPepGen.seqvar.VariantRecord import ALTERNATIVE_SPLICING_TYPES
from moPepGen.circ import CircRNAModel


DNA_ALPHABET = ['A', 'T', 'G', 'C']
ALT_SPLICE_TYPE_TO_RMATS = {
    'Deletion': ['SE', 'RI', 'A3SS', 'A5SS', 'MXE'],
    'Insertion': ['SE', 'RI', 'A3SS', 'A5SS', 'MXE'],
    'Substitution': ['MXE']
}

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

def fake_rmats_record(anno:GenomicAnnotation, genome:DNASeqDict, tx_id:str
        ) -> VariantRecord:
    """ Create an alternative splicing variant """
    var_type = random.choice(ALTERNATIVE_SPLICING_TYPES)
    rmats_type = random.choice(ALT_SPLICE_TYPE_TO_RMATS[var_type])

    err_message = f"Not valid variant type of {rmats_type} {var_type}"

    if var_type == 'Insertion':
        return fake_intron_insertion(anno, genome, tx_id, rmats_type)
    if var_type == 'Deletion':
        return fake_exon_deletion(anno, genome, tx_id, rmats_type)
    if var_type == 'Substitution':
        return fake_mxe_substitution(anno, genome, tx_id)
    raise ValueError(err_message)

def fake_exon_deletion(anno:GenomicAnnotation, genome:DNASeqDict, tx_id:str,
        rmats_type:str) -> VariantRecord:
    """ Create an deletion of an entire or a part of an exon. """
    tx_model = anno.transcripts[tx_id]
    gene_id = tx_model.gene_id
    gene_model = anno.genes[gene_id]
    chrom = gene_model.chrom
    gene_seq = gene_model.get_gene_sequence(genome[chrom])

    # sample an exon
    if gene_model.strand == 1:
        i = random.choice(list(range(1, len(tx_model.exon))))
    else:
        i = random.choice(list(range(len(tx_model.exon) - 1)))

    exon = tx_model.exon[i]

    if rmats_type in ['SE', 'MXE']:
        start_genomic = exon.location.start
        end_genomic = exon.location.end
    elif rmats_type == 'RI':
        start_genomic = random.randint(exon.location.start + 1, exon.location.end - 2)
        end_genomic = random.randint(start_genomic, exon.location.end - 1)
    elif (rmats_type == 'A5SS' and gene_model.strand == 1)\
            or (rmats_type == 'A3SS' and gene_model.strand == -1):
        start_genomic = random.randint(exon.location.start + 1, exon.location.end - 2)
        end_genomic = exon.location.end
    elif (rmats_type == 'A5SS' and gene_model.strand == -1)\
            or (rmats_type == 'A3SS' and gene_model.strand == 1):
        start_genomic = exon.location.start
        end_genomic = random.randint(start_genomic, exon.location.end - 1)
    else:
        raise ValueError(
            "Alternative splicing even could not be recognized with rmats_type"
            f"={rmats_type} and strand={gene_model.strand}"
        )

    start_gene = anno.coordinate_genomic_to_gene(start_genomic, gene_id)
    end_gene = anno.coordinate_genomic_to_gene(end_genomic - 1, gene_id)

    if gene_model.strand == -1:
        start_gene, end_gene = end_gene, start_gene

    end_gene += 1
    ref = str(gene_seq.seq[start_gene])

    genomic_position = f"{chrom}:{start_genomic+1}:{end_genomic}"

    if gene_model.strand == 1:
        insert_position = anno.coordinate_genomic_to_gene(
            tx_model.exon[i - 1].location.end - 1, gene_id
        )
    else:
        insert_position = anno.coordinate_genomic_to_gene(
            tx_model.exon[i + 1].location.start, gene_id
        )

    location = FeatureLocation(
        seqname=gene_id, start=insert_position, end=insert_position + 1
    )

    attrs = {
        'TRANSCRIPT_ID': tx_id,
        'START': start_gene,
        'END': end_gene,
        'GENE_SYMBOL': gene_model.gene_name,
        'GENOMIC_POSITION': genomic_position
    }

    return VariantRecord(
        location=location,
        ref=ref,
        alt='<DEL>',
        _type='Deletion',
        _id=f"{rmats_type}_{start_gene}",
        attrs=attrs
    )

def fake_intron_insertion(anno:GenomicAnnotation, genome:DNASeqDict,
        tx_id:str, rmats_type:str) -> VariantRecord:
    """ Create an insertion of an entire intron (RI) or a part of it (e.g.
    A3SS & A5SS). """
    tx_model = anno.transcripts[tx_id]
    gene_id = tx_model.gene_id
    gene_model = anno.genes[gene_id]
    chrom = gene_model.chrom
    gene_seq = gene_model.get_gene_sequence(genome[chrom])

    # sample an intron
    i = random.choice(list(range(len(tx_model.exon) - 1)))
    intron_start = tx_model.exon[i].location.end
    intron_end = tx_model.exon[i + 1].location.start

    if rmats_type in ['SE', 'MXE']:
        start_genomic = random.randint(intron_start + 1, intron_end - 2)
        end_genomic = random.randint(start_genomic, intron_end - 1)
    elif rmats_type == 'RI':
        start_genomic = intron_start
        end_genomic = intron_end
    elif (rmats_type == 'A3SS' and gene_model.strand == 1)\
            or (rmats_type == 'A5SS' and gene_model.strand == -1):
        start_genomic = random.randint(intron_start + 1, intron_end - 2)
        end_genomic = intron_end
    elif (rmats_type == 'A3SS' and gene_model.strand == -1)\
            or (rmats_type == 'A5SS' and gene_model.strand == 1):
        start_genomic = intron_start
        end_genomic = random.randint(start_genomic, intron_end - 1)
    else:
        raise ValueError(
            "Alternative splicing even could not be recognized with rmats_type"
            f"={rmats_type} and strand={gene_model.strand}"
        )

    start_gene = anno.coordinate_genomic_to_gene(start_genomic, gene_id)
    end_gene = anno.coordinate_genomic_to_gene(end_genomic - 1, gene_id)

    if gene_model.strand == -1:
        start_gene, end_gene = end_gene, start_gene

    end_gene += 1

    if gene_model.strand == 1:
        insert_position = anno.coordinate_genomic_to_gene(
            tx_model.exon[i].location.end - 1, gene_id
        )
    else:
        insert_position = anno.coordinate_genomic_to_gene(
            tx_model.exon[i + 1].location.start, gene_id
        )

    location = FeatureLocation(
        seqname=gene_id, start=insert_position, end=insert_position + 1
    )

    ref = str(gene_seq.seq[insert_position])

    genomic_position = f"{chrom}:{start_genomic+1}:{end_genomic}"

    attrs = {
        'TRANSCRIPT_ID': tx_id,
        'DONOR_START': start_gene,
        'DONOR_END': end_gene,
        'DONOR_GENE_ID': gene_id,
        'GENE_SYMBOL': gene_model.gene_name,
        'GENOMIC_POSITION': genomic_position
    }

    return VariantRecord(
        location=location,
        ref=ref,
        alt='<INS>',
        _type='Insertion',
        _id=f"{rmats_type}_{start_gene}",
        attrs=attrs
    )

def fake_mxe_substitution(anno:GenomicAnnotation, genome:DNASeqDict,
        tx_id:str) -> VariantRecord:
    """ Create a MXE substitution. """
    tx_model = anno.transcripts[tx_id]
    gene_id = tx_model.gene_id
    gene_model = anno.genes[gene_id]
    chrom = gene_model.chrom
    gene_seq = gene_model.get_gene_sequence(genome[chrom])

    # sample an exon
    i = random.choice(list(range(1, len(tx_model.exon) - 1)))
    exon = tx_model.exon[i]
    if random.choice([-1,1]) == 1:
        intron_start = tx_model.exon[i].location.end
        intron_end = tx_model.exon[i+1].location.start
    else:
        intron_start = tx_model.exon[i-1].location.end
        intron_end = tx_model.exon[i].location.start

    first_start_genomic = exon.location.start
    first_end_genomic = exon.location.end

    second_start_genomic = random.randint(intron_start + 1, intron_end - 2)
    second_end_genomic = random.randint(second_start_genomic, intron_end - 1)

    first_start_gene = anno.coordinate_genomic_to_gene(first_start_genomic, gene_id)
    first_end_gene = anno.coordinate_genomic_to_gene(first_end_genomic - 1, gene_id)
    second_start_gene = anno.coordinate_genomic_to_gene(second_start_genomic, gene_id)
    second_end_gene = anno.coordinate_genomic_to_gene(second_end_genomic - 1, gene_id)

    if gene_model.strand == -1:
        first_start_gene, first_end_gene = first_end_gene, first_start_gene
        second_start_gene, second_end_gene = second_end_gene, second_start_gene
    second_end_gene += 1
    first_end_gene += 1

    location = FeatureLocation(
        seqname=gene_id, start=first_start_gene, end=first_end_gene
    )
    ref = str(gene_seq.seq[first_start_gene])
    genomic_position = f"{chrom}-{first_start_gene + 1}:{first_end_gene}" +\
        f"-{second_start_gene + 1}:{second_end_gene}"
    _id=f"MXE_{first_start_gene + 1}-{first_end_gene}" +\
            f"{second_start_gene}-{second_end_gene}"
    attrs = {
        'TRANSCRIPT_ID': tx_id,
        'START': first_start_gene,
        'END': first_end_gene,
        'DONOR_START': second_start_gene,
        'DONOR_END': second_end_gene,
        'DONOR_GENE_ID': gene_id,
        'COORDINATE': 'gene',
        'GENE_SYMBOL': gene_model.gene_name,
        'GENOMIC_POSITION': genomic_position
    }

    return VariantRecord(
        location=location,
        ref=ref,
        alt='<SUB>',
        _type='Substitution',
        _id=_id,
        attrs=attrs
    )

def fake_circ_rna_model(anno:GenomicAnnotation, tx_id:str) -> CircRNAModel:
    """ Create a fake circRNA model for either CIRC or CI. """
    circ_type = random.choice(['CIRC', 'CI'])
    if circ_type == 'CIRC':
        return fake_circ_rna_model_circ(anno, tx_id)
    return fake_circ_rna_model_ci(anno, tx_id)

def fake_circ_rna_model_circ(anno:GenomicAnnotation, tx_id:str
        ) -> CircRNAModel:
    """ Create a fake circRNA model of CIRC """
    tx_model = anno.transcripts[tx_id]
    gene_id = tx_model.gene_id
    chrom = tx_model.transcript.chrom

    n_exons = random.randint(1, len(tx_model.exon) - 1)
    first_exon_ind = random.randint(0, len(tx_model.exon) - n_exons - 1)
    exons = tx_model.exon[first_exon_ind:first_exon_ind + n_exons]

    fragments = []
    for exon in exons:
        fragment = SeqFeature(
            chrom=tx_id, location=exon.location, attributes={}, type='exon'
        )
        fragments.append(fragment)

    exon_indices = [f"E{x}" for x in
        range(first_exon_ind + 1, first_exon_ind + n_exons + 1)]

    _id = f'CIRC-{tx_id}-' + '-'.join(exon_indices)

    if tx_model.transcript.strand == 1:
        start = anno.coordinate_genomic_to_gene(exons[0].location.start, gene_id)
    else:
        start = anno.coordinate_genomic_to_gene(exons[-1].location.end - 1, gene_id)

    genomic_location = f"{chrom}:{start}"

    return CircRNAModel(
        transcript_id=tx_id,
        fragments=fragments,
        intron=[],
        _id=_id,
        gene_id=gene_id,
        gene_name=tx_model.transcript.gene_name,
        genomic_location=genomic_location
    )

def fake_circ_rna_model_ci(anno:GenomicAnnotation, tx_id:str
        ) -> CircRNAModel:
    """ Create a fake circRNA model of CI (circular intron).  """
    tx_model = anno.transcripts[tx_id]
    gene_id = tx_model.gene_id
    chrom = tx_model.transcript.chrom

    exon_ind = random.randint(0, len(tx_model.exon) - 2)
    start = tx_model.exon[exon_ind].location.end
    end = tx_model.exon[exon_ind + 1].location.start - 1

    location = FeatureLocation(start=start, end=end)
    fragment = SeqFeature(chrom=tx_id, location=location, attributes={}, type='intron')
    introns = [fragment]

    _id = f'CI-{tx_id}-I{exon_ind + 1}'

    if tx_model.transcript.strand == 1:
        start_gene = anno.coordinate_genomic_to_gene(introns[0].location.start, gene_id)
    else:
        start_gene = anno.coordinate_genomic_to_gene(introns[-1].location.end - 1, gene_id)

    genomic_location = f"{chrom}:{start_gene}"

    return CircRNAModel(
        transcript_id=tx_id,
        fragments=introns,
        intron=[0],
        _id=_id,
        gene_id=gene_id,
        gene_name=tx_model.transcript.gene_name,
        genomic_location=genomic_location
    )
