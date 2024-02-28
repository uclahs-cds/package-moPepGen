""" Module to generate fake data """
import copy
import random
from typing import Tuple, List, Dict
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.dna.DNASeqRecord import DNASeqRecord
from moPepGen.dna.DNASeqDict import DNASeqDict
from moPepGen.gtf.GenomicAnnotation import GenomicAnnotation, GeneAnnotationModel, \
    TranscriptAnnotationModel, GTFSeqFeature
from moPepGen.seqvar import VariantRecord
from moPepGen.circ import CircRNAModel
from moPepGen import constant


# pylint: disable=R0912, R0915

DNA_ALPHABET = ['A', 'T', 'G', 'C']
ALT_SPLICE_TYPE_TO_RMATS = {
    'Deletion': ['SE', 'RI', 'A3SS', 'A5SS', 'MXE'],
    'Insertion': ['SE', 'RI', 'A3SS', 'A5SS', 'MXE'],
    'Substitution': ['MXE']
}
STOP_CODONS = ['TAA', 'TAG', 'TGA']

def create_random_dna_sequence(size:int) -> str:
    """ Create a randeom DNA sequence in str """
    return ''.join(random.choices(DNA_ALPHABET, k=size))

def fake_variant_record(anno:GenomicAnnotation, genome:DNASeqDict,
        tx_id:str, var_type:str, max_size:int, exonic_only:bool
        )->VariantRecord:
    """ Create a fake VariantRecord object """
    tx_model = anno.transcripts[tx_id]
    gene_id = tx_model.transcript.gene_id
    gene_model = anno.genes[gene_id]
    chrom = gene_model.chrom
    gene_seq = gene_model.get_gene_sequence(genome[chrom])

    if var_type == 'SNV':
        frames_shifted = 0
    else:
        while True:
            frames_shifted = random.randint(-max_size + 1, max_size - 1)
            if frames_shifted != 0:
                break
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
            start = random.randint(tx_start, tx_end - 1 + (frames_shifted - 1))
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
        cds_end_genomic = anno.coordinate_transcript_to_genomic(cds_end_tx - 1, tx_id)
        cds_end_gene = anno.coordinate_genomic_to_gene(cds_end_genomic, donor_gene_id) + 1
    else:
        cds_start_genomic = anno.coordinate_gene_to_genomic(0, donor_gene_id)
        cds_start_gene = anno.coordinate_genomic_to_gene(cds_start_genomic, donor_gene_id)
        cds_end_tx = len(donor_tx_seq) - 1
        cds_end_genomic = anno.coordinate_transcript_to_genomic(cds_end_tx - 1, tx_id)
        cds_end_gene = anno.coordinate_genomic_to_gene(cds_end_genomic, donor_gene_id) + 1
    donor_breakpoint = random.randint(cds_start_gene + 1, cds_end_gene - 1)
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
    cds_start_genomic = anno.coordinate_transcript_to_genomic(0, accepter_tx_id)
    cds_start_gene = anno.coordinate_genomic_to_gene(cds_start_genomic, accepter_gene_id)
    cds_end_tx = len(accepter_tx_seq) - 1
    cds_end_genomic = anno.coordinate_transcript_to_genomic(cds_end_tx - 1, accepter_tx_id)
    cds_end_gene = anno.coordinate_genomic_to_gene(cds_end_genomic, accepter_gene_id) + 1

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
    while True:
        var_type = random.choice(constant.ALTERNATIVE_SPLICING_TYPES)
        if var_type != 'Substitution' or len(anno.transcripts[tx_id].exon) > 3:
            break
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
        if rmats_type in ['A3SS', 'MXE']:
            i = random.choice(list(range(1, len(tx_model.exon) - 1)))
        else:
            i = random.choice(list(range(1, len(tx_model.exon))))
    else:
        if rmats_type in ['A5SS', 'MXE']:
            i = random.choice(list(range(1, len(tx_model.exon) - 1)))
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
        _id=f"{rmats_type}_{start_gene}-{end_gene}",
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
        end_genomic = random.randint(start_genomic + 1, intron_end - 1)
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
        end_genomic = random.randint(start_genomic + 1, intron_end - 2)
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
        _id=f"{rmats_type}_{insert_position}-{start_gene}-{end_gene}",
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

    second_start_genomic = random.randint(intron_start + 1, intron_end - 5)
    second_end_genomic = random.randint(second_start_genomic + 2, intron_end - 1)

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

def fake_circ_rna_model(anno:GenomicAnnotation, tx_id:str, ci_ratio:float=0.2
        ) -> CircRNAModel:
    """ Create a fake circRNA model for either CIRC or CI. """
    weights = [1-ci_ratio, ci_ratio]
    circ_type = random.choices(['CIRC', 'CI'], weights=weights, k=1)[0]
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
        start=anno.coordinate_genomic_to_gene(exon.location.start, gene_id)
        end=anno.coordinate_genomic_to_gene(exon.location.end - 1, gene_id)
        if tx_model.transcript.strand == -1:
            start, end = end, start
        end += 1
        location = FeatureLocation(start=start, end=end)
        fragment = SeqFeature(
            chrom=tx_id, location=location, attributes={}, type='exon'
        )
        fragments.append(fragment)

    fragments.sort()

    exon_indices = [f"E{x}" for x in
        range(first_exon_ind + 1, first_exon_ind + n_exons + 1)]

    _id = f'CIRC-{tx_id}-' + '-'.join(exon_indices)

    if tx_model.transcript.strand == 1:
        start_genomic = anno.coordinate_gene_to_genomic(
            fragments[0].location.start, gene_id
        )
    else:
        start_genomic = anno.coordinate_gene_to_genomic(
            fragments[-1].location.end - 1, gene_id
        )

    genomic_location = f"{chrom}:{start_genomic}"

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
    start = anno.coordinate_genomic_to_gene(tx_model.exon[exon_ind].location.end - 1, gene_id)
    end = anno.coordinate_genomic_to_gene(tx_model.exon[exon_ind + 1].location.start, gene_id)
    if tx_model.transcript.strand == -1:
        start, end = end, start
    end += 1

    location = FeatureLocation(start=start, end=end)
    fragment = SeqFeature(chrom=tx_id, location=location, attributes={}, type='intron')
    introns = [fragment]

    _id = f'CI-{tx_id}-I{exon_ind + 1}'

    if tx_model.transcript.strand == 1:
        start_genomic = introns[0].location.start
    else:
        start_genomic = introns[-1].location.end - 1

    genomic_location = f"{chrom}:{start_genomic}"

    return CircRNAModel(
        transcript_id=tx_id,
        fragments=introns,
        intron=[0],
        _id=_id,
        gene_id=gene_id,
        gene_name=tx_model.transcript.gene_name,
        genomic_location=genomic_location
    )


def fake_transcript_model(n_exons:int, is_coding:bool, is_selenoprotein:bool,
        chrom:str, strand:int, start_pos:int, gene_id:str, transcript_id:str,
        protein_id:str, cds_start_nf:bool, mrna_end_nf:bool, min_exon_size:int,
        max_exon_size:int, min_intron_size:int, max_intron_size:int
        ) -> TranscriptAnnotationModel:
    """ Create a fake transcript model """
    exon_set:List[GTFSeqFeature] = []
    cds_set:List[GTFSeqFeature] = []
    utr_set:List[GTFSeqFeature] = []
    three_utr_set:List[GTFSeqFeature] = []
    five_utr_set:List[GTFSeqFeature] = []
    offset = start_pos
    tx_id = f"FAKET{random.randint(1, 10000):08}"
    attributes = {
        'gene_id': gene_id,
        'transcript_id': transcript_id,
        'protein_id': protein_id,
        'tag': []
    }
    if cds_start_nf:
        attributes['tag'].append('cds_start_NF')
    if mrna_end_nf:
        attributes['tag'].append('mrna_end_NF')
    for i in range(n_exons):
        exon_len = random.randint(min_exon_size, max_exon_size)
        loc = FeatureLocation(start=offset, end=offset + exon_len, seqname=chrom, strand=strand)
        exon = GTFSeqFeature(
            location=loc, type='exon', id=tx_id, attributes=attributes, chrom=chrom
        )
        exon_set.append(exon)
        if is_coding:
            if i == 0:
                if strand == 1:
                    if cds_start_nf:
                        cds_start = offset
                    else:
                        cds_start = offset + random.randint(1, exon_len - 3)
                else:
                    if mrna_end_nf:
                        cds_start = offset
                    else:
                        cds_start = offset + random.randint(1, exon_len - 1)
            else:
                cds_start = offset

            if i == n_exons - 1:
                if strand == -1:
                    if cds_start_nf:
                        cds_end = offset + exon_len
                    else:
                        cds_end = offset + exon_len - random.randint(3, exon_len - 3)
                else:
                    if mrna_end_nf:
                        cds_end = offset + exon_len
                    else:
                        cds_end = offset + exon_len - random.randint(3, exon_len - 3)
            else:
                cds_end = offset + exon_len

            loc = FeatureLocation(start=cds_start, end=cds_end, seqname=chrom, strand=strand)
            cds = GTFSeqFeature(
                location=loc, type='CDS', id=tx_id, attributes=attributes, chrom=chrom
            )
            cds_set.append(cds)

        intron_size = random.randint(min_intron_size, max_intron_size)
        offset += exon_len + intron_size

    # define frame for each cds
    if is_coding:
        cds_len = sum(len(x.location) for x in cds_set)
        if strand == 1:
            for i,cds in enumerate(cds_set):
                if i == 0:
                    cds.frame = cds_len % 3
                else:
                    prev_exon = exon_set[i-1]
                    prev_cds = cds_set[i-1]
                    cds.frame = (3 - (len(prev_exon.location) - prev_cds.frame) % 3) % 3
        else:
            for i,cds in reversed(list(enumerate(cds_set))):
                if i == len(cds_set) - 1:
                    cds.frame = cds_len % 3
                else:
                    prev_exon = exon_set[i+1]
                    prev_cds = cds_set[i+1]
                    cds.frame = (3 - (len(prev_exon.location) - prev_cds.frame) % 3) % 3

    # add selenocysteine
    sec_set = []
    sec_pos_cds_set = []
    if is_selenoprotein:
        n_sec = random.choices([1,2,3], weights=(10,5,1))[0]
        cds_iter = iter(cds_set) if strand == 1 else reversed(cds_set)
        cds_len = 0
        cds_junctions = []
        for cds in cds_iter:
            cds_len += cds.location.end - cds.location.start
            cds_junctions.append(cds_len)
        cds_len -= cds_set[0].frame if strand == 1 else cds_set[-1].frame
        n_codons = int(cds_len/3)
        offset = 1
        for i in range(n_sec):
            if offset >= n_codons - 2:
                break
            while True:
                sec_pos_codon = random.randint(offset, n_codons - 2)
                sec_pos_cds = sec_pos_codon * 3
                if strand == 1:
                    sec_pos_cds += cds_set[0].frame
                    if not any(sec_pos_cds < j <= sec_pos_cds + 3 for j in cds_junctions):
                        break
                else:
                    sec_pos_cds += cds_set[-1].frame
                    if not any(sec_pos_cds <= j < sec_pos_cds + 3 for j in cds_junctions):
                        break
            offset = sec_pos_codon + 1
            sec_pos_cds_set.append(sec_pos_cds)

        sec_pos_cds_iter = iter(sec_pos_cds_set)
        sec_pos_cds = next(sec_pos_cds_iter, None)
        k = 0
        cds_iter = iter(cds_set) if strand == 1 else reversed(cds_set)
        cds = next(cds_iter, None)
        while cds and sec_pos_cds:
            cds_size = cds.location.end - cds.location.start
            if k + cds_size <= sec_pos_cds:
                k += cds_size
                cds = next(cds_iter, None)
                continue

            if strand == 1:
                sec_pos = sec_pos_cds - k + cds.location.start
                loc = FeatureLocation(sec_pos, sec_pos + 3, strand=strand)
            else:
                sec_pos = cds.location.end - (sec_pos_cds - k)
                loc = FeatureLocation(sec_pos - 3, sec_pos, strand=strand)
            sec = GTFSeqFeature(
                location=loc, type='selenocysteine', id=tx_id,
                attributes=attributes, chrom=chrom
            )
            sec_set.append(sec)
            sec_pos_cds = next(sec_pos_cds_iter, None)

    # add UTR
    if is_coding:
        if cds_set[0].location.start != exon_set[0].location.start:
            loc = FeatureLocation(
                start=exon_set[0].location.start, end=cds_set[0].location.start,
                seqname=chrom, strand=strand
            )
            utr = GTFSeqFeature(
                location=loc, type='UTR', id=tx_id, attributes=attributes, chrom=chrom
            )
            utr_set.append(utr)
            if strand == 1:
                five_utr_set.append(utr)
            else:
                three_utr_set.append(utr)

        if cds_set[-1].location.end != exon_set[-1].location.end:
            loc = FeatureLocation(
                start=cds_set[-1].location.end, end=exon_set[-1].location.end,
                seqname=chrom, strand=strand
            )
            utr = GTFSeqFeature(
                location=loc, type='UTR', id=tx_id, attributes=attributes, chrom=chrom
            )
            utr_set.append(utr)
            if strand == 1:
                three_utr_set.append(utr)
            else:
                five_utr_set.append(utr)

    loc = FeatureLocation(
        start=exon_set[0].location.start, end=exon_set[-1].location.end,
        seqname=chrom, strand=strand
    )
    transcript = GTFSeqFeature(
        location=loc, type='transcript', id=tx_id,
        attributes=attributes, chrom=chrom
    )

    return TranscriptAnnotationModel(
        transcript=transcript, cds=cds_set, exon=exon_set, utr=utr_set,
        selenocysteine=sec_set, five_utr=five_utr_set, three_utr=three_utr_set,
        is_protein_coding=is_coding, transcript_id=transcript_id,
        gene_id=gene_id, protein_id=protein_id
    )

def fake_genomic_annotation(n_genes:int, chrom:str, min_exons:int, max_exons:int,
        min_exon_size:int, max_exon_size:int, min_intron_size:int, max_intron_size:int,
        min_intergenic_size:int, max_intergenic_size:int) -> GenomicAnnotation:
    """ Create afke genomic annotation object. """
    anno = GenomicAnnotation()

    offset = 0
    for _ in range(n_genes):
        strand = random.choice((1, -1))
        is_coding = random.choice((True, False))
        n_exons = random.randint(min_exons, max_exons)
        while True:
            x = random.randint(1, 1000)
            gene_id = f"FAKEG{x:08}"
            tx_id = f"FAKET{x:08}"
            protein_id = f"FAKEP{x:08}"
            if gene_id not in anno.genes and tx_id not in anno.transcripts:
                break
        if is_coding:
            cds_start_nf = random.choice((True, False))
            mrna_end_nf = random.choice((True, False))
        else:
            cds_start_nf, mrna_end_nf = False, False

        is_selenoprotein = random.choices([True, False], weights=(19,1))[0] \
            if is_coding else False

        tx_model = fake_transcript_model(
            n_exons=n_exons, is_coding=is_coding, is_selenoprotein=is_selenoprotein,
            chrom=chrom, strand=strand, start_pos=offset, gene_id=gene_id,
            transcript_id=tx_id, protein_id=protein_id, cds_start_nf=cds_start_nf,
            mrna_end_nf=mrna_end_nf, min_exon_size=min_exon_size,
            max_exon_size=max_exon_size, min_intron_size=min_intron_size,
            max_intron_size=max_intron_size
        )
        gene_start = tx_model.transcript.location.start
        gene_end = tx_model.transcript.location.end
        loc = FeatureLocation(start=gene_start, end=gene_end, seqname=chrom)
        gene_model = GeneAnnotationModel(
            location=loc, chrom=chrom, transcripts=[tx_id], type='gene',
            attributes=copy.deepcopy(tx_model.transcript.attributes)
        )
        anno.genes[gene_id] = gene_model
        anno.transcripts[tx_id] = tx_model

        offset += len(tx_model.transcript) \
            + random.randint(min_intergenic_size, max_intergenic_size)
    return anno

def fake_genome(anno:GenomicAnnotation) -> DNASeqDict:
    """ Create a fake genome from an genomic annotation. """
    genome_sizes:Dict[str, int] = {}
    for gene_model in anno.genes.values():
        chrom = gene_model.chrom
        if chrom not in genome_sizes:
            genome_sizes[chrom] = gene_model.location.end
        else:
            genome_sizes[chrom] = max(genome_sizes[chrom], gene_model.location.end)

    genome_nts = {k:random.choices(DNA_ALPHABET, k=v) for k,v in genome_sizes.items()}

    for tx_model in anno.transcripts.values():
        chrom_nts = genome_nts[tx_model.transcript.chrom]
        if tx_model.is_protein_coding:
            attrs = tx_model.transcript.attributes
            # make sure cds start is start codon
            if 'tag' not in attrs or 'cds_start_NF' not in attrs:
                if tx_model.transcript.strand == 1:
                    cds_start = tx_model.cds[0].location.start + tx_model.cds[0].frame
                    chrom_nts[cds_start] = 'A'
                    chrom_nts[cds_start + 1] = 'T'
                    chrom_nts[cds_start + 2] = 'G'
                else:
                    cds_start = tx_model.cds[-1].location.end - tx_model.cds[-1].frame - 1
                    chrom_nts[cds_start] = 'T'
                    chrom_nts[cds_start - 1] = 'A'
                    chrom_nts[cds_start - 2] = 'C'

            # make sure cds end is stop codon
            if 'tag' not in attrs or 'mrna_end_NF' not in attrs['tag']:
                stop_codon = random.choice(STOP_CODONS)
                if tx_model.transcript.strand == 1:
                    cds_end = tx_model.cds[-1].location.end
                    chrom_nts[cds_end] = stop_codon[0]
                    chrom_nts[cds_end + 1] = stop_codon[1]
                    chrom_nts[cds_end + 2] = stop_codon[2]
                else:
                    stop_codon = str(Seq(stop_codon).reverse_complement())
                    cds_end = tx_model.cds[0].location.start
                    chrom_nts[cds_end - 1] = stop_codon[2]
                    chrom_nts[cds_end - 2] = stop_codon[1]
                    chrom_nts[cds_end - 3] = stop_codon[0]

            # make sure no stop codon in cds
            if tx_model.transcript.strand == 1:
                i = tx_model.cds[0].location.start + tx_model.cds[0].frame
                cds_iter = iter(tx_model.cds)
                cds = next(cds_iter, None)
                while True:
                    if not cds:
                        break
                    if i >= cds.location.end:
                        if cds is tx_model.cds[-1]:
                            break
                        cds = next(cds_iter, None)
                        i = cds.location.start

                    if i + 2 > cds.location.end and cds is tx_model.cds[-1]:
                        break

                    codon = []
                    positions = []
                    k = i
                    while len(codon) < 3:
                        if k >= cds.location.end:
                            cds = next(cds_iter, None)
                            k = cds.location.start
                        codon.append(chrom_nts[k])
                        positions.append(k)
                        k += 1

                    while ''.join(codon) in STOP_CODONS:
                        codon = random.choices(DNA_ALPHABET, k=3)
                        for pos, nt in zip(positions, codon):
                            chrom_nts[pos] = nt

                    i = k
            else:
                i = tx_model.cds[-1].location.end - tx_model.cds[-1].frame - 1
                cds_iter = reversed(tx_model.cds)
                cds = next(cds_iter, None)
                while True:
                    if not cds:
                        break
                    if i < cds.location.start:
                        if cds is tx_model.cds[0]:
                            break
                        cds = next(cds_iter, None)
                        i = cds.location.end - 1
                    if i - 2 < cds.location.start and cds is tx_model.cds[0]:
                        break

                    codon = []
                    positions = []
                    k = i
                    while len(codon) < 3:
                        if k < cds.location.start:
                            cds = next(cds_iter, None)
                            k = cds.location.end - 1
                        codon.append(chrom_nts[k])
                        positions.append(k)
                        k -= 1

                    while str(Seq(''.join(codon)).complement()) in STOP_CODONS:
                        codon = random.choices(DNA_ALPHABET, k=3)
                        for pos, nt in zip(positions, codon):
                            chrom_nts[pos] = nt

                    i = k

            # make sure selenocysteine position is TGA
            for sec in tx_model.selenocysteine:
                if tx_model.transcript.strand == 1:
                    chrom_nts[sec.location.start] = 'T'
                    chrom_nts[sec.location.start + 1] = 'G'
                    chrom_nts[sec.location.start + 2] = 'A'
                else:
                    chrom_nts[sec.location.start] = 'T'
                    chrom_nts[sec.location.start + 1] = 'C'
                    chrom_nts[sec.location.start + 2] = 'A'

    genome = DNASeqDict()
    for chrom_id, chrom_nts in genome_nts.items():
        chrom_seq = DNASeqRecord(
            seq=Seq(''.join(chrom_nts)),
            id=chrom_id, name=chrom_id, description=chrom_id
        )
        genome[chrom_id] = chrom_seq

    return genome

def fake_genome_and_annotation(n_genes:int) -> Tuple[DNASeqDict, GenomicAnnotation]:
    """ Create a fake genome and annotation. """
    params = {
        'min_intron_size': 20,
        'max_intron_size': 500,
        'min_exon_size': 10,
        'max_exon_size': 300,
        'min_intergenic_size': 10,
        'max_intergenic_size': 50,
        'min_exons': 3,
        'max_exons': 10
    }
    chrom = 'chrF'

    anno = fake_genomic_annotation(
        n_genes=n_genes, chrom=chrom,
        min_exons=params['min_exons'], max_exons=params['max_exons'],
        min_exon_size=params['min_exon_size'], max_exon_size=params['max_exon_size'],
        min_intron_size=params['min_intron_size'], max_intron_size=params['max_intron_size'],
        min_intergenic_size=params['min_intergenic_size'],
        max_intergenic_size=params['max_intergenic_size']
    )

    genome = fake_genome(anno)

    return genome, anno
