""""""
from typing import List, Dict
import argparse
from Bio import SeqIO, SeqUtils
from moPepGen import svgraph, dna, gtf, aa, seqvar, logger, get_equivalent
from moPepGen.SeqFeature import FeatureLocation


def call_variant_peptides(args:argparse.Namespace) -> None:
    """"""
    variant_files:List[str] = args.input_variants
    genome_fasta:str = args.genome_fasta
    proteome_fasta:str = args.proteome_fasta
    annotation_gtf:str = args.annotation_gtf
    output_fasta:str = args.output_fasta
    verbose:bool = args.verbose
    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    exception = 'trypsin_exception' if rule == 'trypsin' else None

    if verbose:
        logger('moPepGen callPeptide started.')

    annotations = gtf.TranscriptGTFDict()
    annotations.dump_gtf(annotation_gtf)
    if verbose:
        logger('Annotation GTF loaded.')
    
    genome = dna.DNASeqDict()
    genome.dump_fasta(genome_fasta)
    if verbose:
        logger('Genome assembly FASTA loaded.')
    
    proteome = aa.AminoAcidSeqDict()
    proteome.dump_fasta(proteome_fasta)
    if verbose:
        logger('Proteome FASTA loaded.')
    
    variants:Dict[str, List[seqvar.VariantRecord]] = {}
    for file in variant_files:
        with open(file, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                transcript_id = fields[0]
                record = seqvar.VariantRecord(
                    location=FeatureLocation(
                        seqname=transcript_id,
                        start=int(fields[1]),
                        end=int(fields[2])
                    ),
                    ref=fields[3],
                    alt=fields[4],
                    type=fields[5],
                    id=fields[6]
                )
                if transcript_id not in variants:
                    variants[transcript_id] = [record]
                else:
                    variants[transcript_id].append(record)
        logger(f'Variant file {file} loaded.')
    
    for records in variants.values():
        records.sort()
    if verbose:
        logger(f'Variant records sorted.')

    carnonical_peptides = proteome.create_unique_peptide_pool(
        rule=rule, exception=exception, miscleavage=miscleavage, min_mw=min_mw
    )
    if verbose:
        logger('Carnonical peptide pool generated.')

    variant_peptides = set()

    i = 0
    for transcript_id, variant_records in variants.items():
        anno:gtf.TranscriptAnnotationModel = annotations[transcript_id]
        chrom = anno.transcript.location.seqname
        transcript_seq = anno.get_transcript_sequence(genome[chrom])
        protein_seq:aa.AminoAcidSeqRecord = proteome[transcript_id]

        dgraph = svgraph.TranscriptVariantGraph(
            seq=transcript_seq,
            transcript_id=transcript_id
        )
        dgraph.create_variant_graph(variant_records)
        dgraph.fit_into_codons()
        pgraph = dgraph.translate()
    
        pgraph.form_cleavage_graph(rule=rule, exception=exception)
        peptides = pgraph.call_vaiant_peptides(miscleavage=miscleavage)

        for peptide in peptides:
            if SeqUtils.molecular_weight(peptide.seq, 'protein') < min_mw:
                continue
            if str(peptide.seq) in carnonical_peptides:
                continue
            same_peptide = get_equivalent(variant_peptides, peptide)
            if same_peptide:
                new_label = peptide.id
                same_peptide.id += ('||' + new_label)
                same_peptide.name = same_peptide.id
                same_peptide.description = same_peptide.id
            else:
                variant_peptides.add(peptide)

        if verbose:    # for logging
            i += 1
            if i % 500 == 0:
                logger(f'{i} transcripts processed.')
    
    with open(output_fasta, 'w') as handle:
        SeqIO.write(variant_peptides, handle, 'fasta')
    
    if verbose:
        logger('Variant peptide FASTA file written to disk.')


if __name__ == '__main__':
    args = argparse.Namespace
    args.input_variants = [
        'test/files/downsampled_set/CPCG0100_gencode_v34_moPepGen.txt'
    ]
    args.genome_fasta = 'test/files/downsampled_set/gencode_v34_genome_chr22.fasta'
    args.proteome_fasta = 'test/files/downsampled_set/gencode_v34_translations_chr22.fasta'
    args.annotation_gtf = 'test/files/downsampled_set/gencode_v34_chr22.gtf'
    args.output_fasta = 'test/files/downsampled_set/CPCG0100_gencode_v34_moPepGen.fasta'
    args.verbose = True
    args.cleavage_rule = 'trypsin'
    args.miscleavage = 2
    args.min_mw = 500.
    call_variant_peptides(args)