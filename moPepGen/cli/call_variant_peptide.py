""" Module for call variant paptide """
from __future__ import annotations
from typing import List, Dict, Set
import argparse
import pickle
from Bio import SeqUtils
from Bio.SeqIO import FastaIO
from Bio.Seq import Seq
from moPepGen import svgraph, dna, gtf, aa, seqvar, logger, get_equivalent, CircRNA
from moPepGen.SeqFeature import FeatureLocation


def call_variant_peptide(args:argparse.Namespace) -> None:
    """ Main entry point for calling variant peptide """
    variant_files:List[str] = args.input_variant
    circ_rna_bed:str = args.circ_rna_bed
    output_fasta:str = args.output_fasta
    index_dir:str = args.index_dir
    verbose:bool = args.verbose
    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    min_length:int = args.min_length
    max_length:int = args.max_length

    if verbose:
        logger('moPepGen callPeptide started.')

    if index_dir:
        with open(f'{index_dir}/genome.pickle', 'rb') as handle:
            genome = pickle.load(handle)

        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            annotation = pickle.load(handle)

        with open(f"{index_dir}/canonical_peptides.pickle", 'rb') as handle:
            canonical_peptides = pickle.load(handle)
    else:
        genome_fasta:str = args.genome_fasta
        proteome_fasta:str = args.proteome_fasta
        annotation_gtf:str = args.annotation_gtf

        annotation = gtf.GenomicAnnotation()
        annotation.dump_gtf(annotation_gtf)
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

        canonical_peptides = proteome.create_unique_peptide_pool(
            rule=rule, exception=exception, miscleavage=miscleavage,
            min_mw=min_mw
        )
        if verbose:
            logger('canonical peptide pool generated.')

    variants:Dict[str, List[seqvar.VariantRecord]] = {}
    for file in variant_files:
        for record in seqvar.io.parse(file):
            transcript_id = record.location.seqname
            if transcript_id not in variants:
                variants[transcript_id] = [record]
            else:
                variants[transcript_id].append(record)
        if verbose:
            logger(f'Variant file {file} loaded.')

    for records in variants.values():
        records.sort()
    if verbose:
        logger('Variant records sorted.')

    variant_peptides = set()

    i = 0
    for transcript_id in variants:

        try:
            peptides = call_peptide_main(variants, transcript_id, annotation,
                genome, rule, exception, miscleavage)
        except:
            logger(f'Exception raised from {transcript_id}')
            raise

        for peptide in peptides:
            if SeqUtils.molecular_weight(peptide.seq, 'protein') < min_mw:
                continue
            if len(peptide.seq) < min_length or len(peptide.seq) > max_length:
                continue
            if str(peptide.seq) in canonical_peptides:
                continue
            same_peptide = get_equivalent(variant_peptides, peptide)
            if same_peptide:
                same_peptide:Seq
                new_label = peptide.id
                same_peptide.id += ('||' + new_label)
                same_peptide.name = same_peptide.id
                same_peptide.description = same_peptide.id
            else:
                variant_peptides.add(peptide)

        if verbose:    # for logging
            i += 1
            if i % 1000 == 0:
                logger(f'{i} transcripts processed.')

    if circ_rna_bed:
        for circ in CircRNA.parse(circ_rna_bed):
            peptides = call_peptide_circ_rna(circ, annotation, genome,
                variants, rule, exception, miscleavage)

        for peptide in peptides:
            if SeqUtils.molecular_weight(peptide.seq, 'protein') < min_mw:
                continue
            if len(peptide.seq) < min_length or len(peptide.seq) > max_length:
                continue
            if str(peptide.seq) in canonical_peptides:
                continue
            same_peptide = get_equivalent(variant_peptides, peptide)
            if same_peptide:
                new_label = peptide.id
                same_peptide.id += ('||' + new_label)
                same_peptide.name = same_peptide.id
                same_peptide.description = same_peptide.id
            else:
                variant_peptides.add(peptide)

    with open(output_fasta, 'w') as handle:
        writer = FastaIO.FastaWriter(handle, record2title=lambda x: x.description)
        for record in variant_peptides:
            writer.write_record(record)

    if verbose:
        logger('Variant peptide FASTA file written to disk.')

def call_peptide_main(variants:Dict[str, List[seqvar.VariantRecord]],
        transcript_id:str, annotation:gtf.GenomicAnnotation,
        genome:dna.DNASeqDict, rule:str, exception:str, miscleavage:int
        ) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for main variants (except cirRNA). """
    variant_records = variants[transcript_id]
    anno = annotation.transcripts[transcript_id]
    chrom = anno.transcript.location.seqname
    transcript_seq = anno.get_transcript_sequence(genome[chrom])

    dgraph = svgraph.TranscriptVariantGraph(
        seq=transcript_seq,
        _id=transcript_id
    )

    ## Create transcript variant graph
    # dgraph.create_variant_graph(variant_records)
    dgraph.add_null_root()
    variant_iter = iter(variant_records)
    variant = next(variant_iter, None)
    cur = dgraph.root.get_reference_next()
    while variant:
        if cur.seq.locations[0].ref.start > variant.location.start:
            variant = next(variant_iter, None)
            continue

        if cur.seq.locations[-1].ref.end <= variant.location.start:
            cur = cur.get_reference_next()
            continue

        if variant.type == 'Fusion':
            donor_transcript_id = variant.attrs['DONOR_TRANSCRIPT_ID']
            donor_anno = annotation.transcripts[donor_transcript_id]
            donor_chrom = donor_anno.transcript.location.seqname
            donor_seq = donor_anno.get_transcript_sequence(genome[donor_chrom])

            donor_variant_records = []
            if donor_transcript_id in variants:
                for record in variants[donor_transcript_id]:
                    # If the donor sequence has another fusion event, it is not
                    # considered. The chance of two fusion events happened
                    # together should be low.
                    if record.type == 'Fusion':
                        continue
                    if record.location.start > int(variant.attrs['DONOR_POS']):
                        donor_variant_records.append(record)

            cur = dgraph.apply_fusion(cur, variant, donor_seq, donor_variant_records)
            variant = next(variant_iter, None)
            continue

        if variant.type in 'Insertion':
            if variant.attrs['COORDINATE'] == 'gene':
                gene_id = variant.attrs['GENE_ID']
                gene_model = annotation.genes[gene_id]
                chrom = gene_model.chrom
                donor_start = variant.get_donor_start()
                donor_end = variant.get_donor_end()
                gene_seq = gene_model.get_gene_sequence(genome[chrom])
                insert_seq = gene_seq[donor_start:donor_end]
                exclude_type = ['Insertion', 'Deletion', 'Substitution', 'Fusion']
                insert_variants = find_gene_variants(gene_id, annotation,
                    variants, donor_start, donor_end, exclude_type)
                cur = dgraph.apply_insertion(cur, variant, insert_seq, insert_variants)
                variant = next(variant_iter, None)
                continue

        if variant.type == 'Substitution':
            if variant.attrs['COORDINATE'] == 'gene':
                gene_id = variant.attrs['GENE_ID']
                gene_model = annotation.genes[gene_id]
                chrom = gene_model.chrom
                donor_start = variant.get_donor_start()
                donor_end = variant.get_donor_end()
                gene_seq = gene_model.get_gene_sequence(genome[chrom])
                sub_seq = gene_seq[donor_start:donor_end]
                exclude_type = ['Insertion', 'Deletion', 'Substitution', 'Fusion']
                sub_variants = find_gene_variants(gene_id, annotation,
                    variants, donor_start, donor_end, exclude_type)
                cur = dgraph.apply_substitution(cur, variant, sub_seq, sub_variants)
                variant = next(variant_iter, None)
                continue

        if variant.type == 'Deletion':
            variant = seqvar.VariantRecord(
                location=FeatureLocation(
                    seqname=variant.location.seqname,
                    start=variant.location.start,
                    end=int(variant.attrs['END'])
                ),
                ref=variant.ref,
                alt=variant.alt,
                _type=variant.type,
                _id=variant.id,
                attrs=variant.attrs
            )
            variant = next(variant_iter, None)
            continue

        cur = dgraph.apply_variant(cur, variant)
        if len(cur.in_edges) == 0:
            dgraph.root = cur
        variant = next(variant_iter, None)

    dgraph.find_all_orfs()
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()

    pgraph.form_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)


def call_peptide_circ_rna(circ:CircRNA.CircRNAModel,
        annotation:gtf.GenomicAnnotation, genome:dna.DNASeqDict,
        variants:Dict[str,List[seqvar.VariantRecord]], rule:str,
        exception:str, miscleavage:int)-> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides from a given circRNA """
    gene_id = circ.gene_id
    gene_model = annotation.genes[gene_id]
    chrom = gene_model.location.seqname
    gene_seq = gene_model.get_gene_sequence(genome[chrom])
    circ_seq = circ.get_circ_rna_sequence(gene_seq)

    variant_records = set()
    for transcript_id in circ.transcript_ids:
        if transcript_id not in variants:
            continue
        for record in variants[transcript_id]:
            record = annotation.variant_coordinates_to_gene(record, gene_id)
            # Alternative splicing should not be included. Alternative splicing
            # are represented as Insertion, Deletion or Substitution.
            if record.type in ['Insertion', 'Deletion', 'Substitution']:
                continue
            for feature in circ.fragments:
                if feature.location.is_superset(record.location):
                    variant_records.add(record)
                    break
    variant_records = list(variant_records)
    variant_records.sort()

    cgraph = svgraph.CircularVariantGraph(circ_seq, _id=circ.id)

    cgraph.create_variant_graph(variant_records)
    cgraph.align_all_variants()
    tgraph = cgraph.find_all_orfs()
    tgraph.fit_into_codons()
    pgraph = tgraph.translate()
    pgraph.form_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)

def find_gene_variants(gene_id:str, annotation:gtf.GenomicAnnotation,
        variants:Dict[str,seqvar.VariantRecord], start:int, end:int,
        exclude_type:List[str]) -> List[seqvar.VariantRecord]:
    """ Find all unique variants of a gene within a given range. Fusions
    are skipped. """
    records = set()
    for transcript_id in annotation.genes[gene_id].transcripts:
        if transcript_id not in variants:
            continue
        record:seqvar.VariantRecord
        for record in variants[transcript_id]:
            if record.type in exclude_type:
                continue
            record = annotation.variant_coordinates_to_gene(record, gene_id)
            if record.location.start > start and record.location.end < end:
                records.add(record)
    records = list(records)
    records.sort()
    return records

if __name__ == '__main__':
    test_args = argparse.Namespace()
    test_args.input_variant = [
        'test/files/CPCG0100_gencode_aa_indel_ENST00000558600.1.tvf'
    ]
    test_args.index_dir = 'test/files/downsampled_index/ENST00000558600.1'
    test_args.circ_rna_bed = None
    test_args.output_fasta = 'test/files/CPCG0100_gencode_aa_indel_ENST00000558600.1.fasta'
    test_args.verbose = True
    test_args.cleavage_rule = 'trypsin'
    test_args.miscleavage = 2
    test_args.min_mw = 500.
    test_args.min_length = 7
    test_args.max_length = 25
    call_variant_peptide(args=test_args)
