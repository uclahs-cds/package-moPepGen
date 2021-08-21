""" Module for call variant paptide """
from __future__ import annotations
from typing import List, Dict, Set
import argparse
from moPepGen import svgraph, dna, gtf, aa, seqvar, logger, circ
from moPepGen.SeqFeature import FeatureLocation
from .common import add_args_cleavage, add_args_verbose, print_start_message, \
    print_help_if_missing_args, add_args_reference, load_references


# pylint: disable=W0212
def add_subparser_call_peptides(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callPeptides """
    p = subparsers.add_parser(
        name='callPeptide',
        help='Call non-canonical peptides from genomic variants.',
        description="Genomic variant data must be generated by one of the"
        "moPepGen parser. See moPepGen --help"
    )

    p.add_argument(
        '-i', '--input-variant',
        type=str,
        nargs='+',
        help='Path to input variant files. Must be generated by any of the'
        'moPepGen parser. This can be multiple.',
        metavar='',
        required=True
    )
    p.add_argument(
        '-b', '--circ-rna-bed',
        type=str,
        help="Path to the circRNA BED file. Must be generated by moPepGen's"
        "caller. Defaults to None.",
        metavar='',
        default=None,
        dest='circ_rna_bed'
    )
    p.add_argument(
        '-o', '--output-fasta',
        type=str,
        help='Filename for the output FASTA.',
        metavar='',
        required=True
    )

    add_args_reference(p)
    add_args_cleavage(p)
    add_args_verbose(p)

    p.set_defaults(func=call_variant_peptide)
    print_help_if_missing_args(p)

def call_variant_peptide(args:argparse.Namespace) -> None:
    """ Main entry point for calling variant peptide """
    variant_files:List[str] = args.input_variant
    circ_rna_bed:str = args.circ_rna_bed
    output_fasta:str = args.output_fasta
    verbose:bool = args.verbose
    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    min_length:int = args.min_length
    max_length:int = args.max_length

    print_start_message(args)

    genome, annotation, canonical_peptides = load_references(args=args)

    variants = seqvar.VariantRecordPool.load_variants(variant_files,
        annotation, genome, verbose)

    variant_pool = aa.VariantPeptidePool()

    i = 0
    for tx_id in variants.transcriptional:

        try:
            peptides = call_peptide_main(variants, tx_id, annotation,
                genome, rule, exception, miscleavage)
        except:
            logger(f'Exception raised from {tx_id}')
            raise

        for peptide in peptides:
            variant_pool.add_peptide(peptide, canonical_peptides, min_mw,
                min_length, max_length)

        if verbose:
            i += 1
            if i % 1000 == 0:
                logger(f'{i} transcripts processed.')

    if circ_rna_bed:
        for record in circ.io.parse(circ_rna_bed):
            peptides = call_peptide_circ_rna(record, annotation, genome,
                variants, rule, exception, miscleavage)

        for peptide in peptides:
            variant_pool.add_peptide(peptide, canonical_peptides, min_mw,
                min_length, max_length)

    variant_pool.write(output_fasta)

    if verbose:
        logger('Variant peptide FASTA file written to disk.')

def call_peptide_main(variant_pool:seqvar.VariantRecordPool,
        tx_id:str, annotation:gtf.GenomicAnnotation,
        genome:dna.DNASeqDict, rule:str, exception:str, miscleavage:int
        ) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for main variants (except cirRNA). """
    tx_variants = variant_pool.transcriptional[tx_id]
    anno = annotation.transcripts[tx_id]
    chrom = anno.transcript.location.seqname
    transcript_seq = anno.get_transcript_sequence(genome[chrom])
    cds_start_nf = 'tag' in anno.transcript.attributes and \
        'cds_start_NF' in anno.transcript.attributes['tag']

    dgraph = svgraph.TranscriptVariantGraph(
        seq=transcript_seq,
        _id=tx_id,
        cds_start_nf=cds_start_nf
    )

    ## Create transcript variant graph
    # dgraph.create_variant_graph(variant_records)
    dgraph.add_null_root()
    variant_iter = iter(tx_variants)
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
            accepter_transcript_id = variant.attrs['ACCEPTER_TRANSCRIPT_ID']
            accepter_anno = annotation.transcripts[accepter_transcript_id]
            accepter_chrom = accepter_anno.transcript.location.seqname
            accepter_seq = accepter_anno.get_transcript_sequence(genome[accepter_chrom])

            accepter_variant_records = []
            if accepter_transcript_id in variant_pool.transcriptional:
                for record in variant_pool.transcriptional[accepter_transcript_id]:
                    # If the donor sequence has another fusion event, it is not
                    # considered. The chance of two fusion events happened
                    # together should be low.
                    if record.type == 'Fusion':
                        continue
                    if record.location.start > variant.get_accepter_position():
                        accepter_variant_records.append(record)

            cur = dgraph.apply_fusion(cur, variant, accepter_seq, accepter_variant_records)
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
                    variant_pool, donor_start, donor_end, exclude_type)
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
                    variant_pool, donor_start, donor_end, exclude_type)
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

def call_peptide_circ_rna(record:circ.CircRNAModel,
        annotation:gtf.GenomicAnnotation, genome:dna.DNASeqDict,
        variant_pool:seqvar.VariantRecordPool, rule:str,
        exception:str, miscleavage:int)-> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides from a given circRNA """
    gene_id = record.gene_id
    gene_model = annotation.genes[gene_id]
    chrom = gene_model.location.seqname
    gene_seq = gene_model.get_gene_sequence(genome[chrom])
    circ_seq = record.get_circ_rna_sequence(gene_seq)

    variant_records = set()
    for tx_id in record.transcript_ids:
        if tx_id not in variant_pool.transcriptional:
            continue
        for variant in variant_pool.transcriptional[tx_id]:
            variant = annotation.variant_coordinates_to_gene(variant, gene_id)
            # Alternative splicing should not be included. Alternative splicing
            # are represented as Insertion, Deletion or Substitution.
            if variant.type in ['Insertion', 'Deletion', 'Substitution']:
                continue
            for fragment in record.fragments:
                if fragment.location.is_superset(variant.location):
                    variant_records.add(variant)
                    break
    variant_records = list(variant_records)
    variant_records.sort()

    cgraph = svgraph.CircularVariantGraph(
        circ_seq, _id=record.id, circ_record=record
    )

    cgraph.create_variant_graph(variant_records)
    cgraph.align_all_variants()
    tgraph = cgraph.find_all_orfs()
    tgraph.fit_into_codons()
    pgraph = tgraph.translate()
    pgraph.form_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)

def find_gene_variants(gene_id:str, annotation:gtf.GenomicAnnotation,
        variant_pool:seqvar.VariantRecordPool, start:int, end:int,
        exclude_type:List[str]) -> List[seqvar.VariantRecord]:
    """ Find all unique variants of a gene within a given range. Fusions
    are skipped. """
    records = set()
    for transcript_id in annotation.genes[gene_id].transcripts:
        if transcript_id not in variant_pool.transcriptional:
            continue
        record:seqvar.VariantRecord
        for record in variant_pool.transcriptional[transcript_id]:
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
        'test/files/vep/CPCG0100_gencode_aa_indel_ENST00000515757.5.tvf'
    ]
    test_args.index_dir = 'test/files/downsampled_index/ENST00000515757.5'
    test_args.circ_rna_bed = None
    test_args.output_fasta = 'test/files/vep/CPCG0100_gencode_aa_indel_ENST00000515757.5.fasta'
    test_args.verbose = True
    test_args.cleavage_rule = 'trypsin'
    test_args.miscleavage = 2
    test_args.min_mw = 500.
    test_args.min_length = 7
    test_args.max_length = 25
    call_variant_peptide(args=test_args)
