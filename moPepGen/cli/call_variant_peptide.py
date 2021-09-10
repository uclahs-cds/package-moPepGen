""" Module for calling variant paptide """
from __future__ import annotations
import argparse
from typing import List, Set, TYPE_CHECKING
from pathlib import Path
from moPepGen import svgraph, aa, seqvar, logger, circ
from .common import add_args_cleavage, add_args_verbose, print_start_message, \
    print_help_if_missing_args, add_args_reference, load_references


if TYPE_CHECKING:
    from moPepGen import dna, gtf

# pylint: disable=W0212
def add_subparser_call_variant(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callPeptides """
    p = subparsers.add_parser(
        name='callVariant',
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
        type=Path,
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

    genome, annotation, _, canonical_peptides = load_references(args=args)

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
        tx_id:str, anno:gtf.GenomicAnnotation,
        genome:dna.DNASeqDict, rule:str, exception:str, miscleavage:int
        ) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for main variants (except cirRNA). """
    tx_variants = variant_pool.transcriptional[tx_id]
    tx_model = anno.transcripts[tx_id]
    chrom = tx_model.transcript.location.seqname
    transcript_seq = tx_model.get_transcript_sequence(genome[chrom])
    cds_start_nf = 'tag' in tx_model.transcript.attributes and \
        'cds_start_NF' in tx_model.transcript.attributes['tag']

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
            cur = dgraph.apply_fusion(cur, variant, variant_pool, genome, anno)
            variant = next(variant_iter, None)
            continue

        if variant.type in 'Insertion':
            cur = dgraph.apply_insertion(cur, variant, variant_pool, genome,
                anno)
            variant = next(variant_iter, None)
            continue

        if variant.type == 'Substitution':
            cur = dgraph.apply_substitution(cur, variant, variant_pool,
                genome, anno)
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

    # Alternative splicing should not be included. Alternative splicing
    # represented as Insertion, Deletion or Substitution.
    exclusion_variant_types = ['Insertion', 'Deletion', 'Substitution']

    variant_records = variant_pool.filter_variants(gene_id, annotation, genome,
        exclusion_variant_types, intron=False, segments=record.fragments)

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

if __name__ == '__main__':
    _args = argparse.Namespace()
    _args.command = 'callPeptides'
    _args.input_variant = [
        'test/files/vep/CPCG0100_gencode_aa_snv_ENST00000588049.5.gvf'
    ]
    _args.index_dir = 'test/files/gencode_34_index'
    _args.genome_fasta = None
    _args.annotation_gtf = None
    _args.proteome_fasta = None
    _args.circ_rna_bed = None
    _args.output_fasta = 'test/files/vep/CPCG0100_gencode_aa_snv_ENST00000588049.5.fasta'
    _args.verbose = True
    _args.cleavage_rule = 'trypsin'
    _args.miscleavage = 2
    _args.min_mw = 500.
    _args.min_length = 7
    _args.max_length = 25
    call_variant_peptide(args=_args)
