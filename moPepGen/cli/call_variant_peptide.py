""" Module for calling variant paptide """
from __future__ import annotations
import argparse
from typing import List, Set, TYPE_CHECKING
from pathlib import Path
from moPepGen import svgraph, aa, seqvar, logger, circ
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.seqvar import GVFMetadata
from moPepGen.cli.common import add_args_cleavage, add_args_verbose, \
    print_start_message, print_help_if_missing_args, add_args_reference, \
    load_references


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
    output_fasta:str = args.output_fasta
    verbose:bool = args.verbose
    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    min_length:int = args.min_length
    max_length:int = args.max_length

    print_start_message(args)

    genome, anno, _, canonical_peptides = load_references(args=args)

    variant_pool = seqvar.VariantRecordPool()
    circ_rna_pool = []

    for file in variant_files:
        with open(file, 'rt') as handle:
            metadata = GVFMetadata.parse(handle)
            if metadata.is_circ_rna():
                for record in circ.io.parse(handle):
                    circ_rna_pool.append(record)
            else:
                variant_pool.load_variants(handle, anno, genome)

        if verbose:
            logger(f'Variant file {file} loaded.')

    variant_pool.sort()
    if verbose:
        logger('Variant records sorted.')

    variant_peptides = aa.VariantPeptidePool()

    i = 0
    for tx_id in variant_pool.transcriptional:

        try:
            peptides = call_peptide_main(variant_pool, tx_id, anno,
                genome, rule, exception, miscleavage)
        except:
            logger(f'Exception raised from {tx_id}')
            raise

        for peptide in peptides:
            variant_peptides.add_peptide(peptide, canonical_peptides, min_mw,
                min_length, max_length)

        if verbose:
            i += 1
            if i % 1000 == 0:
                logger(f'{i} transcripts processed.')

    for circ_rna in circ_rna_pool:
        peptides = call_peptide_circ_rna(circ_rna, anno, genome,
            variant_pool, rule, exception, miscleavage)

        for peptide in peptides:
            variant_peptides.add_peptide(peptide, canonical_peptides, min_mw,
                min_length, max_length)
    if circ_rna_pool:
        if verbose:
            logger('circRNA processed')

    variant_peptides.write(output_fasta)

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

    start = transcript_seq.orf.start
    start_codon = FeatureLocation(start=start, end=start+3)

    has_start_altering = any(x.location.overlaps(start_codon) for x in
        tx_variants)

    dgraph = svgraph.ThreeFrameTVG(
        seq=transcript_seq,
        _id=tx_id,
        cds_start_nf=cds_start_nf,
        has_known_orf=tx_model.is_protein_coding(),
        has_start_altering=has_start_altering and not cds_start_nf
    )
    dgraph.init_three_frames()
    dgraph.create_variant_graph(tx_variants, variant_pool, genome, anno)
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()

    pgraph.create_cleavage_graph(rule=rule, exception=exception)
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

    cgraph = svgraph.ThreeFrameCVG(
        circ_seq, _id=record.id, circ_record=record
    )
    cgraph.init_three_frames()
    cgraph.create_variant_circ_graph(variant_records)
    cgraph.extend_loop()
    cgraph.truncate_three_frames()
    cgraph.fit_into_codons()
    pgraph = cgraph.translate()
    pgraph.create_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)
