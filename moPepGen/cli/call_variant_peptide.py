""" `callVariant` is the core of moPepGen. It takes multiple GVF files, generated
by any moPepGen parser, and calls variant peptides caused by genomic variants
using a graph-based algorithm. For any transcript, it creates a three-frame
transcript variant graph by incorporating all variants from any sources (SNV,
INDEL, fusion, alternative splicing, RNA editing, and circRNA). The transcript
variant graph is then translated into a peptide variant graph, followed by
converting to a cleavage graph based on the enzymatic cleavage rule. The
variant peptide graph is than used to call for variant peptides that contains
at least one variant, and do not present in the canonical peptide pool.
 """
from __future__ import annotations
import argparse
from typing import List, Set, TYPE_CHECKING
from pathlib import Path
from moPepGen import svgraph, aa, seqvar, logger, circ
from moPepGen.circ.CircRNA import CircRNAModel
from moPepGen.seqvar import GVFMetadata
from moPepGen.cli.common import add_args_cleavage, add_args_quiet, \
    print_start_message, print_help_if_missing_args, add_args_reference, \
    load_references


if TYPE_CHECKING:
    from moPepGen import dna, gtf

# pylint: disable=W0212
def add_subparser_call_variant(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callPeptides """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='callVariant',
        help='Call non-canonical peptides from genomic variants.',
        description="Genomic variant data must be generated by one of the"
        "moPepGen parser. See moPepGen --help",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        '-i', '--input-variant',
        type=str,
        nargs='+',
        help='Path to input variant files. Must be generated by any of the'
        'moPepGen parser. This can be multiple.',
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '-o', '--output-fasta',
        type=Path,
        help='Filename for the output FASTA.',
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '--max-variants-per-node',
        type=int,
        help='Maximal number of variants per node. This argument can be useful'
        ' when there are local regions that are heavily mutated. When creating'
        ' the cleavage graph, nodes containing variants larger than this value'
        ' are skipped. Setting to -1 will avoid checking for this.',
        default=-1,
        metavar='<number>'
    )
    p.add_argument(
        '--verbose-level',
        type=int,
        help='Level of verbose for logging.',
        default=1,
        metavar='<number>'
    )
    add_args_reference(p)
    add_args_cleavage(p)
    add_args_quiet(p)

    p.set_defaults(func=call_variant_peptide)
    print_help_if_missing_args(p)
    return p

def call_variant_peptide(args:argparse.Namespace) -> None:
    """ Main entry point for calling variant peptide """
    variant_files:List[str] = args.input_variant
    output_fasta:str = args.output_fasta
    quiet:bool = args.quiet
    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    min_length:int = args.min_length
    max_length:int = args.max_length
    max_variants_per_node:int = args.max_variants_per_node
    verbose = args.verbose_level
    if quiet is True:
        verbose = 0

    print_start_message(args)

    genome, anno, _, canonical_peptides = load_references(args=args)

    variant_pool = seqvar.VariantRecordPool()
    circ_rna_pool:List[CircRNAModel] = []

    for file in variant_files:
        with open(file, 'rt') as handle:
            metadata = GVFMetadata.parse(handle)
            if metadata.is_circ_rna():
                for record in circ.io.parse(handle):
                    circ_rna_pool.append(record)
            else:
                variant_pool.load_variants(handle, anno, genome)

        if verbose >= 1:
            logger(f'Variant file {file} loaded.')

    for tx_model in anno.transcripts.values():
        tx_model.remove_cached_seq()

    variant_pool.sort()
    if verbose > 1:
        logger('Variant records sorted.')

    variant_peptides = aa.VariantPeptidePool()

    i = 0
    for tx_id in variant_pool.transcriptional:
        if verbose >= 2:
            logger(tx_id)
        try:
            peptides = call_peptide_main(
                variant_pool=variant_pool, tx_id=tx_id, anno=anno,
                genome=genome, rule=rule, exception=exception,
                miscleavage=miscleavage,
                max_variants_per_node=max_variants_per_node
            )
        except:
            logger(f'Exception raised from {tx_id}')
            raise

        for peptide in peptides:
            variant_peptides.add_peptide(peptide, canonical_peptides, min_mw,
                min_length, max_length)

        if verbose >= 1:
            i += 1
            if i % 1000 == 0:
                logger(f'{i} transcripts processed.')

    for variant in variant_pool.fusion:
        if verbose >= 2:
            logger(variant.id)
        try:
            peptides = call_peptide_fusion(
                variant=variant, variant_pool=variant_pool, anno=anno,
                genome=genome, rule=rule, exception=exception,
                miscleavage=miscleavage,
                max_variants_per_node=max_variants_per_node
            )
        except:
            logger(f"Exception raised from fusion {variant.id}")
            raise

        for peptide in peptides:
            variant_peptides.add_peptide(peptide, canonical_peptides, min_mw,
                min_length, max_length)

    if variant_pool.fusion and verbose >= 1:
        logger('Fusion processed.')

    for circ_rna in circ_rna_pool:
        try:
            peptides = call_peptide_circ_rna(
                record=circ_rna, annotation=anno, genome=genome,
                variant_pool=variant_pool, rule=rule, exception=exception,
                miscleavage=miscleavage,
                max_variants_per_node=max_variants_per_node
            )
        except:
            logger(f"Exception raised from {circ_rna.id}")
            raise

        for peptide in peptides:
            variant_peptides.add_peptide(peptide, canonical_peptides, min_mw,
                min_length, max_length)

    if circ_rna_pool and verbose >= 1:
        logger('circRNA processed')

    variant_peptides.write(output_fasta)

    if verbose >= 1:
        logger('Variant peptide FASTA file written to disk.')

def call_peptide_main(variant_pool:seqvar.VariantRecordPool,
        tx_id:str, anno:gtf.GenomicAnnotation,
        genome:dna.DNASeqDict, rule:str, exception:str, miscleavage:int,
        max_variants_per_node:int) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for main variants (except cirRNA). """
    tx_variants = variant_pool.transcriptional[tx_id]
    tx_model = anno.transcripts[tx_id]
    chrom = tx_model.transcript.location.seqname
    transcript_seq = tx_model.get_transcript_sequence(genome[chrom])

    dgraph = svgraph.ThreeFrameTVG(
        seq=transcript_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        max_variants_per_node=max_variants_per_node
    )
    dgraph.init_three_frames()
    dgraph.create_variant_graph(tx_variants, variant_pool, genome, anno)
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()

    pgraph.create_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)

def call_peptide_fusion(variant:seqvar.VariantRecord,
    variant_pool:seqvar.VariantRecordPool, anno:gtf.GenomicAnnotation,
    genome:dna.DNASeqDict, rule:str, exception:str, miscleavage:int,
    max_variants_per_node:int) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for fusion """
    tx_id = variant.location.seqname
    tx_model = anno.transcripts[tx_id]
    chrom = tx_model.transcript.chrom
    tx_seq = tx_model.get_transcript_sequence(genome[chrom])
    tx_seq = tx_seq[:variant.location.end]

    if tx_id in variant_pool.transcriptional:
        tx_variants = [x for x in variant_pool.transcriptional[tx_id]
            if x.location.end < variant.location.end]
    else:
        tx_variants = []

    tx_variants.append(variant)

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=tx_model.is_protein_coding,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        max_variants_per_node=max_variants_per_node
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
        exception:str, miscleavage:int, max_variants_per_node:int
        )-> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides from a given circRNA """
    gene_id = record.gene_id
    gene_model = annotation.genes[gene_id]
    chrom = gene_model.location.seqname
    gene_seq = gene_model.get_gene_sequence(genome[chrom])
    circ_seq = record.get_circ_rna_sequence(gene_seq)

    # Alternative splicing should not be included. Alternative splicing
    # represented as Insertion, Deletion or Substitution.
    exclusion_variant_types = ['Insertion', 'Deletion', 'Substitution']

    # TODO(Trevor): we should leave out alternative splicing and only keep
    # SNV and INDEL. Please point it out if you see this!
    variant_records = variant_pool.filter_variants(
        gene_id=gene_id, anno=annotation, exclude_type=exclusion_variant_types,
        intron=False, segments=record.fragments
    )

    cgraph = svgraph.ThreeFrameCVG(
        circ_seq, _id=record.id, circ_record=record,
        max_variants_per_node=max_variants_per_node
    )
    cgraph.init_three_frames()
    cgraph.create_variant_circ_graph(variant_records)
    cgraph.extend_loop()
    cgraph.truncate_three_frames()
    cgraph.fit_into_codons()
    pgraph = cgraph.translate()
    pgraph.create_cleavage_graph(rule=rule, exception=exception)
    return pgraph.call_variant_peptides(miscleavage=miscleavage)
