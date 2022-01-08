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
from moPepGen.cli.common import add_args_cleavage, add_args_quiet, \
    print_start_message, print_help_if_missing_args, add_args_reference, \
    load_references


if TYPE_CHECKING:
    from moPepGen import dna, gtf
    from moPepGen.circ.CircRNA import CircRNAModel

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
        type=Path,
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
        '--noncanonical-transcripts',
        action='store_true',
        help='Process only noncanonical transcripts of fusion transcripts and'
        'circRNA. Canonical transcripts are skipped.'
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

class VariantPeptideCaller():
    """ Helper class to call variant peptides """
    def __init__(self, args:argparse.Namespace):
        """ constructor """
        self.args = args
        self.variant_files:List[str] = args.input_variant
        self.output_fasta:str = args.output_fasta
        self.quiet:bool = args.quiet
        self.rule:str = args.cleavage_rule
        self.miscleavage:int = int(args.miscleavage)
        self.min_mw:float = float(args.min_mw)
        self.exception = 'trypsin_exception' if self.rule == 'trypsin' else None
        self.min_length:int = args.min_length
        self.max_length:int = args.max_length
        self.max_variants_per_node:int = args.max_variants_per_node
        self.verbose = args.verbose_level
        if self.quiet is True:
            self.verbose = 0
        self.genome:dna.DNASeqDict = None
        self.anno:gtf.GenomicAnnotation = None
        self.canonical_peptides:Set[str] = None
        self.variant_pool:seqvar.VariantRecordPoolOnDisk = None
        self.circ_rna_pool:List[CircRNAModel] = []
        self.variant_peptides = aa.VariantPeptidePool()

    def load_reference(self):
        """ load reference genome, annotation, and canonical peptides """
        self.genome, self.anno, _, self.canonical_peptides = load_references(self.args)

    def remove_cached_sequences(self):
        """ remove chaced transcript sequences to reduce memory usage. """
        for tx_model in self.anno.transcripts.values():
            tx_model.remove_cached_seq()

    def create_in_disk_vairant_pool(self):
        """ Create in disk variant pool """
        self.variant_pool = seqvar.VariantRecordPoolOnDisk(
            pointers=None, gvf_files=self.variant_files,
            anno=self.anno, genome=self.genome
        )

    def call_variant_peptides(self):
        """ call variant peptides """
        with seqvar.VariantRecordPoolOnDisk(pointers=None, gvf_files=self.variant_files,
                anno=self.anno, genome=self.genome) as pool:
            tx_rank = self.anno.get_transcirpt_rank()
            sorted_key = sorted(pool.pointers.keys(), key=lambda x:tx_rank[x])
            i = 0
            for key in sorted_key:
                if self.verbose >= 2:
                    logger(key)
                series = pool[key]

                if series.transcriptional:
                    self.call_variants_main(key, series.transcriptional, pool)

                for fusion_variant in series.fusion:
                    self.call_variants_fusion(fusion_variant, pool)

                for circ_model in series.circ_rna:
                    self.call_variants_circ_rna(circ_model, pool)

                self.anno.remove_cached_tx_seq()

                if self.verbose >= 1:
                    i += 1
                    if i % 1000 == 0:
                        logger(f'{i} transcripts processed.')

    def call_variants_main(self, tx_id:str, tx_variants:List[seqvar.VariantRecord],
            pool:seqvar.VariantRecordPoolOnDisk):
        """ main variant peptide caller """
        try:
            peptides = call_peptide_main(
                tx_id=tx_id, tx_variants=tx_variants, variant_pool=pool,
                anno=self.anno, genome=self.genome, rule=self.rule,
                exception=self.exception, miscleavage=self.miscleavage,
                max_variants_per_node=self.max_variants_per_node
            )
        except:
            logger(f'Exception raised from {tx_id}')
            raise

        for peptide in peptides:
            self.variant_peptides.add_peptide(
                peptide=peptide, canonical_peptides=self.canonical_peptides,
                min_mw=self.min_mw, min_length=self.min_length,
                max_length=self.max_length
            )

    def call_variants_fusion(self, variant:seqvar.VariantRecord,
            pool:seqvar.VariantRecordPoolOnDisk):
        """ call variant peptides from fusion transcripts """
        if self.verbose >= 2:
            logger(variant.id)
        try:
            peptides = call_peptide_fusion(
                variant=variant, variant_pool=pool,
                anno=self.anno, genome=self.genome, rule=self.rule,
                exception=self.exception, miscleavage=self.miscleavage,
                max_variants_per_node=self.max_variants_per_node
            )
        except:
            logger(f"Exception raised from fusion {variant.id}")
            raise

        for peptide in peptides:
            self.variant_peptides.add_peptide(
                peptide=peptide, canonical_peptides=self.canonical_peptides,
                min_mw=self.min_mw, min_length=self.min_length,
                max_length=self.max_length
            )

    def call_variants_circ_rna(self, circ_model:CircRNAModel,
            pool:seqvar.VariantRecordPoolOnDisk):
        """ call variant peptides from circRNA """
        try:
            peptides = call_peptide_circ_rna(
                record=circ_model, annotation=self.anno, genome=self.genome,
                variant_pool=pool, rule=self.rule,
                exception=self.exception, miscleavage=self.miscleavage,
                max_variants_per_node=self.max_variants_per_node
            )
        except:
            logger(f"Exception raised from {circ_model.id}")
            raise

        for peptide in peptides:
            self.variant_peptides.add_peptide(
                peptide=peptide, canonical_peptides=self.canonical_peptides,
                min_mw=self.min_mw,min_length=self.min_length,
                max_length=self.max_length
            )
        if self.circ_rna_pool and self.verbose >= 1:
            logger('circRNA processed.')

    def write_fasta(self):
        """ write variant peptides to fasta. """
        self.variant_peptides.write(self.output_fasta)
        if self.verbose >= 1:
            logger('Variant peptide FASTA file written to disk.')

def call_variant_peptide(args:argparse.Namespace) -> None:
    """ Main entry point for calling variant peptide """
    caller = VariantPeptideCaller(args)

    print_start_message(args)
    caller.load_reference()
    caller.call_variant_peptides()
    caller.write_fasta()

def call_peptide_main(tx_id:str, tx_variants:List[seqvar.VariantRecord],
        variant_pool:seqvar.VariantRecordPoolOnDisk, anno:gtf.GenomicAnnotation,
        genome:dna.DNASeqDict, rule:str, exception:str, miscleavage:int,
        max_variants_per_node:int) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for main variants (except cirRNA). """
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
    variant_pool:seqvar.VariantRecordPoolOnDisk, anno:gtf.GenomicAnnotation,
    genome:dna.DNASeqDict, rule:str, exception:str, miscleavage:int,
    max_variants_per_node:int) -> Set[aa.AminoAcidSeqRecord]:
    """ Call variant peptides for fusion """
    tx_id = variant.location.seqname
    tx_model = anno.transcripts[tx_id]
    chrom = tx_model.transcript.chrom
    tx_seq = tx_model.get_transcript_sequence(genome[chrom])
    tx_seq = tx_seq[:variant.location.end]

    if tx_id in variant_pool:
        tx_variants = [x for x in variant_pool[tx_id].transcriptional
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

    variant_records = variant_pool.filter_variants(
        gene_id=gene_id, exclude_type=exclusion_variant_types,
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
