""" `generateIndex` takes the reference genome FASTA, annotation GTF, and the
translated proteome FASTA file, converts them to the moPepGen objects,
serializes them and saves into disk. The outputted index files also contain the
canonical peptide pool. The index files can then be used in any moPepGen
command. It is recommended to run `generateIndex` before any analysis using
moPepGen to avoid processing the reference files repeatedly and save massive
time. """
from __future__ import annotations
from typing import TYPE_CHECKING
import argparse
from pathlib import Path
import pickle
import lzma
from moPepGen import dna, aa, gtf, logger
from moPepGen.gtf import GtfIO
from moPepGen.cli import common


if TYPE_CHECKING:
    from moPepGen.gtf.GTFPointer import GenePointer, TranscriptPointer

# pylint: disable=W0212
def add_subparser_generate_index(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen generateIndex """
    p = subparsers.add_parser(
        name='generateIndex',
        help='Generate genome and proteome index files for moPepGen',
        description='Generate genome and proteome index files for moPepGen'
        ' parsers and peptide caller.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-o', '--output-dir',
        type=Path,
        help='Ouput directory for index files.',
        metavar='<file>',
        dest='output_dir',
        required=True
    )
    common.add_args_reference(p, index=False)
    common.add_args_cleavage(p)
    common.add_args_quiet(p)
    p.set_defaults(func=generate_index)
    common.print_help_if_missing_args(p)
    return p


def index_gtf(file:Path, source:str=None, proteome:aa.AminoAcidSeqDict=None,
        invalid_protein_as_noncoding:bool=True):
    """"""
    anno = gtf.GenomicAnnotationOnDisk()
    with open(file, 'rb') as handle:
        anno.generate_index(handle, source)

    if proteome:
        anno.check_protein_coding(proteome, invalid_protein_as_noncoding)

    gene_idx_file, tx_idx_file = anno.get_index_files(file)

    with open(gene_idx_file, 'wt') as handle:
        handle.write(f"# source={anno.source}\n")
        for gene in anno.genes.keys():
            gene_pointer:GenePointer = anno.genes.get_pointer(gene)
            handle.write(gene_pointer.to_line() + '\n')

    with open(tx_idx_file, 'wt') as handle:
        handle.write(f"# source={anno.source}\n")
        for tx in anno.transcripts.keys():
            tx_pointer:TranscriptPointer = anno.transcripts.get_pointer(tx)
            handle.write(tx_pointer.to_line() + '\n')

def generate_index(args:argparse.Namespace):
    """ Generate  """
    path_genome:Path = args.genome_fasta
    path_gtf:Path = args.annotation_gtf
    parth_proteome:Path = args.proteome_fasta

    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    min_length:int = int(args.min_length)
    max_length:int = int(args.max_length)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    invalid_protein_as_noncoding:bool = args.invalid_protein_as_noncoding
    quiet:bool = args.quiet

    output_dir:Path = args.output_dir
    output_genome = output_dir/"genome.pkl"
    output_proteome = output_dir/"proteome.pkl"
    output_anno = output_dir/"annotation.dat"
    output_peptides = output_dir/"canonical_peptides.pkl"
    output_coding_tx = output_dir/"coding_transcripts.pkl"

    common.print_start_message(args)

    output_dir.mkdir(exist_ok=True)

    genome = dna.DNASeqDict()
    genome.dump_fasta(path_genome)
    if not quiet:
        logger('Genome FASTA loaded')
    with open(output_genome, 'wb') as handle:
        pickle.dump(genome, handle)
    if not quiet:
        logger('Genome FASTA saved to disk.')
    del genome

    anno = gtf.GenomicAnnotation()
    anno.dump_gtf(path_gtf, source=args.reference_source)
    if not quiet:
        logger('Genome annotation GTF loaded.')

    proteome = aa.AminoAcidSeqDict()
    proteome.dump_fasta(parth_proteome, source=args.reference_source)

    with lzma.open(output_anno, 'wt') as handle:
        GtfIO.write(handle, anno)
    if not quiet:
        logger('Genome annotation GTF saved to disk.')

    if not quiet:
        logger('Proteome FASTA loaded.')
    with open(output_proteome, 'wb') as handle:
        pickle.dump(proteome, handle)
    if not quiet:
        logger('Proteome FASTA saved to disk.')

    # canoincal peptide pool
    canonical_peptides = proteome.create_unique_peptide_pool(
        anno=anno, rule=rule, exception=exception, miscleavage=miscleavage,
        min_mw=min_mw, min_length = min_length, max_length = max_length
    )
    if not quiet:
        logger('canonical peptide pool generated.')
    with open(output_peptides, 'wb') as handle:
        pickle.dump(canonical_peptides, handle)
    if not quiet:
        logger('canonical peptide pool saved to disk.')

    # create list of coding transcripts
    coding_tx = {tx_id for tx_id, tx_model in anno.transcripts.items()
        if tx_model.is_protein_coding}
    with open(output_coding_tx, 'wb') as handle:
        pickle.dump(coding_tx, handle)
