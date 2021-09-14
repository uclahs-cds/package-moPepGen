""" Module for calling noncoding peptide """
from __future__ import annotations
import copy
from typing import TYPE_CHECKING, Set, List, Tuple, IO
from pathlib import Path
import pkg_resources
from Bio.SeqIO import FastaIO
from moPepGen import svgraph, aa, logger
from moPepGen.err import ReferenceSeqnameNotFoundError, warning
from moPepGen.cli.common import add_args_cleavage, add_args_verbose, add_args_reference, \
    print_start_message, print_help_if_missing_args, load_references


if TYPE_CHECKING:
    import argparse
    from moPepGen.gtf import TranscriptAnnotationModel
    from moPepGen.dna import DNASeqDict

# pylint: disable=W0212
def add_subparser_call_noncoding(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callNoncoding """
    p = subparsers.add_parser(
        name='callNoncoding',
        help='Call non-canonical peptides from noncoding transcripts.'
    )
    p.add_argument(
        '-t', '--min-tx-length',
        type=int,
        help='Minimal transcript length.',
        metavar='',
        default=21
    )
    p.add_argument(
        '-i', '--inclusion-biotypes',
        type=Path,
        help='Inclusion biotype list.',
        metavar='',
        default=None
    )
    p.add_argument(
        '-e', '--exclusion-biotypes',
        type=Path,
        help='Exclusion biotype list.',
        metavar='',
        default=None
    )
    p.add_argument(
        '-o', '--output-prefix',
        type=Path,
        help='File prefix for the output FASTA.',
        metavar='',
        required=True
    )

    add_args_reference(p)
    add_args_cleavage(p)
    add_args_verbose(p)

    p.set_defaults(func=call_noncoding_peptide)
    print_help_if_missing_args(p)

def call_noncoding_peptide(args:argparse.Namespace) -> None:
    """ Main entry poitn for calling noncoding peptide """
    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    exception = 'trypsin_exception' if rule == 'trypsin' else None
    min_length:int = args.min_length
    max_length:int = args.max_length
    peptide_fasta = f"{args.output_prefix}_peptide.fasta"
    orf_fasta = f"{args.output_prefix}_orf.fasta"

    print_start_message(args)

    genome, anno, proteome, canonical_peptides = load_references(
        args=args, load_proteome=True
    )

    inclusion_biotypes = []
    if args.inclusion_biotypes:
        with open(args.inclusion_biotypes, 'rt') as handle:
            for line in handle:
                inclusion_biotypes.append(line.rstrip())

    exclusion_path = args.exclusion_biotypes
    if not exclusion_path:
        exclusion_path = pkg_resources.resource_filename(
            'moPepGen', 'data/gencode_hs_exclusion_list.txt'
        )

    exclusion_biotypes = []
    if exclusion_path:
        with open(exclusion_path, 'rt') as handle:
            for line in handle:
                exclusion_biotypes.append(line.rstrip())

    noncanonical_pool = aa.VariantPeptidePool()
    orf_pool = []

    i = 0
    for tx_id, tx_model in anno.transcripts.items():
        if inclusion_biotypes and \
                tx_model.transcript.biotype not in inclusion_biotypes:
            continue
        if exclusion_biotypes and \
                tx_model.transcript.biotype in exclusion_biotypes:
            continue
        if tx_id in proteome:
            continue
        if tx_model.transcript_len() < args.min_tx_length:
            continue

        try:
            peptides, orfs = call_noncoding_peptide_main(tx_id, tx_model, genome,
                rule, exception, miscleavage)

            if not orfs:
                continue

            orf_pool.extend(orfs)

            for peptide in peptides:
                noncanonical_pool.add_peptide(peptide, canonical_peptides,
                    min_mw, min_length, max_length)
        except ReferenceSeqnameNotFoundError as e:
            if not ReferenceSeqnameNotFoundError.raised:
                warning(e.args[0] + ' Make sure your GTF and FASTA files match.')
                ReferenceSeqnameNotFoundError.mute()
        except:
            logger(f'Exception raised from {tx_id}')
            raise

        if args.verbose:
            i += 1
            if i % 5000 == 0:
                logger(f'{i} transcripts processed.')

    noncanonical_pool.write(peptide_fasta)
    with open(orf_fasta, 'w') as handle:
        write_orf(orf_pool, handle)

    if args.verbose:
        logger('Noncanonical peptide FASTA file written to disk.')


def call_noncoding_peptide_main(tx_id:str, tx_model:TranscriptAnnotationModel,
        genome:DNASeqDict, rule:str, exception:str, miscleavage:int
        ) -> Tuple[Set[aa.AminoAcidSeqRecord],List[aa.AminoAcidSeqRecord]]:
    """ Call noncoding peptides """
    chrom = tx_model.transcript.location.seqname
    try:
        contig_seq = genome[chrom]
    except KeyError as e:
        raise ReferenceSeqnameNotFoundError(chrom) from e

    tx_seq = tx_model.get_transcript_sequence(contig_seq)

    dgraph = svgraph.TranscriptVariantGraph(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=True
    )
    dgraph.add_null_root()
    dgraph.find_all_orfs()
    if not dgraph.root.out_edges:
        return None, None
    pgraph = dgraph.translate()
    orfs = get_orf_sequences(pgraph)
    pgraph.form_cleavage_graph(rule=rule, exception=exception)
    peptides = pgraph.call_variant_peptides(
        miscleavage=miscleavage,
        check_variants=False,
        check_orf=True
    )
    return peptides, orfs

def get_orf_sequences(pgraph:svgraph.PeptideVariantGraph) -> List[aa.AminoAcidSeqRecord]:
    """ Get the full ORF sequences """
    seqs = []
    if not pgraph.orf_id_map:
        pgraph.create_orf_id_map()
    orf_id_map = pgraph.orf_id_map
    for node in pgraph.root.out_nodes:
        orf_start = node.orf[0]
        orf_end = orf_start + len(node.seq.seq) * 3
        orf_id = orf_id_map[orf_start]
        seqname = f"{pgraph.id}|{orf_id}|{orf_start}-{orf_end}"
        seq = copy.copy(node.seq)
        seq.id = seqname
        seq.name = seqname
        seq.description = seqname
        seqs.append(seq)
    return seqs

def write_orf(orfs:List[aa.AminoAcidSeqRecord], handle:IO):
    """ Write ORF sequences """
    record2title = lambda x: x.description
    writer = FastaIO.FastaWriter(handle, record2title=record2title)
    for record in orfs:
        writer.write_record(record)
