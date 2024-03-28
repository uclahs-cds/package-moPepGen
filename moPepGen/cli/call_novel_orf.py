""" `callNovelORF` calls noncanonical peptide sequences from novel ORFs.
It finds all start codons of any novel ORF gene. """
from __future__ import annotations
import argparse
from typing import TYPE_CHECKING, Set, List, Tuple, IO, Dict
from pathlib import Path
from Bio.SeqIO import FastaIO
from moPepGen import params, svgraph, aa, get_logger, VARIANT_PEPTIDE_SOURCE_DELIMITER
from moPepGen.dna.DNASeqRecord import DNASeqRecordWithCoordinates
from moPepGen.err import ReferenceSeqnameNotFoundError
from moPepGen.cli import common


if TYPE_CHECKING:
    from Bio.Seq import Seq
    from moPepGen.gtf import TranscriptAnnotationModel
    from moPepGen.dna import DNASeqDict

OUTPUT_FILE_FORMATS = ['.fa', '.fasta']

# pylint: disable=W0212
def add_subparser_call_novel_orf(subparsers:argparse._SubParsersAction):
    """ CLI for moPepGen callNovelORF """
    p:argparse.ArgumentParser = subparsers.add_parser(
        name='callNovelORF',
        help='Call non-canonical peptides from novel ORFs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        '-o', '--output-path',
        type=Path,
        help='Output path to the novel ORF peptide FASTA.'
        f" Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=True
    )
    p.add_argument(
        '--output-orf',
        type=Path,
        help='Output path to the FASTA file with novel ORF sequences.'
        f" Valid formats: {OUTPUT_FILE_FORMATS}",
        metavar='<file>',
        required=False,
        default=None
    )
    p.add_argument(
        '--min-tx-length',
        type=int,
        help='Minimal transcript length.',
        metavar='<number>',
        default=21
    )
    p.add_argument(
        '--orf-assignment',
        type=str,
        help='Defines how ORF assignment should be done. The last ORF upstream'
        ' to the peptide is used for `max` and the first (most upstream) one is'
        ' used for `min`',
        choices=['max', 'min'],
        default='max',
        metavar='<choice>'
    )
    p.add_argument(
        '--coding-novel-orf',
        action='store_true',
        help='Include coding transcripts to find alternative ORFs.'
    )
    p.add_argument(
        '--w2f-reassignment',
        action='store_true',
        help='Include peptides with W > F (Tryptophan to Phenylalanine) '
        'reassignment.'
    )
    p.add_argument(
        '--inclusion-biotypes',
        type=Path,
        help='Inclusion biotype list.',
        metavar='<file>',
        default=None
    )
    p.add_argument(
        '--exclusion-biotypes',
        type=Path,
        help='Exclusion biotype list.',
        metavar='<file>',
        default=None
    )

    common.add_args_reference(p)
    common.add_args_cleavage(p)
    common.add_args_debug_level(p)

    p.set_defaults(func=call_novel_orf_peptide)
    common.print_help_if_missing_args(p)
    return p

def call_novel_orf_peptide(args:argparse.Namespace) -> None:
    """ Main entrypoint for calling novel ORF peptide """
    logger = get_logger()
    common.validate_file_format(
        args.output_path, OUTPUT_FILE_FORMATS, check_writable=True
    )
    if args.output_orf:
        common.validate_file_format(
            args.output_orf, OUTPUT_FILE_FORMATS, check_writable=True
        )

    cleavage_params = params.CleavageParams(
        enzyme=args.cleavage_rule,
        exception=args.cleavage_exception,
        miscleavage=int(args.miscleavage),
        min_mw=float(args.min_mw),
        min_length=args.min_length,
        max_length=args.max_length
    )

    common.print_start_message(args)

    genome, anno, proteome, canonical_peptides = common.load_references(
        args=args, load_proteome=True, cleavage_params=cleavage_params
    )

    inclusion_biotypes, exclusion_biotypes = common.load_inclusion_exclusion_biotypes(args)

    novel_orf_peptide_pool = aa.VariantPeptidePool()
    orf_pool = []

    i = 0
    for tx_id in anno.transcripts:
        tx_model = anno.transcripts[tx_id]
        if tx_model.is_protein_coding:
            if not args.coding_novel_orf:
                pass
        else:
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
            peptides, orfs = call_noncoding_peptide_main(
                tx_id=tx_id, tx_model=tx_model, genome=genome,
                canonical_peptides=canonical_peptides,
                cleavage_params=cleavage_params,
                orf_assignment=args.orf_assignment,
                w2f_reassignment=args.w2f_reassignment
            )

            if not orfs:
                continue

            orf_pool.extend(orfs)

            for peptide in peptides:
                novel_orf_peptide_pool.add_peptide(peptide, canonical_peptides,
                    cleavage_params)
        except ReferenceSeqnameNotFoundError as e:
            if not ReferenceSeqnameNotFoundError.raised:
                logger.warning('%s: Make sure your GTF and FASTA files match.', e.args[0])
                ReferenceSeqnameNotFoundError.mute()
        except:
            logger.error('Exception raised from %s', tx_id)
            raise

        i += 1
        if i % 5000 == 0:
            logger.info('%i transcripts processed.', i)

    novel_orf_peptide_pool.write(args.output_path)
    if args.output_orf:
        with open(args.output_orf, 'w') as handle:
            write_orf(orf_pool, handle)

    logger.info('Noncanonical peptide FASTA file written to disk.')


def call_noncoding_peptide_main(tx_id:str, tx_model:TranscriptAnnotationModel,
        genome:DNASeqDict, canonical_peptides:Set[str],
        cleavage_params:params.CleavageParams, orf_assignment:str,
        w2f_reassignment:bool
        ) -> Tuple[Set[aa.AminoAcidSeqRecord],List[aa.AminoAcidSeqRecord]]:
    """ Call novel ORF peptides """
    chrom = tx_model.transcript.location.seqname
    try:
        contig_seq = genome[chrom]
    except KeyError as e:
        raise ReferenceSeqnameNotFoundError(chrom) from e

    tx_seq = tx_model.get_transcript_sequence(contig_seq)

    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=True,
        has_known_orf=False,
        cleavage_params=cleavage_params,
        gene_id=tx_model.gene_id,
        coordinate_feature_type='transcript',
        coordinate_feature_id=tx_id
    )
    dgraph.init_three_frames()
    pgraph = dgraph.translate()
    pgraph.create_cleavage_graph()
    peptide_anno = pgraph.call_variant_peptides(
        check_variants=False,
        check_orf=True,
        denylist=canonical_peptides,
        orf_assignment=orf_assignment,
        w2f=w2f_reassignment,
        check_external_variants=False
    )
    peptide_map:Dict[Seq, Set[str]] = {}
    for seq, annotated_labels in peptide_anno.items():
        for label in annotated_labels:
            if seq in peptide_map:
                peptide_map[seq].add(label.label)
            else:
                peptide_map[seq] = {label.label}
    peptides = set()
    for seq, labels in peptide_map.items():
        label = VARIANT_PEPTIDE_SOURCE_DELIMITER.join(labels)
        peptides.add(
            aa.AminoAcidSeqRecord(
                seq=seq,
                description=label,
                name=label
            )
        )
    orfs = get_orf_sequences(
        pgraph=pgraph,
        tx_id=tx_id,
        gene_id=tx_model.gene_id,
        tx_seq=tx_seq,
        exclude_canonical_orf=False
    )
    return peptides, orfs

def get_orf_sequences(pgraph:svgraph.PeptideVariantGraph, tx_id:str, gene_id:str,
        tx_seq:DNASeqRecordWithCoordinates, exclude_canonical_orf:bool
        ) -> List[aa.AminoAcidSeqRecord]:
    """ Get the full ORF sequences """
    seqs = []
    translate_seqs = [tx_seq[i:].translate() for i in range(3)]
    for orf, orf_id in pgraph.orf_id_map.items():
        orf_start = orf[0]
        if exclude_canonical_orf and tx_seq.orf and tx_seq.orf.start == orf_start:
            continue
        seq_start = int(orf_start / 3)
        reading_frame_index = orf_start % 3
        #seq_start = int((orf_start - node.orf[0]) / 3)
        translate_seq = translate_seqs[reading_frame_index]
        seq_len = translate_seq.seq[seq_start:].find('*')
        if seq_len == -1:
            seq_len = len(translate_seq.seq) - seq_start
        orf_end = orf_start + seq_len * 3
        seq_end = seq_start + seq_len
        seqname = f"{tx_id}|{gene_id}|{orf_id}|{orf_start}-{orf_end}"
        seq = translate_seq[seq_start:seq_end]
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
