"""
Worker functions for calling novel ORF peptides.

This module contains the core logic for identifying peptides from novel open
reading frames (ORFs) in non-coding or coding transcripts, decoupled from CLI
concerns.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Set, List, Tuple, Dict, IO

from Bio.Seq import Seq
from Bio.SeqIO import FastaIO

from moPepGen import svgraph, aa, VARIANT_PEPTIDE_SOURCE_DELIMITER
from moPepGen.err import ReferenceSeqnameNotFoundError

if TYPE_CHECKING:
    from moPepGen.gtf import TranscriptAnnotationModel
    from moPepGen.dna import DNASeqDict, DNASeqRecordWithCoordinates
    from moPepGen.params import CodonTableInfo, CleavageParams


def call_novel_orf_for_transcript(
    tx_id: str,
    tx_model: TranscriptAnnotationModel,
    genome: DNASeqDict,
    canonical_peptides: Set[str],
    codon_table: CodonTableInfo,
    cleavage_params: CleavageParams,
    orf_assignment: str,
    w2f_reassignment: bool
) -> Tuple[Set[aa.AminoAcidSeqRecord], List[aa.AminoAcidSeqRecord]]:
    """
    Call novel ORF peptides for a single transcript.

    This function builds a three-frame transcript variant graph without a known
    ORF, translates it, creates a cleavage graph, and calls peptides from all
    discovered ORFs that are not in the canonical peptide set.

    Args:
        tx_id: Transcript ID
        tx_model: Transcript annotation model
        genome: Genome sequence dictionary
        canonical_peptides: Set of canonical peptides to exclude
        codon_table: Codon table information for the chromosome
        cleavage_params: Parameters for peptide cleavage
        orf_assignment: How to assign ORFs ('max' or 'min')
        w2f_reassignment: Include W>F (Trp to Phe) reassignment peptides

    Returns:
        Tuple of (peptide records, ORF sequences)

    Raises:
        ReferenceSeqnameNotFoundError: If chromosome not found in genome
    """
    chrom = tx_model.transcript.location.seqname
    try:
        contig_seq = genome[chrom]
    except KeyError as e:
        raise ReferenceSeqnameNotFoundError(chrom) from e

    tx_seq = tx_model.get_transcript_sequence(contig_seq)

    # Build three-frame TVG without known ORF
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

    # Translate to PVG
    pgraph = dgraph.translate(
        table=codon_table.codon_table,
        start_codons=codon_table.start_codons
    )

    # Create cleavage graph and call peptides with ORF checking
    pgraph.create_cleavage_graph()
    peptide_anno = pgraph.call_variant_peptides(
        check_variants=False,
        check_orf=True,
        denylist=canonical_peptides,
        orf_assignment=orf_assignment,
        w2f=w2f_reassignment,
        check_external_variants=False
    )

    # Aggregate labels for peptides
    peptide_map: Dict[Seq, Set[str]] = {}
    for seq, annotated_labels in peptide_anno.items():
        for label in annotated_labels:
            if seq in peptide_map:
                peptide_map[seq].add(label.label)
            else:
                peptide_map[seq] = {label.label}

    # Create AminoAcidSeqRecord objects
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

    # Extract ORF sequences
    orfs = get_orf_sequences(
        pgraph=pgraph,
        tx_id=tx_id,
        gene_id=tx_model.gene_id,
        tx_seq=tx_seq,
        exclude_canonical_orf=False,
        codon_table=codon_table.codon_table,
        start_codons=codon_table.start_codons
    )

    return peptides, orfs


def get_orf_sequences(
    pgraph: svgraph.PeptideVariantGraph,
    tx_id: str,
    gene_id: str,
    tx_seq: DNASeqRecordWithCoordinates,
    exclude_canonical_orf: bool,
    codon_table: str,
    start_codons: List[str]
) -> List[aa.AminoAcidSeqRecord]:
    """
    Extract full ORF sequences from the peptide variant graph.

    Args:
        pgraph: Peptide variant graph with discovered ORFs
        tx_id: Transcript ID
        gene_id: Gene ID
        tx_seq: Transcript sequence with coordinates
        exclude_canonical_orf: Whether to exclude the canonical ORF
        codon_table: Codon table name for translation
        start_codons: List of valid start codons

    Returns:
        List of ORF amino acid sequences
    """
    seqs = []
    translate_seqs = [tx_seq[i:].translate(table=codon_table) for i in range(3)]

    for orf, orf_id in pgraph.orf_id_map.items():
        orf_start = orf[0]

        # Skip canonical ORF if requested
        if exclude_canonical_orf and tx_seq.orf and tx_seq.orf.start == orf_start:
            continue

        # Validate start codon
        assert tx_seq.seq[orf_start:orf_start+3] in start_codons

        # Calculate sequence coordinates
        seq_start = int(orf_start / 3)
        reading_frame_index = orf_start % 3
        translate_seq = translate_seqs[reading_frame_index]

        # Find stop codon
        seq_len = translate_seq.seq[seq_start:].find('*')
        if seq_len == -1:
            seq_len = len(translate_seq.seq) - seq_start

        orf_end = orf_start + seq_len * 3
        seq_end = seq_start + seq_len

        # Create sequence record
        seqname = f"{tx_id}|{gene_id}|{orf_id}|{orf_start}-{orf_end}"
        seq = translate_seq[seq_start:seq_end]

        # Ensure first amino acid is methionine
        if seq.seq[0] != 'M':
            seq.seq = Seq('M') + seq.seq[1:]

        seq.id = seqname
        seq.name = seqname
        seq.description = seqname
        seqs.append(seq)

    return seqs


def write_orf_fasta(orfs: List[aa.AminoAcidSeqRecord], handle: IO) -> None:
    """
    Write ORF sequences to a FASTA file.

    Args:
        orfs: List of ORF amino acid sequence records
        handle: Open file handle for writing
    """
    record2title = lambda x: x.description
    writer = FastaIO.FastaWriter(handle, record2title=record2title)
    for record in orfs:
        writer.write_record(record)
