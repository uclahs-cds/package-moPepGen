"""
Worker functions for calling alternative translation peptides.

This module contains the core logic for identifying peptides with alternative
translation events (W>F reassignment, selenocysteine termination) from coding
transcripts, decoupled from CLI concerns.
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Dict, Set

from moPepGen import svgraph, aa, VARIANT_PEPTIDE_SOURCE_DELIMITER

if TYPE_CHECKING:
    from Bio.Seq import Seq
    from moPepGen.gtf import TranscriptAnnotationModel, GenomicAnnotation
    from moPepGen.dna import DNASeqDict
    from moPepGen.params import CodonTableInfo, CleavageParams


def call_alt_translation_for_transcript(
    tx_id: str,
    tx_model: TranscriptAnnotationModel,
    genome: DNASeqDict,
    anno: GenomicAnnotation,
    codon_table: CodonTableInfo,
    cleavage_params: CleavageParams,
    w2f_reassignment: bool,
    sec_truncation: bool
) -> Set[aa.AminoAcidSeqRecord]:
    """
    Call alternative translation peptides for a single transcript.

    This function builds a three-frame transcript variant graph, translates it,
    creates a cleavage graph, and calls peptides with alternative translation
    events (W>F reassignment or Sec truncation).

    Args:
        tx_id: Transcript ID
        tx_model: Transcript annotation model
        genome: Genome sequence dictionary
        anno: Genomic annotation
        codon_table: Codon table information for the chromosome
        cleavage_params: Parameters for peptide cleavage
        w2f_reassignment: Include W>F (Trp to Phe) reassignment peptides
        sec_truncation: Include selenocysteine termination peptides

    Returns:
        Set of alternative translation peptide records
    """
    chrom = tx_model.transcript.chrom
    tx_seq = tx_model.get_transcript_sequence(genome[chrom])

    # Build three-frame TVG
    dgraph = svgraph.ThreeFrameTVG(
        seq=tx_seq,
        _id=tx_id,
        cds_start_nf=tx_model.is_cds_start_nf(),
        has_known_orf=True,
        mrna_end_nf=tx_model.is_mrna_end_nf(),
        cleavage_params=cleavage_params,
        coordinate_feature_type='transcript',
        coordinate_feature_id=tx_id
    )

    # Gather Sec variants from annotation
    dgraph.gather_sect_variants(anno)
    dgraph.init_three_frames()

    # Translate to PVG
    pgraph = dgraph.translate(
        table=codon_table.codon_table,
        start_codons=codon_table.start_codons
    )

    # Create cleavage graph and call peptides
    pgraph.create_cleavage_graph()
    peptide_anno = pgraph.call_variant_peptides(
        check_variants=True,
        truncate_sec=sec_truncation,
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

    return peptides
