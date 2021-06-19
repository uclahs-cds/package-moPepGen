""" Module for call variant paptide """
from typing import List, Dict
import argparse
import pickle
from Bio import SeqIO, SeqUtils
from moPepGen import svgraph, dna, gtf, aa, seqvar, logger, get_equivalent


def call_variant_peptide(args:argparse.Namespace) -> None:
    """ Main entry point for calling variant peptide """
    variant_files:List[str] = args.input_variant
    output_fasta:str = args.output_fasta
    index_dir:str = args.index_dir
    verbose:bool = args.verbose
    rule:str = args.cleavage_rule
    miscleavage:int = int(args.miscleavage)
    min_mw:float = float(args.min_mw)
    exception = 'trypsin_exception' if rule == 'trypsin' else None

    if verbose:
        logger('moPepGen callPeptide started.')

    if index_dir:
        with open(f'{index_dir}/genome.pickle', 'rb') as handle:
            genome = pickle.load(handle)

        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            annotation = pickle.load(handle)

        with open(f"{index_dir}/canonical_peptides.pickle", 'rb') as handle:
            canonical_peptides = pickle.load(handle)
    else:
        genome_fasta:str = args.genome_fasta
        proteome_fasta:str = args.proteome_fasta
        annotation_gtf:str = args.annotation_gtf

        annotation = gtf.TranscriptGTFDict()
        annotation.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

        genome = dna.DNASeqDict()
        genome.dump_fasta(genome_fasta)
        if verbose:
            logger('Genome assembly FASTA loaded.')

        proteome = aa.AminoAcidSeqDict()
        proteome.dump_fasta(proteome_fasta)
        if verbose:
            logger('Proteome FASTA loaded.')

        canonical_peptides = proteome.create_unique_peptide_pool(
            rule=rule, exception=exception, miscleavage=miscleavage,
            min_mw=min_mw
        )
        if verbose:
            logger('canonical peptide pool generated.')

    variants:Dict[str, List[seqvar.VariantRecord]] = {}
    for file in variant_files:
        for record in seqvar.io.parse(file):
            transcript_id = record.location.seqname
            if transcript_id not in variants:
                variants[transcript_id] = [record]
            else:
                variants[transcript_id].append(record)
        logger(f'Variant file {file} loaded.')

    for records in variants.values():
        records.sort()
    if verbose:
        logger('Variant records sorted.')

    variant_peptides = set()

    i = 0
    for transcript_id, variant_records in variants.items():
        anno:gtf.TranscriptAnnotationModel = annotation[transcript_id]
        chrom = anno.transcript.location.seqname
        transcript_seq = anno.get_transcript_sequence(genome[chrom])

        dgraph = svgraph.TranscriptVariantGraph(
            seq=transcript_seq,
            transcript_id=transcript_id
        )

        ## Create transcript variant graph
        # dgraph.create_variant_graph(variant_records)
        dgraph.add_null_root()
        variant_iter = iter(variant_records)
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
                # TODO: inplement fusion
                donor_transcript_id = variant.attrs['DONOR_TRANSCRIPT_ID']
                donor_anno = annotation[donor_transcript_id]
                donor_chrom = donor_anno.transcript.location.seqname
                donor_seq = anno.get_transcript_sequence(genome[donor_chrom])
                donor_variant_records = variants[donor_transcript_id] \
                    if donor_transcript_id in variants else []
                cur = dgraph.apply_fusion(cur, variant, donor_seq, donor_variant_records)
                variant = next(variant_iter, None)
                continue

            cur = dgraph.apply_variant(cur, variant)
            if len(cur.in_edges) == 0:
                dgraph.root = cur
            variant = next(variant_iter, None)

        dgraph.fit_into_codons()
        pgraph = dgraph.translate()

        pgraph.form_cleavage_graph(rule=rule, exception=exception)
        peptides = pgraph.call_vaiant_peptides(miscleavage=miscleavage)

        for peptide in peptides:
            if SeqUtils.molecular_weight(peptide.seq, 'protein') < min_mw:
                continue
            if str(peptide.seq) in canonical_peptides:
                continue
            same_peptide = get_equivalent(variant_peptides, peptide)
            if same_peptide:
                new_label = peptide.id
                same_peptide.id += ('||' + new_label)
                same_peptide.name = same_peptide.id
                same_peptide.description = same_peptide.id
            else:
                variant_peptides.add(peptide)

        if verbose:    # for logging
            i += 1
            if i % 1000 == 0:
                logger(f'{i} transcripts processed.')

    with open(output_fasta, 'w') as handle:
        SeqIO.write(variant_peptides, handle, 'fasta')

    if verbose:
        logger('Variant peptide FASTA file written to disk.')


if __name__ == '__main__':
    test_args = argparse.Namespace()
    test_args.input_variant = [
        'test/files/vep_test_files/vep.tvf',
        'test/files/fusion/fusion.tvf'
    ]
    test_args.index_dir = 'test/files/vep_test_files/index'
    test_args.output_fasta = 'test/files/vep_test_files/vep_moPepGen.fasta'
    test_args.verbose = True
    test_args.cleavage_rule = 'trypsin'
    test_args.miscleavage = 2
    test_args.min_mw = 500.
    call_variant_peptide(args=test_args)
