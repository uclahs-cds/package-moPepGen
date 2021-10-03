""" Profile the moPepGen callPeptide process """
import pathlib
from Bio import SeqIO
from moPepGen import dna, aa, svgraph, seqvar
from moPepGen.gtf import GenomicAnnotation
from moPepGen.SeqFeature import FeatureLocation, MatchedLocation


def run_task():
    """ run task """
    transcript_id = 'ENST00000422127.5'
    file_dir = f'{pathlib.Path(__file__).parent.absolute()}/files/{transcript_id}'

    # gtf
    anno = GenomicAnnotation()
    anno.dump_gtf(f'{file_dir}/annotation.gtf')

    # transcript
    transcript_seq = None
    for transcript_seq in SeqIO.parse(f'{file_dir}/transcript.fasta', 'fasta'):
        transcript_seq.__class__ = dna.DNASeqRecordWithCoordinates
        location = MatchedLocation(
            query=FeatureLocation(start=0, end=len(transcript_seq)),
            ref=FeatureLocation(start=0, end=len(transcript_seq))
        )
        transcript_seq.locations = [location]
        transcript_seq.orf = FeatureLocation(
            start=anno.transcripts[transcript_id].get_cds_start_index(),
            end=anno.transcripts[transcript_id].get_cds_end_index()
        )

    # protein
    proteins = aa.AminoAcidSeqDict()
    proteins.dump_fasta(f'{file_dir}/translate.fasta')

    # vep
    variant_pool = seqvar.VariantRecordPool()

    variant_file = f'{file_dir}/vep_moPepGen.txt'

    with open(variant_file, 'rt') as handle:
        seqvar.GVFMetadata.parse(handle)
        variant_pool.load_variants(handle, anno, None)

    tx_variants = variant_pool.transcriptional[transcript_id]

    dgraph = svgraph.ThreeFrameTVG(
        seq=transcript_seq,
        _id=transcript_id,
        has_known_orf=True
    )
    dgraph.create_variant_graph(tx_variants, variant_pool, None, anno)
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()
    pgraph.create_cleavage_graph('trypsin')
    peptides = pgraph.call_variant_peptides()
    print(peptides)

if __name__ == '__main__':
    run_task()
