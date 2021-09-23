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
    gtf = GenomicAnnotation()
    gtf.dump_gtf(f'{file_dir}/annotation.gtf')

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
            start=gtf.transcripts[transcript_id].get_cds_start_index(),
            end=gtf.transcripts[transcript_id].get_cds_end_index()
        )

    # protein
    proteins = aa.AminoAcidSeqDict()
    proteins.dump_fasta(f'{file_dir}/translate.fasta')

    # vep
    variants = {}

    variant_file = f'{file_dir}/vep_moPepGen.txt'

    with open(variant_file, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            transcript_id = fields[0]
            record = seqvar.VariantRecord(
                location=FeatureLocation(
                    seqname=transcript_id,
                    start=int(fields[1]),
                    end=int(fields[2])
                ),
                ref=fields[3],
                alt=fields[4],
                _type=fields[5],
                _id=fields[6]
            )
            if transcript_id not in variants:
                variants[transcript_id] = [record]
            else:
                variants[transcript_id].append(record)

    dgraph = svgraph.TranscriptVariantGraph(
        seq=transcript_seq,
        _id=transcript_id
    )
    dgraph.create_variant_graph(variants[transcript_id])
    dgraph.fit_into_codons()
    pgraph = dgraph.translate()
    pgraph.form_cleavage_graph('trypsin')
    peptides = pgraph.call_variant_peptides()
    print(peptides)

if __name__ == '__main__':
    run_task()
