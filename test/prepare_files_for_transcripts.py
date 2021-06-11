""" Prepare files for profiling """
import argparse
from pathlib import Path
import subprocess as sp
from Bio import SeqIO


VEP_SNP= '/hot/projects/cpcgene/noncanonical_peptides/Mutation/gencodev34_gr'+\
    'ch38/VEP/germline/filtered_snv/CPCG0100.gencode.aa.tsv'
VEP_INDEL = '/hot/projects/cpcgene/noncanonical_peptides/Mutation/gencodev34'+\
    '_grch38/VEP/germline/filtered_indel/CPCG0100.gencode.aa.tsv'
TRANSCRIPT_FASTA = '/hot/ref/reference/GRCh38-EBI-GENCODE34/gencode.v34.tran'+\
    'scripts.fa'
TRANSLATE_FASTA = '/hot/ref/reference/GRCh38-EBI-GENCODE34/gencode.v34.pc_tr'+\
    'anslations.fa'
GTF = '/hot/ref/reference/GRCh38-EBI-GENCODE34/gencode.v34.chr_patch_hapl_sc'+\
    'aff.annotation.gtf'

def subset_vep(transcript_id, dirname):
    """ Subset a VEP output file by transcript ID. """
    cmd = f'''
    cat {VEP_SNP} | grep {transcript_id} > {dirname}/vep_snp.txt
    cat {VEP_INDEL} | grep {transcript_id} > {dirname}/vep_indel.txt
    '''
    sp.run(cmd, shell=True, check=False)

def subset_gtf(transcript_id, dirname):
    """ Subset a GTF by transcript ID """
    cmd = f'''
    cat {GTF} | grep {transcript_id} > {dirname}/annotation.gtf
    '''
    sp.run(cmd, shell=True, check=False)

def subset_transcript(transcript_id, dirname):
    """ Create a FASTA file with only one transcript. """
    record = None
    for record in SeqIO.parse(TRANSCRIPT_FASTA, 'fasta'):
        if record.id.split('|')[0] == transcript_id:
            break
    SeqIO.write(record, f'{dirname}/transcript.fasta', 'fasta')

def subset_translate(transcript_id, dirname):
    """ Create a FASTA file with only one protein """
    record = None
    for record in SeqIO.parse(TRANSLATE_FASTA, 'fasta'):
        if record.id.split('|')[1] == transcript_id:
            break
    SeqIO.write(record, f'{dirname}/translate.fasta', 'fasta')

def main():
    """ Main entry """
    args = parse_args()
    work_dir = Path(__file__).parent.absolute()
    dirname = f"{work_dir}/files/{args.transcript_id}"
    Path(dirname).mkdir(parents=True, exist_ok=True)
    subset_vep(args.transcript_id, dirname)
    subset_gtf(args.transcript_id, dirname)
    subset_transcript(args.transcript_id, dirname)
    subset_translate(args.transcript_id, dirname)


def parse_args():
    """ parse args """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-t', '--transcript-id',
        type=str,
        help='transcript ID'
    )
    return parser.parse_args()

if __name__ == '__main__':
    main()
