r""" Downsample reference FASTA and GTF

command line usage:
    python -m test.downsample_reference --help

    python -m test.downsample_reference \
        --genome-fasta path/to/genome.fasta \
        --annotation-gtf path/to/annotation.gtf \
        --protein-fasta path/to/translate.fasta \
        --tx-list ENST0001 ENST0002 \
        --output-dir path/to/downsampled_index

    python -m test.downsample_reference \
        --genome-fasta path/to/genome.fasta \
        --annotation-gtf path/to/annotation.gtf \
        --protein-fasta path/to/translate.fasta \
        --gene-list ENSG0001 ENSG0002 \
        --output-dir path/to/downsampled_index
"""
import argparse
from typing import List, Tuple, Iterable, Dict
from pathlib import Path
from Bio import SeqIO
from moPepGen import gtf, dna, aa
from moPepGen.gtf import GtfIO
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from moPepGen.cli.common import add_args_cleavage, add_args_reference, \
    print_help_if_missing_args


# pylint: disable=W0212
def add_subparser_downsample_reference(subparsers:argparse._SubParsersAction):
    """ Parse args """
    parser:argparse.ArgumentParser = subparsers.add_parser(
        name='downsampleReference',
        help='Create downsampleed reference files to contain only specified'
        ' transcripts.'
    )
    parser.add_argument(
        '--gene-list',
        type=str,
        help='Gene list to filter.',
        nargs='*',
        metavar=''
    )
    parser.add_argument(
        '--tx-list',
        type=str,
        help='Transcript list to filter.',
        nargs='*',
        metavar=''
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Output directory',
        metavar=''
    )
    add_args_reference(parser, index=False)
    add_args_cleavage(parser)
    parser.set_defaults(func=downsample_reference)
    print_help_if_missing_args(parser)
    return parser

GeneTranscriptModel = Tuple[gtf.GeneAnnotationModel, Dict[str, \
    gtf.TranscriptAnnotationModel]]

def parse_gtf(path:Path) -> Iterable[GeneTranscriptModel]:
    """ Parse GTF """
    it = GtfIO.parse(path)
    gene:gtf.GeneAnnotationModel = None
    transcripts:Dict[str, gtf.TranscriptAnnotationModel] = {}

    while True:
        record:GTFSeqFeature = next(it, None)
        if not record:
            if gene:
                yield gene, transcripts
            return
        if record.type.lower() == 'gene':
            if gene:
                yield gene, transcripts
            gene = record
            gene.__class__ = gtf.GeneAnnotationModel
            gene.transcripts = []
            transcripts = {}
            continue
        feature = record.type.lower()
        if feature not in gtf.GTF_FEATURE_TYPES:
            continue
        transcript_id = record.transcript_id
        if transcript_id not in transcripts:
            transcripts[transcript_id] = gtf.TranscriptAnnotationModel()
        transcripts[transcript_id].add_record(feature, record)
        if transcript_id not in transcripts:
            gene.transcripts.append(transcript_id)


def downsample_gtf(path:Path, gene_list=None, tx_list=None
        ) -> gtf.GeneAnnotationModel:
    """ Downsample a GTF file """
    anno = gtf.GenomicAnnotation()

    for model in parse_gtf(path):
        if gene_list:
            gene_id = model[0].gene_id
            if gene_id in gene_list:
                anno.genes[gene_id] = model[0]
                anno.transcripts.update(model[1])
        if tx_list:
            gene_wanted = False
            tx_ids = []
            for key,val in model[1].items():
                if key in tx_list:
                    gene_wanted = True
                    anno.transcripts[key] = val
                    tx_ids.append(key)
            if gene_wanted:
                gene_id = model[0].gene_id
                anno.genes[gene_id] = model[0]
                anno.genes[gene_id].transcripts = tx_ids
    for transcript in anno.transcripts.values():
        transcript.exon.sort()
        transcript.cds.sort()
    return anno


def get_gene_sequences(path:Path, anno:gtf.GenomicAnnotation) -> dna.DNASeqDict:
    """ Get gene sequences """
    gene_seqs = dna.DNASeqDict()
    chrom_genes:Dict[str, List[gtf.GeneAnnotationModel]] = {}
    for gene_model in anno.genes.values():
        chrom = gene_model.location.seqname
        if chrom not in chrom_genes:
            chrom_genes[chrom] = []
        chrom_genes[chrom].append(gene_model)

    for record in SeqIO.parse(path, 'fasta'):
        record.__class__ = dna.DNASeqRecord
        record.id = record.id.split(' ')[0]
        record.id = record.id.split('|')[0]
        if record.id not in chrom_genes:
            continue
        for gene_model in chrom_genes[record.id]:
            gene_id = gene_model.gene_id
            gene_seqs[gene_id] = gene_model.get_gene_sequence(record)
    return gene_seqs

def downsample_proteins(path:Path, anno:gtf.GenomicAnnotation
        ) -> aa.AminoAcidSeqDict:
    """ Downsample proteins """
    proteins = aa.AminoAcidSeqDict()
    tx_ids = list(anno.transcripts.keys())

    count = 0
    source = None
    infered = set()

    for record in SeqIO.parse(path, 'fasta'):
        record.__class__ = aa.AminoAcidSeqRecord
        record:aa.AminoAcidSeqRecord
        if count > 100 and not source:
            source = infered.pop()
        if not source:
            count += 1
            infered.add(record.infer_ids())
        else:
            record.infer_ids(style=source)

        if record.transcript_id in tx_ids:
            proteins[record.transcript_id] = record
    return proteins

def shift_seq_feature(feature:SeqFeature, offset:int, seqname:str='chr1'):
    """ Shift SeqFeature with a given offset """
    feature.chrom = seqname
    start = feature.location.start + offset
    end = feature.location.end + offset
    strand = feature.location.strand
    location = FeatureLocation(seqname=seqname, start=start, end=end, strand=strand)
    feature.location = location

def shift_reference(gene_seqs:dna.DNASeqDict, anno:gtf.GenomicAnnotation
        ) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotation]:
    """ Shift reference """
    seqname = 'chr1'
    genome = dna.DNASeqDict({seqname: None})
    for gene_id, gene_seq in gene_seqs.items():
        gene_seq:dna.DNASeqRecord
        strand = anno.genes[gene_id].location.strand
        seq = gene_seq
        if strand == -1:
            seq = seq.reverse_complement()
        if not genome[seqname]:
            genome[seqname] = seq
        else:
            genome[seqname] += seq
    genome[seqname].id = seqname

    gene_start = 0
    for gene in anno.genes.values():
        shift_offset = gene_start - gene.location.start
        gene_start = gene_start + gene.location.end - gene.location.start
        shift_seq_feature(gene, shift_offset, seqname)

        for tx_id in gene.transcripts:
            tx_model = anno.transcripts[tx_id]
            shift_seq_feature(tx_model.transcript, shift_offset, seqname)
            for exon in tx_model.exon:
                shift_seq_feature(exon, shift_offset, seqname)
            for cds in tx_model.cds:
                shift_seq_feature(cds, shift_offset, seqname)
    return genome, anno

def downsample_reference(args:argparse.Namespace):
    """ Downsample reference FASTA and GTF """
    genome_fasta:Path = args.genome_fasta
    annotation_gtf:Path = args.annotation_gtf
    protein_fasta:Path = args.protein_fasta
    gene_list:List[str] = args.gene_list
    tx_list:List[str] = args.tx_list
    output_dir:Path = args.output_dir

    output_dir.mkdir(exist_ok=True)

    anno = downsample_gtf(annotation_gtf, gene_list, tx_list)
    gene_seqs = get_gene_sequences(genome_fasta, anno)
    genome, anno = shift_reference(gene_seqs, anno)
    proteins = downsample_proteins(protein_fasta, anno)

    GtfIO.write(output_dir/'annotation.gtf', anno)

    with open(output_dir/'genome.fasta', 'wt') as handle:
        writer = SeqIO.FastaIO.FastaWriter(handle, record2title=lambda x: x.id)
        for record in genome.values():
            writer.write_record(record)

    with open(output_dir/'proteome.fasta', 'wt') as handle:
        SeqIO.write(proteins.values(), handle, 'fasta')