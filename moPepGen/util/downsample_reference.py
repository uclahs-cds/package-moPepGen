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
import copy
from typing import List, Tuple, Iterable, Dict
from pathlib import Path
import math
from Bio import SeqIO
from moPepGen import gtf, dna, aa
from moPepGen.gtf import GtfIO
from moPepGen.SeqFeature import FeatureLocation, SeqFeature
from moPepGen.gtf.GTFSeqFeature import GTFSeqFeature
from moPepGen.cli.common import add_args_cleavage, add_args_reference, \
    print_help_if_missing_args, add_args_debug_level


# pylint: disable=W0212
def parse_args(subparsers:argparse._SubParsersAction):
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
        metavar='<value>'
    )
    parser.add_argument(
        '--tx-list',
        type=str,
        help='Transcript list to filter.',
        nargs='*',
        metavar='<value>'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Output directory',
        metavar='<file>'
    )
    parser.add_argument(
        '--translate-noncoding',
        type=str,
        choices=['true', 'false'],
        help='Create translate sequences for noncoding at any possible ORF.',
        default='true'
    )
    add_args_reference(parser, index=False)
    add_args_cleavage(parser)
    add_args_debug_level(parser)
    parser.set_defaults(func=main)
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
            gene.transcripts.append(transcript_id)
        transcripts[transcript_id].add_record(feature, record)

def downsample_gtf(path:Path, gene_list=None, tx_list=None) -> gtf.GeneAnnotationModel:
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
        transcript.sort_records()
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
            seq = gene_model.get_gene_sequence(record)
            seq.__class__ = dna.DNASeqRecord
            gene_seqs[gene_id] = seq
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
    for gene_id in gene_seqs.keys():
        gene = anno.genes[gene_id]
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
            for utr in tx_model.utr:
                shift_seq_feature(utr, shift_offset, seqname)
            for sec in tx_model.selenocysteine:
                shift_seq_feature(sec, shift_offset, seqname)
    return genome, anno

def get_noncoding_translate(tx_id:str, anno:gtf.GenomicAnnotation,
        genome:dna.DNASeqDict) -> Dict[str, aa.AminoAcidSeqRecord]:
    """ Translate all possible ORF of a noncoding transcript """
    tx_model = anno.transcripts[tx_id]
    protein_id = tx_id.replace('ENST', 'ENSP')
    gene_id = tx_model.transcript.gene_id
    chrom = tx_model.transcript.chrom
    tx_seq = tx_model.get_transcript_sequence(genome[chrom])
    start_positions = tx_seq.find_all_start_codons()

    translates = {}
    for start in start_positions:
        end = start + math.floor((len(tx_seq) - start) / 3) * 3
        aa_seq = tx_seq[start:end].translate(to_stop=True)
        end = start + len(aa_seq) * 3
        orf = f"ORF{start}:{end}"
        alt_protein_id = f"{protein_id}-{orf}"
        alt_tx_id = f"{tx_id}-{orf}"
        description = f"{alt_protein_id}|{alt_tx_id}|{gene_id}|-"
        aa_seq.id = alt_protein_id
        aa_seq.protein_id = alt_protein_id
        aa_seq.transcript_id = alt_tx_id
        aa_seq.gene_id = gene_id
        aa_seq.description = description
        translates[alt_tx_id] = aa_seq
    return translates

def create_dummy_tx_model(tx_model:gtf.TranscriptAnnotationModel, tx_id:str,
        protein_id:str):
    """ Create dummy transcript model """
    location = FeatureLocation(
        start=tx_model.transcript.location.start,
        end=tx_model.transcript.location.end,
        strand=tx_model.transcript.location.strand
    )
    attrs = copy.deepcopy(tx_model.transcript.attributes)
    attrs['protein_id'] = protein_id
    attrs['transcript_id'] = tx_id
    transcript = GTFSeqFeature(
        location=location, chrom=tx_model.transcript.chrom,
        attributes=attrs, type='transcript'
    )
    model = gtf.TranscriptAnnotationModel(transcript)
    return model

def subset_genome(genome:dna.DNASeqDict, anno:gtf.GenomicAnnotation,
        start:int, end:int) -> Tuple[dna.DNASeqDict, gtf.GenomicAnnotation]:
    """ Subset genome """
    if anno.genes > 1:
        raise ValueError(
            "Don't know how to subset the genome when there are multiple genes."
        )
    gene_id = list(anno.genes.keys())[0]
    gene_model = anno.genes[gene_id]
    if start > gene_model.location.end or end <= gene_model.location.start:
        raise ValueError("Subset location is invalid.")

    for tx_id in list(anno.transcripts.keys()):
        tx_model = anno.transcripts[tx_id]
        if tx_model.transcript.location.end < start:
            anno.transcripts.pop(tx_id)
            continue
        original_start = tx_model.transcript.location.start
        original_end = tx_model.transcript.location.end
        tx_start = max(original_start - start, 0)
        tx_end = min(original_end - start, end - start)
        tx_model.transcript.location = FeatureLocation(
            start=tx_start, end=tx_end, seqname=tx_model.transcript.seqname
        )

        tx_model.exon = subset_and_filter(tx_model.exon, start, end)
        tx_model.cds = subset_and_filter(tx_model.cds, start, end, True)
        tx_model.utr = subset_and_filter(tx_model.utr, start, end)

    seqname = list(genome.keys())[0]
    if len(genome[seqname]) < start or len(genome[seqname]) < end:
        raise ValueError("Don't know how to subset genome like this.")

    genome[seqname] = genome[seqname][start:end]

def subset_and_filter(features:List[SeqFeature], start:int, end:int,
        cds:bool=False) -> List[SeqFeature]:
    """ Subset and filter SeqFeatures """
    new_features = []
    for i, feature in enumerate(features):
        original_start = feature.location.start
        original_end = feature.location.end
        if original_start > end or start > original_end:
            continue
        new_start = max(original_start - start, 0)
        new_end = min(original_end - start, end - start)
        if cds:
            if i == 0 and feature.strand == 1:
                if new_start == 0:
                    new_start += (3 - (start - original_start) % 3) % 3
            if i == len(features) and feature.strand == -1:
                if new_end == end - start:
                    new_end -= (3 - (original_end - end) % 3) % 3
        location = FeatureLocation(
            start=new_start, end=new_end, seqname=feature.seqname,
            strand=feature.strand
        )
        new_feature = SeqFeature(
            chrom=feature.chrom,
            attributes=feature.attributes,
            location=location
        )
        new_features.append(new_feature)
    return new_features


def main(args:argparse.Namespace):
    """ Downsample reference FASTA and GTF """
    genome_fasta:Path = args.genome_fasta
    annotation_gtf:Path = args.annotation_gtf
    protein_fasta:Path = args.proteome_fasta
    gene_list:List[str] = args.gene_list
    tx_list:List[str] = args.tx_list
    output_dir:Path = args.output_dir
    translate_noncoding:bool = args.translate_noncoding == 'true'

    output_dir.mkdir(exist_ok=True)

    anno = downsample_gtf(annotation_gtf, gene_list, tx_list)
    gene_seqs = get_gene_sequences(genome_fasta, anno)
    genome, anno = shift_reference(gene_seqs, anno)
    proteins = downsample_proteins(protein_fasta, anno)

    if translate_noncoding:
        for tx_id in copy.copy(list(anno.transcripts.keys())):
            if tx_id not in proteins:
                aa_seqs = get_noncoding_translate(tx_id, anno, genome)
                proteins.update(aa_seqs)

    with open(output_dir/'annotation.gtf', 'wt') as handle:
        GtfIO.write(handle, anno)

    with open(output_dir/'genome.fasta', 'wt') as handle:
        writer = SeqIO.FastaIO.FastaWriter(handle, record2title=lambda x: x.id)
        for record in genome.values():
            writer.write_record(record)

    with open(output_dir/'proteome.fasta', 'wt') as handle:
        writer = SeqIO.FastaIO.FastaWriter(handle, record2title=lambda x: x.description)
        for record in proteins.values():
            writer.write_record(record)
