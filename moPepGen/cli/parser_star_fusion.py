""" Module for STAR-Fusion parser """
from typing import List, Dict
import pathlib
import argparse
import pickle
from moPepGen import logger, gtf, seqvar, parser, dna


def parse_star_fusion(args:argparse.Namespace) -> None:
    """ Parse the STAR-Fusion's output and save it in TVF format. """
    # unpack args
    fusion = args.fusion
    index_dir:str = args.index_dir
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.tvf'
    verbose = args.verbose

    if verbose:
        logger('moPepGen parseREDItools started.')

    if index_dir:
        with open(f'{index_dir}/genome.pickle', 'rb') as handle:
            genome:dna.DNASeqDict = pickle.load(handle)

        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            anno:gtf.GenomicAnnotation = pickle.load(handle)

        if verbose:
            logger('Indexed genome and annotation loaded.')

    else:
        genome_fasta:str = args.genome_fasta
        annotation_gtf:str = args.annotation_gtf

        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

        genome = dna.DNASeqDict()
        genome.dump_fasta(genome_fasta)
        if verbose:
            logger('Genome assembly FASTA loaded.')

    anno2:Dict[str, Dict[str, gtf.TranscriptAnnotationModel]] = {}
    val:gtf.TranscriptAnnotationModel
    for key, val in anno.transcripts.items():
        gene_id = val.transcript.attributes['gene_id']
        if gene_id not in anno2:
            anno2[gene_id] = {}
        anno2[gene_id][key] = val

    variants:List[seqvar.VariantRecord] = []

    for record in parser.STARFusionParser.parse(fusion):
        var_records = record.convert_to_variant_records(anno2, genome)
        variants.extend(var_records)

    if verbose:
        logger(f'STAR-Fusion output {fusion} loaded.')

    variants.sort()

    if verbose:
        logger('Variants sorted.')

    if index_dir:
        reference_index = pathlib.Path(index_dir).absolute()
        genome_fasta = None
        annotation_gtf = None
    else:
        reference_index = None
        genome_fasta = pathlib.Path(genome_fasta).absolute()
        annotation_gtf = pathlib.Path(annotation_gtf).absolute()

    metadata = seqvar.TVFMetadata(
        parser='parseREDITools',
        reference_index=reference_index,
        genome_fasta=genome_fasta,
        annotation_gtf=annotation_gtf
    )

    seqvar.io.write(variants, output_path, metadata)
