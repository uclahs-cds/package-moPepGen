""" VEP2VariantPeptides module """
from typing import Dict, List
import argparse
from moPepGen.vep import VepIO
from moPepGen import gtf, dna, vep, seqvar, logger


_VALID_BIOTYPES = ['protein_coding']

def parse_vep(args:argparse.Namespace) -> None:
    """ Main entry point for the VEP parser. """
    # unpack args
    vep_files:List[str] = args.vep_txt
    genome_fasta:str = args.genome_fasta
    annotation_gtf:str = args.annotation_gtf
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '_moPepGen.txt'
    verbose = args.verbose

    if verbose:
        logger('moPepGen parseVEP started.')

    anno = gtf.TranscriptGTFDict()
    anno.dump_gtf(annotation_gtf)
    if verbose:
        logger('Annotation GTF loaded.')
    
    genome = dna.DNASeqDict()
    genome.dump_fasta(genome_fasta)
    if verbose:
        logger('Genome assembly FASTA loaded.')
    
    vep_records:Dict[str, List[seqvar.VariantRecord]] = {}
    
    for vep_file in vep_files:
        for record in VepIO.parse(vep_file):
            transcript_id = record.feature

            biotype = anno[transcript_id].transcript.biotype

            if biotype not in _VALID_BIOTYPES:
                continue

            if transcript_id not in vep_records.keys():
                vep_records[transcript_id] = []
            
            chrom_seqname = record.location.split(':')[0]
            transcript_seq = anno[transcript_id]\
                .get_transcript_sequence(genome[chrom_seqname])

            record = record.convert_to_variant_record(transcript_seq)
            
            vep_records[transcript_id].append(record)
            
        if verbose:
            logger(f'VEP file {vep_file} loaded.')

    for records in vep_records.values():
        records.sort()
    
    if verbose:
        logger('VEP sorting done.')

    with open(output_path, 'w') as handle:
        headers = ['transcript_id', 'start', 'end', 'ref', 'alt', 'type', 'id']
        handle.write('#' + '\t'.join(headers) + '\n')
        for transcript_id, records in vep_records.items():
            record:seqvar.VariantRecord
            for record in records:
                line = [transcript_id, str(int(record.location.start)),
                    str(int(record.location.end)), str(record.ref),
                    str(record.alt), record.type, record.id]
                handle.write('\t'.join(line) + '\n')
    
    if verbose:
        logger('Variant info written to disk.')
        