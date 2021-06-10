""""""
from typing import Dict, List
import argparse
import pickle
from moPepGen import logger, gtf, seqvar, parser


def parse_reditools(args:argparse.Namespace) -> None:
    """"""
    # unpack args
    table_file = args.table_file
    transcript_id_column = args.transcript_id_column
    index_dir:str = args.index_dir
    output_prefix:str = args.output_prefix
    output_path = output_prefix + '.mop'
    verbose = args.verbose

    if verbose:
        logger('moPepGen parseREDITools started.')

    if index_dir:
        with open(f'{index_dir}/annotation.pickle', 'rb') as handle:
            anno = pickle.load(handle)

        if verbose:
            logger('Index annotation loaded')

    else:
        annotation_gtf:str = args.annotation_gtf

        anno = gtf.TranscriptGTFDict()
        anno.dump_gtf(annotation_gtf)
        if verbose:
            logger('Annotation GTF loaded.')

    variants:Dict[str, List[seqvar.VariantRecord]] = {}

    for record in parser.REDItoolsParser.parse(table_file, transcript_id_column):
        _vars = record.convert_to_variant_records(anno)
        for variant in _vars:
            transcript_id = variant.location.seqname
            if transcript_id not in variants:
                variants[transcript_id] = []
            variants[transcript_id].append(variant)

    if verbose:
        logger(f'REDITools table {table_file} loaded.')

    for records in variants.values():
        records.sort()

    if verbose:
        logger('Variants sorted.')

    with open(output_path, 'w') as handle:
        mode = 'w'
        for records in variants.values():
            seqvar.io.write(records, handle, mode)
            if mode == 'w':
                mode = 'a'

    if verbose:
        logger("Variants written to disk.")


if __name__ == '__main__':
    test_args = argparse.Namespace()
    test_args.table_file = 'test/files/reditools/CPT0208690010_merged_chr22.txt'
    test_args.transcript_id_column = 16
    test_args.index_dir = 'test/files/downsampled_set/gencode_v36_index'
    test_args.output_prefix = 'test/files/reditools/CPT0208690010_merged_chr22.mop'
    test_args.verbose = True
    parse_reditools(test_args)
