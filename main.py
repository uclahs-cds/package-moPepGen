""" Module for mkdocs-macro """
import argparse
from contextlib import redirect_stdout
import io
from moPepGen import cli


COMMAND_MAPPER = {
    'generateIndex': cli.add_subparser_generate_index,
    'parseVEP': cli.add_subparser_parse_vep,
    'parseREDItools': cli.add_subparser_parse_reditools,
    'parseSTARFusion': cli.add_subparser_parse_star_fusion,
    'parseFusionCatcher': cli.add_subparser_parse_fusion_catcher,
    'parseArriba': cli.add_subparser_parse_arriba,
    'parseRMATS': cli.add_subparser_parse_rmats,
    'parseCIRCexplorer': cli.add_subparser_parse_circexplorer,
    'callVariant': cli.add_subparser_call_variant,
    'callNoncoding': cli.add_subparser_call_noncoding,
    'filterFasta': cli.add_subparser_filter_fasta,
    'splitFasta': cli.add_subparser_split_fasta,
    'mergeFasta': cli.add_subparser_merge_fasta,
    'encodeFasta': cli.add_subparser_encode_fasta,
    'decoyFasta': cli.add_subparser_decoy_fasta
}

def define_env(env):
    "Hook function"
    # pylint: disable=W0612
    # pylint: disable=W0212

    @env.macro
    def get_arg_data(command:str):
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        add_parser = COMMAND_MAPPER[command]
        p = add_parser(subparsers)
        return p._actions

    @env.macro
    def get_arg_usage(command:str):
        stream = io.StringIO()
        with redirect_stdout(stream):
            parser = argparse.ArgumentParser()
            subparsers = parser.add_subparsers()
            add_parser = COMMAND_MAPPER[command]
            p = add_parser(subparsers)
            p.print_help()
        return stream.getvalue()
