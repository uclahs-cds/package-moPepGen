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
    'parseRMATS': cli.add_subparser_parse_rmats,
    'parseCIRCexplorer': cli.add_subparser_parse_circexplorer,
    'callVariant': cli.add_subparser_call_variant,
    'callNoncoding': cli.add_subparser_call_noncoding,
    'filterFasta': cli.add_subparser_filter_fasta,
    'splitDatabase': cli.add_subparser_split_database
}

def define_env(env):
    "Hook function"
    # pylint: disable=W0612

    @env.macro
    def get_arg_data(command:str):
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        add_parser = COMMAND_MAPPER[command]
        p = add_parser(subparsers)
        return p._actions
