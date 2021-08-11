""" Module for command line interface """
from .generate_index import add_subparser_generate_index, generate_index
from .call_variant_peptide import add_subparser_call_peptides, call_variant_peptide
from .parse_vep import add_subparser_parse_vep, parse_vep
from .parse_reditools import add_subparser_parse_reditools, parse_reditools
from .parse_star_fusion import add_subparser_parse_star_fusion, parse_star_fusion
from .parse_fusion_catcher import add_subparser_parse_fusion_catcher, parse_fusion_catcher
from .parse_rmats import add_subparser_parse_rmats, parse_rmats
