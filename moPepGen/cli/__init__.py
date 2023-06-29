""" Module for command line interface """
from .generate_index import add_subparser_generate_index, generate_index, index_gtf
from .call_variant_peptide import add_subparser_call_variant, call_variant_peptide
from .call_noncoding_peptide import add_subparser_call_noncoding, call_noncoding_peptide
from .call_alt_translation import add_subparser_call_alt_translation, call_alt_translation
from .parse_vep import add_subparser_parse_vep, parse_vep
from .parse_reditools import add_subparser_parse_reditools, parse_reditools
from .parse_star_fusion import add_subparser_parse_star_fusion, parse_star_fusion
from .parse_fusion_catcher import add_subparser_parse_fusion_catcher, parse_fusion_catcher
from .parse_arriba import add_subparser_parse_arriba, parse_arriba
from .parse_rmats import add_subparser_parse_rmats, parse_rmats
from .parse_circexplorer import add_subparser_parse_circexplorer, parse_circexplorer
from .split_fasta import add_subparser_split_fasta, split_fasta
from .filter_fasta import add_subparser_filter_fasta, filter_fasta
from .index_gvf import add_subparser_index_gvf, index_gvf
from .merge_fasta import add_subparser_merge_fasta, merge_fasta
from .encode_fasta import add_subparser_encode_fasta, encode_fasta
from .decoy_fasta import add_subparser_decoy_fasta, decoy_fasta
from .summarize_fasta import add_subparser_summarize_fasta, summarize_fasta
