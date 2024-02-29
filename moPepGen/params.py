""" This module defined classes in order to group certain parameters together. """
from __future__ import annotations
from typing import TYPE_CHECKING, Set

from moPepGen import aa


if TYPE_CHECKING:
    from moPepGen import dna, gtf


class CleavageParams():
    """ Cleavage related parameters.

    ## Attributes:
        - enzyme (str): Enzyme name.
        - exception (str): Enzymatic cleavage exception.
        - miscleavage (int): Number of cleavages to allow per non-canonical peptide.
        - min_mw (int): The minimal molecular weight of the non-canonical peptides.
        - min_length (int): The minimal length of non-canonical peptides, inclusive.
        - max_length (int): The maximum length of non-canonical peptides, inclusive.
        - max_variants_per_node (int): Maximal number of variants per node. This
            can be useful when there are local regions that are heavily mutated.
            When creating the cleavage graph, nodes containing variants larger
            than this value are skipped. Setting to -1 will avoid checking for
            this.
        - additional_variants_per_misc (int): Additional variants allowed for
            every miscleavage. This argument is used together with
            max_variants_per_node to handle hypermutated regions. Setting to -1
            will avoid checking for this.
        - min_nodes_to_collapse (int): When making the cleavage graph, the minimal
            number of nodes to trigger pop collapse.
        - naa_to_collapse (int): The number of bases used for pop collapse.
    """
    def __init__(self, enzyme:str=None, exception:str=None, miscleavage:int=2,
            min_mw:int=500, min_length:int=7, max_length:int=25,
            max_variants_per_node:int=7, additional_variants_per_misc:int=2,
            min_nodes_to_collapse:int=30, naa_to_collapse:int=5):
        """ constructor """
        self.enzyme = enzyme
        self.exception = exception
        self.miscleavage = miscleavage
        self.min_mw = min_mw
        self.min_length = min_length
        self.max_length = max_length
        self.max_variants_per_node = max_variants_per_node
        self.additional_variants_per_misc = additional_variants_per_misc
        self.min_nodes_to_collapse = min_nodes_to_collapse
        self.naa_to_collapse = naa_to_collapse
        if self.exception == 'auto':
            if enzyme == 'trypsin':
                self.exception = 'trypsin_exception'
            else:
                self.exception = None

    def jsonfy(self, graph_params:bool=False):
        """ jsonfy """
        data = {
            'enzyme': self.enzyme,
            'exception': self.exception,
            'miscleavage': self.miscleavage,
            'min_mw': self.min_mw,
            'min_length': self.min_length,
            'max_length': self.max_length
        }
        if graph_params:
            data.update({
                'max_variants_per_node': self.max_variants_per_node,
                'additional_variants_per_misc': self.additional_variants_per_misc,
                'min_nodes_to_collapse': self.min_nodes_to_collapse,
                'naa_to_collapse': self.naa_to_collapse
            })
        return data

class ReferenceData():
    """ Reference related parameters

    ## Attributes
        - genome (dna.DNASeqDict)
        - anno (gtf.GeneAnnotationModel)
        - canonical_peptides (Set[str])
        - proteome (aa.AminoAcidSeqDict)
    """
    def __init__(self, genome:dna.DNASeqDict, anno:gtf.GenomicAnnotation,
            canonical_peptides:Set[str], proteome:aa.AminoAcidSeqDict=None):
        """ constructor """
        self.genome = genome
        self.anno = anno
        self.canonical_peptides = canonical_peptides
        self.proteome = proteome  or aa.AminoAcidSeqDict()
