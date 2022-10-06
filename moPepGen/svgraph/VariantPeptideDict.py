""" Module for variant peptide dict """
from __future__ import annotations
from collections import deque
import copy
from typing import Deque, Dict, Iterable, List, Set, Tuple, TYPE_CHECKING
from Bio.Seq import Seq
from Bio import SeqUtils
from moPepGen import aa, get_equivalent
from moPepGen.svgraph.PVGNode import PVGNode
from moPepGen.aa import VariantPeptideIdentifier as vpi


if TYPE_CHECKING:
    from moPepGen.seqvar.VariantRecord import VariantRecord
    from moPepGen import params


class MiscleavedNodeSeries():
    """ Helper class when calling for miscleavage peptides. The nodes contained
    by this class can be joined to make a miscleavage peptide. """
    def __init__(self, nodes:List[PVGNode], additional_variants:Set[VariantRecord]):
        """ """
        self.nodes = nodes
        self.additional_variants = additional_variants

class MiscleavedNodes():
    """ Helper class for looking for peptides with miscleavages. This class
    defines the collection of nodes in sequence that starts from the same node,
    with up to X allowed miscleavages that makes the miscleaved peptide
    sequences.

    Attributes:
        - `data` (Deque[List[PVGNode]]): The collection of nodes that make all
          miscleaved sequences.
        - `cleavage_params` (CleavageParams): Cleavage related parameters.
        - `orf` (Tuple[int,int]): The ORF start and end locations.
        - `tx_id` (str): Transcript ID.
        - `leading_node` (PVGNode): The start node that the miscleaved peptides
          are called from. This node must present in the PVG graph.
    """
    def __init__(self, data:Deque[MiscleavedNodeSeries],
            cleavage_params:params.CleavageParams,
            orf:Tuple[int,int]=None, tx_id:str=None,
            leading_node:PVGNode=None):
        """ constructor """
        self.data = data
        self.orf = orf
        self.tx_id = tx_id
        self.cleavage_params = cleavage_params
        self.leading_node = leading_node

    @staticmethod
    def find_miscleaved_nodes(node:PVGNode, orf:Tuple[int,int],
            cleavage_params:params.CleavageParams, tx_id:str,
            leading_node:PVGNode) -> MiscleavedNodes:
        """ Find all miscleaved nodes.

        node vs leading_node:
            - When calling this function from a node within a PVG graph, a copy
              of the node should be first created. This is because if the ORF
              start site is located in the node, the peptide sequence should be
              called form the ORF start site, rather than the beginning of the
              sequence. So then a copy should be created from this node and
              only keep the sequence after the ORF start.
            - The `leading_node` must tbe the original copy of the node within
              the graph, that the `node` is copied from. This is used to
              retrieve INDELs from the downstream node.

        Args:
            - `node` (PVGNode): The node that miscleaved peptide sequences
              starts from it should be called.
            - `orf` (Tuple[int, int]): The ORF start and end locations.
            - `cleavage_params` (CleavageParams): Cleavage related parameters.
            - `tx_id` (str): Transcript ID.
            - `leading_node` (PVGNode): The start node that the miscleaved
              peptides are called from. This node must present in the PVG graph.
        """
        if orf[0] is None:
            raise ValueError('orf is None')
        queue = deque([[node]])
        nodes = MiscleavedNodes(
            data=deque([]),
            cleavage_params=cleavage_params,
            orf=orf,
            tx_id=tx_id,
            leading_node=leading_node
        )

        if not node.cpop_collapsed:
            additional_variants = leading_node.get_downstream_stop_altering_variants()
            series = MiscleavedNodeSeries([node], additional_variants)
            nodes.data.append(series)

        while queue:
            cur_batch = queue.pop()
            cur_node = cur_batch[-1]
            batch_vars = set()

            for _node in cur_batch:
                for var in _node.variants:
                    batch_vars.add(var.variant)

            n_cleavages = len([x for x in cur_batch if not x.cpop_collapsed]) - 1

            if n_cleavages >= cleavage_params.miscleavage:
                continue

            # This is done to reduce the complexity when the node has too
            # many out nodes, and its out nodes also have too many out nodes.
            if cleavage_params.additional_variants_per_misc == -1:
                allowed_n_vars = float('Inf')
            else:
                allowed_n_vars = cleavage_params.max_variants_per_node
                if n_cleavages > 0:
                    allowed_n_vars += n_cleavages * cleavage_params.additional_variants_per_misc

            for _node in cur_node.out_nodes:
                if not _node.out_nodes:
                    continue
                if _node.truncated:
                    continue
                new_batch = copy.copy(cur_batch)

                is_stop = _node.seq.seq == '*'
                if is_stop:
                    continue

                new_batch.append(_node)

                if not _node.cpop_collapsed:
                    additional_variants = _node.get_downstream_stop_altering_variants()

                    cur_vars = copy.copy(batch_vars)
                    for var in _node.variants:
                        cur_vars.add(var.variant)

                    if len(cur_vars) + len(additional_variants) > allowed_n_vars:
                        continue

                    series = MiscleavedNodeSeries(copy.copy(new_batch), additional_variants)
                    nodes.data.append(series)

                    if n_cleavages + 1 == cleavage_params.miscleavage:
                        continue
                    if n_cleavages + 1 > cleavage_params.miscleavage:
                        raise ValueError('Something just went wrong')
                queue.appendleft(new_batch)
        return nodes

    def is_valid_seq(self, seq:Seq, blacklist:Set[str]) -> bool:
        """ Check whether the seq is valid """
        min_length = self.cleavage_params.min_length
        max_length = self.cleavage_params.max_length
        min_mw = self.cleavage_params.min_mw

        return seq not in blacklist \
            and min_length <= len(seq) <= max_length \
            and 'X' not in seq \
            and SeqUtils.molecular_weight(seq, 'protein') >= min_mw

    def join_miscleaved_peptides(self, check_variants:bool,
            additional_variants:List[VariantRecord], blacklist:Set[str],
            is_start_codon:bool=False
            ) -> Iterable[Tuple[aa.AminoAcidSeqRecord, VariantPeptideMetadata]]:
        """ join miscleaved peptides and update the peptide pool.

        Args:
            - `check_variants` (bool): When true, only peptides that carries at
              least 1 variant are kept. And when false, all unique peptides are
              reported (e.g. noncoding).
            - `additional_variants` (List[VariantRecord]): Additional variants,
              e.g., start gain, cleavage gain, stop lost.
            - `blacklist` (Set[str]): Peptide sequences that should be excluded.
            - `is_start_codon` (bool): Whether the node contains start codon.
        """
        for series in self.data:
            queue = series.nodes
            metadata = VariantPeptideMetadata(orf=self.orf)
            seq:str = None

            variants:Set[VariantRecord] = set()

            for i, node in enumerate(queue):
                other = str(node.seq.seq)
                if seq is None:
                    seq = other
                else:
                    seq = seq + other
                if check_variants:
                    for variant in node.variants:
                        variants.add(variant.variant)
                    if i < len(queue) - 1:
                        _node = self.leading_node if i == 0 else node
                        indels = queue[i + 1].upstream_indel_map.get(_node)
                        if indels:
                            for variant in indels:
                                variants.add(variant)

            if check_variants:
                variants.update(additional_variants)
                variants.update(series.additional_variants)

            cleavage_gain_down = queue[-1].get_cleavage_gain_from_downstream()
            variants.update(cleavage_gain_down)

            if any(v.is_circ_rna() for v in variants)\
                    and any(n.is_missing_any_variant(variants) for n in queue):
                continue

            if not seq:
                continue

            if check_variants and not variants:
                continue

            seq = Seq(seq)

            is_valid = self.is_valid_seq(seq, blacklist)
            is_valid_start = is_start_codon and seq.startswith('M') and\
                self.is_valid_seq(seq[1:], blacklist)

            if not (is_valid or is_valid_start):
                continue

            metadata.is_pure_circ_rna = len(variants) == 1 and \
                list(variants)[0].is_circ_rna()

            label = vpi.create_variant_peptide_id(self.tx_id, variants, None)
            metadata.label = label

            if is_valid:
                aa_seq = aa.AminoAcidSeqRecord(seq=seq)
                yield aa_seq, metadata

            if is_valid_start:
                aa_seq = aa.AminoAcidSeqRecord(seq=seq[1:])
                yield aa_seq, metadata

def update_peptide_pool(seq:aa.AminoAcidSeqRecord,
        peptide_pool:Set[aa.AminoAcidSeqDict],
        label_counter:Set[aa.AminoAcidSeqRecord], label:str,
        update_label:bool=True) -> None:
    """ Add a peptide sequence to a peptide pool if not already exists,
    otherwise update the same peptide's description """
    if 'X' in seq.seq:
        return
    if update_label:
        if label not in label_counter:
            label_counter[label] = 0
        label_counter[label] += 1
        seq.description = label + '|' + str(label_counter[label])
        seq.id = seq.description
        seq.name = seq.description
        same_peptide = get_equivalent(peptide_pool, seq)
        if same_peptide:
            same_peptide:aa.AminoAcidSeqRecord
            same_peptide.description += ' ' + seq.description
            same_peptide.id = same_peptide.description
            same_peptide.name = same_peptide.description
        else:
            peptide_pool.add(seq)
    else:
        if label not in label_counter:
            label_counter[label] = 0
        label_counter[label] += 1
        seq.description = label + '|'\
            + str(label_counter[label])
        peptide_pool.add(seq)

class VariantPeptideMetadata():
    """ Variant peptide metadata """
    def __init__(self, label:str=None, orf:Tuple[int,int]=None,
            is_pure_circ_rna:bool=False):
        """  """
        self.label = label
        self.orf = orf
        self.is_pure_circ_rna = is_pure_circ_rna

T = Dict[aa.AminoAcidSeqRecord, Set[VariantPeptideMetadata]]

class VariantPeptideDict():
    """ Variant peptide pool as dict.

    Attributes:
        - `tx_id` (str): Transcript ID.
        - `peptides` (Dict[AminoAcidSeqRecord, Set[VariantPeptideMetadata]]):
          The peptide data pool, with keys being the AminoAcidRecord, and
          values being set of VariantPeptdieMetadata (variants and ORF).
        - `labels` (Dict[str,int]): Label counter, as a dict with key being the
          variant peptide label, and values being the number of occurrence
          of this label.
        - `global_variant` (VariantRecord): A variant that should be applied to
          every node. This is useful for PVG derived from a circRNA, and this
          argument is the circRNA variant.
    """
    def __init__(self, tx_id:str, peptides:T=None, labels:Dict[str,int]=None,
            global_variant:VariantRecord=None):
        """ constructor """
        self.tx_id = tx_id
        self.peptides = peptides or {}
        self.labels = labels or {}
        self.global_variant = global_variant

    def add_miscleaved_sequences(self, node:PVGNode, orf:List[int, int],
            cleavage_params:params.CleavageParams,
            check_variants:bool, is_start_codon:bool,
            additional_variants:List[VariantRecord], blacklist:Set[str],
            leading_node:PVGNode=None):
        """ Add amino acid sequences starting from the given node, with number
        of miscleavages no more than a given number. The sequences being added
        are the sequence of the current node, and plus n of downstream nodes,
        with n ranges from 1 to miscleavages.

        Args:
            - `node` (PVGNode): The start node.
            - `orf` (List[int, int]): The start and end position of ORF.
            - `cleavage_params` (CleavageParams): Cleavage related parameters.
            - `check_variants` (bool): Whether to check variants.
            - `is_start_codon` (bool): Whether the peptide starts with the start
              codon (M).
            - `additional_variants` (List[VariantRecord]): Additional variants,
              e.g., start gain, cleavage gain, stop lost.
            - `blacklist` (Set[str]): Peptide sequences that should be excluded.
            - `leading_node` (PVGNode): The start node that the miscleaved
              peptides are called from. This node must present in the PVG graph.
        """
        if leading_node is None:
            leading_node = node
        miscleaved_nodes = MiscleavedNodes.find_miscleaved_nodes(
            node=node, orf=orf, cleavage_params=cleavage_params,
            tx_id=self.tx_id, leading_node=leading_node
        )
        if self.global_variant and self.global_variant not in additional_variants:
            additional_variants.append(self.global_variant)
        seqs = miscleaved_nodes.join_miscleaved_peptides(
            check_variants=check_variants,
            additional_variants=additional_variants,
            blacklist=blacklist,
            is_start_codon=is_start_codon
        )
        for seq, metadata in seqs:
            if 'X' in seq.seq:
                continue
            if '*' in seq:
                raise ValueError('Invalid amino acid symbol found in the sequence.')
            val = self.peptides.setdefault(seq, set())
            val.add(metadata)

    def get_peptide_sequences(self, keep_all_occurrence:bool=True,
            orf_id_map:Dict[Tuple[int,int],str]=None
            ) -> Set[aa.AminoAcidSeqRecord]:
        """ Get the peptide sequence and check for miscleavages.

        Args:
            - `keep_all_occurrence` (bool): Whether to keep all occurance of
              the peptide within the graph/transcript.
            - `orf_id_map` (Dict[int, str]): An ID map of ORF IDs from ORF
              start position.
        """
        peptide_pool:Set[aa.AminoAcidSeqRecord] = set()

        while self.peptides:
            seq, metadatas = self.peptides.popitem()
            if '*' in seq:
                raise ValueError('Invalid amino acid symbol found in the sequence.')
            labels = []
            metadatas = list(metadatas)
            metadatas.sort(key=lambda x: x.orf[0])
            had_pure_circ_rna = False
            for metadata in metadatas:
                if metadata.is_pure_circ_rna:
                    if had_pure_circ_rna:
                        continue
                    had_pure_circ_rna = True
                orf_id = None
                if orf_id_map:
                    orf_id = orf_id_map[metadata.orf]

                label = metadata.label
                if orf_id:
                    label += f"|{orf_id}"

                if label in self.labels:
                    self.labels[label] += 1
                else:
                    self.labels[label] = 1
                label += f"|{self.labels[label]}"
                labels.append(label)

                if not keep_all_occurrence:
                    break
            seq.description = " ".join(labels)
            seq.id = seq.description
            seq.name = seq.description

            same_peptide = get_equivalent(peptide_pool, seq)
            if same_peptide:
                same_peptide:aa.AminoAcidSeqRecord
                same_peptide.description += ' ' + seq.description
                same_peptide.id = same_peptide.description
                same_peptide.name = same_peptide.description
            else:
                peptide_pool.add(seq)

        return peptide_pool
