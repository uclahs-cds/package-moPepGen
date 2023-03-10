""" Module for variant peptide dict """
from __future__ import annotations
from collections import deque
import copy
from typing import Deque, Dict, Iterable, List, Set, Tuple, TYPE_CHECKING
from Bio.Seq import Seq
from Bio import SeqUtils
from moPepGen import aa, circ, get_equivalent, seqvar
from moPepGen.SeqFeature import FeatureLocation
from moPepGen.svgraph.PVGNode import PVGNode
from moPepGen.svgraph.PVGOrf import PVGOrf
from moPepGen.aa import VariantPeptideIdentifier as vpi
from moPepGen.svgraph.SubgraphTree import SubgraphTree


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
            orfs:List[PVGOrf]=None, tx_id:str=None, gene_id:str=None,
            leading_node:PVGNode=None, subgraphs:SubgraphTree=None):
        """ constructor """
        self.data = data
        self.orfs = orfs or []
        self.tx_id = tx_id
        self.gene_id = gene_id
        self.cleavage_params = cleavage_params
        self.leading_node = leading_node
        self.subgraphs = subgraphs

    @staticmethod
    def find_miscleaved_nodes(node:PVGNode, orfs:List[PVGOrf],
            cleavage_params:params.CleavageParams, tx_id:str, gene_id:str,
            leading_node:PVGNode, subgraphs:SubgraphTree, is_circ_rna:bool
            ) -> MiscleavedNodes:
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
            - `orfs` (List[PVGPOrf]): The ORF start and end locations.
            - `cleavage_params` (CleavageParams): Cleavage related parameters.
            - `tx_id` (str): Transcript ID.
            - `leading_node` (PVGNode): The start node that the miscleaved
              peptides are called from. This node must present in the PVG graph.
        """
        if not orfs:
            raise ValueError('ORFs are empty')
        queue = deque([[node]])
        nodes = MiscleavedNodes(
            data=deque([]),
            cleavage_params=cleavage_params,
            orfs=orfs,
            tx_id=tx_id,
            gene_id=gene_id,
            leading_node=leading_node,
            subgraphs=subgraphs
        )

        if not (node.cpop_collapsed or node.truncated):
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
                allowed_n_vars += (n_cleavages + 1) * cleavage_params.additional_variants_per_misc

            for _node in cur_node.out_nodes:
                if is_circ_rna and _node.is_hybrid_node(subgraphs):
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
                    cur_vars_plus_additional = copy.copy(cur_vars)
                    cur_vars_plus_additional.update(additional_variants)
                    if len(cur_vars_plus_additional) > allowed_n_vars:
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
            is_start_codon:bool=False, circ_rna:circ.CircRNAModel=None,
            truncate_sec:bool=False
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
            metadata = VariantPeptideMetadata()
            seq:str = None
            variants:Set[VariantRecord] = set()
            in_seq_variants:Set[VariantRecord] = set()

            selenocysteins = []

            for i, node in enumerate(queue):
                other = str(node.seq.seq)
                if seq is None:
                    seq = other
                    if truncate_sec:
                        selenocysteins = copy.copy(node.selenocysteines)
                else:
                    if truncate_sec:
                        selenocysteins += [x.shift(len(seq)) for x in node.selenocysteines]
                    seq = seq + other

                if check_variants:
                    for variant in node.variants:
                        if variant.is_silent:
                            continue
                        if variant.upstream_cleavage_altering:
                            continue
                        variants.add(variant.variant)
                        if not variant.variant.is_circ_rna():
                            in_seq_variants.add(variant.variant)
                    if i < len(queue) - 1:
                        _node = self.leading_node if i == 0 else node
                        indels = queue[i + 1].upstream_indel_map.get(_node)
                        if indels:
                            for variant in indels:
                                variants.add(variant)

            valid_orf = None
            for orf in self.orfs:
                if any(v.is_circ_rna() for v in variants) \
                        or any(v.is_circ_rna() for v in orf.start_gain):
                    if orf.is_valid_orf_to_misc_nodes(queue, self.subgraphs, circ_rna):
                        if any(n.is_missing_any_variant(in_seq_variants) for n in queue):
                            continue
                        metadata.orf = tuple(orf.orf)
                        valid_orf = orf
                        break
                else:
                    metadata.orf = tuple(orf.orf)
                    valid_orf = orf
                    break

            if valid_orf is None:
                continue

            cleavage_gain_down = queue[-1].get_cleavage_gain_from_downstream()
            variants.update(cleavage_gain_down)

            if check_variants:
                variants.update(valid_orf.start_gain)
                variants.update(additional_variants)
                variants.update(series.additional_variants)

            if any(v.is_circ_rna() for v in variants) \
                    and queue[-1].any_unaccounted_downstream_cleavage_or_stop_altering(
                        {x for x in variants if not x.is_circ_rna()}):
                continue

            if not seq:
                continue

            if check_variants and not (variants or selenocysteins):
                continue

            metadata.is_pure_circ_rna = len(variants) == 1 and \
                list(variants)[0].is_circ_rna()

            seqs = self.translational_modification(Seq(seq), metadata, blacklist,
                variants, is_start_codon, selenocysteins, check_variants)
            for seq, metadata in seqs:
                yield seq, metadata

    def translational_modification(self, seq:Seq, metadata:VariantPeptideMetadata,
            blacklist:Set[str], variants:List[VariantRecord], is_start_codon:bool,
            selenocysteins:List[seqvar.VariantRecordWithCoordinate],
            check_variants:bool
            ) -> Iterable[Tuple[str,VariantPeptideMetadata]]:
        """ Apply any modification that could happen during translation. The
        kinds of modifications that could happen are:
        1. Leading Methionine truncation.
        2. Selenocysteine truncation.
        """
        if variants or not check_variants:
            is_valid = self.is_valid_seq(seq, blacklist)
            is_valid_start = is_start_codon and seq.startswith('M') and\
                self.is_valid_seq(seq[1:], blacklist)

            if is_valid or is_valid_start:
                cur_metadata = copy.copy(metadata)
                label = vpi.create_variant_peptide_id(
                    transcript_id=self.tx_id, variants=variants, orf_id=None,
                    gene_id=self.gene_id
                )
                cur_metadata.label = label

                if is_valid:
                    aa_seq = aa.AminoAcidSeqRecord(seq=seq)
                    yield aa_seq, cur_metadata

                if is_valid_start:
                    aa_seq = aa.AminoAcidSeqRecord(seq=seq[1:])
                    yield aa_seq, cur_metadata

        for sec in selenocysteins:
            seq_mod = seq[:sec.location.start]
            is_valid = self.is_valid_seq(seq_mod, blacklist)
            is_valid_start = is_start_codon and seq_mod.startswith('M') and\
                self.is_valid_seq(seq_mod[1:], blacklist)

            if is_valid or is_valid_start:
                cur_metadata = copy.copy(metadata)
                cur_variants = [v for v in variants if v.location.end
                    <= sec.variant.location.start]
                cur_variants.append(sec.variant)
                label = vpi.create_variant_peptide_id(
                    transcript_id=self.tx_id, variants=cur_variants, orf_id=None,
                    gene_id=self.gene_id
                )
                cur_metadata.label = label

                if is_valid:
                    aa_seq = aa.AminoAcidSeqRecord(seq=seq_mod)
                    yield aa_seq, cur_metadata

                if is_valid_start:
                    aa_seq = aa.AminoAcidSeqRecord(seq=seq_mod[1:])
                    yield aa_seq, cur_metadata

    @staticmethod
    def any_overlaps_and_all_missing_variant(nodes:Iterable[PVGNode], variant:VariantRecord):
        """ Whether any node overlaps with a given variant and all nodes are
        missing it. Should be only used for circRNA. """
        any_overlap = False
        all_missing = True
        for node in nodes:
            locs = []
            if node.seq.locations:
                for loc in node.seq.locations:
                    tx_loc = FeatureLocation(
                        start=int(loc.ref.start) * 3 + node.reading_frame_index,
                        end=int(loc.ref.end) * 3 + node.reading_frame_index
                    )
                    locs.append(tx_loc)
            if node.variants:
                locs += [x.variant.location for x in node.variants
                    if not x.variant.is_circ_rna()]
            node_variants = {x.variant for x in node.variants}
            if not any(loc.overlaps(variant.location) for loc in locs):
                continue
            any_overlap = True
            for loc in locs:
                if loc.overlaps(variant.location) \
                        and any(variant == v for v in node_variants):
                    all_missing = False
        return any_overlap and all_missing


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
            global_variant:VariantRecord=None, gene_id:str=None,
            truncate_sec:bool=False):
        """ constructor """
        self.tx_id = tx_id
        self.peptides = peptides or {}
        self.labels = labels or {}
        self.global_variant = global_variant
        self.gene_id = gene_id
        self.truncate_sec = truncate_sec

    def add_miscleaved_sequences(self, node:PVGNode, orfs:List[PVGOrf],
            cleavage_params:params.CleavageParams,
            check_variants:bool, is_start_codon:bool,
            additional_variants:List[VariantRecord], blacklist:Set[str],
            leading_node:PVGNode=None, subgraphs:SubgraphTree=None,
            circ_rna:circ.CircRNAModel=None):
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
            node=node, orfs=orfs, cleavage_params=cleavage_params,
            tx_id=self.tx_id, gene_id=self.gene_id, leading_node=leading_node,
            subgraphs=subgraphs, is_circ_rna=circ_rna is not None
        )
        if self.global_variant and self.global_variant not in additional_variants:
            additional_variants.append(self.global_variant)
        seqs = miscleaved_nodes.join_miscleaved_peptides(
            check_variants=check_variants,
            additional_variants=additional_variants,
            blacklist=blacklist,
            is_start_codon=is_start_codon,
            circ_rna=circ_rna,
            truncate_sec=self.truncate_sec
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
            unique_labels = set()
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

                if label in unique_labels:
                    continue
                unique_labels.add(label)

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
