""""""
from __future__ import annotations
import copy
from moPepGen.aa.AminoAcidSeqRecord import AminoAcidSeqRecord
from typing import List, Set, Tuple
from collections import deque
import itertools
import copy
from Bio.Seq import Seq
from moPepGen import get_equivalent
from moPepGen import aa
from moPepGen import vep
from moPepGen.svgraph.VariantRecordWithCoordinate import VariantRecordWithCoordinate


class PeptideNode():
    """"""
    def __init__(self, seq:aa.AminoAcidSeqRecord,
            variants:List[VariantRecordWithCoordinate]=None,
            in_nodes:Set[PeptideNode]=None,
            out_nodes:Set[PeptideNode]=None,
            frameshifts:Set[vep.VEPVariantRecord]=None):
        """"""
        self.seq = seq
        self.variants = variants
        self.in_nodes = set() if in_nodes is None else in_nodes
        self.out_nodes = set() if out_nodes is None else out_nodes
        self.frameshifts = set() if frameshifts is None else frameshifts
    
    def add_out_edge(self, node:PeptideNode):
        self.out_nodes.add(node)
        node.in_nodes.add(self)
    
    def remove_out_edge(self, node:PeptideNode):
        try:
            self.out_nodes.remove(node)
        except KeyError:
            pass
        
        try:
            node.in_nodes.remove(self)
        except KeyError:
            pass

        return
        
    
    def find_reference_next(self) -> PeptideNode:
        if not self.out_nodes:
            return None
        for node in self.out_nodes:
            if not node.variants:
                return node
        raise ValueError('None of the out nodes is reference.')
    
    def find_reference_prev(self) -> PeptideNode:
        if not self.out_nodes:
            return None
        for node in self.in_nodes:
            if not node.variants:
                return node
        raise ValueError('None of the in nodes is reference.')
    
    def split_node(self, index:int) -> PeptideNode:
        """"""
        left_seq = self.seq[:index]
        right_seq = self.seq[index:]
        
        left_variants = []
        right_variants = []
        for variant in self.variants:
            if variant.location.start < index:
                left_variants.append(variant)
            if variant.location.end > index:
                right_variants.append(variant.shift(-index))
        self.seq = left_seq
        self.variants = left_variants
        
        new_node = PeptideNode(seq=right_seq, variants=right_variants)
        new_node.frameshifts = copy.copy(self.frameshifts)

        while self.out_nodes:
            node = self.out_nodes.pop()
            self.remove_out_edge(node)
            new_node.add_out_edge(node)
        self.add_out_edge(new_node)
        return new_node
    
    
class PeptideVariantGraph():
    """"""
    def __init__(self, root:PeptideNode):
        """"""
        self.root = root
        self.stop = PeptideNode(aa.AminoAcidSeqRecord(Seq('*')))
    
    def add_stop(self, node:PeptideNode):
        node.add_out_edge(self.stop)
    
    def next_is_stop(self, node:PeptideNode):
        return self.stop in node.out_nodes
    
    def cleave_if_posible(self, node:PeptideNode, rule:str,
            exception:str=None, return_first:bool=False) -> PeptideNode:
        """"""
        first_node = node
        site = node.seq.find_first_enzymatic_cleave_site(rule=rule,
            exception=exception)
        while site > -1:
            node = node.split_node(site)
            site = node.seq.find_first_enzymatic_cleave_site(rule=rule,
                exception=exception)
        return first_node if return_first else node

    def expand_alignment_backword(self, node:PeptideNode, rule:str,
            exception:str=None) -> Tuple[PeptideNode, Set[PeptideNode]]:
        """"""
        branches = set()
        to_return = self.stop
        if len(node.in_nodes) > 1:
            raise ValueError('Inbond node has multiple')
        
        upstream = next(iter(node.in_nodes))
        
        while node.out_nodes:
            cur = node.out_nodes.pop()
            node.remove_out_edge(cur)
            if cur is self.stop:
                cur = PeptideNode(
                    seq=node.seq,
                    variants=node.variants,
                    frameshifts=copy.copy(node.frameshifts)
                )
                self.add_stop(cur)
            else:
                cur.seq = node.seq + cur.seq

                if not cur.variants:
                    to_return = cur

                variants = []
                for variant in node.variants:
                    variants.append(variant)
                for variant in cur.variants:
                    variants.append(variant.shift(len(cur.seq)))
                cur.variants = variants
                
                # only branches out if there is new frameshifts
                for frameshift in cur.frameshifts:
                    if frameshift not in node.frameshifts:
                        branches.add(cur)
            upstream.add_out_edge(cur)
            last = self.cleave_if_posible(cur, rule=rule, exception=exception)
            if to_return is cur:
                to_return = last
            
        
        upstream.remove_out_edge(node)
        return to_return, branches
    
    def expand_alignment_forward(self, node:PeptideNode, rule:str,
            exception:str=None) -> Tuple[PeptideNode, Set[PeptideNode]]:
        """"""
        if len(node.out_nodes) > 1:
            raise ValueError('Outbond node has multiple')
        
        downstream = next(iter(node.out_nodes))

        while node.in_nodes:
            cur = node.in_nodes.pop()
            cur.remove_out_edge(node)
            cur.seq = cur.seq + node.seq
            cur.remove_out_edge(node)
            cur.add_out_edge(downstream)
            self.cleave_if_posible(cur, rule=rule, exception=exception)
        
        node.remove_out_edge(downstream)
        return downstream
    
    def join_merge_alignments(self, node:PeptideNode, rule:str, 
            exception:str=None) -> Tuple[PeptideNode, Set[PeptideNode]]:
        """"""
        branches = set()
        to_return = self.stop
        prev_ref = node.find_reference_prev()
        if len(prev_ref.in_nodes) > 1:
            raise ValueError('upstream has multiple')
        upstream = next(iter(prev_ref.in_nodes))
        primary_nodes = copy.copy(node.in_nodes)
        secondary_nodes = copy.copy(node.out_nodes)
        for first, second in itertools.product(primary_nodes, secondary_nodes):
            variants = copy.copy(first.variants)
            for variant in second.variants:
                variants.append(variant.shift(len(first.seq) + len(node.seq)))
            
            if second is self.stop:
                seq = first.seq + node.seq
            else:
                seq = first.seq + node.seq + second.seq

            frameshifts = copy.copy(first.frameshifts)
            frameshifts_before = len(frameshifts)
            frameshifts.update(second.frameshifts)
            frameshifts_after = len(frameshifts)
            
            new_node = PeptideNode(
                seq=seq,
                variants=variants,
                frameshifts=frameshifts
            )
            upstream.add_out_edge(new_node)

            if frameshifts_before != frameshifts_after:
                branches.add()
            
            if second is self.stop:
                self.add_stop(new_node)
            else:
                for downstream in second.out_nodes:
                    new_node.add_out_edge(downstream)
            
            new_node = self.cleave_if_posible(new_node, rule=rule,
                exception=exception)

            if not new_node.variants:
                to_return = new_node
            
        
        for first in primary_nodes:
            upstream.remove_out_edge(first)
        
        for second in secondary_nodes:
            downstreams = copy.copy(second.out_nodes)
            for downstream in downstreams:
                second.remove_out_edge(downstream)
        
        if to_return is not self.stop:
            to_return = to_return.find_reference_next()
        return to_return, branches
    
    def cleave_merge_alignments(self, node:PeptideNode, site:int, rule:str,
            exception:str=None) -> Tuple[PeptideNode, Set[PeptideNode]]:
        """"""
        to_return = self.stop
        branches = set()
        left = node
        right = node.split_node(site)
        primary_nodes = copy.copy(left.in_nodes)
        secondary_nodes = copy.copy(right.out_nodes)

        primary_nodes2 = set()
        while primary_nodes:
            first = primary_nodes.pop()
            first.seq = first.seq + left.seq
            while first.out_nodes:
                out_node = first.out_nodes.pop()
                first.remove_out_edge(out_node)
            last = self.cleave_if_posible(first, rule=rule, 
                exception=exception)
            primary_nodes2.add(last)
        primary_nodes = primary_nodes2
        
        for second in secondary_nodes:
            second.seq = right.seq + second.seq
            
            variants = []
            for variant in right.variants:
                variants.append(variant)
            for variant in second.variants:
                variants.append(variant.shift(len(right.seq)))
            second.variants = variants

            while second.in_nodes:
                in_node = second.in_nodes.pop()
                in_node.remove_out_edge(second)
            
            # only branch out when there is new frameshifts
            for frameshift in second.frameshifts:
                if frameshift not in right.frameshifts:
                    branches.add(second)
            
            last = self.cleave_if_posible(first, rule=rule,
                exception=exception)
            if not variants:
                to_return = last
            
        for first, second in itertools.product(primary_nodes, secondary_nodes):
            first.add_out_edge(second)
        
        if to_return is not self.stop:
            to_return = to_return.find_reference_next()
        return to_return, branches
    

    def to_cleavage_graph(self, rule:str, exception:str=None):
        """"""
        queue = deque([self.root])
        branch = set()

        while queue:
            cur = queue.pop()

            # each is a branch
            if cur.seq is None:
                for out_node in cur.out_nodes:
                    queue.appendleft(out_node)
                continue
            
            site = cur.seq.find_first_enzymatic_cleave_site(rule=rule,
                exception=exception)
            while site > -1:
                cur = cur.split_node(site)
                site = cur.seq.find_first_enzymatic_cleave_site(rule=rule,
                    exception=exception)
            
            if self.next_is_stop(cur):
                if len(cur.out_nodes) > 1:
                    # return branch
                    cur, branches = self.expand_alignment_backword(cur,
                        rule=rule, exception=exception)
                    for branch in branches:
                        queue.appendleft(branch)
                continue
            
            # expand the variant alignments to the left
            # returns the new reference node in the alignment.
            # cur is NOT stop
            cur, branches = self.expand_alignment_backword(cur,
                        rule=rule, exception=exception)
            for branch in branches:
                queue.appendleft(branch)

            if self.next_is_stop(cur):
                # leave it to the next iteration to handle
                queue.append(cur)
                continue

            cur = cur.find_reference_next()

            sites = cur.seq.find_all_enzymatic_cleave_sites(rule=rule,
                exception=exception)

            branches = set()
            if len(sites) == 0:
                if self.next_is_stop(cur):
                    self.expand_alignment_forward(cur,
                        rule=rule, exception=exception)
                    continue
                cur, branches = self.join_merge_alignments(cur,
                    rule=rule, exception=exception)
            elif len(sites) == 1:
                cur, branches = self.cleave_merge_alignments(cur, sites[0],
                    rule=rule, exception=exception)
            else:
                cur.split_node(sites[0])
                cur = self.expand_alignment_forward(cur,
                    rule=rule, exception=exception)
            for branch in branches:
                queue.appendleft(branch)
            if cur is not self.stop:
                queue.appendleft(cur)

    def call_vaiant_peptides(self, miscleavage:int=2, min_mw:float=500.):
        """"""
        queue:deque[PeptideNode] = deque([self.root])
        visited:Set[PeptideNode] = set()
        variant_peptides:Set[AminoAcidSeqRecord] = set()
        while queue:
            cur = queue.pop()
            if cur.seq is None:
                for out_node in cur.out_nodes:
                    queue.appendleft(out_node)
                continue
            
            # skip the node if already visited
            visited_len_before = len(visited)
            visited.add(cur)
            visited_len_after = len(visited)
            if visited_len_before == visited_len_after:
                continue
            
            # add out nodes to the queue
            for out_node in cur.out_nodes:
                if out_node is not self.stop:
                    queue.appendleft(out_node)
            
            node_paths:deque[deque[PeptideNode]] = deque([deque([cur])])
            i = 0
            while i < miscleavage and node_paths:
                i += 1
                next_batch:deque[deque[PeptideNode]] = deque()
                while node_paths:
                    nodes = node_paths.pop()

                    variants = []
                    seq = None
                    for node in nodes:
                        if seq is None:
                            seq = node.seq
                        else:
                            seq = seq + node.seq
                        if node.variants:
                            variants.extend(node.variants)
                    if not variants:
                        continue
                    variant_label = ''
                    for variant in variants:
                        variant_label += ('|' + str(variant))
                    same_peptide = get_equivalent(variant_peptides, seq)
                    if same_peptide:
                        same_peptide.id += variant_label
                        same_peptide.name = same_peptide.id
                        same_peptide.description = same_peptide.id
                    else:
                        seq.id = seq.transcript_id + variant_label
                        seq.name = seq.id
                        seq.description = seq.id
                        variant_peptides.add(seq)
                    
                    if i + 1 < miscleavage:
                    # create a new batch of node paths for the next iteration
                        last_node = nodes[-1]
                        for node in last_node.out_nodes:
                            if node is self.stop:
                                continue
                            new_list = copy.copy(nodes)
                            new_list.append(node)
                            next_batch.appendleft(new_list)

                node_paths = next_batch
        
        return variant_peptides