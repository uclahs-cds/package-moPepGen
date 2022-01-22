""" SubgraphMapper """
from __future__ import annotations
from typing import Dict, Tuple
import uuid
from moPepGen.SeqFeature import FeatureLocation


class SubgraphLocation():
    """ ParentGraphLocation """
    def __init__(self, _id:str, level:int, parent_id:str,  location:FeatureLocation):
        """ constructur """
        self.id = _id
        self.level = level
        self.parent_id = parent_id
        self.location = location

class SubgraphTree():
    """ This SubgraphTree class defines the relationship of subgraphes and
    parents/children in a tree-like data structure. """
    def __init__(self, data:Dict[str, SubgraphLocation]=None,
            root:SubgraphLocation=None):
        """ Constructor """
        self.data = data or {}
        self.root = root

    def __getitem__(self, key:str) -> SubgraphLocation:
        """ Get item """
        return self.data[key]

    def add_root(self, _id:str):
        """ Add the main graph """
        subgraph = SubgraphLocation(
            _id=_id, level=0, parent_id=None, location=None
        )
        self.data[_id] = subgraph
        self.root = subgraph

    def add_subgraph(self, child_id:str, parent_id:str, level:int,
            start:int, end:int):
        """ Add subgraph """
        location = FeatureLocation(
            seqname=parent_id, start=start, end=end
        )
        subgraph = SubgraphLocation(
            _id=child_id, level=level, parent_id=parent_id,
            location=location
        )
        self.data[child_id] = subgraph

    @staticmethod
    def generate_subgraph_id() -> str:
        """ Generate a random subgraph ID """
        return str(uuid.uuid4())

    def find_compatible_parents(self, first:str, second:str
            ) -> Tuple[SubgraphLocation, SubgraphLocation]:
        subgraph1 = self[first.subgraph_id]
        subgraph2 = self[second.subgraph_id]

        while subgraph1.id != subgraph2.id and \
                subgraph1.parent_id != subgraph2.parent_id and \
                subgraph1.id !=  subgraph2.parent_id and \
                subgraph1.parent_id != subgraph2.id:
            level1 = subgraph1.level
            level2 = subgraph2.level
            if level1 >= level2:
                subgraph1 = self.subgraphs[subgraph1.parent_id]
            if level2 >= level1:
                subgraph2 = self.subgraphs[subgraph2.parent_id]

        return subgraph1, subgraph2

