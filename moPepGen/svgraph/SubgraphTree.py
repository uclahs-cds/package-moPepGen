""" SubgraphMapper """
from __future__ import annotations
from typing import Dict
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
    def __init__(self, data:Dict[str, SubgraphLocation]=None):
        """ Constructor """
        self.data = data or {}

    def __getitem__(self, key:str) -> SubgraphLocation:
        """ Get item """
        return self.data[key]

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
