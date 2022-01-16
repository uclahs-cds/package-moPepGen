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

class SubgraphMapper():
    """ SubgraphMapper """
    def __init__(self, data:Dict[str, SubgraphLocation]=None):
        """ """
        self.data = data or {}

    def __getitem__(self, key:str) -> SubgraphLocation:
        """ """
        return self.data[key]

    def add_subgraph(self, child_id:str, parent_id:str, level:int,
            location:FeatureLocation):
        """ """
        subgraph = SubgraphLocation(
            _id=child_id, level=level, parent_id=parent_id,
            location=location
        )
        self.data[child_id] = subgraph

    @staticmethod
    def get_subgraph_id() -> str:
        """ Generate a random subgraph ID """
        return str(uuid.uuid4())
