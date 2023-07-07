""" """
from __future__ import annotations
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from moPepGen.annotation.models import Base


class AnnotationDB():
    """ """
    def __init__(self, db_string:str):
        """ Constructor """
        self.engine = create_engine(db_string)
        self.session = Session(self.engine)
        Base.metadata.create_all(self.engine)

    def __del__(self):
        self.session.close()
