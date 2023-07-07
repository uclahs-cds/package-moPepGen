""" pytest configuration """
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm.session import Session
from moPepGen.annotation.models import Base


@pytest.fixture
def session() -> Session:
    """ Fixture to create a in-memory sqlite database """
    engine = create_engine('sqlite://')
    with Session(engine) as session:
        Base.metadata.create_all(engine)
        return session
