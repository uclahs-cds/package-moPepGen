""" """
from __future__ import annotations
from typing import List, Tuple
import pytest
from sqlalchemy.orm.session import Session
from moPepGen.annotation.models import Plex, Sample, Run, Peptide


@pytest.mark.parametrize(
    "data,n",
    [(
        [('1', 'TMT')],
        1
    ),(
        [('1', 'TMT8plex'), ('2', 'TMT8plex')],
        2
    )]
)
def test__plex__rows_added_successfully(
        session: Session, data:List[Tuple[str,str]], n:int):
    """ Test rows are added successfully to the Plex table """
    for plex_id, experiment_type in data:
        plex = Plex(plex_id=plex_id, experiment_type=experiment_type)
        session.add(plex)
    session.commit()
    assert session.query(Plex).count() == n


@pytest.mark.parametrize(
    'plex_data,sample_data,n',
    [(
        [('1', 'TMT8plex')],
        [('UCLA0001', 1), ('UCLA0002', 1)],
        2
    ), (
        [('1', 'TMT8plex'), ('2', 'TMT8plex')],
        [('UCLA0001', 1), ('UCLA0002', 1), ('UCLA0003', 2), ('UCLA0004', 2)],
        4
    )]
)
def test__sample__rows_added_successfully(session: Session,
        plex_data:List[Tuple[str,str]], sample_data:List[Tuple[str, int]],
        n:int):
    """ Test rows are added successfully to the Sample table """
    for plex_id, experiment_type in plex_data:
        plex = Plex(plex_id=plex_id, experiment_type=experiment_type)
        session.add(plex)
    session.commit()

    for sample_id, _plex_id in sample_data:
        sample = Sample(sample_id=sample_id, _plex_id=_plex_id)
        session.add(sample)
    session.commit()

    assert session.query(Sample).count() == n

@pytest.mark.parametrize(
    'plex_data,run_data,n',
    [(
        [('1', 'TMT8plex')],
        [('Coding', 1), ('Noncoding', 1)],
        2
    ), (
        [('1', 'TMT8plex'), ('2', 'TMT8plex')],
        [('Coding', 1), ('Noncidng', 1), ('Coding', 2), ('Noncoding', 2)],
        4
    )]
)
def test__run__rows_added_successfully(session: Session,
        plex_data:List[Tuple[str,str]], run_data:List[Tuple[str, int]],
        n:int):
    """ Test rows are added successfully to the Run table """
    for plex_id, experiment_type in plex_data:
        plex = Plex(plex_id=plex_id, experiment_type=experiment_type)
        session.add(plex)
    session.commit()

    for db_id, _plex_id in run_data:
        run = Run(db_id=db_id, _plex_id=_plex_id)
        session.add(run)
    session.commit()

    assert session.query(Run).count() == n

@pytest.mark.parametrize(
    'run_data,peptide_data,n',
    [(
        [('Coding', 1)],
        [('AAAAAR', 'abcde-12345', 1), ('THTHTHK', 'cdefg-54321', 1)],
        2
    ), (
        [('Coding', 1), ('Noncidng', 1)],
        [('AAAAAR', 'abcde-12345', 1), ('THTHTHK', 'cdefg-54321', 2)],
        2
    )]
)
def test__peptide__rows_added_successfully(session: Session,
        run_data:List[Tuple[str, int]], peptide_data:List[Tuple[str,str,int]],
        n:int):
    """ Test rows are added successfully to the Peptide table """
    for db_id, _plex_id in run_data:
        run = Run(db_id=db_id, _plex_id=_plex_id)
        session.add(run)
    session.commit()

    for sequence, uuid, _run_id in peptide_data:
        peptide = Peptide(sequence=sequence, uuid=uuid, _run_id=_run_id)
        session.add(peptide)
    session.commit()

    assert session.query(Peptide).count() == n

    peptide:Peptide = session.query(Peptide).filter(Peptide._peptide_id == 1).first()
    assert peptide.run
