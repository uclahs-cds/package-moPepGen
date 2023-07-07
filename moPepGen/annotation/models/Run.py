""" Run """
from __future__ import annotations
from typing import List, TYPE_CHECKING
from sqlalchemy import Integer, Unicode, ForeignKey
from sqlalchemy.orm import relationship, Mapped, mapped_column
from moPepGen.annotation.models import Base

if TYPE_CHECKING:
    from moPepGen.annotation.models import Plex, Peptide


class Run(Base):
    __tablename__ = 'run'

    _run_id: Mapped[int] = mapped_column(
        '_run_id', Integer, primary_key=True
    )
    _plex_id: Mapped[int] = mapped_column(
        ForeignKey("plex._plex_id"),
        nullable=False,
        index=True
    )
    db_id: Mapped[str] = mapped_column(
        'db_id', Unicode, nullable=False
    )

    plex: Mapped[Plex] = relationship(back_populates='runs', uselist=False)

    peptides: Mapped[List[Peptide]] = relationship(back_populates='run', uselist=True)

    def __repr__(self) -> str:
        """ repr """
        return "<Run " +\
            f"_run_id={self._run_id} " +\
            f"_plex_id={self._plex_id} " +\
            f"db_id={self.db_id} " +\
            ">"
