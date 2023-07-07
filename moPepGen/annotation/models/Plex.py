""" Run """
from __future__ import annotations
from typing import TYPE_CHECKING, List
from sqlalchemy import Integer, Unicode
from sqlalchemy.orm import relationship, Mapped, mapped_column
from moPepGen.annotation.models import Base

if TYPE_CHECKING:
    from moPepGen.annotation.models import Sample, Run


class Plex(Base):
    __tablename__ = 'plex'

    _plex_id: Mapped[int] = mapped_column(
        '_plex_id', Integer, primary_key=True
    )
    plex_id: Mapped[str] = mapped_column(
        'plex_id', Unicode, nullable=False
    )
    experiment_type: Mapped[str] = mapped_column(
        'experiment_type', Unicode, nullable=False
    )

    samples: Mapped[List[Sample]] = relationship(back_populates='plex', uselist=True)
    runs: Mapped[List[Run]] = relationship(back_populates='plex', uselist=True)

    def __repr__(self) -> str:
        """ repr """
        return "<Plex" +\
            f" _plex_id={self._plex_id} " +\
            f"plex_id={self.plex_id} " +\
            f"experiment_type={self.experiment_type}" +\
            ">"
