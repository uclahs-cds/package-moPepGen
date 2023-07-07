""" Sample """
from __future__ import annotations
from typing import TYPE_CHECKING
from sqlalchemy import Integer, Unicode, ForeignKey
from sqlalchemy.orm import relationship, Mapped, mapped_column
from moPepGen.annotation.models import Base

if TYPE_CHECKING:
    from moPepGen.annotation.models import Plex


class Sample(Base):
    __tablename__ = 'sample'

    _sample_id: Mapped[int] = mapped_column(
        '_sample_id', Integer, primary_key=True
    )
    sample_id: Mapped[str] = mapped_column(
        'sample_id', Unicode, nullable=False
    )
    _plex_id: Mapped[int] = mapped_column(
        ForeignKey("plex._plex_id"),
        nullable=False,
        index=True
    )

    plex: Mapped[Plex] = relationship(back_populates='samples', uselist=False)

    def __repr__(self) -> str:
        """ repr """
        return f"<Sample: _sample_id={self._sample_id}, _plex_id={self._plex_id}>"
