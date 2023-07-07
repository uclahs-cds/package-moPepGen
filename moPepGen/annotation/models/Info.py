""" Info """
from __future__ import annotations
from typing import TYPE_CHECKING, List
from sqlalchemy import Integer, Unicode, ForeignKey
from sqlalchemy.orm import relationship, Mapped, mapped_column
from moPepGen.annotation.models import Base

if TYPE_CHECKING:
    from moPepGen.annotation.models import FastaHeader


class Info(Base):
    __tablename__ = 'info'

    _info_id: Mapped[int] = mapped_column(
        '_info_id', Integer, primary_key=True
    )

    info: Mapped[str] = mapped_column(
        'info', Unicode, nullable=False
    )

    fasta_headers: Mapped[List[FastaHeader]] = relationship(
        secondary='fasta_header_to_info', back_populates='infos', uselist=True
    )

    def __repr__(self) -> str:
        """ repr """
        return "<Info " +\
            f"_info_id={self._info_id} " +\
            f"info={self.info} " +\
            ">"
