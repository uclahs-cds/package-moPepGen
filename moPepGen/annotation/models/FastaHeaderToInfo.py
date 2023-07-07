""" FastaHeaderToInfo """
from __future__ import annotations
from typing import TYPE_CHECKING
from sqlalchemy import Integer, Unicode, ForeignKey
from sqlalchemy.orm import relationship, Mapped, mapped_column
from moPepGen.annotation.models import Base

if TYPE_CHECKING:
    from moPepGen.annotation.models import FastaHeader, Info


class FastaHeaderToInfo(Base):
    __tablename__ = 'fasta_header_to_info'

    _fasta_header_to_info_id: Mapped[int] = mapped_column(
        '_fasta_header_to_info_id', Integer, primary_key=True
    )
    _fasta_header_id: Mapped[int] = mapped_column(
        ForeignKey("fasta_header._fasta_header_id"),
        nullable=False,
        index=True
    )
    _info_id: Mapped[int] = mapped_column(
        ForeignKey("info._info_id"),
        nullable=False,
        index=True
    )

    # FastaHeader: Mapped[FastaHeader] = relationship(back_populates='fasta_header_to_info', uselist=False)
    # Info: Mapped[Info] = relationship(back_populates='fasta_header_to_info', uselist=False)

    def __repr__(self) -> str:
        """ repr """
        return "<FastaHeaderToInfo " +\
            f"_fasta_header_to_info_id={self._fasta_header_to_info_id} " +\
            f"_fasta_header_id={self._fasta_header_id} " +\
            f"_info_id={self._info_id} " +\
            ">"
