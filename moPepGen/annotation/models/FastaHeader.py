""" FastaHeader """
from __future__ import annotations
from typing import TYPE_CHECKING, List
from sqlalchemy import Integer, Unicode, ForeignKey
from sqlalchemy.orm import relationship, Mapped, mapped_column
from moPepGen.annotation.models import Base

if TYPE_CHECKING:
    from moPepGen.annotation.models import Peptide, Info


class FastaHeader(Base):
    __tablename__ = 'fasta_header'

    _fasta_header_id: Mapped[int] = mapped_column(
        '_fasta_header_id', Integer, primary_key=True
    )
    header: Mapped[str] = mapped_column(
        'header', Unicode, nullable=False
    )
    _peptide_id: Mapped[int] = mapped_column(
        ForeignKey("peptide._peptide_id"),
        nullable=False,
        index=True
    )

    peptide: Mapped[Peptide] = relationship(back_populates='fasta_headers', uselist=False)

    infos: Mapped[List[Info]] = relationship(
        secondary='fasta_header_to_info', back_populates='fasta_headers', uselist=True
    )

    def __repr__(self) -> str:
        """ repr """
        return "<FastaHeader " +\
            f"_fasta_header_id={self._fasta_header_id} " +\
            f"_peptide_id={self._peptide_id} " +\
            f"header={self.header} " +\
            ">"
