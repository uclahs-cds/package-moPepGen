""" Peptide """
from __future__ import annotations
from typing import TYPE_CHECKING
from sqlalchemy import Integer, Unicode, ForeignKey
from sqlalchemy.orm import relationship, Mapped, mapped_column
from moPepGen.annotation.models import Base

if TYPE_CHECKING:
    from moPepGen.annotation.models import Run


class Peptide(Base):
    __tablename__ = 'peptide'

    _peptide_id: Mapped[int] = mapped_column(
        '_peptide_id', Integer, primary_key=True
    )
    _run_id: Mapped[int] = mapped_column(
        ForeignKey("run._run_id"),
        nullable=False,
        index=True
    )
    sequence: Mapped[str] = mapped_column(
        'sequence', Unicode, nullable=False
    )
    uuid: Mapped[str] = mapped_column(
        'uuid', Unicode, nullable=False
    )

    run: Mapped[Run] = relationship(back_populates='peptides', uselist=False)

    def __repr__(self) -> str:
        """ repr """
        return f"<Peptide: " +\
            f"_peptide_id={self._peptide_id}, " +\
            f"_run_id={self._run_id}, " +\
            f"sequence={self.sequence}>" +\
            f"uuid={self.uuid}"
