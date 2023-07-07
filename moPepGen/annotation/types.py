from __future__ import annotations
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    # pylint: disable=unsupported-binary-operation
    from datetime import date, datetime
    from decimal import Decimal
    from typing import Optional

    from sqlalchemy import (
        Boolean,
        Column,
        Date,
        DateTime,
        Integer,
        Numeric,
        String,
        Unicode,
    )

    IntegerColumn = Union[Column[Integer], int]
    NumericColumn = Union[Column[Numeric], Decimal]
    StringColumn = Union[Column[Unicode], Column[String], str]
    BooleanColumn = Union[Column[Boolean], bool]
    DateColumn = Union[Column[Date], date]
    DateTimeColumn = Union[Column[DateTime], datetime]

    OptionalIntegerColumn = Optional[IntegerColumn]
    OptionalNumericColumn = Optional[NumericColumn]
    OptionalStringColumn = Optional[StringColumn]
    OptionalBooleanColumn = Optional[BooleanColumn]
    OptionalDateColumn = Optional[DateColumn]
    OptionalDateTimeColumn = Optional[DateTimeColumn]
