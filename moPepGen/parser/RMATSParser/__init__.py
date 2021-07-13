""" Module for rMATS parsers """
from typing import Iterable
from .RMATSRecord import RMATSRecord
from .SERecord import SERecord
from .A5SSRecord import A5SSRecord
from .A3SSRecord import A3SSRecord
from .MXERecord import MXERecord
from .RIRecord import RIRecord

def parse(path:str, event_type:str) -> Iterable[RMATSRecord]:
    """ parse """
    with open(path, 'rt') as handle:
        # first line is header
        next(handle, None)
        for line in handle:
            if event_type == 'SE':
                yield SERecord.readline(line)
            if event_type == 'A5SS':
                yield A5SSRecord.readline(line)
            if event_type == 'A3SS':
                yield A3SSRecord.readline(line)
            if event_type == 'MXE':
                yield MXERecord.readline(line)
            if event_type == 'RI':
                yield RIRecord.readline(line)
