""" ExPASy's PeptideCutter rules. Adopted from pyteomics at:
https://github.com/levitsky/pyteomics """
from typing import Dict, Tuple


EXPASY_RULES: Dict[str, str] = {
    'arg-c': r'R',
    'asp-n': r'\w(?=D)',
    'bnps-skatole': r'W',
    'caspase 1': r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2': r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3': r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4': r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5': r'(?<=[LW]EH)D',
    'caspase 6': r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7': r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8': r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9': r'(?<=LEH)D',
    'caspase 10': r'(?<=IEA)D',
    'chymotrypsin high specificity': r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity': r'([FLY](?=[^P]))|(W(?=[^MP]))|'
    r'(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain': r'R',
    'cnbr': r'M',
    'enterokinase': r'(?<=[DE]{3})K',
    'factor xa': r'(?<=[AFGILTVM][DE]G)R',
    'formic acid': r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b': r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc': r'K',
    'lysn': r'\w(?=K)',
    'ntcb': r'\w(?=C)',
    'pepsin ph1.3':  r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
    r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'pepsin ph2.0': r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
    r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k': r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin': r'[^DE](?=[AFILMV])',
    'thrombin': r'((?<=G)R(?=G))|'
    r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin': r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    'trypsin_exception': r'((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|'
    r'((?<=R)R(?=[HR]))',
}

EXPASY_RULES2: Dict[str, str] = {
    'arg-c': r'R',
    'asp-n': r'\wD',
    'bnps-skatole': r'W',
    'caspase 1': r'[FWYL]\w[HAT]D[^PEDQKR]',
    'caspase 2': r'DVAD[^PEDQKR]',
    'caspase 3': r'DMQD[^PEDQKR]',
    'caspase 4': r'LEVD[^PEDQKR]',
    'caspase 5': r'[LW]EHD',
    'caspase 6': r'VE[HI]D[^PEDQKR]',
    'caspase 7': r'DEVD[^PEDQKR]',
    'caspase 8': r'[IL]ETD[^PEDQKR]',
    'caspase 9': r'LEHD',
    'caspase 10': r'IEAD',
    'chymotrypsin high specificity': r'([FY][^P])|(W[^MP])',
    'chymotrypsin low specificity': r'([FLY][^P])|(W[^MP])|'
    r'(M[^PY])|(H[^DMPW])',
    'clostripain': r'R',
    'cnbr': r'M',
    'enterokinase': r'[DE]{3}K',
    'factor xa': r'[AFGILTVM][DE]GR',
    'formic acid': r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b': r'IEPD',
    'hydroxylamine': r'NG',
    'iodosobenzoic acid': r'W',
    'lysc': r'K',
    'lysn': r'\wK',
    'ntcb': r'\wC',
    'pepsin ph1.3':  r'([^HKR][^P][^R][FL][^P])|([^HKR][^P][FL]\w[^P])',
    'pepsin ph2.0': r'([^HKR][^P][^R][FLWY][^P])|([^HKR][^P][FLWY]\w[^P])',
    'proline endopeptidase': r'[HKR]P[^P]',
    'proteinase k': r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'[^E]E',
    'thermolysin': r'[^DE][AFILMV]',
    'thrombin': r'(GRG)|([AFGILTVM][AFGILTVWA]PR[^DE][^DE])',
    'trypsin': r'([KR][^P])|(WKP)|(MRP)',
    'trypsin_exception': r'([CD]KD)|(CK[HY])|(CRK)|(RR[HR])',
}

EXPASY_RULES_WINGS_SIZE: Dict[str, Tuple[int, int]] = {
    'arg-c': (1, 0),
    'asp-n': (0, 1),
    'bnps-skatole': (1, 0),
    'caspase 1': (4, 1),
    'caspase 2': (2, 1),
    'caspase 3': (2, 1),
    'caspase 4': (2, 1),
    'caspase 5': (2, 0),
    'caspase 6': (3, 1),
    'caspase 7': (2, 1),
    'caspase 8': (4, 1),
    'caspase 9': (4, 0),
    'caspase 10': (4, 0),
    'chymotrypsin high specificity': (1, 1),
    'chymotrypsin low specificity': (1, 1),
    'clostripain': (1, 0),
    'cnbr': (1, 0 ),
    'enterokinase': (4, 0),
    'factor xa': (4, 0),
    'formic acid': (1, 0),
    'glutamyl endopeptidase': (1, 0),
    'granzyme b': (4, 0),
    'hydroxylamine': (1, 1),
    'iodosobenzoic acid': (1, 0),
    'lysc': (1, 0),
    'lysn': (1, 1),
    'ntcb': (0, 1),
    'pepsin ph1.3': (3, 2),
    'pepsin ph2.0': (3, 2),
    'proline endopeptidase': (2, 1),
    'proteinase k': (1, 0),
    'staphylococcal peptidase i': (2, 0),
    'thermolysin': (1, 1),
    'thrombin': (4, 2),
    'trypsin': (2, 1),
    'trypsin_exception': (2, 1),
}
