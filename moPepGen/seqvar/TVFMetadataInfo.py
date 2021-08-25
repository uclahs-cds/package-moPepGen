""" TVF Metadata Info """
from typing import Dict, Union


TVF_METADATA_INFO:Dict[str,Dict[str,Union[int,str]]] = {
    'Base': {
        'TRANSCRIPT_ID': {
            'Number': 1,
            'Type': 'String',
            'Description': 'Transcript ID'
        },
        'GENE_SYMBOL': {
            'Number': 1,
            'Type': 'String',
            'Description': 'Gene Symbol'
        },
        'GENOMIC_POSITION': {
            'Number': 1,
            'Type': 'String',
            'Description': 'Genomic Position'
        }
    },
    'Fusion': {
        'ACCEPTER_GENE_ID': {
            'Number': 1,
            'Type': 'String',
            'Description': "3' Accepter Transcript's Gene ID"
        },
        'ACCEPTER_TRANSCRIPT_ID': {
            'Number': 1,
            'Type': 'String',
            'Description': "3' Accepter Transcript's Transcript ID"
        },
        'ACCEPTER_POSITION': {
            'Number': 1,
            'Type': 'Integer',
            'Description': "Position of the break point of the 3' accepter transcript"
        }
    },
    'Insertion': {
        'DONOR_START': {
            'Number': 1,
            'Type': 'Integer',
            'Description': 'Donor Start Position'
        },
        'DONOR_END': {
            'Number': 1,
            'Type': 'Integer',
            'Description': 'Donor End Position'
        }
    },
    'Deletion': {
        'START': {
            'Number': 1,
            'Type': 'Integer',
            'Description': 'Start Position'
        },
        'END': {
            'Number': 1,
            'Type': 'Integer',
            'Description': 'End Position'
        }
    },
    'Substitution': {
        'START': {
            'Number': 1,
            'Type': 'Integer',
            'Description': 'Start Position'
        },
        'END': {
            'Number': 1,
            'Type': 'Integer',
            'Description': 'End Position'
        },
        'DONOR_START': {
            'Number': 1,
            'Type': 'Integer',
            'Description': 'Donor Start Position'
        },
        'DONOR_END': {
            'Number': 1,
            'Type': 'Integer',
            'Description': 'Donor End Position'
        },
        'COORDINATE': {
            'Number': 1,
            'Type': 'String',
            'Description': 'Coordinate for Insertion or Substitution'
        }
    }
}
