import unittest
from Bio.Seq import Seq
from moPepGen.SeqFeature import FeatureLocation
from moPepGen import aa, dna, svgraph, vep


def create_graph(seq, variants) -> svgraph.TranscriptVariantGraph:
    location = dna.MatchedLocation(
        query=FeatureLocation(start=0, end=len(seq)),
        ref=FeatureLocation(start=0, end=len(seq))
    )
    seq = dna.DNASeqRecordWithCoordinates(
        seq=Seq(seq),
        locations=[location],
        orf=FeatureLocation(start=0, end=len(seq))
    )
    graph = svgraph.TranscriptVariantGraph(seq, '', None)
    records = []
    for start, end, ref, alt, consequences in variants:
        location = FeatureLocation(start=start, end=end)
        records.append(vep.VEPVariantRecord(
            variant=vep.VariantRecord(location=location, ref=ref, alt=alt),
            consequences=consequences
        ))
    graph.create_variant_graph(records)
    return graph

class TestPeptideVariantGraph(unittest.TestCase):
    def test_merge_bubbles(self):
        #    V-P
        #   /
        #  M-VC-PRNLKI-Y-GEV
        #   \  /      \ /
        #    DC        F
        seq = 'ATGGTCTGCCCTAGGAACCTCAAAATCTACGGGGAGGTGTGA'
        variants = [
            (4, 7, 'TCT', 'T', ['deletion', 'frameshifting']),
            (4, 5, 'T', 'A', ['missense']),
            (28, 29, 'A', 'T', ['missense'])
        ]
        dgraph = create_graph(seq, variants)
        dgraph.fit_into_codons()
        pgraph = dgraph.translate('ESNT0001', 'ESNP0001', 'ESNG0001')
        pgraph.to_cleavage_graph('trypsin')
        peptides = pgraph.call_vaiant_peptides()
        return
    

if __name__ == '__main__':
    unittest.main()