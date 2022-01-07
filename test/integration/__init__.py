""" Integration Test """
from pathlib import Path
import shutil
import unittest
from moPepGen import gtf, seqvar


class TestCaseIntegration(unittest.TestCase):
    """ Test cases for the command line interface """

    def __init__(self, *args, **kwargs):
        """"""
        super().__init__(*args, **kwargs)
        self.base_dir = Path(__file__).parent.parent.absolute()
        self.work_dir = self.base_dir/'work_dir'
        self.data_dir = self.base_dir/'files'

    def setUp(self):
        """ set up working directory """
        super().setUp()
        shutil.rmtree(self.work_dir, ignore_errors=True)
        self.work_dir.mkdir(parents=False, exist_ok=True)

    def tearDown(self):
        """ remove working files """
        super().tearDown()
        shutil.rmtree(self.work_dir, ignore_errors=True)


    def assert_gvf_order(self, gvf_file:Path, anno_gtf:Path) -> None:
        """ """
        anno = gtf.GenomicAnnotation()
        anno.dump_gtf(anno_gtf)

        with open(gvf_file, 'rt') as handle:
            variants = list(seqvar.io.parse(handle))

        gene_ids = []
        for variant in variants:
            gene_id = variant.location.seqname
            if gene_id not in gene_ids:
                gene_ids.append(gene_id)

        genes_rank = anno.get_genes_rank()
        sorted_gene_ids = sorted(gene_ids, key=lambda x:genes_rank[x])
        self.assertEqual(gene_ids, sorted_gene_ids)
