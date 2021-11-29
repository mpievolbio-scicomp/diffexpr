#!/usr/bin/env python
import os
import pandas as pd 
import numpy as np
from diffexpr.py_deseq import py_DESeq2
import warnings
import unittest
warnings.filterwarnings("ignore")
test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'

class DeseqTest(unittest.TestCase):
       
    def setUp(self):
        
        # Set up for htseq input.
        self.htseq_dir = "/home/grotec/mnt/wallace/mnt/micropop/Data/NGS/MPGC/project4293/proc/htseq/chromosome/"
        self.htseq_sample_table = pd.read_csv("/home/grotec/mnt/wallace/mnt/micropop/Data/NGS/MPGC/project4293/ref/4293_sampleTable.csv", 
                                              sep="\t",
                                              header=0,
                                             )
        self.htseq_sample_table.index = self.htseq_sample_table["fileName"]
       
        # Setup for dataframe input.
        self.df = pd.read_csv(test_data_path + '/ercc.tsv', sep='\t')

        self.sample_df = pd.DataFrame({'samplename': self.df.columns}) \
            .query('samplename != "id"')\
            .assign(sample = lambda d: d.samplename.str.extract('([AB])_', expand=False)) \
            .assign(batch = lambda d: d.samplename.str.extract('_([123])', expand=False)) 
        self.sample_df.index = self.sample_df.samplename
        
        
    def test_construction(self):
        """ Test constructing a deseq object with count data and samplesheet."""

        dds = py_DESeq2(count_matrix = self.df,
                   design_matrix = self.sample_df,
                   design_formula = '~ batch + sample',
                   gene_column = 'id')
        
        self.assertIsInstance(dds, py_DESeq2)
            
    @unittest.expectedFailure
    def test_get_count_data(self):
        raise NotImplementedError()
        df = run_deseq._get_count_data(self.test_datadir)
        
    def test_construction_htseq(self):
        """ Test construction with htseq data. """
        
        dds = py_DESeq2(
                   htseq_dir = self.htseq_dir,
                   design_matrix = self.htseq_sample_table,
                   design_formula = '~ condition',
        )
        
        self.assertIsInstance(dds, py_DESeq2)
        
    def test_run_deseq(self):

        df = pd.read_csv(test_data_path + '/ercc.tsv', sep='\t')
        """
            id     A_1     A_2     A_3     B_1     B_2     B_3
        0  ERCC-00002  111461  106261  107547  333944  199252  186947
        1  ERCC-00003    6735    5387    5265   13937    8584    8596
        2  ERCC-00004   17673   13983   15462    5065    3222    3353
        3  ERCC-00009    4669    4431    4211    6939    4155    3647
        4  ERCC-00012       0       2       0       0       0       0
        """


        sample_df = pd.DataFrame({'samplename': df.columns}) \
            .query('samplename != "id"')\
            .assign(sample = lambda d: d.samplename.str.extract('([AB])_', expand=False)) \
            .assign(batch = lambda d: d.samplename.str.extract('_([123])', expand=False)) 
        sample_df.index = sample_df.samplename

        dds = py_DESeq2(count_matrix = df,
                   design_matrix = sample_df,
                   design_formula = '~ batch + sample',
                   gene_column = 'id')

        dds.run_deseq() 
        dds.get_deseq_result()
        res = dds.deseq_result 
        self.assertEqual(res.query('padj < 0.05').shape, (35,7))

        dds.get_deseq_result(contrast = ['sample','B','A'])
        res = dds.deseq_result 
        self.assertEqual(res.query('padj < 0.05').shape, (35,7))

        lfc_res = dds.lfcShrink(coef=4, method='apeglm')
        self.assertEqual(lfc_res[lfc_res.padj < 0.05].shape[0], 35)

        norm_df = dds.normalized_count()
        res.to_csv(test_data_path + '/py_deseq.tsv', index=False, sep='\t')


    def test_result(self):
        os.chdir(os.path.dirname(test_data_path))
        os.system('Rscript deseq.R')

        py = pd.read_table(test_data_path + '/py_deseq.tsv') 
        R = pd.read_table(test_data_path + '/R_deseq.tsv') 

        for col in py.columns:
            if py.columns.dtype == 'float64':
                self.assertTrue(np.all(np.isclose(py[col].fillna(0), R[col].fillna(0))))

        
if __name__ == "__main__":
    unittest.main()