###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import shutil
import tempfile
import unittest
import pandas as pd

from scie2g.csv import Csv, Epi2GeneException


class TestClass(unittest.TestCase):

    @classmethod
    def setup_class(self):
        local = True
        # Create a base object since it will be the same for all the tests
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))

        self.data_dir = os.path.join(THIS_DIR, 'data/')
        if local:
            self.tmp_dir = os.path.join(THIS_DIR, 'data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='loc2gene_tmp_')
        # Setup the default data for each of the tests
        self.h3k27ac = os.path.join(self.data_dir, 'test_H3K27ac.bed')
        self.h3k27me3 = os.path.join(self.data_dir, 'test_H3K27me3.bed')
        self.dmrseq = os.path.join(self.data_dir, 'test_dmrseq.csv')
        self.methyl = os.path.join(self.data_dir, 'test_methyl.csv')
        self.methyl_overlaps = os.path.join(self.data_dir, 'test_methyl_overlaps.csv')

        self.mm10_annot = os.path.join(self.data_dir, 'mmusculus_gene_ensembl-GRCm38.p6.csv')
        self.hg38_annot = os.path.join(self.data_dir, 'hsapiens_gene_ensembl-GRCh38.p13.csv')

    @classmethod
    def teardown_class(self):
        shutil.rmtree(self.tmp_dir)


class TestCsv(TestClass):

    def test_csv_limma(self):
        self.setup_class()
        file_name = 'data/Illunima_CPTAC_Infiltrating_duct_carcinoma__NOS_DMCLimma.csv'
        f = pd.read_csv(file_name)
        f['chrs'] = [c.replace('chr', '') for c in f['chrom'].values]
        f.to_csv('data/Illunima_CPTAC_Infiltrating_duct_carcinoma__NOS_DMCLimma_2.csv', index=False)

        mapper = Csv('data/Illunima_CPTAC_Infiltrating_duct_carcinoma__NOS_DMCLimma_2.csv',
                     'chrs', 'start', 'start', 't',
                     ['direction', 'logFC', 'AveExpr', 'P.Value', 'adj.P.Val', 'UCSC_RefGene_Name', 'Strand'],
                     overlap_method='in_promoter',
                     direction_aware=True)
        mapper.set_annotation_from_file('data/NCBIRefSeq_GRCh38_hg38.2.csv')
        mapper.assign_locations_to_genes()
        mapper.save_loc_to_csv(f'{file_name.replace(".csv", "_annot.csv")}')

    def test_csv_dmrseq(self):
        self.setup_class()
        """ Tests the generic function of the bed data """
        f = Csv(self.dmrseq,  'seqnames', 'start', 'end', 'stat', ['index.start', 'index.end', 'qval'])
        # Add the gene annot
        f.set_annotation_from_file(self.hg38_annot)
        # Now we can run the assign values
        f.assign_locations_to_genes()
        f.save_loc_to_csv(f'{self.data_dir}test_dmrseq_output.csv')

    def test_csv_methylkit(self):
        self.setup_class()
        """ Tests the generic function of the csv data """
        f = Csv(self.methyl,  'chr', 'start', 'end', 'meth.diff', ['pvalue', 'qvalue', 'genes'])
        # Add the gene annot
        f.set_annotation_from_file(self.hg38_annot)
        # Now we can run the assign values
        f.assign_locations_to_genes()
        f.save_loc_to_csv(f'{self.data_dir}test_methyl_output.csv')
        f.convert_to_bed(f.loc_df, f'{self.data_dir}test_methyl_output.bed', 'MethylOutput')
        # Check that we got the correct (had an overlap and one that should be missing because it is > 500 bp from
        # the TSS)
        gene_names = f.loc_df['genes'].values
        i = 0
        count_nlgn4y = 0
        for gene in f.loc_df['external_gene_name'].values:
            if 'AC00' not in gene and '-AS2' not in gene:  # Ommit HOXA-AS2 since there is a difference in the annotation for this gene
                assert gene in gene_names[i]  # Predefined gene list!
                if 'NLGN4Y' == gene:
                    count_nlgn4y += 1
            i += 1
        # We should only have 2 of these
        assert count_nlgn4y == 2

    def test_csv_gene_overlaps(self):
        self.setup_class()
        """ This test ensures that we don't miss genes or events that overlaps with the genes """
        f = Csv(self.methyl_overlaps,  'chr', 'start', 'end', 'meth.diff', ['pvalue', 'qvalue', 'description', 'genes'])
        # Add the gene annot
        f.set_annotation_from_file(self.hg38_annot)
        # Now we can run the assign values
        f.assign_locations_to_genes()
        f.save_loc_to_csv(f'{self.data_dir}test_methyl_overlaps_output.csv')
        # Read in the output and check the genes were all gotten (ommit AC00 transcripts)
        loc_df = pd.read_csv(f'{self.data_dir}test_methyl_overlaps_output.csv')
        gene_names = loc_df['genes'].values
        i = 0
        for gene in loc_df['external_gene_name'].values:
            if 'AC00' not in gene and '-AS2' not in gene:  # Ommit HOXA-AS2 since there is a difference in the annotation for this gene
                assert gene in gene_names[i]  # Predefined gene list!

    def test_csv_arg_parse_err(self):
        self.setup_class()
        # Test fails when we are trying to use a string not in the header
        with self.assertRaises(Epi2GeneException) as context:
            f = Csv(self.methyl_overlaps, 'chr', 'startNotReal', 'end', 'meth.diff',
                    ['pvalue', 'qvalue', 'description', 'genes'])
            print(context)

        # Test fails when we are trying to use a string not in the header as the extra info
        with self.assertRaises(Epi2GeneException) as context:
            f = Csv(self.methyl_overlaps, 'chr', 'start', 'end', 'meth.diff',
                    ['pvalue', 'qvalue', 'description', 'startNotReal'])
            print(context)

        # Test fails when we are trying to use a string rather than a integer value for the start
        with self.assertRaises(Epi2GeneException) as context:
            f = Csv(self.methyl_overlaps, 'chr', 'genes', 'end', 'meth.diff',
                    ['pvalue', 'qvalue', 'description', 'genes'])
            print(context)

        # Test fails when we are trying to use a string rather than a integer value for the end
        with self.assertRaises(Epi2GeneException) as context:
            f = Csv(self.methyl_overlaps, 'chr', 'start', 'genes', 'meth.diff',
                    ['pvalue', 'qvalue', 'description', 'genes'])
            print(context)

        # Test if we use an incorrect separator
        with self.assertRaises(Epi2GeneException) as context:
            f = Csv(self.methyl_overlaps, 'chr', 'start', 'genes', 'meth.diff',
                    ['pvalue', 'qvalue', 'description', 'genes'], sep='\t')
            print(context)
