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
import pandas as pd
import unittest

from scie2g.bed import Bed, Epi2GeneException
from multiprocessing.dummy import Pool as ThreadPool


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

        self.mm10_annot = os.path.join(self.data_dir, 'mmusculus_gene_ensembl-GRCm38.p6.csv')
        self.hg38_annot = os.path.join(self.data_dir, 'hsapiens_gene_ensembl-GRCh38.p13.csv')

    @classmethod
    def teardown_class(self):
        shutil.rmtree(self.tmp_dir)


class TestBed(TestClass):

    def test_bed(self):
        self.setup_class()
        """ Tests the generic function of the bed data """
        bed = Bed(self.h3k27ac)
        # Add the gene annot
        bed.set_annotation_from_file(self.mm10_annot)
        # Now we can run the assign values
        bed.assign_locations_to_genes()
        bed.save_loc_to_csv(f'{self.data_dir}test_bed_h3k27ac_output.csv')
        # Re-read it in and confirm the genes were found (These are the last column in the original file)
        loc_df = pd.read_csv(self.h3k27ac, sep=' ', header=None)
        genes = loc_df[loc_df.columns[-1]].values
        # Compare this to found gene names
        annotated_genes = bed.loc_df['external_gene_name'].values
        found_genes = []
        for g in genes:
            if g in annotated_genes:
                found_genes.append(g)
        assert len(found_genes) == len(genes)
        assert len(annotated_genes) > len(genes)

    def test_bed_annot(self):
        self.setup_class()
        """ Tests the generic function of the bed data """
        bed = Bed('data/test_altmapping_H3K27ac.bed')
        # Add the gene annot
        bed.set_annotation_from_bed_file('data/GCF_000001635.27_GRCm39.bed')
        # Now we can run the assign values
        bed.assign_locations_to_genes()
        bed.save_loc_to_csv(f'{self.data_dir}test_altmapping_H3K27ac_output.csv')
        found_loc_df = bed.loc_df
        print(found_loc_df.head())
        # Re-read it in and confirm the genes were found (These are the last column in the original file)

    def test_bed_unsorted_annot(self):
        self.setup_class()
        # Test it on an unsorted list as well
        """ Tests the generic function of the bed data """
        bed = Bed(self.h3k27ac, overlap_method='overlaps', buffer_after_tss=2500,
                  buffer_before_tss=2500, buffer_gene_overlap=500,
                  gene_start=3, gene_end=4, gene_chr=2,
                  gene_direction=5, gene_name=0,
                  chr_idx=0, start_idx=1,
                  end_idx=2, peak_value=4, header_extra='"0","1","2","3","4"')
        # Add the gene annot
        self.mm10_annot = os.path.join(self.data_dir, 'mmusculus_gene_ensembl-GRCm38.p6_unsorted.csv')
        bed.set_annotation_from_file(self.mm10_annot)
        # Now we can run the assign values
        bed.assign_locations_to_genes()
        bed.save_loc_to_csv(f'{self.data_dir}test_bed_h3k27ac_output.csv')
        # Re-read it in and confirm the genes were found (These are the last column in the original file)
        loc_df = pd.read_csv(self.h3k27ac, sep=' ', header=None)
        genes = loc_df[loc_df.columns[-1]].values
        # Compare this to found gene names
        annotated_genes = bed.loc_df['external_gene_name'].values
        found_genes = []
        for g in genes:
            if g in annotated_genes:
                found_genes.append(g)
        assert len(found_genes) == len(genes)
        assert len(annotated_genes) > len(genes)

    def test_bed_arg_parse_err(self):
        self.setup_class()
        # Test raises an exception when we pass a value that isn't within the range
        with self.assertRaises(Epi2GeneException) as context:
            Bed(self.h3k27ac, overlap_method='overlaps', buffer_after_tss=2500,
                  buffer_before_tss=2500, buffer_gene_overlap=500,
                  gene_start=3, gene_end=4, gene_chr=2,
                  gene_direction=5, gene_name=0,
                  chr_idx=0, start_idx=1,
                  end_idx=2, peak_value=10, header_extra='"0","1","2","3","10"')
            print(context)

        # Test raises an exception when we pass a non numeric start
        with self.assertRaises(Epi2GeneException) as context:
            Bed(self.h3k27ac, overlap_method='overlaps', buffer_after_tss=2500,
                  buffer_before_tss=2500, buffer_gene_overlap=500,
                  gene_start=3, gene_end=4, gene_chr=2,
                  gene_direction=5, gene_name=0,
                  chr_idx=0, start_idx=0,
                  end_idx=2, peak_value=4, header_extra='"0","1","2","3","4"')
            print(context)

        # Test raises an exception when we pass a non numeric end
        with self.assertRaises(Epi2GeneException) as context:
            Bed(self.h3k27ac, overlap_method='overlaps', buffer_after_tss=2500,
                  buffer_before_tss=2500, buffer_gene_overlap=500,
                  gene_start=3, gene_end=4, gene_chr=2,
                  gene_direction=5, gene_name=0,
                  chr_idx=0, start_idx=1,
                  end_idx=0, peak_value=4, header_extra='"0","1","2","3","4"')
            print(context)

        # Test we have a header that is out of range throws an exception
        with self.assertRaises(Epi2GeneException) as context:
            Bed(self.h3k27ac, overlap_method='overlaps', buffer_after_tss=2500,
                  buffer_before_tss=2500, buffer_gene_overlap=500,
                  gene_start=3, gene_end=4, gene_chr=2,
                  gene_direction=5, gene_name=0,
                  chr_idx=0, start_idx=1,
                  end_idx=0, peak_value=4, header_extra='"0","1","2","3","10"')
            print(context)

        # Test if we use an incorrect separator
        with self.assertRaises(Epi2GeneException) as context:
            Bed(self.h3k27ac, overlap_method='overlaps', buffer_after_tss=2500,
                  buffer_before_tss=2500, buffer_gene_overlap=500,
                  gene_start=3, gene_end=4, gene_chr=2,
                  gene_direction=5, gene_name=0,
                  chr_idx=0, start_idx=1,
                  end_idx=0, peak_value=4, header_extra='"0","1","2","3","6"', sep=',')
            print(context)


    def test_chromhmm_bed(self):
        self.setup_class()
        header = "chr,start,end,label,annotation"
        bed = Bed(os.path.join(self.data_dir, 'e16.5_forebrain_15_segments.bed'),
                  header=None, overlap_method='in_promoter',
                  output_bed_file=os.path.join(self.data_dir, 'e16.5_forebrain_15_segments_selectedPeaks.bed'),
                  buffer_after_tss=0, buffer_before_tss=10, buffer_gene_overlap=0,
                  gene_start=3, gene_end=4, gene_chr=2, gene_direction=5, gene_name=0,
                  chr_idx=0, start_idx=1, end_idx=2, peak_value=4, header_extra="3", sep='\t'
                  )
        # Add the gene annot
        bed.set_annotation_from_file(self.mm10_annot)
        # Now we can run the assign values
        bed.assign_locations_to_genes()
        bed.save_loc_to_csv(f'{self.data_dir}test_e16.5_forebrain_15_segments.csv')

    # def test_parallel(self):
    #     data_dir = 'encode/sorted_bed/'
    #     files = os.listdir(data_dir)
    #     files_to_run = []
    #     for filename in files:
    #         files_to_run.append(filename)
    #     # Run in paralell
    #     pool = ThreadPool(12)
    #     results = pool.map(run_bed, files_to_run)