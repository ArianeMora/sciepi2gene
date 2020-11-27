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
import pandas as pd
import pytest
import os
import shutil
import tempfile
import unittest

from scie2g import Epi2Gene, Epi2GeneException


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
            self.tmp_dir = tempfile.mkdtemp(prefix='scie2g_tmp_')
        # Setup the default data for each of the tests
        self.h3k27ac = os.path.join(self.data_dir, 'test_H3K27ac.bed')
        self.h3k27me3 = os.path.join(self.data_dir, 'test_H3K27me3.bed')
        self.dmrseq = os.path.join(self.data_dir, 'test_dmrseq.csv')

        self.mm10_annot = os.path.join(self.data_dir, 'mmusculus_gene_ensembl-GRCm38.p6.csv')
        self.hg38_annot = os.path.join(self.data_dir, 'hsapiens_gene_ensembl-GRCh38.p13.csv')

    @classmethod
    def teardown_class(self):
        shutil.rmtree(self.tmp_dir)


class TestEpi2Gene(TestClass):

    def test_generate_gene_info(self):
        self.setup_class()
        l2g = Epi2Gene('', [])
        l2g.set_annotation_using_biomart('ENSEMBL_MART_ENSEMBL', 'mmusculus_gene_ensembl',
                                         {'ensembl_gene_id': 'ENSMUSG00000029844,ENSMUSG00000032446,'
                                                             'ENSMUSG00000020875,ENSMUSG00000038210'})
        assert l2g.gene_annot_values[0][0] == 'ENSMUSG00000020875'
        assert len(l2g.gene_annot_df) == 4
        # Check that we can save this
        l2g.save_annotation(self.tmp_dir)
        files = os.listdir(self.tmp_dir)
        found = False
        for f in files:
            if 'mmusculus_gene_ensembl' in f:
                found = True
        assert found

    def test_get_start_end(self):
        l2g = Epi2Gene('', [])
        start, end = 0, 1
        start1, end1 = l2g.get_start_end(start, end)
        assert start == start1
        # Check when we have it flipped
        start2, end2 = l2g.get_start_end(end, start)
        assert start2 == start

    def test_in_gene_promoter(self):
        """ in_promotor(self, gene_start: int, gene_end: int, gene_direction: int, loc_start: int, loc_end: int) """
        l2g = Epi2Gene('', [])
        # Set overlaps so it's easy to work out
        l2g.buffer_gene_overlap = 10
        l2g.buffer_before_tss = 5
        # Should overlap
        was_in = l2g.in_promotor(0, 40, 1, 0, 3)
        assert was_in
        # Check on the edge of the boundry
        was_in = l2g.in_promotor(10, 40, 1, 4, 5)
        assert was_in
        # Check just out the boundry
        was_in = l2g.in_promotor(10, 40, 1, 3, 4)
        assert not was_in
        # Do the same for reverse transcribed
        was_in = l2g.in_promotor(10, 40, -1, 10, 20)
        assert not was_in  # This one we are looking at the wrong side (would need to be within dist of the 40)
        was_in = l2g.in_promotor(10, 40, -1, 39, 41)
        assert was_in
        was_in = l2g.in_promotor(10, 40, -1, 46, 50)  # Just out
        assert not was_in
        was_in = l2g.in_promotor(10, 40, -1, 45, 50)  # Just in
        assert was_in
        was_in = l2g.in_promotor(10, 40, -1, 30, 35)  # in on GB side
        assert was_in
        was_in = l2g.in_promotor(10, 40, -1, 25, 30)  # in on GB side
        assert was_in
        was_in = l2g.in_promotor(10, 40, -1, 25, 29)  # in on GB side
        assert not was_in

    def test_overlaps_gene(self):
        """overlaps_gene(self, gene_start: int, gene_end: int, loc_start: int, loc_end: int)"""
        # Similar to above but allows for overlap of the whole gene
        l2g = Epi2Gene('', [])
        # Set overlaps so it's easy to work out
        l2g.buffer_after_tss = 10
        l2g.buffer_before_tss = 5

        was_in = l2g.overlaps_gene(10, 40, 1, 10, 20)
        assert was_in
        was_in = l2g.overlaps_gene(10, 40, 1, 4, 20)
        assert was_in
        was_in = l2g.overlaps_gene(10, 40, 1, 0, 5)
        assert was_in
        was_in = l2g.overlaps_gene(10, 40, 1, 0, 4)
        assert not was_in
        was_in = l2g.overlaps_gene(10, 40, 1, 40, 45)
        assert was_in
        was_in = l2g.overlaps_gene(10, 40, 1, 50, 51)
        assert was_in
        was_in = l2g.overlaps_gene(10, 40, 1, 51, 52)
        assert not was_in
        # Check for reverse gene
        was_in = l2g.overlaps_gene(10, 40, -1, 50, 52)
        assert not was_in
        was_in = l2g.overlaps_gene(10, 40, -1, 10, 20)
        assert was_in
        was_in = l2g.overlaps_gene(10, 40, -1, 0, 1)
        assert was_in
        was_in = l2g.overlaps_gene(20, 40, -1, 0, 1)
        assert not was_in

    def test_overlaps(self):
        l2g = Epi2Gene('', [])
        # Check a specific one
        was_in = l2g.overlaps(45060061, 47301371, 1, 45060358, 45060620)
        assert was_in

        # Set overlaps so it's easy to work out
        l2g.buffer_after_tss = 10
        l2g.buffer_before_tss = 5
        l2g.buffer_gene_overlap = 5

        # Check for in promoter & then check overlaps body
        l2g.overlap_method = 'in_promoter'
        was_in = l2g.overlaps(10, 40, 1, 4, 5)
        assert was_in
        # Check just out the boundry
        was_in = l2g.overlaps(10, 40, 1, 3, 4)
        assert not was_in

        # Do the same for reverse transcribed
        was_in = l2g.overlaps(10, 40, -1, 10, 20)
        assert not was_in  # This one we are looking at the wrong side (would need to be within dist of the 40)
        # Check for overlaps
        l2g.overlap_method = 'overlaps'
        was_in = l2g.overlaps(10, 40, -1, 10, 20)  # This one wasn't in the promoter but was overlapping the gene!
        assert was_in
        was_in = l2g.overlaps(10, 40, -1, 10, 20)
        assert was_in
        was_in = l2g.overlaps(10, 40, -1, 0, 1)
        assert was_in
        was_in = l2g.overlaps(20, 40, -1, 0, 1)
        assert not was_in

    def test_assign_loc_value(self):
        l2g = Epi2Gene('', [])
        l2g.cur_loc_idx, l2g.cur_chr, l2g.cur_loc_start, l2g.cur_loc_end, l2g.cur_gene_idx = 1, 2, 3, 4, 5
        l2g.location_to_gene_dict = {1: []}
        l2g.gene_to_location_dict = {5: []}
        l2g.assign_loc_value()
        assert l2g.gene_to_location_dict[5][0] == 1
        assert l2g.location_to_gene_dict[1][0] == 5

    def test_update_loc_value(self):
        l2g = Epi2Gene('', [])
        l2g.cur_loc_idx, l2g.cur_chr, l2g.cur_loc_start, l2g.cur_loc_end, l2g.cur_gene_idx = 1, 2, 3, 4, 5
        args = {'A': 'he', 'B': 'she'}
        l2g.rows_with_genes = []
        l2g.location_to_gene_dict = {1: []}
        l2g.gene_to_location_dict = {5: []}
        l2g.update_loc_value(args)
        assert l2g.rows_with_genes[0] == [1, 2, 3, 4, 5, 'he', 'she']
        assert l2g.gene_to_location_dict[5][0] == 1
        assert l2g.location_to_gene_dict[1][0] == 5

    def test_check_chr(self):
        l2g = Epi2Gene('', [])
        with pytest.raises(Epi2GeneException):
            l2g.check_chr('chr1', '2')
        with pytest.raises(Epi2GeneException):
            l2g.check_chr('chr3', '1')
        with pytest.raises(Epi2GeneException):
            l2g.check_chr('chr1', 1)
        l2g.gene_annot_df = pd.DataFrame()
        l2g.gene_annot_df['chromosome_name'] = ['1', '2']
        assert not l2g.check_chr('chr2', '1')

    def test_get_columns_in_gene_info(self):
        l2g = Epi2Gene('', [])
        with pytest.raises(Epi2GeneException):
            l2g.get_columns_in_gene_info()
        # Add gene annot info
        l2g.set_annotation_using_biomart('ENSEMBL_MART_ENSEMBL', 'mmusculus_gene_ensembl',
                                     {'ensembl_gene_id': 'ENSMUSG00000029844,ENSMUSG00000032446,'
                                                         'ENSMUSG00000020875,ENSMUSG00000038210'})
        columns = l2g.get_columns_in_gene_info()
        assert columns == ['ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position',
                           'end_position', 'strand']
