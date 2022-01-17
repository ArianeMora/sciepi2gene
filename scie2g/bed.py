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

from collections import defaultdict
import pandas as pd
from tqdm import tqdm
import os

from scie2g import Epi2Gene, Epi2GeneException


class Bed(Epi2Gene):

    def __init__(self, filename: str, header=None, overlap_method='overlaps', output_bed_file=None,
                 buffer_after_tss=500, buffer_before_tss=2500, buffer_gene_overlap=500,
                 gene_start=3, gene_end=4, gene_chr=2, gene_direction=5, gene_name=0,
                 chr_idx=0, start_idx=1, end_idx=2, peak_value=6, header_extra="8,9", sep='\t'):
        super().__init__(filename, header, overlap_method=overlap_method,
                         buffer_after_tss=buffer_after_tss,
                         buffer_before_tss=buffer_before_tss,
                         buffer_gene_overlap=buffer_gene_overlap,
                         gene_start=gene_start, gene_end=gene_end, gene_direction=gene_direction, gene_chr=gene_chr,
                         gene_name=gene_name)
        self.filename = filename
        self.location_to_gene_dict, self.loc_idxs_np, self.gene_to_location_dict = defaultdict(list), None,\
                                                                                   defaultdict(list)
        self.rows = []
        self.rows_with_genes = []
        self.hdr_gene_idx = 1
        self.header = ['peak_idx', 'gene_idx'] + header.split(',') if header else ['peak_idx', 'gene_idx',
                                                                                    'chr', 'start', 'end',
                                                                                   'peak_value'] + header_extra.split(',')
        self.header.append('width')  # This is automatically added for every peak
        self.output_bed_file = None if not output_bed_file else open(output_bed_file, 'a+')
        self.cur_line = None  # Keep track of the current line in the file
        self.chr_idx, self.start_idx, self.end_idx, self.peak_value = chr_idx, start_idx, end_idx, peak_value
        self.hdr_idx = [chr_idx, start_idx, end_idx, peak_value] + [int(h.strip().replace('"', '')) for h
                                                                    in header_extra.split(',')]
        self.sep = sep
        # Check parameters
        if not self.check_args():
            raise Epi2GeneException('Parsing arguments failed. Please read detailed error message printed to STDOUT.')

    def check_args(self):
        # First check the file exists
        if not os.path.isfile(self.filename):
            self.u.err_p(['check_args: CSV: your filename does not exist!', self.filename, '\nExiting... '])
            return False
        # Check supplied header is in the columns
        cols = []
        first_row = None
        with open(self.filename) as f:
            cols = f.readline().split(self.sep)
            cols = [c.replace('"', '').strip() for c in cols]

            first_row = f.readline().split(self.sep)
            first_row = [c.replace('"', '').strip() for c in first_row]
        if len(cols) == 1:
            self.u.warn_p(['check_args: Warning! The program will most likely fail ---> we were only able to read '
                           'in one column'
                           ' from your file: \n', cols, '\n Check that the separator is correct: sep:', self.sep])

        # Check all the columns in the header index can be accessed (i.e. doesn't exceed the params in the bed file
        for h in self.hdr_idx:
            try:
                value = first_row[int(h)]
            except:
                self.u.err_p(['check_args: Header index had a value that was inaccessible: ',
                              h, '\nYour first row is: ', first_row, 'and your header indexs were', self.hdr_idx,
                              '\nExiting... '])
                return False
        try:
            value = first_row[self.peak_value]
        except:
            self.u.err_p(['check_args: Peak value index passed was out of the range for you first row: ',
                          self.peak_value, '\nYour first row is: ', first_row, '\nExiting... '])
            return False

        # Check the first row has numeric data in the places where it needs to have
        try:
            start = int(first_row[self.start_idx])
        except:
            self.u.err_p(['check_args: Your first row in your start column was not numeric: ',
                          first_row[self.start_idx], '\nExiting... '])
            return False
        try:
            end = int(first_row[self.end_idx])
        except:
            self.u.err_p(['check_args: Your first row in your end column was not numeric: ',
                          first_row[self.end_idx], '\nExiting... '])
            return False
        return True

    def _assign_values(self):
        CHROM, START, END, SIGNAL = self.chr_idx, self.start_idx, self.end_idx, self.peak_value
        # Here we have overridden the default method and we're interested in keeping the value
        # that has the largest median value in a single condition
        self.cur_loc_idx, self.cur_gene_idx, count = 0, 0, 0
        first = True
        num_genes = len(self.gene_annot_values)
        with open(self.filename, 'r+') as bedfile:
            for bed_idx, line in enumerate(tqdm(bedfile)):
                self.cur_line = line
                line = line.split('\t')
                loc_chr = line[CHROM]
                if first:
                    self.check_chr(loc_chr, self.gene_annot_values[0][self.gene_chr])
                    first = False
                self.cur_loc_idx = bed_idx
                if self.cur_gene_idx >= num_genes:
                    self.u.warn_p(['read_bed: \t Current gene index greater than our length of chrom starts. '
                                   'Returning. \n Filename: \t', self.filename])
                    break

                # Get the information for this location
                loc_start, loc_end, loc_signal = int(line[START]),  int(line[END]), line[SIGNAL]
                gene_chr, gene_start, gene_end, gene_direction = self.get_current_gene_params(loc_chr, loc_start,
                                                                                              loc_end)
                loc_width = loc_end - loc_start
                self.cur_loc_start, self.cur_loc_end = loc_start, loc_end
                # Create the argument object that will be used.
                loc_args = {'peak_idx': self.cur_loc_idx, 'gene_idx': self.cur_gene_idx}
                for h in self.hdr_idx:
                    loc_args[self.header[len(loc_args)]] = line[int(h)].strip()   # Create the arguments based on the header
                loc_args['width'] = loc_width  # Add in the peak width
                if loc_chr == gene_chr:
                    # Run loop and check if we have a gene match for this location.
                    self.check_for_gene_match(gene_chr, gene_start, gene_end, gene_direction, loc_start,
                                              loc_end, loc_args)

        # Create dataframe based on rows and columns
        self.df = pd.DataFrame(self.rows, columns=self.header)
        # Close the file
        if self.output_bed_file:
            self.output_bed_file.close()

    def update_loc_value(self, loc_args: dict):
        # First we assign that this gene mapped to this location
        self.assign_loc_value()
        tmp_row = []
        loc_args['gene_idx'] = self.cur_gene_idx  # Update the gene in-case we have changed this
        for h in self.header:
            tmp_row.append(loc_args[h])
        self.rows_with_genes.append(tmp_row)
        # If we have an output file to write (which is just the filtered bed file) then write that
        if self.output_bed_file:
            self.output_bed_file.write(self.cur_line)

