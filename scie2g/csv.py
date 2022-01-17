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

"""
The aim of this class is to convert epigenetic events stored in a csv to genes.
"""

from collections import defaultdict
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

from scie2g import Epi2Gene, Epi2GeneException


class Csv(Epi2Gene):

    def __init__(self, filename: str, chr_str: str, start: str, end: str, value: str, header_extra: None,
                 overlap_method='in_promoter',
                 buffer_after_tss=500,
                 buffer_before_tss=2500,
                 buffer_gene_overlap=500,
                 gene_start=3, gene_end=4, gene_chr=2,
                 gene_direction=5, gene_name=0, sep=','
                 ):
        self.chr_str, self.start_str, self.end_str, self.value_str = chr_str, start, end, value
        header = ['idx', self.chr_str, self.start_str, self.end_str, 'gene_idx', value]
        self.header_extra = header_extra or []
        if header_extra:
            header += header_extra
        super().__init__(filename, header, overlap_method=overlap_method,
                         buffer_after_tss=buffer_after_tss,
                         buffer_before_tss=buffer_before_tss,
                         buffer_gene_overlap=buffer_gene_overlap,
                         gene_start=gene_start, gene_end=gene_end, gene_direction=gene_direction, gene_chr=gene_chr,
                         gene_name=gene_name
                         )
        self.filename = filename
        # Set to only look for an in promoter region
        self.location_to_gene_dict, self.loc_idxs_np, self.gene_to_location_dict = defaultdict(list), \
                                                                                   None, \
                                                                                   defaultdict(list)
        self.rows = []
        self.rows_with_genes = []
        self.header = header
        self.sep = sep
        # Check parameters
        try:
            if not self.check_args():
                raise Epi2GeneException(
                    'Parsing arguments failed. Please read detailed error message printed to STDOUT.')
        except:
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

        for h in self.header:
            if h not in cols and h != 'idx' and h not in ["location_chr", "location_start", "location_end", 'gene_idx']:
                self.u.err_p(['check_args: Error: you included header columns that were not in your columns!\n', h,
                              f'\nPossible issues include: Is your separator correct ({self.sep}), does your file have '
                              f'an extra header (if so delete this line - we expect the header to be the first row). '
                              f'Terminating program, please fix your issue and retry.\nFYI Columns in your dataset:',
                              cols, '\nHeader passed: ', self.header])
                raise Epi2GeneException
        # Check the first row has numeric data in the places where it needs to have
        col_dict = {}
        for i, c in enumerate(cols):
            col_dict[c] = i

        # Now we want to check that start and end are both ints
        try:
            start = int(first_row[col_dict[self.start_str]])
        except:
            self.u.err_p(['check_args: Your first row in your start column was not numeric: \n',
                          first_row[col_dict[self.start_str]], '\nExiting... '])
            return False
        try:
            end = int(first_row[col_dict[self.end_str]])
        except:
            self.u.err_p(['check_args: Your first row in your end column was not numeric: ',
                          first_row[col_dict[self.end_str]], '\nExiting... '])
            return False
        return True

    def convert_to_bed(self, df: pd.DataFrame, filename: str, track_name: str, pvalue_str=None, name_str=None):
        chrs, starts, ends, values = df[self.chr_str].values, df[self.start_str].values, df[self.end_str].values, \
                                                 df[self.value_str].values
        if pvalue_str:
            vals = df[pvalue_str].values
            min_value = 0
            if len(np.nonzero(vals)) > 1:
                min_value = np.min(vals[np.nonzero(vals)])
            vals = -1 * np.log10(vals + min_value)
        else:
            vals = values
        if name_str:    # Override using the name
            values = df[name_str].values
        i = 0
        with open(filename, 'w+') as f:
            f.write(f'track name="{track_name}" description="{track_name}" visibility=2 itemRgb="On"\n')
            for row in df.values:
                # Assign the row values to variables
                c = '0,0,255'
                # chr start end name score
                f.write(f'{chrs[i]}    {starts[i]}    {ends[i]}    {values[i]}    {vals[i]}    {starts[i]}    '
                        f'{ends[i]}    {c}\n')
                i += 1

    def format_df(self, df):
        """ Format the csv & sort for efficiency """
        # Also ensure the start and ends are integers
        convert_dict = {self.start_str: int,
                        self.end_str: int,
                        self.chr_str: str}
        df = df.astype(convert_dict)

        # Lastly, we want to at chr to all the seqnames so that they are in ncbi format
        if 'chr' not in str(df[self.chr_str].values[0]):
            df['chr'] = 'chr' + df[self.chr_str].astype(str)
        else:
            df['chr'] = df[self.chr_str].values

        df = df.sort_values(['chr', self.start_str], ascending=[True, True])
        return df

    def _assign_values(self):
        # Since the output file is in a simple csv format we can just read in using pandas
        df = pd.read_csv(self.filename, sep=self.sep)

        df = self.format_df(df)

        self.check_chr(df['chr'][0], self.gene_annot_values[0][self.gene_chr])

        # Now we are ready to iterate through and annotate our DMRs to genes
        num_genes = len(self.gene_annot_values)

        chrs, starts, ends, values = df['chr'].values, df[self.start_str].values, df[self.end_str].values, \
                                                 df[self.value_str].values
        header_idxs = {}
        i = 0
        for c in df.columns:
            if c in self.header_extra:
                header_idxs[c] = i
            i += 1

        rows = df.values
        for loc_idx in tqdm(range(len(df))):
            # Assign the row values to variables
            loc_chr, loc_start, loc_end, loc_value = chrs[loc_idx], starts[loc_idx], ends[loc_idx] + 1, values[loc_idx]
            loc_start, loc_end = int(loc_start), int(loc_end)
            self.cur_loc_idx = loc_idx

            if self.cur_gene_idx >= num_genes:
                self.u.warn_p(['_assign_values: \t Current gene index greater than our length of chrom starts. '
                               'Returning. \n Filename: \t', self.filename])
                break
            # Add in the other header information the user has requested
            loc_args = {'idx': loc_idx, self.chr_str: loc_chr, self.start_str: loc_start,
                        self.end_str: loc_end, 'gene_idx': None, self.value_str: loc_value}
            for h in self.header_extra:
                loc_args[h] = rows[loc_idx][header_idxs[h]]

            gene_chr, gene_start, gene_end, gene_direction = self.get_current_gene_params(loc_chr, loc_start, loc_end)

            # Update row values
            self.cur_loc_start, self.cur_loc_end = loc_start, loc_end

            if loc_chr == gene_chr:
                self.check_for_gene_match(gene_chr, gene_start, gene_end, gene_direction, loc_start,
                                          loc_end, loc_args)
            self.cur_loc_idx = loc_idx

        # Create dataframe based on rows and columns
        self.df = pd.DataFrame(self.rows, columns=self.header)

    def update_loc_value(self, loc_args: dict):
        # First we assign that this gene mapped to this location
        self.assign_loc_value()
        tmp_row = []
        loc_args['gene_idx'] = self.cur_gene_idx  # Update the gene in-case we have changed this
        for h in self.header:
            tmp_row.append(loc_args[h])
        self.rows_with_genes.append(tmp_row)
