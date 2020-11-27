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
import numpy as np
import pandas as pd
from typing import Tuple
from tqdm import tqdm

from sciutil import SciUtil, SciException
from scibiomart import SciBiomartApi

# Errors
errors = {'GENE_ANNOT_ERR': 'Err: assign_locations_to_genes, You have not initialised a gene information object yet.'
                            '\nSee set_annotation_from_file if you already have an annotation file or '
                            'set_annotation_using_biomart to make a new one (uses biomart).\nDetails to '
                            ' make an annotation file can be found in the package scibiomart.'}


class Epi2GeneException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class Epi2Gene:

    def __init__(self, filename: str, header: list, overlap_method='in_promoter', buffer_after_tss=500,
                 buffer_before_tss=2500,
                 buffer_gene_overlap=500, gene_start=3, gene_end=4, gene_chr=2,
                 gene_direction=5, gene_name=0, gene_id_type=None, output_dir='.', sciutil=None,
                 hdr_gene_idx=4):

        self.u = SciUtil() if sciutil is None else sciutil
        # Settings for choosing the overlap
        self.overlap_method, self.buffer_gene_overlap = overlap_method, buffer_gene_overlap
        self.buffer_after_tss, self.buffer_before_tss = buffer_after_tss, buffer_before_tss
        # Settings for the biomart
        self.gene_start, self.gene_end, self.gene_chr = gene_start, gene_end, gene_chr
        self.gene_direction, self.gene_name, self.gene_id_type = gene_direction, gene_name, gene_id_type
        self.output_dir, self.filename = output_dir, filename
        self.location_to_gene_dict, self.gene_to_location_dict, self.df = None, None, None
        self.cur_gene_idx, self.cur_loc_idx, self.cur_loc_start, self.cur_loc_end, self.cur_chr = 0, 0, 0, 0, None
        self.rows_with_genes, self.header, self.loc_df = [], header, None
        self.hdr_gene_idx, self.biomart, self.gene_info_df = hdr_gene_idx, None, None
        self.gene_annot_df, self.gene_annot_values = pd.DataFrame(), []
        self.num_genes = 0
        # set the gene information of throw warning that they need to generate this

    def assign_locations_to_genes(self):
        """ Wrapper for the main _assign values method, here we just perform some generic tests & setups """
        if len(self.gene_annot_df) < 1:
            self.u.err_p([errors.get('GENE_ANNOT_ERR')])
            return
        self.loc_idxs_np = np.full(len(self.gene_annot_df), -1)
        # Run assignment
        self._assign_values()

    def _assign_values(self):
        self.u.warn_p(['Warn: _assign_values not performed. Please use the correct wrapper for your file type.'
                       '\nDMRseq, Generic, or Bed.'])
        return
    """
    -----------------------------------------------------------------
    Generation of gene data.
    -----------------------------------------------------------------
    """
    def set_annotation_from_file(self, gene_annotation_file):
        # Assume the file is the correct format (i.e. from scibiomart)
        self.gene_annot_df = pd.read_csv(gene_annotation_file)
        # Ensure it is sorted
        self.biomart = SciBiomartApi()
        self.gene_annot_df = self.biomart.sort_df_on_starts(self.gene_annot_df)
        # Gene information is just all the values from our annot df
        self.gene_annot_values = self.gene_annot_df.values
        self.num_genes = len(self.gene_annot_values)

    def set_annotation_from_bed_file(self, gene_annotation_file):
        # Assume the file is the correct format (i.e. from scibiomart)
        self.gene_annot_df = pd.read_csv(gene_annotation_file, sep='\t')
        self.gene_annot_df.columns = ['external_gene_name', 'chr', 'start', 'end', 'direction', 'sample_type']
        # Ensure it is sorted
        self.biomart = SciBiomartApi()
        self.gene_annot_df = self.biomart.sort_df_on_starts(self.gene_annot_df)
        # Gene information is just all the values from our annot df
        self.gene_annot_values = self.gene_annot_df.values
        self.num_genes = len(self.gene_annot_values)

    def set_annotation_using_biomart(self, mart: str, dataset: str, filter_dict=None):
        self.biomart = SciBiomartApi()
        self.biomart.set_mart(mart)
        self.biomart.set_dataset(dataset)
        # Finally lets make our regions of interest
        self.gene_annot_df = self.biomart.run_default(filter_dict)
        # Sort the values
        self.gene_annot_df = self.biomart.sort_df_on_starts(self.gene_annot_df)

        # Gene information is just all the values from our annot df
        self.gene_annot_values = self.gene_annot_df.values
        self.num_genes = len(self.gene_annot_values)

    def save_annotation(self, output_dir=None):
        output_dir = output_dir or self.output_dir
        self.biomart.save_as_csv(self.gene_annot_df, output_dir)

    """
    -----------------------------------------------------------------
    Functions for checking overlaps.
    -----------------------------------------------------------------
    """
    @staticmethod
    def get_start_end(start: int, end: int) -> Tuple[int, int]:
        """
        Returns the start and end of a gene. Here we want to switch these around if we have a reverse
        transcribed gene.

        Parameters
        ----------
        start:      int: the starting index of a gene
        end:        int: the ending index of a gene

        Returns
        -------
        start, end
        """
        if start < end:
            return start, end
        return end, start

    def in_promotor(self, gene_start: int, gene_end: int, gene_direction: int, loc_start: int, loc_end: int) -> bool:
        """
        Check if our loc is within the buffer size of the promotor region of the gene of interest.
        Here we check first which direction the gene is transcribing then we check if location/region is
        less than gene TSS + buffer and gene end is greater than this value.

        Parameters
        ----------
        self
        gene_start:         int: gene start
        gene_end:           int: gene end
        gene_direction:     int: [-1, 1] direction of gene
        loc_start:          int: location of the gene start
        loc_end:            int: location of the gene end

        Returns
        -------
        Boolean
        """
        # i.e. forward transcribing gene
        if gene_direction > 0:
            if (loc_start <= (gene_start + self.buffer_gene_overlap)) and \
                    (loc_end >= (gene_start - self.buffer_before_tss)):
                return True
        elif (loc_start <= gene_end + self.buffer_before_tss) and (loc_end >= gene_end - self.buffer_gene_overlap):
            return True
        return False

    def overlaps_gene(self, gene_start: int, gene_end: int, gene_direction: int, loc_start: int, loc_end: int) -> bool:
        " Flip the starts and ends if reversed --> Tests for overlap with gene body. "
        gene_start, gene_end = self.get_start_end(gene_start, gene_end)
        if gene_direction > 0:
            if (loc_start <= gene_end + self.buffer_after_tss) and loc_end >= gene_start - self.buffer_before_tss:
                return True
        else:
            if (loc_start <= gene_end + self.buffer_before_tss) and loc_end >= gene_start - self.buffer_after_tss:
                return True
        return False

    def overlaps(self, gene_start: int, gene_end: int, gene_direction: int, loc_start_i: int, loc_end_i: int) -> bool:
        """
        Tests if the location overlaps the gene. This has two methods, first checking if it overlaps in
        the promotor, second if it overlaps the whole gene body.

        Parameters
        ----------
        gene_start:     int: annotated gene start
        gene_end:       int: annotated gene end
        gene_direction: int: either -1 or 1 depending if it is reverse transcribed
        loc_start_i:    int: start value of the location
        loc_end_i:      int: end value of the location (i.e. CpG island)

        Returns
        -------
        Boolean as to whether it overlaps.
        """
        if self.overlap_method == 'in_promoter':
            return self.in_promotor(gene_start, gene_end, gene_direction, loc_start_i, loc_end_i)
        elif self.overlap_method == 'overlaps':
            return self.overlaps_gene(gene_start, gene_end, gene_direction, loc_start_i, loc_end_i)
        else:
            msg = self.u.msg.msg_arg_err("overlaps", "self.overlap_method", self.overlap_method,
                                         ['in_promoter', 'overlaps'])
            self.u.err_p([msg])
            raise Epi2GeneException(msg)

    """
    -----------------------------------------------------------------
    Functions for running in the loop.
    -----------------------------------------------------------------
    """
    def update_loc_value(self, loc_args: dict):
        """
        Update the location's value with the new score & value (if this is better than the stored one).

        Here we first assign the location value (to keep track of all locations assigned to genes).
        Second, we look test if the new score is better than the previous score and if so we replace the
        CpG island we have stored.

        Parameters
        ----------
        loc_args:      An ordered dictionary with the arguments stored that one wants added ot the row

        Returns
        -------
        stored_score, stored_value
        """
        self.assign_loc_value()
        tmp_row = [self.cur_loc_idx, self.cur_chr, self.cur_loc_start, self.cur_loc_end, self.cur_gene_idx]
        for arg in loc_args:
            tmp_row.append(loc_args[arg])
        self.rows_with_genes.append(tmp_row)

    def check_for_gene_match(self, gene_chr: str, gene_start: int, gene_end: int, gene_direction: int,
                             loc_start: int, loc_end: int, loc_args: dict) -> None:
        """
        This is the general search loop that we go through each time we are trying to assign a location to a gene.
        """
        self.cur_chr = gene_chr

        # Check if we have locs around or accross the starting area
        # First check if this loc is overlapping the gene, if it is, then
        if self.overlaps(gene_start, gene_end, gene_direction, loc_start, loc_end):
            # Add this to our rows
            # Before we move onto the next loc, we want to see if there are any other genes that
            # meet the same criteria
            self.update_loc_value(loc_args)
            # We've already assigned this gene, lets move up 1
            prev_gene_idx = self.cur_gene_idx
            self.cur_gene_idx += 1
            self.find_genes_in_loc(loc_start, loc_end, loc_args)
            self.cur_gene_idx = prev_gene_idx
        else:
            prev_gene_idx = self.cur_gene_idx
            self.cur_gene_idx += 1
            self.find_genes_in_loc(loc_start, loc_end, loc_args)
            # Again we may need to re visit this gene
            self.cur_gene_idx = prev_gene_idx

    def assign_loc_value(self):
        """
        Creates a dictionary using the location indicies as the keys and the gene index as the value
        Returns
        -------

        """
        # Add that gene to the location dictionary
        self.location_to_gene_dict[self.cur_loc_idx].append(self.cur_gene_idx)
        # Assign the location to the gene
        self.gene_to_location_dict[self.cur_gene_idx].append(self.cur_loc_idx)

    def get_current_gene_params(self, loc_chr: str,  loc_start_i: int, loc_end_i: int) -> Tuple[str, int, int, int]:
        """
        This is run in each of the main file loops. Here we check whether the current gene chr matches the current
        location chr. If not, we need to iterate through the gene file until we match the location chr. Since we
        assume that the location file is sorted this means the gene index always has to catch up.

        Parameters
        ----------
        loc_chr
        gene_chr

        Returns
        -------

        """
        gene_chr = self.gene_annot_values[self.cur_gene_idx][self.gene_chr]
        if not loc_chr == gene_chr:
            # We need to iterate through until we get to the next chromosone
            for i in range(0, len(self.gene_annot_values)):
                gene_chr = self.get_gene_chr(i)
                if gene_chr == loc_chr:
                    self.cur_chr = loc_chr
                    self.cur_gene_idx = i
                    break
        # We also need to see if this is before the current location if so, we have passed the TSS and need to move
        # onto the next gene
        for i in range(self.cur_gene_idx, self.num_genes):
            gene_end = self.gene_annot_values[i][self.gene_end]
            gene_direction = self.gene_annot_values[i][self.gene_direction]
            gene_start = self.gene_annot_values[i][self.gene_start]

            if self.overlaps(gene_start, gene_end, gene_direction, loc_start_i, loc_end_i):
                # Update the gene index
                self.cur_gene_idx = i
                break
                # Check if we have passed our loc (i.e. both the loc start and end are before our gene
                # We need to do this in a way where we are aware of the reverse transcribed genes
            if self.overlap_method == 'in_promoter':
                if (gene_direction > 0 and (loc_end_i < (gene_start - self.buffer_before_tss))) or (gene_direction < 0 and
                                                                                                  (loc_end_i < (gene_end -
                                                                                                   self.buffer_gene_overlap))):
                    self.cur_gene_idx = i
                    # We want to see if there are any locs for this gene
                    break
            elif self.overlap_method == 'overlaps':
                if gene_direction > 0 and (loc_end_i < (gene_start - self.buffer_before_tss)):
                    # Here our gene has overshot the mark, so we go back and see if there are any new locations in
                    # this gene
                    self.cur_gene_idx = i
                    break
                elif gene_direction < 0 and (loc_end_i < (gene_start - self.buffer_after_tss)):
                    # Here our gene has overshot the mark, so we go back and see if there are any new locations in
                    # this gene
                    self.cur_gene_idx = i
                    # We want to see if there are any locs for this gene
                    break
        # Set current chr
        gene_chr = self.gene_annot_values[self.cur_gene_idx][self.gene_chr]
        gene_end = self.gene_annot_values[self.cur_gene_idx][self.gene_end]
        gene_direction = self.gene_annot_values[self.cur_gene_idx][self.gene_direction]
        gene_start = self.gene_annot_values[self.cur_gene_idx][self.gene_start]
        return gene_chr, gene_start, gene_end, gene_direction

    def check_chr(self, loc_chr, gene_chr):
        """
        Checks if the chr format is the same for each of the methods. i.e. either both don't have chr or both have chr.
        Parameters
        ----------
        loc_chr
        gene_chr

        Returns
        -------

        """
        try:
            if 'chr' in loc_chr:
                # Since we're assuming that they have used the scibiomart package (which doesn't have chrs by default)
                # we need to add those
                self.gene_annot_df['chromosome_name'] = 'chr' + self.gene_annot_df['chromosome_name'].astype(str)
                self.gene_annot_values = self.gene_annot_df.values
                msg = f'Warning: Your input file used different chr conventions to ensembl: {loc_chr} vs {gene_chr} ' \
                      f'\nWe have added chr prefix to the ensembl chrs.' \
                      f'\nfile: {self.filename}'
                self.u.warn_p([msg])
        except Exception as e:
            # Catch any exceptions and raise an error
            msg = f'Issue with your chromosone conventions: {loc_chr} vs {gene_chr} ' \
                  f'\nfile: {self.filename}'
            self.u.err_p([msg])
            raise Epi2GeneException(msg)

    def find_genes_in_loc(self, loc_start_i: int, loc_end_i: int, loc_args: dict) -> None:
        """
        Searchers the region for any genes. This is done so that if we have a broad location e.g.
        H3K27me3 that might have a peak accross many genes we want to indicate that these genes are associated
        with the region.

        Parameters
        ----------
        loc_start_i:        int: start index of the location
        loc_end_i:          int: end index of the location
        loc_args:           dict: arguments for the location

        Returns
        -------
        None
        """
        for i in range(self.cur_gene_idx, self.num_genes):
            gene_start = self.get_gene_start(i)
            gene_chr = self.get_gene_chr(i)
            gene_end = self.get_gene_end(i)
            gene_direction = self.get_gene_direction(i)
            # We have moved onto the next chromosone, so return.
            if gene_chr != self.cur_chr:
                self.cur_gene_idx = i
                break
            if self.overlaps(gene_start, gene_end, gene_direction, loc_start_i, loc_end_i):
                # Update the gene index
                self.cur_gene_idx = i
                self.update_loc_value(loc_args)
                # Check if we have passed our loc (i.e. both the loc start and end are before our gene
                # We need to do this in a way where we are aware of the reverse transcribed genes
            if self.overlap_method == 'in_promoter':
                if (gene_direction > 0 and (loc_end_i < (gene_start - self.buffer_before_tss))) or (gene_direction < 0 and
                                                                                                  (loc_end_i < (gene_end -
                                                                                                   self.buffer_before_tss))):
                    self.cur_gene_idx = i
                    # We want to see if there are any locs for this gene
                    break
            elif self.overlap_method == 'overlaps':
                if gene_direction > 0 and (loc_end_i < (gene_start - self.buffer_before_tss)):
                    # Here our gene has overshot the mark, so we go back and see if there are any new locations in
                    # this gene
                    self.cur_gene_idx = i
                    break
                elif gene_direction < 0 and (loc_end_i < (gene_start - self.buffer_after_tss)):
                    # Here our gene has overshot the mark, so we go back and see if there are any new locations in
                    # this gene
                    self.cur_gene_idx = i
                    # We want to see if there are any locs for this gene
                    break

    """
    -----------------------------------------------------------------
    Functions for saving.
    -----------------------------------------------------------------
    """
    def get_columns_in_gene_info(self):
        """
        Gets the columns in the gene information dictionary
        Returns
        -------

        """
        if len(self.gene_annot_df) < 1:
            msg = errors.get('GENE_ANNOT_ERR')
            self.u.err_p([msg])
            raise Epi2GeneException(msg)

        return list(self.gene_annot_df.columns)

    def get_gene_info_as_df(self, columns=None, keep_unassigned=True) -> pd.DataFrame:
        """
        Convert the gene information info a dataframe that can be interrogated.

        Parameters
        ----------
        columns:    list: the columns we want to save, if None all are kept.

        Returns
        -------
        DataFrame
        """
        columns = columns if columns is not None else self.get_columns_in_gene_info()
        gene_info_df = pd.DataFrame()
        print("Running assign_loc_info_to_gene_list")
        for c in columns:
            gene_info_df[c] = self.gene_annot_df[c].values

        new_rows = []
        new_columns = columns + list(self.df.columns)
        # Assign gene values to the locations that had genes assigned.
        print("Running assign_gene_info_to_loc_df")
        df_location_values = self.df.values
        num_location_values = len(self.df.values[0])
        df_gene_values = gene_info_df.values
        for i in tqdm(range(self.num_genes)):
            values = self.gene_to_location_dict.get(i)
            if values is not None:
                # For each of the genes assigned we want to make a new row
                for location_idx in values:
                    vals = list(df_gene_values[i])
                    # Get the values associated with that location.
                    loc_info = df_location_values[location_idx]
                    for loc_value in loc_info:
                        vals.append(loc_value)
                    new_rows.append(vals)
            elif keep_unassigned:
                vals = list(df_gene_values[i])
                for j in range(num_location_values):
                    vals.append(None)
                new_rows.append(vals)
        # Now we want to add these rows to the dataframe
        c_idx = 0
        new_df = pd.DataFrame()
        for c in new_columns:
            values = []
            for v in new_rows:
                values.append(v[c_idx])
            if c_idx < len(columns):
                new_df[f'{c}'] = values
            else:
                new_df[f'loc_{c}'] = values
            c_idx += 1

        # Check if we want to drop the unassigned rows
        if not keep_unassigned:
            new_df = new_df.dropna()

        self.gene_info_df = new_df
        return gene_info_df

    def save_gene_info_to_csv(self, filename: str, dropnull=False) -> None:
        """
        Save the gene information to the csv file.

        Parameters
        ----------
        filename:       str: filename to save it to
        dropnull:       Bool: if true we remove the null rows

        Returns
        -------
        None saves df to csv
        """
        self.gene_info_df = self.gene_info_df if self.gene_info_df is not None else self.get_gene_info_as_df()
        if dropnull:
            df = self.gene_info_df.dropna()
        else:
            df = self.gene_info_df
        self.u.save_df(df, filename)

    def save_loc_to_csv(self, filename: str, keep_unassigned=False) -> None:
        """
        Save the information of the location data
        Parameters
        ----------
        filename

        Returns
        -------

        """
        self.loc_df = self.loc_df if self.loc_df is not None else \
            self.assign_gene_info_to_loc_df(self.get_columns_in_gene_info(), keep_unassigned=keep_unassigned)
        self.u.save_df(self.loc_df, filename)

    def assign_gene_info_to_loc_df(self, gene_info_columns: list, col_prefix='gene_',
                                   keep_unassigned=False) -> pd.DataFrame:
        """
        Assign the gene information parameters to the dataframe (make a new dataframe).

        Parameters
        ----------
        gene_info_columns:            The columns from the gene_info that we want (id, chr, gc, start, end, direction, go_terms, ncbi)
        keep_unassigned:    Keep locations that weren't assigned to a gene
        col_prefix

        Returns
        -------
        pd.DataFrame:       Copy of the original dataframe with the new columns added
        """
        if self.loc_df is not None:
            return self.loc_df

        new_df = pd.DataFrame(self.rows_with_genes, columns=self.header)
        # Copy over the elements from the gene info df based on the index
        for g in range(0, len(gene_info_columns)):
            values = []
            for r in self.rows_with_genes:
                values.append(self.gene_annot_values[r[self.hdr_gene_idx]][g])
            new_df[gene_info_columns[g]] = values

        # Check if we want to drop the unassigned rows
        if not keep_unassigned:
            new_df = new_df.dropna()

        self.loc_df = new_df
        return new_df


    """
    -----------------------------------------------------------------
    Simple getters.
    -----------------------------------------------------------------
    """
    def get_gene_start(self, idx):
        return self.gene_annot_values[idx][self.gene_start]

    def get_gene_end(self, idx):
        return self.gene_annot_values[idx][self.gene_end]

    def get_gene_chr(self, idx):
        return self.gene_annot_values[idx][self.gene_chr]

    def get_gene_direction(self, idx):
        return self.gene_annot_values[idx][self.gene_direction]

    def get_gene_name(self, idx):
        return self.gene_annot_values[idx][self.gene_name]

    def get_gene_value_by_key(self, idx, key):
        return self.gene_annot_values[idx].get(key)
