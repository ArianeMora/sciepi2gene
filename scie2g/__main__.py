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

import argparse
import os
import sys

from sciutil import SciUtil

from scie2g import __version__
from scie2g import Bed, Csv


def print_help():
    lines = ['-h Print help information.']
    print('\n'.join(lines))


def run(args):
    if args.t == 'd':
        if not args.value:
            u = SciUtil()
            u.warn_p(['WARNING: You did not pass a column name for the value! Please use --value for the column '
                      'you would like to use as your value in your output file.\n Returning ...'])
            return
        c = Csv(args.l2g, chr_str=args.chr, start=args.start, end=args.end, value=args.value, header_extra=args.hdr,
                overlap_method=args.m, buffer_after_tss=args.downflank,
                buffer_before_tss=args.upflank, buffer_gene_overlap=args.overlap,
                gene_start=args.gstart, gene_end=args.gend, gene_chr=args.gchr,
                gene_direction=args.gdir, gene_name=args.gname
                )
        c.set_annotation_from_file(args.a)
        c.assign_locations_to_genes()  # Now we can run the assign values
        c.save_loc_to_csv(args.o)
        if args.b:
            c.convert_to_bed(c.loc_df, args.b, args.b)
    elif args.t == 'b':
        bed = Bed(args.l2g, overlap_method=args.m, buffer_after_tss=args.downflank,
                  buffer_before_tss=args.upflank, buffer_gene_overlap=args.overlap,
                  gene_start=args.gstart, gene_end=args.gend, gene_chr=args.gchr,
                  gene_direction=args.gdir, gene_name=args.gname, chr_idx=args.chridx, start_idx=args.startidx,
                  end_idx=args.endidx, peak_value=args.valueidx, header_extra=args.hdridx
                  )
        # Add the gene annot
        bed.set_annotation_from_file(args.a)
        # Now we can run the assign values
        bed.assign_locations_to_genes()
        bed.save_loc_to_csv(args.o)


def gen_parser():
    parser = argparse.ArgumentParser(description='scie2g')
    parser.add_argument('--a', type=str, help='Annotation with the gene locations')
    parser.add_argument('--o', type=str, default='l2g_outputfile.csv', help='Output file (csv)')
    parser.add_argument('--b', type=str, default='l2g_outputfile.bed', help='Output file (bed)')
    parser.add_argument('--l2g', type=str, help='Input file to run scie2g on')
    parser.add_argument('--t', type=str, default='b', help='The input file type: d=CSV, b=Bed')
    parser.add_argument('--upflank', type=int, default=2500, help='Maximum distance upstream from TSS'
                                                                  ' (default = 2500) for overlaps and in_promoter')
    parser.add_argument('--downflank', type=int, default=500, help='Maximum distance downstream from gene end '
                                                                   '(default = 500) only used in overlaps')
    parser.add_argument('--overlap', type=int, default=500, help='Overlap with gene body (default = 500) used in'
                                                                 ' in_promoter')
    parser.add_argument('--m', type=str, default='in_promoter', help='Overlap method'
                                                                     ' (overlaps or in_promoter <- default).')

    parser.add_argument('--chr', type=str, default="chr", help='CSV only: name of your chromosone column')
    parser.add_argument('--start', type=str, default="start", help='CSV only: name of your start column')
    parser.add_argument('--end', type=str, default="end", help='CSV only: name of your end column')
    parser.add_argument('--value', type=str, default=None, help='CSV only: name of your value column')
    parser.add_argument('--hdr', type=str, default="", help='CSV only: comma separated list of other columns you '
                                                            'want to include in the output e.g "stat,pvalue"')

    parser.add_argument('--chridx', type=int, default=0, help='BED only: index of your chromosone column')
    parser.add_argument('--startidx', type=int, default=1, help='BED only: index of your start column')
    parser.add_argument('--endidx', type=int, default=2, help='BED only: index of your end column')
    parser.add_argument('--valueidx', type=int, default=7, help='BED only: index of your value column')
    parser.add_argument('--hdridx', type=str, default="0,1,2,3,6,8", help='BED only: comma separated list of indexs')
    parser.add_argument('--hdrlbl', type=str, default='"chr","start","end","peak_name","signal","qvalue"',
                        help='BED only: comma separated list of header in human readable format as output '
                             'to your csv file.')

    parser.add_argument('--gchr', type=int, default=2, help='Position in annotation file that your chr annotation is.')
    parser.add_argument('--gstart', type=int, default=3, help='Position in annotation file that your start is.')
    parser.add_argument('--gend', type=int, default=4, help='Position in annotation file that your end is.')
    parser.add_argument('--gdir', type=int, default=5, help='Position in annotation file that your gene direction is.')
    parser.add_argument('--gname', type=int, default=0, help='Position in annotation file that gene name is.')

    return parser


def main(args=None):
    parser = gen_parser()
    u = SciUtil()
    if args:
        sys.argv = args
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'scie2g v{__version__}')
        sys.exit(0)
    else:
        print(f'scie2g v{__version__}')
        args = parser.parse_args(args)
        # Validate the input arguments.
        if not os.path.isfile(args.a):
            u.err_p([f'The annotation file could not be located, file passed: {args.a}'])
            sys.exit(1)
        if not os.path.isfile(args.l2g):
            u.err_p([f'The input file could not be located, file passed: {args.l2g}'])
            sys.exit(1)
        if args.t != 'b' and args.t != 'd':
            u.err_p([f'The file type passed is not supported: {args.t}, '
                     f'filetype must be "b" for bed or "d" for dmrseq.'])
            sys.exit(1)
        # Otherwise we have need successful so we can run the program
        u.dp(['Running scie2g on input file: ', args.l2g,
              '\nWith annotation file: ', args.a,
              '\nSaving to output file: ', args.o,
              '\nOverlap method:', args.m,
              '\nUpstream flank: ', args.upflank,
              '\nDownstream flank:', args.downflank,
              '\nGene overlap: ', args.overlap])
        u.warn_p(['Assuming your annotation file and your input file are SORTED!'])
        # RUN!
        run(args)
    # Done - no errors.
    sys.exit(0)


if __name__ == "__main__":
    main()
    # ----------- Example below -----------------------
    """
    root_dir = '../scie2g/'
    main(["--a", f'{root_dir}data/hsapiens_gene_ensembl-GRCh38.p13.csv',
          "--o", f'{root_dir}output_file.csv',
          "--l2g", f'{root_dir}tests/data/test_H3K27ac.bed',
          "--t", "b",
          "--upflank", "3000"])
          
    root_dir = '../'
    main(["--a", f'{root_dir}data/hsapiens_gene_ensembl-GRCh38.p13.csv', # mmusculus_gene_ensembl-GRCm38.p6.csv', #
      "--o", f'{root_dir}output_file_2.csv',
      "--l2g", f'{root_dir}tests/data/test_dmrseq.csv',
      "--t", "d",
      "--upflank", "20", "--chr", "seqnames", "--value", "stat"])
    """
