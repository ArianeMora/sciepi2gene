# Examples

## Current data types implemented

#### 1) CSV:  
This is just a generic csv data input. You need to input the header information (i.e. column name for chr, start, end, 
value). Does not require to be sorted (we sort it).  

#### 2) Bed:  
Bed files are able to be mapped using this tool as well. There are examples in the tests. Here we 
keep both the width and the signal for a peak and assign that to a gene if it overlaps the gene body by default. We 
enforce a pre-specified ordering (e.g. following the protocol used in Encode). i.e. [CHROM, START, END, ID, SCORE, SIGNAL]

**Note bed files are assumed to be sorted (used bedtools to sort if they aren't already sorted).**

### Workflow:

1) Ensure your input file is sorted by the event start (i.e. peak start for bed files)  
2) Initialise your class (i.e. Bed, DMRseq...)   
        2.a) Here you choose the overlap methods i.e.:  
                 *overlap_method='in_promoter' or 'overlaps'* --> overlaps entire gene e.g. for H3K36me3  
                 *buffer_after_tss=500* --> buffer after gene stop when using overlaps  
                 *buffer_before_tss=2500* --> buffer before TSS in overlaps and in_promoter  
                 *buffer_gene_overlap=500* --> overlap with the gene when using in_promoter  
3) Add an annotation file - this has the gene identifiers and the locations of the gene (start, end, direction)  
        3.a) see examples here: https://github.com/ArianeMora/sciepi2ene/data/  
        3.b) function called: set_annotation_using_biomart or set_annotation_from_file  
4) Assign locations to genes  
        4.a) assign_locations_to_genes()  
5) Save output file  
        5.a) save_loc_to_csv(filename)
        5.b) convert_to_bed(dataframe, output_filename, label_for_track, value_column, name/id_column)

## Bed

#### Example Bed input
```
chr1	3669319	3669933	Peak_39246	23	.	2.95584	8.15439	6.38865	445
chr1	3670510	3670815	Peak_111700	14	.	2.17341	4.00125	2.48284	45
chr1	3672158	3672625	Peak_23195	30	.	3.47746	11.49358	9.59823	247

```

#### Bed usage

```
from scie2g.bed import Bed

# Create an annotation file
bed = Bed(input_bed_file.bed, overlap_method='in_promoter',
                   output_bed_file=f'selected_peaks/filename.bed')
# This is orgamism specific see scibiomart for more details
bed.set_annotation_using_biomart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')
# Save the annotation so next time you don't have to do this
bed.save_annotation(path_to_file + 'hg38')
# bed.set_annotation_from_file(filename) --> if you have already run the previous line 

# Now we can annotate the bed file to the genes
bed.assign_locations_to_genes()
bed.save_loc_to_csv(f'path_to_save_to/filename.csv')

```
The above code snippet will annotate peaks to genes if they fall in the promoter region. It will output a csv file
and also the selected peaks as a bed file.

If you only want the largest peak to be assigned to a gene, you can use pandas to filter the CSV to groupby the gene 
name and take the max value.

#### Example Bed output
```
peak_idx,chr,start,end,gene_idx,padj,signal,width,id
3,chr1,3649114,3649487,10,3.02989,3.1721,373,Gm19938
3,chr1,3649114,3649487,11,3.02989,3.1721,373,Xkr4
4,chr1,3668534,3668732,11,2.46958,3.08673,198,Xkr4
```

## 3) CSV
An example of a CSV file is the output from DMRseq. Here we want to keep track of the region ID. 
We also need to keep the stat and the pvalue. One peak may map to multiple genes 
however we leave this in there you can easily sumarise this using R or Python.

#### Example csv input
```
"","seqnames","start","end","width","strand","L","area","beta","stat","pval","qval","index.start","index.end","index.width"
"1","16",8855737,8856501,765,"*",39,33,2.9,40,0.005,0.01,3303,3303,39
"2","10",11430660,11430937,278,"*",22,14,1.5,102,0.001,0.02,2284,2284,22
"3","19",1872607,1872632,266,"*",15,14,2.3,59,0.005,0.01,3766,3766,15
```

#### CSV usage
```
from scie2g.csv import CSV
# seqnames = column name for chr
# start = column name for start
# end = column name for end
# stat = column name for value
# [] = the list contains extra features we want in the output
f = Csv('csv_filename',  'seqnames', 'start', 'end', 'stat', ['index.start', 'index.end', 'qval'])
# Add the gene annot
f.set_annotation_from_file('annotation_file_name')
# Now we can run the assign values
f.assign_locations_to_genes()
f.save_loc_to_csv(output_filename)
# Save as bed for viewing in IGV --> loc_df contains the dataframe with the locations
f.convert_to_bed(f.loc_df, 'output_bedfilename', 'tracklabel', 'qval', 'external_gene_name') # Uses the gene name to annotate the peak
```

#### Example CSV output
```
region_idx,chr,start,end,gene_idx,padj,stat,num_cpgs,ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position,strand
1,chr16,8855737,8856501,20576,0.01,39.0,40.0,ENSG00000260276,AC022167.2,chr16,8848105,8860456,1
1,chr16,8855737,8856501,20579,0.01,39.0,40.0,ENSG00000153048,CARHSP1,chr16,8852942,8869012,-1
3,chr19,1872607,1872632,27126,0.01,15.0,59.0,ENSG00000267232,AC012615.5,chr19,1875016,1875992,1
```

## Multiprocessing

Here is an example using the python multiprocessing functions, we imagine that there is a list of files that you want
to process
```
from scie2g.bed import Bed

from multiprocessing.dummy import Pool as ThreadPool

def run_bed(filename):
     try:
         data_dir = 'directory containing bed files'
         output_dir = 'directory for output'
         bed = Bed(f'{data_dir}{filename}', overlap_method='in_promoter',
                   output_bed_file=f'selected_peaks/{filename}')
         # Use your annotation file i.e. for humans or mouse
         bed.set_annotation_from_file('supps/entrez-mmusculus_gene_ensembl-GRCm38.p6')
         # Now we can run the assign values
         bed.assign_locations_to_genes()
         bed.save_loc_to_csv(f'{output_dir}{filename[:-3]}csv',  keep_unassigned=True)
     except:
         print(filename)

files_to_process = [file1, file2, file...]

# Run in paralell
pool = ThreadPool(10)
results = pool.map(run_bed, files_to_process)
```


