CLI
===

CLI instructions for running spi-epi2gene.


Examples
--------

The examples cover some simple case where we map 1) a bed file with genes that overlap (e.g. H3K36me3) and
2) a bed file with peaks we would expect to land in the promoter (H3K27ac). Lastly, we show how a DMRseq file would be
annotated to genes.


Genes are assigned to a peak if the peak overlaps with 2500 (default) upstream of the TSS or 500 (default) base pairs on the gene body. The peak
could overlap with any part of this region and it will be assigned to the gene. Peaks can be assigned to multiple genes i.e.
if we have a broad peak then it could be assigned to many genes (i.e. H3K27me3).


Annotate gene body overlapping peaks
------------------------------------
Here we read in the file `test_H3K27me3.bed` and annotate it to peaks that fall in the promoter region of genes
annotated in `hsapiens_gene_ensembl-GRCh38.p13.csv` (generated using `sci-biomart <https://arianemora.github.io/scibiomart/>`_).
Genes are assigned to a peak if the peak overlaps with 2500 upstream of the TSS or 500 base pairs from the gene end.

.. code-block:: bash

    scie2g --a data/hsapiens_gene_ensembl-GRCh38.p13.csv --o data/output_file.csv --l2g data/test_H3K36me3.bed --t b --upflank 3000 --downflank 500 --m overlaps


Annotate promoter region peaks
------------------------------
Here we read in the file `test_H3K27me3.bed` and annotate it to peaks that fall in the promoter region of genes
annotated in `hsapiens_gene_ensembl-GRCh38.p13.csv` (generated using `sci-biomart <https://arianemora.github.io/scibiomart/>`_).
Genes are assigned to a peak if the peak overlaps with 2500 upstream of the TSS or 500 base pairs on the gene body. The peak
could overlap with any part of this region and it will be assigned to the gene.

.. code-block:: bash

    scie2g --a data/hsapiens_gene_ensembl-GRCh38.p13.csv --o data/output_file.csv --l2g data/test_H3K27me3.bed --t b --upflank 2500 --overlap 500 --m in_promoter

Annotate DMRseq regions (CSV) to genes
--------------------------------------

Here we have had to override the column: **'chr'** with 'seqnames' seen with the tag --chr and the **'value'** term, with 'stat'
.. code-block:: bash

    scie2g --a data/hsapiens_gene_ensembl-GRCh38.p13.csv --o data/output_file.csv --l2g data/test_dmrseq.csv --t d --upflank 2500 --m overlaps --chr seqnames --value stat

Annotate MethylKit DMCs (CSV) to genes
--------------------------------------

    scie2g --a data/hsapiens_gene_ensembl-GRCh38.p13.csv --o data/output_file.csv --l2g data/test_H3K27me3.bed --t b --upflank 2500 --overlap 500 --m in_promoter

.. code-block:: bash

    scie2g --a data/hsapiens_gene_ensembl-GRCh38.p13.csv --o data/output_file.csv --l2g data/test_methyl.csv --t d --upflank 2500 --m overlaps --value meth.diff


Arguments
---------

.. argparse::
   :module: scie2g
   :func: gen_parser
   :prog: scie2g
   :nodefaultconst:
