# sciepi2gene
[![codecov.io](https://codecov.io/github/ArianeMora/sciepi2gene/coverage.svg?branch=master)](https://codecov.io/github/ArianeMora/sciepi2gene?branch=master)
[![PyPI](https://img.shields.io/pypi/v/sciepi2gene)](https://pypi.org/project/sciepi2gene/)



Sci-epi2gene maps events annotated to a genome location to nearby genes - i.e. peaks from histone modification data
ChIP-seq experiemnts stored as bed data, or DNA methylation data in csv format (e.g. output from DMRseq or methylKit).

The user provides a SORTED gene annotation file with start, end, and direction for each gene (we recommend using
[sci-biomart](https://github.com/ArianeMora/scibiomart), see examples for detail.  

The user then selects how to annotate, i.e. whether it is in the promoter region, or overlaps the gene body. Finally,
the parameters for overlap on each side are chosen.

It is available under the [GNU General Public License (Version 3) ](https://www.gnu.org/licenses/gpl-3.0.en.html).

This package is a wrapper that allows various epigenetic data types to be annotated to genes.

I also wanted to have different upper flanking and lower flanking distances that took into account the directionality of the strand
and also an easy output csv file that can be filtered and used in downstream analyses. This is why I keep all features
that fall within the annotation region of a gene (example below):

The overlapping methods are as follows:
    1) overlaps: this means does ANY part of the peak/feature overlap the gene body + some buffer before the TSS and some buffer on the non-TSS side
    2) promoter: does ANY part of the peak/feature overlap with the TSS of the gene taking into account buffers on either side of the TSS.

.. image:: _static/example_overlaps.png
   :width: 600

As you can see from the above screenshot using IGV, the input peaks are in purple, and the green are the output
peaks as annotated to genes. The function *convert_to_bed* converts the output csv to bed files for viewing. This example
shows that a peak/feature can be annotated to multiple genes. Peaks/features outside of the regions of genes (e.g.
the first peak) are dropped from the output.

We show this example in the notebook (see examples folder), where we use [IGV](https://github.com/igvteam/igv-jupyter#igvjs-jupyter-extension)
to view the tracks (see image below).

.. image:: _static/igv_jupyter.png
   :width: 600
   
Lastly, there are sometimes differences between annotations (i.e. the TSS on your annotation in IGV may differ to the
annotation you input to sciepi2gene), naturally, how your genes/features are annotated depends on the input file so if you see differences check this first!

Please post questions and issues related to sci-epi2gene on the `Issues <https://github.com/ArianeMora/sciepi2gene/issues>`_  section of the GitHub repository.

