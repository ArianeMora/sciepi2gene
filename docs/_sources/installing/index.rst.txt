.. _installing:

Installing sci-epi2gene
=======================

Note: It is strongly recommended to create a `virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_
for each version of sci-loc2gene installed.

.. code-block:: python

    pip install scie2g


sci-epi2gene reference data
---------------------------

For ease of use we provide reference files for mouse and human, however, we recommend following the examples
section which shows how to build these.

.. code-block:: bash

    wget https://github.com/ArianeMora/sciepi2gene/data/hsapiens_gene_ensembl-GRCh38.p13.csv
    wget https://github.com/ArianeMora/sciepi2gene/data/mmusculus_gene_ensembl-GRCm38.p6.csv

Examples using scibiomart to build references:

*Human*

.. code-block:: bash

    scibiomart --m ENSEMBL_MART_ENSEMBL --d hsapiens_gene_ensembl --a "ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position,strand" --o "humanSorted" --s t


*Mouse*

.. code-block:: bash

    scibiomart --m ENSEMBL_MART_ENSEMBL --d mmusculus_gene_ensembl --a "ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position,strand" --o "mm10Sorted" --s t


Below are supported versions of sci-epi2gene:

.. list-table::
   :widths: 10 10
   :header-rows: 1

   * - sci-epi2gene Release
     - Version
   * - 1
     - 1.0.0

