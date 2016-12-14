Microarray data
---------------

Expression arrays 
~~~~~~~~~~~~~~~~~

Microarray normalisation 
^^^^^^^^^^^^^^^^^^^^^^^^

When investigating differential gene expression using microarrays, it’s
often the case that the expression levels of genes that should not
change given different conditions (e.g. housekeeping genes) report an
expression ratio other than 1. This can be caused by a variety of
reasons, for instance: variation caused by differential labelling
efficiency of the two fluorescent dyes used or different amounts of
starting mRNA. You can read more about this here_.

Normalisation is a process that eliminates such variations in order to
allow users to observe the actual biological differences in gene
expression levels. On Genestack, we have three different Microarray
Normalisation apps - one for each of the three commonly used chips:
Affymetrix, Agilent, and GenePix.

Affymetrix Microarray Normalisation.

Best used: for Affymetrix microarray data

Action: Normalisation of microarray data

Normalisation performed using using RMA.(Robust Multi-array Average).

The normalised microarray output can be assessed using the Microarrays
Quality Control application to detect and remove potential outliers.
Normalised microarrays that are of good quality can then be processed
for downstream processing such as Dose Response
Analysis.

RMA normalisation is based on the affy R package developed by Gautier L,
Cope L, Bolstad BM and Irizarry RA (2004), distributed under the GNU
Lesser General Public License (LGPL) version 2.0 or later license.

Agilent Microarray Normalisation

Best used: for Agilent microarray data

Action: Normalisation of microarray data

For 1-channel Agilent microarrays, various procedures for background
correction (e.g. "subtract", "half", "minimum", "normexp"), and
between-array normalisation (e.g. "quantile", "scale"), can be applied.

For 2-channel Agilent microarrays, procedures for within-array
normalisation (e.g. "loess", "median") can also be applied.

The normalised microarray output can be assessed using the Microarrays
Quality Control application to detect and remove potential outliers.
Normalised microarrays that are of good quality, can then be used for
downstream processing, such as Dose Response Analysis.

Normalisation procedures for Agilent are based on the limma R package
developed by Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth
GK (2015), distributed under GNU General Public License (GPL) version
2.0 or later.

Link:

-  https://www.bioconductor.org/packages/3.3/bioc/html/affy.html
-  https://www.bioconductor.org/packages/3.3/bioc/html/limma.html

References:

-  `Gautier L, Cope L, Bolstad BM and Irizarry RA (2004). “affy—analysis
   of Affymetrix GeneChip data at the probe
   level.” <https://www.google.com/url?q=http://www.ncbi.nlm.nih.gov/pubmed/14960456&sa=D&ust=1480960532113000&usg=AFQjCNFtiT-91LNqFrgRk1EAgdkumx1r9A>`__`Bioinformatics <https://www.google.com/url?q=http://www.ncbi.nlm.nih.gov/pubmed/14960456&sa=D&ust=1480960532114000&usg=AFQjCNGv5JUcbSDpTnNCyxm5J-sW76IfVQ>`__`, <https://www.google.com/url?q=http://www.ncbi.nlm.nih.gov/pubmed/14960456&sa=D&ust=1480960532114000&usg=AFQjCNGv5JUcbSDpTnNCyxm5J-sW76IfVQ>`__`20 <https://www.google.com/url?q=http://www.ncbi.nlm.nih.gov/pubmed/14960456&sa=D&ust=1480960532114000&usg=AFQjCNGv5JUcbSDpTnNCyxm5J-sW76IfVQ>`__`(3),
   pp.
   307–315. <https://www.google.com/url?q=http://www.ncbi.nlm.nih.gov/pubmed/14960456&sa=D&ust=1480960532114000&usg=AFQjCNGv5JUcbSDpTnNCyxm5J-sW76IfVQ>`__
-  `Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK
   (2015). “limma powers differential expression analyses for
   RNA-sequencing and microarray
   studies.” <https://www.google.com/url?q=http://europepmc.org/abstract/MED/25605792&sa=D&ust=1480960532115000&usg=AFQjCNFI070iNpVwXIIiQQNFN2Yq6-aqMA>`__`Nucleic
   Acids
   Research <https://www.google.com/url?q=http://europepmc.org/abstract/MED/25605792&sa=D&ust=1480960532115000&usg=AFQjCNFI070iNpVwXIIiQQNFN2Yq6-aqMA>`__`, <https://www.google.com/url?q=http://europepmc.org/abstract/MED/25605792&sa=D&ust=1480960532115000&usg=AFQjCNFI070iNpVwXIIiQQNFN2Yq6-aqMA>`__`43 <https://www.google.com/url?q=http://europepmc.org/abstract/MED/25605792&sa=D&ust=1480960532116000&usg=AFQjCNHHiEdnVuGvopb20Ndrx5PTDAUPkg>`__`(7),
   pp.
   e47. <https://www.google.com/url?q=http://europepmc.org/abstract/MED/25605792&sa=D&ust=1480960532116000&usg=AFQjCNHHiEdnVuGvopb20Ndrx5PTDAUPkg>`__

GenePix Microarray Normalisation

Best used: for GenePix microarray data

Action: Normalisation of GenePix microarray data

For GenePix microarrays, quantile between-array normalisation is
performed and various procedures for background correction (e.g.
"subtract", "half", "minimum", "normexp") can be applied.

The normalised microarrays output can be assessed using the Microarrays
Quality Control application to detect and remove potential outliers.
Normalised microarrays that are of good quality, can then be processed
for downstream processing such as Dose Response
Analysis.

Microarray QC 
~~~~~~~~~~~~~

Expression navigator for expression analysis 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compound dose response analysis 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Good for: Microarray data analysis

Input: Microarray data

Action: Describes the expression profiles of the significant genes as a
function of the dose.

Output:

Further apps to use:

This application performs dose response analysis on microarray data. It
identifies differentially expressed genes using the Bioconductor package
`limma <https://www.google.com/url?q=https://www.bioconductor.org/packages/release/bioc/html/limma.html&sa=D&ust=1480960532122000&usg=AFQjCNG3x6jMZtVXNPxYzOvN8LfePE4Upw>`__.
It then fits various regression models (linear, quadratic and power) to
describe the expression profiles of the significant genes as a function
of the dose.

The results are then reported in an interactive table. For each gene, an
optimal model is suggested based on the Akaike Information Criterion
(AIC), and the benchmark dose (BMD) is computed for that model. The
benchmark dose is estimated based on the method described in the
`Benchmark Dose
Software <https://www.google.com/url?q=http://www2.epa.gov/bmds/benchmark-dose-software-bmds-user-manual&sa=D&ust=1480960532124000&usg=AFQjCNHr41OQZN2zUawMuYPhF2n_To5Okg>`__` (BMDS)
user
manual: <https://www.google.com/url?q=http://www2.epa.gov/bmds/benchmark-dose-software-bmds-user-manual&sa=D&ust=1480960532124000&usg=AFQjCNHr41OQZN2zUawMuYPhF2n_To5Okg>`__

Let m(d) be the expected gene expression at dose d. The BMD then
satisfies the following equation: \|m(BMD)-m(0)\| = 1.349σ0 . In this
formula, σ0 is the standard deviation of the response at dose 0, which
we approximate by the sample standard deviation of the model residuals.

Link:
 `https://www.bioconductor.org/packages/release/bioc/html/limma.html <https://www.google.com/url?q=https://www.bioconductor.org/packages/release/bioc/html/limma.html&sa=D&ust=1480960532126000&usg=AFQjCNE45BYdIt9VkhhDfo8s1lWXV9K96A>`__ (Bioconductor
package limma).

Methylation arrays 
~~~~~~~~~~~~~~~~~~

Methylation array normalisation (coming soon) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Methylation array QC (coming soon) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Expression navigator for methylation arrays (coming soon) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. _here: http://www.mrc-lmb.cam.ac.uk/genomes/madanm/microarray/chapter-final.pdf