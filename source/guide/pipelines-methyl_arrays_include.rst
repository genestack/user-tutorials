Methylation arrays
~~~~~~~~~~~~~~~~~~

DNA methylation arrays are a widely-used tool to assess genome-wide DNA methylation.

Microarrays normalisation
+++++++++++++++++++++++++

**Action**: to perform normalisation of methylation microarray assays.

For methylation microarrays, normalisation can be performed with either "subsetQuantileWithinArray"
or "quantile" method, and in addition, "genomeStudio" background correction may be applied.

.. image:: images/meth-array-normalization-options.png
   :align: center

Further, the quality of normalised microarrays can be checked using the **Microarray QC Report**
application to detect and remove potential outliers. Normalised microarrays that are of good quality
may then be used in **Differential methylation analysis** and **Differential regions methylation analysis**.

The application is based on the minfi_ Bioconductor package.

.. _minfi: https://academic.oup.com/bioinformatics/article/30/10/1363/267584/Minfi-a-flexible-and-comprehensive-Bioconductor

Methylation array QC
++++++++++++++++++++

Quality assessment of microarray data is a crucial step in microarray analysis pipeline,
as it allows us to detect and exclude low-quality assays from the further analysis.

A single array with both red and green channels is used to estimate methylation for each
individual sample. Then, for each CpG locus, both methylated and unmethylated signal
intensities are measured.

Currently, there are two main methods used to estimate DNA methylation level:
*Beta-value* and *M-value*. The Beta-value is the ratio of the methylated probe intensity against the
overall intensity (sum of the methylated and the unmethylated probe intensities) plus a constant (Du P. et al.,
2010). Therefore, Beta-value can reflect roughly the percentage of methylated sites.

.. image:: images/qc-betaval.png

The M-value is the log2 ratio of the intensities of the methylated probe versus the unmethylated probe
(Du P. et al., 2010):

.. image:: images/qc-mval.png

**Action**: to assess quality of methylation microarray assays.

The Methylation array QC application allows the user to export files containing methylation and
unmethylation values, as well as the Beta-values, the M-values and Log median intensity values.

Additionally, you can download and explore "Copy number values" file with
the sum of the methylated and unmethylated signals.

Methylation array QC application provides various types of **quality control plots**.
Let's explore QC-report for the Infinium 450K microarrays:

**1) Log median intensity plot**

The scatterplot represents the log median of the signal intensities in both methylated and unmethylated channels
for each array. In general, samples of good quality cluster together,
while "bad" samples tend to separate, typically with lower median intensities.

.. image:: images/log-median-intensities.png
   :scale: 70 %
   :align: center

**2) Beta-values of the assays are represented by two plots:**

- *Beta density* plot represents the Beta-value densities of samples

.. image:: images/qc-beta-density.png
   :scale: 75 %
   :align: center

- *Beta density bean* plot also shows methylation the Beta-values.

.. image:: images/qc-beta-density-bean.png
   :scale: 75 %
   :align: center

**3) Control probes plots:**

The Infinium 450K arrays have several internal control probes helping to track
the quality on different stages of assay preparation (based on Illumina's `Infinium® HD Assay Methylation Protocol Guide`_):

.. _Infinium® HD Assay Methylation Protocol Guide: https://support.illumina.com/downloads/infinium_hd_methylation_assay_protocol_guide_(15019519_b).html

**Sample-independent controls**

Several sample-independent controls allow the monitoring different steps of
the of microarray assay preparation and include:

- *Staining control strip*, which estimate the efficiency of the staining step
  for both the red and green channels. They are independent of the hybridization
  and extension steps.

.. image:: images/qc-staining.png
   :scale: 75 %
   :align: center

- *Extension control strip*, which tests efficiency of single-base extension
  of the probes that incorporates labeled nucleotides. Both red (A and T,
  labeled with dinitrophenyl) and green (C and G labeled with biotin) channels
  are considered.

.. image:: images/qc-extension.png
   :scale: 75 %
   :align: center

- *Hybridization control strip*, which estimates entire performance of the
  microarray assay.

This kind of controls uses synthetic targets that are complementary to the array probe sequence.
Extension of the target provides signal.
The higher concentration of the target is used, the higher signal intensity will be registered.

.. image:: images/qc-hybridisation.png
   :scale: 75 %
   :align: center

- *Target removal control strip*, which tests whether all targets are removed
  after extension step. During extension reaction the sequence on the array is
  used as template to extend the control oligos. The probe sequences, however,
  are not extendable. The signal is expected to be low in comparison to the
  signal of hybridization control.

.. image:: images/qc-target-removal.png
   :scale: 75 %
   :align: center

**Sample-dependent controls**

A number of sample-dependent controls are provided to assess quality across samples.

- Bisulfite-conversion controls

To estimate methylation of DNA, the 450k assay probe preparation involves
bisulfite conversion of DNA when all unmethylated cytosines are converted
to uracils, while methylated cytosines are remains as they are.

*Bisulphite conversion I control strip*

This control uses Infinium I assay chemistry. There are two types of probes in this control:
bisulphite-converted and bisulphite-unconverted ones.
If the bisulphite conversion was successful, the converted
probes matches the converted DNA, and are extended. If the
sample has some unconverted DNA, the unconverted probes get extended.

.. image:: images/qc-bis-conversion-I.png
   :scale: 75 %
   :align: center

*Bisulphite conversion II control strip*

This control uses the Infinium I chemistry technology. If the bisulphite conversion
went well, the adenin base is added, generating signal in the red channel.
If there is some unconverted DNA, the guanine base is incorporated, resulting to
signal in the green channel.

.. image:: images/qc-bis-conversion-II.png
   :scale: 75 %
   :align: center

- Specificity controls, which monitor potential non-specific primer extension.

*Specificity I control strip* is used to assess allele-specific extention for the Infinium I chemistry assays.

.. image:: images/qc-specificity-I.png
   :scale: 75 %
   :align: center

*Specificity II control strip* allows to estimate specificity of extension for Infinium II assay
and test whether there is any nonspecific methylation signal detected over unmethylated background.

.. image:: images/qc-specificity-II.png
   :scale: 75 %
   :align: center

All the QC-plots shown on the application page may be downloaded in PDF format (see *Minfi PDF Report*).

Finally, based on the QC-results you can exclude particular samples as outliers,
remove them, and re-normalize the rest of the assays together. To do so, click *Sample list* and
select those samples that pass QC-check, then click **Remove outliers and re-normalise button**.

.. image:: images/QC-sample-list.png
   :scale: 75 %
   :align: center

Then, if you are happy with quality of re-normalized arrays, you can proceed to the following
step - **Differential Methylation Analysis**.

The "Methylation array QC" application is based on the minfi_ and the shinyMethyl_ Bioconductor packages.

.. _minfi: https://academic.oup.com/bioinformatics/article/30/10/1363/267584/Minfi-a-flexible-and-comprehensive-Bioconductor
.. _shinyMethyl: https://f1000research.com/articles/3-175/v2

Test differential methylation
+++++++++++++++++++++++++++++

.. Maybe rename the app as "Test differential methylation in CpG sites" or "Analysis of DMRs"?

**Action:** to identify differential methylation in single CpG sites ('a differentially
methylated positions (DMP)') across groups of normalized microarray assays using linear models.
Currently, 450k and EPIC Illumina's Methylation arrays are supported.

The input data for this application is Infinium Methylation Normalization file obtained with
the "Infinium Methylation Normalization” application.

.. Wrong file type! See the ticket https://trac.genestack.com/ticket/8099
.. As a result, the application generates Differential Expression Statistics file that you can further explore
.. with the Methylation Navigator for Sites.

The analysis includes annotating data when the application determines genomic position of the methylated
loci and its location relatively to various genomic features. Differential methylation analysis application
supports custom Methylation Array Annotation that you can upload with Import Data application.

The application computes differential
methylation statistics for each CpG site for the selected group compared to
the average of the other groups. Besides, you can assess differential methylation
for each group compared to a control one.

The application has the following options:

1. **"Group samples by"** option allows to group assays for comparison automatically:
the application helps you to group your samples according to experimental
factor indicated in metainfo for the microarray assays such as disease, tissue or treatment, etc.
(default: None)

2. **Control group** option allows to consider one of the created groups as a control one. In this  case
the application performs differential methylation analysis for each CpG site in the group against the control one.
(default: No control group)

.. image:: images/diff-meth-options.png
   :scale: 80 %
   :align: center

If you specify one or more confounding factors, you can identify differentially methylated sites
between tested groups of samples, while taking into account potential confounders, such as sex,
age, laboratory, etc. In this case the detected methylation changes are only caused by the factor
of interest, for example treatment, while any possible effects of confounding factors are excluded.
As confounding factors must be chosen according to metainfo keys common to all samples, remember
to specify the relevant information for all the samples.


Explore the output with interactive **Methylation Navigator**.

The application is based on the minfi_, limma_ Bioconductor packages.


Test differential regions methylation
+++++++++++++++++++++++++++++++++++++

**Action:** to determine and analyse contiguous regions which are differentially
methylated across groups of normalized microarray assays. Currently, 450k and EPIC Illumina's
Methylation arrays are supported.

As an input the application takes "Infinium Methylation Normalization" file with normalised microarray assays and returns
Differential Methylation Statistics file that you can further explore
with the Methylation Navigator.
Differential methylation analysis application supports custom methylation chip annotations
that you can upload with Import Data application.

The application has the following options:

1. **"Group samples by"** option allows to automatically group assays according to an experimental
factor indicated in metainfo for the selected microarray assays such as disease, tissue or treatment, etc.
(default: None)

2. **Control group** option allows to consider one of the created groups as a control one. In this  case
the application performs differential methylation analysis for each region in the group against the control one.
(default: No control group)

.. image:: images/diff-meth-options.png

If you specify one or more confounding factors, you can identify differentially methylated
regions between tested groups of samples, while taking into account potential confounders,
such as sex, age, laboratory, etc. In this case the detected methylation changes are only
caused by the factor of interest, for example treatment, while any possible effects of
confounding factors are excluded. As confounding factors must be chosen according to metainfo
keys common to all samples, remember to specify the relevant information for all the samples.


The Test Differential Regions Methylation application is based on the minfi_ and DMRcate_ packages.

Explore the output with interactive **Methylation Navigator**.

.. _minfi: https://academic.oup.com/bioinformatics/article/30/10/1363/267584/Minfi-a-flexible-and-comprehensive-Bioconductor
.. _limma: https://www.bioconductor.org/packages/3.3/bioc/html/limma.html
.. _DMRcate: https://bioconductor.org/packages/release/bioc/html/DMRcate.html


Methylation navigator for sites
+++++++++++++++++++++++++++++++

.. REDO pictures on tutorial's files (GSF21398704, i-dev).

**Action**: to view, sort and filter the results of analysis of differential methylation positions (DMPs).

.. image:: images/MN-sites.png

The Methylation Navigator page contains four sections:

1. The **Groups Information** section summarise the information on the created groups of samples to be tested.

2. The **Top Differentially Methylated Sites** table lists all the detected sites that are
differentially methylated in the selected group compared to either the average of the other groups
or a control group (if it is set).

.. image:: images/MN-top-sites.png
   :scale: 80 %
   :align: center

.. NEED TO FIX A REFERENCE NOTE ON THE APP PAGE [?]: we can also compare EACH individual
.. group to a set CONTROL one!

For each DMP (differentially methylated position) or DMR (differentially methylated region),
its Delta Beta, Average Beta, P-value, and FDR are shown.

Click probe ID to get more information about the probe:

.. More detailed description

.. image:: images/MN-sites-annotation.png

You can filter  by maximum acceptable false discovery rate (FDR),
up or down regulation, minimum log fold change (LogFC), and minimum log counts per million (LogCPM).

You can reduce the list of DMPs by filtering the data in the table based on the following criteria:

- *Max FDR* (maximum acceptable false discovery rate) — only shows sites with FDR below the set threshold.
- *Methylation All/ Down/ Up* — to show all sites or just those that are hypo- or hypermethylated.
- *Min Delta Beta* — delta Beta represents the difference between the Beta values in the groups being compared; this filter can be used to get only sites with absolute Delta Beta value of at least this threshold.
- *Min Average Beta* — only shows sites with average Beta value of at least this threshold.

.. image:: images/MN-sites-filter.png
   :align: center

Sort the list of probes by clicking the arrows next to the name of the statistical metrics in the table headers.

.. image:: images/MN-sites-sort.png
   :scale: 80 %
   :align: center

3. **A boxplot of methylation levels**

Each color corresponds to an individual probe you selected; each circle represents an assay belonging to
the tested group. Each boxplot represents the distribution of a methylation in a given group.
The y-axis shows Beta values, while the x-axis shows probe IDs.

.. image:: images/MN-sites-boxplot.png
   :scale: 80 %
   :align: center

4. The bottom-right section contains **a search box** that allows you to explore the results for a particular
probe. Start typing a probe ID and select the probe of interest in the appeared drop-down
list of possible variants.

.. image:: images/MN-sites-search.png
   :scale: 80 %
   :align: center

You can further export either the complete table of differential methylation analysis for all the groups
or the list of values for the specific comparison in TSV format. See **Export Data (for all comparisons, as .tsv)**
and **Download filtered data for current comparison as .tsv** options, respectively.

.. image:: images/MN-sites-export.png
   :align: center


Methylation navigator for regions
+++++++++++++++++++++++++++++++++

.. REDO pictures on tutorial's files (GSF21398704, i-dev).

**Action**: to view, sort and filter the results of analysis of differential methylation regions (DMRs).

.. image:: images/MN-regions.png

The Methylation Navigator page contains the following sections:

1. The **Groups Information** section summarise the information on the created groups of samples to be tested.

.. image:: images/MN-regions-group-info.png

2. The **Top Differentially Methylated Regions** table shows all the detected regions that are
differentially methylated in the selected group compared to either the average of the other
groups or a control group (if it is set).

.. image:: images/MN-top-regions.png

You can further reduce the list of identified DMRs and exclude those regions that do not meet set
filtering criteria. The following filters can be applied:

- *Max FDR* (maximum acceptable Stouffer-transformed false discovery rate) — the FDR is statistical certainty that the given region is differentially methylated. This filter only shows regions with the FDR values below the set threshold. Learn more about Stouffer-test from the paper by `Kim S.C. (2013).`_
- *Methylation* (Down/All/Up) — shows all regions or only hypo- or hypermethylated ones.
- *Min BetaFC* (minimum mean beta fold change within the region) — for every DNA region, each probe has its Beta value, which is defined as relative methylation of the region (B1, B2 etc.). BetaFC, in this case, can be defined as mean Beta fold change; apply the filter to show only regions having BetaFC below the threshold.
- *Min significant CPG sites count* — minimum number of CpG sites inside the genomic region.

.. Suggestion: rename the filter 'Min BetaFC' to 'Min mean BetaFC'
.. Suggestion: rename the filter 'Min significant CPG sites count” to “Min CPG sites count”'

.. _Kim S.C. (2013).: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3653960/

.. add description for the option

.. image:: images/MN-regions-filters.png

You can also sort the list of identified DMRs by clicking the arrows next to the name of
the statistical metrics in the table.

.. image:: images/MN-regions-sort.png

Finally, you can export both the complete table of top differential methylated regions
for all the groups (**Export Data (for all comparisons, as .tsv)**) and the list of
regions with associated statistics for the one comparison in TSV format
(**Download filtered data for current comparison as .tsv**).

.. image:: images/MN-sites-export.png