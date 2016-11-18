Testing Differential Isoform Expression on Genestack Platform
*************************************************************

Detecting differential isoform (i.e., transcript) expression across different biological conditions
from RNA-Seq data is one of the common transcriptome analysis tasks. And the purpose of this
tutorial is to show you how to perform this analysis on `Genestack`_  platform.

|isoforms|

`Alternative splicing`_ of messenger RNA, when different combinations of
exons can be spliced together, produces a wide variety of RNA transcript
isoforms encoding different protein products. Changes in protein sequence may
influence its function, that is why alternative splicing is not only
a source of huge protein diversity but it can also significantly affect various cellular processes.
Undoubtedly, alternative splicing is a major mechanism for gene regulation that is
of great importance to understand.


Setting up an RNA-Seq experiment
--------------------------------
Our analysis will be based on RNA-Seq data coming from `Trapnell et al. 2012`_,
one of the public experiments we have on the platform.
However, you can upload your own RNA-Seq data using `Data Import`_ application
or search through all `public experiments`_ we have on the platform and
choose a suitable one.

Here is some information about this `experiment`_ opened in Metainfo Viewer:

|Metainfo_viewer|

Shortly, the authors conducted differential isoform expression analysis.
They performed RNA interference (RNAi)-mediated knockdown of HOXA1 in
human lung fibroblasts, where HOXA1 was depleted using a pool of
HOXA1-targeting siRNAs. And they compared these HOXA1-depleted
fibroblasts against cells treated with a pool of scrambled siRNAs that
do not target a specific gene. Then, the authors isolated total RNA in
biological triplicate and sequenced RNA on Illumina HiSeq 2000 and
Illumina MiSeq platforms.

As a result, they identified differentially expressed transcripts and genes specific
to cell cycle progression and related to apoptosis induction. And
in this tutorial we will try to reproduce their results.


Building an Isoform-level Differential Expression Analysis pipeline
-------------------------------------------------------------------

Below is a simple dataflow to analyze one of our RNA-Seq samples.
Later we will show you how it is easy to build the same pipeline for the other 11 samples.

|Dataflow_up_1|

The dataflow consists of several steps:

#. Quality control and preprocessing of RNA-Seq reads
#. Mapping raw reads onto reference genome
#. Quality control of mapped reads
#. Calculate FPKM coverage for each isoform
#. Differential isoform expression analysis


Quality control and preprocessing of raw reads
----------------------------------------------
Before mapping raw reads and calculation of isoform abundance, you may
be interested in improving the reads quality. We offer you
various preprocess applications to improve quality  of
your raw sequencing data.

Start with one sample and try to run, for example, Trim Adaptors and Contaminants app.
You will see that each app suggests you to add next analytical
step or to use relevant viewers. 

Generally speaking, quality control checks should be performed on different steps of
the analysis to monitor quality of your data and, therefore, ensure reliability and
accuracy of the results.

Quality control of raw reads read allows to determine sequence quality, GC content
distribution, the presence of adaptors, overrepresented sequences and duplicated reads.

Here we will check quality of raw sequence data with FastQC Report application
based on widely-used FastQC tool.

Look, what quality statistics you can view using FastQC Report app:

**Qualities per base** plot shows the range of quality scores for each
position on the reads.

|Qualities per base|

**GC content distribution** plot shows GC%
(x-axis) and read frequencies (y-axis). In a random library you can
expect a roughly normal distribution of GC content, as in our case. An
unusually shaped distribution could indicate a contaminated library or
some other kinds of biased subset.

|GC content distribution|

You can find more statistics in output Raw Reads QC Reports. We run QC
on all the data in the experiment and collected reports in folder `Raw
Reads QC reports for Trapnell et al. (2012)`_.




Mapping RNA-Seq reads onto reference genome
-------------------------------------------

On the next step, we will use Spliced Mapping app to map RNA-Seq reads
onto the reference genome and discover transcript splice sites. By
default, the app will identify both known and novel alternative splicing
variants, will align reads with no more than 2 mismatches and report
both unique and multiple mappings. Change options and default values,
clicking on “Edit parameters” button. 

You can find all Mapped Reads files in folder `Mapped Reads files for Trapnell et al. (2012)`_.
If you open them in `Genome Browser`_,
you can find out that HOXA1 gene is really non-transcribed for
HOXA1 knockdown data:

|GB_HOXA1|

Quality control of mapped reads
-------------------------------

This is an optional step. There are some apps developed for simple
quality control of your mapped reads. In this tutorial, let’s create QC
report for each Mapped Reads file using Mapped Reads QC Report app and
then open all reports in `Multiple QC plotter app`_:

|Mapped_QC_plotter|

All Mapped Reads QC reports are publicly available and stored in
folder  `Mapped Reads QC reports for Trapnell et al.
(2012)`_.

Calculate FPKM coverage for each isoform
----------------------------------------

We will run Quantify FPKM Coverage in Isoforms app to calculate isoform
abundance. The app takes Mapped Reads file and calculates FPKM
(Fragments Per Kilobase of exon per Million fragments mapped) values for
each transcript.

So, now we have pipeline for one sample. Is it possible to run dataflow
for multiple samples? Yes, it is possible. Open FPKM isoform counts file
in `Data Flow Editor`_
app, for sequencing assay make “Add files” action, choosing the left 11
raw sequencing files and click on blue “Create files” button. Look, we
built 12 pipelines! The last thing is just to “start initialization” for
all of them.

|DF|

We calculate FPKM coverage in all samples and collected result files in
folder  `FPKM isoforms counts for Trapnell et al. (2012)`_.

Differential isoform expression analysis
----------------------------------------

The final step is to perform differential isoform expression
analysis between two groups of samples corresponding to different
conditions. In our case, it is anti-HOXA1 siRNA and scrambled control
fibroblasts.

In `File Browser`_, we choose these 6 Data files with FPKM isoforms counts (let’s consider
only MiSeq data) and click on Test Differential Isoform Expression in
Analyse section. To run the app we need to assign samples to groups. We
can do it manually or apply auto-grouping. Just click, for example on
“GEO transfection” header in the table and the app suggests you to
create two groups according to “HOXA1 knockdown” and “Scramble siRNA”
transfection conditions:

|Diff_iso|

So, we agree and do “Group samples automatically”. Below, you see some
correction parameters you can apply for analysis. We will use default
values. And finally let’s create our file and run the analysis clicking
on “start initialization” in “Other Actions”. We created two
Differential Expression Statistics files (for data from two sequencing
platforms – MiSeq and HiSeq) and put them in folder  `Differential
Isoform Expression Analysis for Trapnell et al. (2012)`_.

When the analysis will be complete, look at the Top Differentially
Expressed Isoforms table. On HiSeq data, more than 800 differentially
expressed isoforms (460 up-regulated and 410 down-regulated) were
identified:

|HiSeq_DIEA|

For selected transcripts, you can see Count Graph with normalised FPKM
counts across samples. This allows you to observe how a gene’s
expression level varies within groups. Look, for example, at first two
down-regulated transcripts for HOXA1 knockdown group:

|graph|

Our results are consistent with paper results. We also found that the
loss of *HOXA1* results in significant expression level changes for
different transcripts encoded by genes which play important role in cell
development.

You can find all tutorial files in folder `Testing Differential Isoform Expression on Genestack Platform`_ and
look at all results we got for each analytical step.

This is the end of this tutorial. We hope you found it useful and that you are now ready to
make the most out of our platform.
If you have any questions and comments, feel free to email us at feedback@genestack.com or
visit our forum_. Also we invite you to follow us on Twitter `@genestack <https://twitter.com/genestack>`__.

.. |isoforms| image:: images/isoforms.png
.. |Metainfo_viewer| image:: images/Metainfo_viewer.png
.. |Dataflow_up_1| image:: images/Dataflow_up_1.png
.. |Qualities per base| image:: images/Qualities-per-base.png
.. |GC content distribution| image:: images/GC-content-distribution.png
.. |GB_HOXA1| image:: images/GB_HOXA1.png
.. |Mapped_QC_plotter| image:: images/Mapped_QC_plotter.png
.. |DF| image:: images/DF.png
.. |Diff_iso| image:: images/Diff_iso.png
.. |HiSeq_DIEA| image:: images/HiSeq_DIEA.png
.. |graph| image:: images/graph.png
.. _Genestack: https://platform.genestack.org/
.. _Alternative splicing: http://en.wikipedia.org/wiki/Alternative_splicing
.. _Data Import: https://platform.genestack.org/endpoint/application/run/genestack/uploader
.. _public experiments: https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF070886&action=viewFile&page=1
.. _Trapnell et al. 2012: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37703
.. _experiment: https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF080230&action=viewFile
.. _Raw Reads QC reports for Trapnell et al. (2012): https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF1018515&action=viewFile&page=1
.. _Mapped Reads files for Trapnell et al. (2012): https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF1018519&action=viewFile&page=1
.. _Genome Browser: https://platform.genestack.org/endpoint/application/run/genestack/genomeBrowser?a=GSF1018248&action=viewFile
.. _Multiple QC plotter app: https://platform.genestack.org/endpoint/application/run/genestack/multiple-qc-plotter?a=GSF1018535&action=viewFile
.. _Mapped Reads QC reports for Trapnell et al. (2012): https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF1018518&action=viewFile&page=1
.. _Data Flow Editor: https://platform.genestack.org/endpoint/application/run/genestack/datafloweditor?a=GSF3725699&action=viewFile
.. _FPKM isoforms counts for Trapnell et al. (2012): https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF1018517&action=viewFile&page=1
.. _File Browser: https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF1018517&action=viewFile&page=1
.. _Differential Isoform Expression Analysis for Trapnell et al. (2012): https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF1018516&action=viewFile&page=1
.. _Testing Differential Isoform Expression on Genestack Platform: https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF123346&action=viewFile
.. _forum: http://forum.genestack.org/