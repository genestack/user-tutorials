WGS QC raw reads
****************

Poorly identified bases, low-quality sequences and contaminants (such as
adaptors) in the raw sequencing data can affect downstream analysis
leading to erroneous results and conclusions. Before starting the WGS
analysis, we will check the initial data quality and decide how to
improve the downstream analysis by a variety of preprocessing options.  
FastQC Report app is based on FastQC tool and produces
several statistics characterising the raw data quality: Phred score
distribution per base and per sequence, GC content distribution, per
base sequence content, read length distribution, and sequence
duplication level.   We will compute quality control statistics with
FastQC Report app for both assays from the experiment by using the "Raw
Reads Quality Control" public data flow. In order to do so, select the
assays and go to "Run data flow on selection..." using the context menu.
Then, from the list of suggested variants choose "Raw Reads Quality
Control" to open it in Data Flow Runner app. |DF list| The app page
presents the quality control part of the pipeline in a graphical form.
To generate QC-reports click on the "Run Data Flow" button and then on
"Start initialization now". |FastQCReport_DF| |Start initializ_DF| If
you don't want to generate QC-reports now, click "Delay initialization
till later" button. When you ready to start the analysis,  you can go
back to the Data Flow Runner page, click "2 files" link and select
"Multiple QC Report". |Start init (FastQC)| The calculations can
be started directly from the Multiple QC Report app page by clicking
"(re)-start computation if possible". |Multiple QC| Follow the progress
of your tasks in Task Manager. |TaskManager| When the computations are
finished, QC reports for both sequencing runs will appear on the
Multiple QC app page. Explore reports for each individual assay in
FastQC Report app by clicking on the app or file name in the Task
Manager.  Alternatively, go to "Created files" folder and look
for a folder containing the files created for "Raw Reads Quality
Control" data flow.   To describe raw reads quality control statistics
we will use the reports from our tutorial folder previously prepared by
our team.   To start with, we will open both of them in  `Multiple QC
Report <https://platform.genestack.org/endpoint/application/run/genestack/multiple-qc-plotter?a=GSF1001533&action=viewFile>`__ app
that interactively represents QC statistics for several raw assays at
once. |Screenshot 2016-02-19 20.53.39|   You can select samples of
interest, for example ones that are suitable for further analysis, and
put them in the separate folder by click on the "New folder with
selection" button. |Screenshot 2016-02-19 20.48.12| For paired reads
the quality control report contains statistics such as total nucleotide
count, GC content, number of reads, and number of distinct reads. Using
the Multiple QC Report app you can sort assays using QC-keys mentioned
above and metainfo-keys, such as "method" or "organism". Now when we
have the general impression of quality of raw reads we can go deeper and
get a more detailed statistics using  `FastQC
report <https://platform.genestack.org/endpoint/application/run/genestack/fastqc-report?a=GSF971377&action=viewFile>`__
for each individual sequencing run.   FastQC report contains several
quality control metrics outlined below:  

-  *Basic statistics* of raw data, for example the total number of
   reads processed, and GC content;

 

-  *Sequence length distribution* describing * * the distribution of
   fragment sizes in the analysed sequencing assay;

 

-  *Per sequence GC content* plot displaying the GC content across the
   whole length of each individual read;

 

-  *Per base sequence quality* plots depicting the range of quality
   scores for each base at each position in the analysed  sequencing
   assay;

 

-  *Per sequence quality scores* plot allowing the detection of poor
   quality sequences in the total sequences;

 

-  *Per base sequence content* plots representing the relative number of
   A, C, T, and G for each position in the tested sample;

 

-  *Sequence duplication level* plots representing the proportion of
   non-unique sequences which could be present in the library;

 

-  *Overrepresented sequences* providing the information on sequences
   that make up more than 0.1% of the total, and may either have a high
   biological significance or indicate contamination of the library. 

Table located on the left side of the page informs us which reports
raise concerns or report failures. In this case it is the  *Per base
sequence content*,  *Sequence duplication
levels * and  *Overrepresented sequences* metrics.   Raw data for both
sequencing runs failed the  *per base sequence content * metric.
Ideally, in a random library we would see four parallel lines
representing the relative base composition. Fluctuations at the
beginning of reads in the tested sample may be caused by adapter
sequences or other contaminations of the library. |Per base sequence
content (Run1)| The warning reported for the  *sequence
duplication * metric for the first sequencing run indicates that the
number of non-unique sequences in the assay has reached more than 20% of
the total. The average duplication levels for read mates are 1.50x and
1.48x.  *Sequence duplication* plot represents the relative number of
sequences having different duplication levels, and for  WGS
experiments, generally characterised by even coverage, this graph should
quickly drop to zero. Duplicates could correspond to PCR amplification
bias generated during library preparation or reading the same
sequence several times.   

|Seq duplication run1|

 Lastly, according to
the reports, the first sequencing run compared to the second one
contains some over-represented sequences — sequences that are highly
duplicated in a sample. In total, the app identified 1,052,139
sequences consisting of 'N'-bases.   The mentioned issues could be fixed
by performing appropriate preprocessing of the raw data. In this case,
we will trim low quality bases at the read ends and remove adaptors and
contaminants. Moreover, we will filter reads by quality score, so that
in further analysis we will only consider reads with high quality
(average Q≥20) score. Despite differences in the raw data quality, we
will apply the same preprocessing steps to both samples. It should be
stressed that after any applied preprocessing step you can check its
influence on the quality of raw reads using the FastQC app.   Now that
we have checked the quality of sequencing assays and decided on
the appropriate preprocessing steps, it is time to create the pipeline
for genetic variants analysis of WGS data from the raw data
preprocessing to the genetic variants annotation and filtering.

.. |DF list| image:: https://genestack.com/wp-content/uploads/2015/11/DF-list.png
.. |FastQCReport_DF| image:: https://genestack.com/wp-content/uploads/2015/12/FastQCReport_DF.png
.. |Start initializ_DF| image:: https://genestack.com/wp-content/uploads/2015/12/Start-initializ_DF.png
.. |Start init (FastQC)| image:: https://genestack.com/wp-content/uploads/2015/12/Start-init-FastQC.png
.. |Multiple QC| image:: https://genestack.com/wp-content/uploads/2015/12/Multiple-QC.png
.. |TaskManager| image:: https://genestack.com/wp-content/uploads/2015/12/TaskManager.png
.. |Screenshot 2016-02-19 20.53.39| image:: https://genestack.com/wp-content/uploads/2016/02/Screenshot-2016-02-19-20.53.39.png
.. |Screenshot 2016-02-19 20.48.12| image:: https://genestack.com/wp-content/uploads/2016/02/Screenshot-2016-02-19-20.48.12.png
.. |Per base sequence content (Run1)| image:: https://genestack.com/wp-content/uploads/2015/11/Per-base-sequence-content-Run1.png
.. |Seq duplication run1| image:: https://genestack.com/wp-content/uploads/2015/11/Seq-duplication-run1.png
