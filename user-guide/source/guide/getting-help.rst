Getting Help
============

Tutorials 
---------

In our tutorials we take you through examples based on public
experiments, re-analysing the data to demonstrate the features of
particular applications. All files used (from raw sequencing reads and reference
genomes to results and data visualisations) are stored on the platform
in the "Tutorials" folder.

We have prepared tutorials on Differential Gene Expression, Methylation
Profiling, Whole Exome Sequencing, and Whole Genome Sequencing and me
plan to add more in the future.

You can find all of our tutorials here: https://genestack.com/tutorial/

Getting in touch with Genestack 
-------------------------------

There are various ways to get in touch with our team.

*Chatra*: a chat window is available at the bottom of every page on the website.
During our working hours, at least one of us will be always available to help.
Outside of working hours, you can leave us a message using the chat
window and we will get back to you as soon as we can.

*Forum*: when you have a question about the platform, try posting
it on our forum. You can find it on forum.genestack.org. Our community
will definitely be keen to help you, and if not, our team regularly checks the forum and will answer any unanswered questions.

*Email*: you can email us at support@genestack.com 

FAQ
---

Where do I find data that was shared with me?
*********************************************

If the files were linked into the group folder, you will find them there.
You can access the group folder from the file browser, under the "Shared with me" section. 
Otherwise, the files can be found via search, if you know their name or accession.

How do I reuse a data flow?
***************************

Open a data flow you would like to run in the "Run Dataflow" application. On the
application page you can set input files and additional files (e.g. reference genome)
that are required for analysis.

What is the difference between Data Flow Runner and Data Flow Editor?
*********************************************************************

Data Flow Editor is used to create data flow templates, e.g. selecting
source files.

When you want to use the data flow to run your analysis, on the Data
Flow Editor page you can click on "Run Data Flow" button, which will
take you to Data Flow Runner. Here you can not only edit source files
and parameters, but also start initialization of your files.

How can I initialize multiple files at once?
********************************************

You can use the File Initializer application for this.
Select the files you want to initialize in File Browser, right-click on them
and go to Manage > File Initializer. From there, click the button "Initialize all"*[]:

How do I create a data flow?
****************************

To create a data flow, select the data you
wish to analyse and choose the first application you wish to use in your
analysis. On the application page, using the "add step" button, add the rest of
the desired steps. Once you are done, click on the name of the file (or
files) at the top of the page, go to Manage, and click on Create New
Data Flow. Your new data flow can be found in the Created Files folder

If you do not want to create a data flow from scratch, but rather re-use
the same analysis pipeline used to create a file, click on the name of
that file, go to Manage, and select Create New Data Flow.

Selecting File Provenance instead of Create New Data Flow will show you
the pipeline (in the form of a data flow) that was used to create this
file. Read more about data flows in this tutorial_.

What is the difference between BWA and Bowtie2?
***********************************************

The biggest differences between the two aligners are:

- the way of accepting or rejecting an alignment;

BWA counts the number of mismatches between the read and the
corresponding genomic position. Bowtie2 aligner uses a quality threshold bases on the probability of the
occurrence of the read sequence given an alignment location.

- accepting colorspace data;

BWA tool does not support data in colorspace data, while Bowtie2 is able to align such files.

How does Genestack process paired-end reads?
********************************************

There are three types of raw sequencing reads that our platform supports:

-  single-end (1 file locally, 1 file in Genestack);
-  paired-end (2 files locally, 1 file in Genestack);
-  paired-with-unpaired (3 or 4 files locally, 2 files in Genestack).

During import, Genestack recognises them and imports them in
their respective format-free form. If the platform
cannot recognise the files automatically, you can allocate the files
manually.

What is the difference between an experiment and a folder?
**********************************************************

Experiments are a special kind of folder, which can only contain
assays, e.g. "raw" experimental data.

What is the difference between masked and unmasked reference genomes?
*********************************************************************

In general, when a genomes is "masked" it means that all repeats and low
complexity regions of your reference genome (detected
by RepeatMasker_ tool)
are hidden away and replaced with "N"s, so that they will not be aligned
to.

We do not recommend using a masked genome, as it always  results in a
 loss of information. Masking can never be 100% accurate, and can lead
to an increase in the number of falsely mapped reads. If you would like to
perform filtering, it is better to do it after the mapping step.

In *soft-masked* genomes, repeated and low complexity regions are still
present, but they have been replaced with lowercased versions of their
nucleic base.

*Unmasked* genomes contain all repeats and low complexity regions
without any changes.

How can I rename a file?
************************

In the File Browser, click on the file name and select the "Rename" option in the context menu.

I've created a file. Where can I find it?
*****************************************

All the files that you create within Genestack go to your "Created files" folder.
Files that you import to Genestack go into your "Imported files" folder.
Both folders are accessible from the dashboard and the file browser.

.. _tutorial: https://genestack.com/tutorial/reproducing-your-work-with-data-flows/
.. _Getting Started: https://genestack.com/blog/2016/01/06/getting-started/
.. _RepeatMasker: http://www.repeatmasker.org/&sa=D&ust=1480960532173000&usg=AFQjCNE4ktR5xI4yZEvRi94d-Tc1QkJnvA
