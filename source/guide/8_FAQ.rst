FAQ
====

Where do I find data shared with me?

If they have been linked, you can find them in the corresponding
group subfolder folder within the “Shared with me” folder. Otherwise,
they can be found using search.

Where do I find the data flows I have created?

How do I reuse a data flow?

Why are my tasks failing?

What’s the difference between Data Flow Runner and Data Flow Editor?

Data Flow Editor is used to create data flow templates: e.g. selecting
source files.

When you want to use the data flow to run your analysis, on the Data
Flow Editor page you can click on “Run Data Flow” button, which will
take you to Data Flow Runner. Here you can not only edit source files
and parameters, but also start initialization of your files.

How do I initialize the files?

How do I create a data flow? To create a data flow, select the data you
wish to analyse and choose the first app you wish to use in your
analysis. On the app page, using the “add step” button, add the rest of
the desired steps. Once you are done, click on the name of the file (or
files) at the top of the page, go to Manage, and click on Create New
Data Flow. Your new data flow can be found in the Created Files folder

If you don’t want to create a data flow from scratch, but rather re-use
the same analysis pipeline used to create a file, click on the name of
that file, go to Manage, and select Create New Data Flow.

Selecting File Provenance instead of Create New Data Flow will show you
the pipeline (in the form of a data flow) that was used to create this
file.

Read more about data flows in this tutorial_:


What’s the difference between BWA and Bowtie2?

The biggest differences between the two aligners are:

1)The way of accepting or rejecting an alignment.

BWA: counts the number of mismatches between the read and the
corresponding genomic position.

Bowtie2: uses a quality threshold bases on the probability of the
occurrence of the read sequence given an alignment location.

2)Accepting colorspace data:

BWA: No.

Bowtie2: Yes.

How does Genestack process paired-end reads?

There are three types of raw reads that our platform supports:

-  single-end (1 file locally, 1 file in Genestack)
-  paired-end (2 files locally, 1 file in Genestack)
-  paired-with-unpaired (3 or 4 files locally, 2 files in Genestack)

During import, Genestack recognises these types and imports them in
their respective format-free form.as 1 or 2 files. If the platform
cannot recognise the files automatically, you can allocate the files
manually.

What’s the difference between an experiment and a folder?

The main difference between a folder and an experiment is in the content
specificity.

Folders work just like folders on your computer and can contain various
biological objects\ :sup:``[bd] <#cmnt56>`__`\ : your assays, processed
files, output reports, etc.

Experiments contain only sequencing assays, but provide additional
features for storing experimental details, e.g. attached pdfs with
experiment notes. When you upload raw reads onto Genestack, they will
automatically be imported as one experiment. On the experiment page, you
can click on “View details” to read the summary of the experiment, get
more information about overall design and experiment type, contributors,
and find links to public databases.

All public experiments available on our platform are provided as
experiments, not folders. Remember you can share both experiments and
folders with other Genestack users using Groups. You can read more about
this in our `“Getting Started”`_ tutorial.

What’s the difference between masked and unmasked reference genomes?

In general, when a genomes is “masked” it means that all repeats and low
complexity regions of your reference genome (detected
by `RepeatMasker`_ tool)
are hidden away and replaced with “N”s, so that they will not be aligned
to.

We do not recommend using a masked genome, as it always  results in a
 loss of information. Masking can never be 100% accurate, and can lead
to an increase in the number of falsely mapped reads. If you’d like to
perform filtering, it’s better to do it after the mapping step.

In “soft-masked” genomes, repeated and low complexity regions are still
present, but they have been replaced with lowercased versions of their
nucleic base.

“Unmasked” genomes contain all repeats and low complexity regions
without any changes.

How do I change the name of the file?

How do I add files to the briefcase?

Where can I find the import templates I created?


.. _tutorial: https://genestack.com/tutorial/reproducing-your-work-with-data-flows/
.. _“Getting Started”: https://genestack.com/blog/2016/01/06/getting-started/
.. _RepeatMasker: http://www.repeatmasker.org/&sa=D&ust=1480960532173000&usg=AFQjCNE4ktR5xI4yZEvRi94d-Tc1QkJnvA

