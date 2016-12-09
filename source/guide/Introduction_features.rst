We work hard to make Genestack Platform the most intuitive and
user-friendly bioinformatics software. That being said, bioinformatics
is always going to be challenging at times and anyone doing any type of
NGS data analysis has hundreds of questions running through their head.

In this guide we will answer any questions you might have, help you
troubleshoot problems and eventually teach you how to make the most out
of our platform.

We wrote this guide for everyone: whether it’s your first or hundredth
time doing bioinformatics you will find something useful here. We will
cover general aspects of bioinformatics and talk about the specifics of
our platform. We recognize that some of you may be bioinformatics
experts and hence won’t need to read about the basics of bioinformatics.
This is why we have colour-coded this guide: general aspects
of bioinformatics are labeled green.

Before we begin, we’d like to remind you that the community edition of
Genestack Platform is meant to be used by academic and non-commercial
users. The Enterprise edition is designed for our commercial users who
will benefit from features such as a broad range of customization
options, tech support or running the platform on or off the cloud, to
name a few things. If you have any further questions, please get in
touch with our team at contact@genestack.com

Introduction to Bioinformatics
==============================

Bioinformatics can be defined as the application of computational
techniques to make sense of and organize large-scale biological
information.

With costs of sequencing dropping below $1000 per genome, enormous
amounts of sequencing data are being produced each day. Soon enough, it
won’t just be trained bioinformaticians who will need to analyse -omics
data daily. Doctors, lab technicians, researchers, PhD students,
postdocs: anyone working with personalised medicine, genetics-related
research, drug design or biotechnology will have to be able to perform
at least basic -omics analysis types.

Luckily, during the last five years, a range of bioinformatics platforms
emerged on the market promising to bridge the gap between bioinformatics
and non-bioinformaticians. Genestack is one of these platforms and our
community edition is freely available online for academic and
non-commercial users. Our motto, “Do Bioinformatics Faster”, is visible
in all we do: the platform helps the user  shorten the path from idea to
results by automating routine tasks. The system also suggests which
applications can be used to analyse the data at each step and our
interactive visual apps make interpreting the results easy.  

Though you don’t need to know how to code to use the platform, a basic
level of understanding of the sequencing data analysis process and the
tools involved is crucial for correct data analysis and results
interpretation.

In this guide, we will talk about the basic concepts of sequencing
analysis, the steps involved and how to interpret your data. Look out
for the green paragraphs in this guide - they will introduce you to the
basic concepts of NGS data analysis.

If you want more information on the concepts and history of
bioinformatics, we’d like to refer you to this awesome article_ by N.M.
Luscombe, D. Greenbaum and M. Gerstein (2001)

Genestack Features
==================

Here are a few of the concepts at the core of Genestack:

Rich metadata system
--------------------

Genestack is a data-centred platform. We believe your data should be
harmonized, searchable and well-managed, so that you can use it over and
over again. This is why we have invested so much of our time into our
data and metadata management system. When using the platform you will
notice excel-like import templates that allow users to specify required
metainfo fields; a context-sensitive metadata editor with ontology
autosuggest support; rapid access to private, shared and public
data; and many, many more useful tools and tricks that leverage our deep
metadata integration. 

Interactive data analysis and visualisation
-------------------------------------------

On Genestack you will find a range of graphical, interactive apps that
help users better understand their data. This ranges from our FastQC
report app that visualises the quality of raw or preprocessed sequencing
reads to our Variant Explorer app, performing real-time interactive
variant filtering by type, impact, quality, frequency, etc. We will
discuss all of the apps in greater detail here.

Format-free files
-----------------

This is one of the very core ideas of our platform and something that
saves our users a lot of time. When you upload your files into
Genestack, they “lose” their original format and become meaningful
biological objects (e.g. “sequencing assays” or “reference genomes”).
Genestack deals with any formatting-issues allowing users to spend more
time on actual data analysis. After you are done with your analysis you
can export the data in the format of your choice.

**NOTE: Formatting Hell**

You might be wondering why we made our platform format-free and why this
is such a big deal. In the current landscape of bioinformatics there
seems to be a never-ending number of formats your data might be saved
in. There are a few prominent formats used in next generation sequencing
analysis, like FastQ, BAM and VCF. But very often new programs come with
new formats. Bioinformaticians say they spend almost 80% of their time
worrying about data grooming and file reformatting and only 20% on
actual data analysis. Now, that’s insane, isn’t it? On Genestack you
don’t have to worry about formats at all – our OS takes care of all the
routine tasks so that you can focus on your work.

Reproducibility
---------------

Every file on Genestack “remembers” how it was made: all applications,
tool versions, parameters, and other settings are saved in the File
Provenance for each of your files. When you want to repeat an identical
analysis on a different set of data, you can do this easily using data
flows. Genestack platform hosts multiple tool versions at any given
time in case you want to reproduce past results.

Data flows and Delayed Initialization
-------------------------------------

Data flows on Genestack are visual representations of bioinformatics
pipelines. When you build a data flow (or use an existing one), the
computations do not start automatically. Instead, all the participating
apps create files. For instance, when you create a pipeline consisting
of three preprocessing apps (e.g. Trim Low Quality Bases, Trim Adaptors
and Contaminants and Filter by Quality Score) and use “Sequencing Assay”
as an input file – three files will be created: “Trimmed Sequencing
Assay”, “Trimmed Trimmed Sequencing Assay” and “Filtered Trimmed Trimmed
Sequencing Assay” (those names may sound funny, but are helpful when
tracking what apps you used already. Remember you can change the file
names!).

**NOTE: How can you create a data flow?**
Forum_ `Tutorial`__

When an application creates a new file, it specifies what should happen
when it is initialised: a script, a download, indexing, computation.
However, before the initialization has begun, you can change the
parameters, replace source files, or add more files to the pipeline.

In practice this  means that uninitialised files are cheap and quick to
create, can be reconfigured, used as inputs to applications to create
other files, and later computed all at once. Remember – you always need
to initialize your files to view your results (same rule goes for data
flows). Once you do, any further changes become impossible.

**NOTE: How do you initialize your files?**
Forum1_  Forum2_ `Tutorial`__

Public experiments collection from databases and archives
---------------------------------------------------------

Genestack platform is preloaded with millions of publicly available
experiments from major repositories like ArrayExpress, GEO, SRA and ENA,
as well as numerous reference genomes for multiple organisms from
Ensembl and UCSC. In practice, this means that the platform can serve as
a data repository, that allows users to work both on private and public
data seamlessly.

.. _article: https://www.google.com/url?q=https://www.ebi.ac.uk/luscombe/docs/imia_review.pdf&sa=D&ust=1480960531653000&usg=AFQjCNFUGLBg9Y8pGX_C7QUt__SuRovLEw
.. _Forum: http://forum.genestack.org/t/creating-new-pipelines-on-genestack/26
.. _Forum1: http://forum.genestack.org/t/initializing-only-1-process-from-the-data-flow/27
.. _Forum2: http://forum.genestack.org/t/how-to-map-or-pre-process-several-raw-reads-files-at-once/28
__ https://genestack.com/blog/2016/06/21/building-pipelines-reproducibility/#buildingapipeline
__ https://genestack.com/blog/2016/06/21/building-pipelines-reproducibility/#fileinitialization
