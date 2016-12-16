
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

.. note:: **Formatting Hell**

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

.. _Forum: http://forum.genestack.org/t/creating-new-pipelines-on-genestack/26
.. _Forum1: http://forum.genestack.org/t/initializing-only-1-process-from-the-data-flow/27
.. _Forum2: http://forum.genestack.org/t/how-to-map-or-pre-process-several-raw-reads-files-at-once/28
__ https://genestack.com/blog/2016/06/21/building-pipelines-reproducibility/#buildingapipeline
__ https://genestack.com/blog/2016/06/21/building-pipelines-reproducibility/#fileinitialization
