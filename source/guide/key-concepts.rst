Genestack key concepts
======================

Genestack is a platform to manage, analyse and visualize bioinformatics data, with a
focus on next-generation sequencing (NGS) and microarray data. Genestack was built around
several key concepts which are described below.

Rich metadata system
--------------------

Bioinformatics data always comes with complex metadata, which is essential for search
and analysis.

On Genestack, each file has metadata directly attached to it, which can store biological
information as well as technical information (file owner, application parameters, etc.).
Metadata values can support different formats (text, numeric data, dates, measurements with units, etc.),
and they can be viewed and edited on Genestack, as well as easily imported from an Excel or CSV spreadsheet.

Genestack files can be searched and filtered by metadata. Additionally, controlled vocabularies and ontologies can be
used to validate metadata, to ensure harmonisation of metadata and improve the quality of the search
results. Finally, metadata templates can be defined by users to ensure that all the collaborators of a group are using
consistent metadata conventions.

Public data
-----------

Genestack is pre-loaded with millions of publicly available
datasets from major repositories like ArrayExpress_,
GEO_, SRA_, and ENA_, as well as numerous reference genomes for multiple organisms from
Ensembl_ and UCSC_. In practice, this means that the platform can serve as
a data repository, that allows users to work both on private and public
data seamlessly. Our Data Browser application can help you to easily find not only
datasets you are interested in but also analysis results performed on these data.

Interactive data analysis and visualisation
-------------------------------------------

On Genestack, you will find many applications to process bioinformatics data, but you will also find
a range of graphical, interactive applications that help users better understand their data.
This ranges from our FastQC Report application
that visualises the quality of raw or preprocessed sequencing
reads, to our Variant Explorer application, performing real-time interactive
genomic variant filtering by type, impact, quality, frequency, etc.
All the Genestack applications are described in the section  :ref:`pipelines-and-applications-label`.

Format-free files
-----------------

When doing bioinformatics, a lot of effort is often spent trying to get your files
in the right format, to work with the right software.

In Genestack, we address this issue by dealing with high-level biological objects
instead of traditional file formats. For instance, this means that when you upload
some sequencing reads to Genestack, whether they are in FASTQ, SRA, paired-end or not,
compressed or not, they will appear in Genestack as "Raw Reads".
And subsequently, you will only be able to use these files with applications that
know how to work with Raw Reads; if any format conversion needs to be done, they
are done automatically and internally, and you will not have to worry about them.

Therefore, we say that files in Genestack are "format-free" in the sense that, as a
user of the system, you do not have to think about file format compatibility and
conversion inside Genestack, because Genestack does that for you.

.. note:: **Formatting Hell**

          You might be wondering why we made our platform format-free and why this
          is such a big deal. In the current landscape of bioinformatics there
          seems to be a never-ending number of formats your data might be saved
          in. There are a few prominent formats used in NGS data analysis, like
          FASTQ, BAM and VCF. But very often new programs come with
          new formats. Bioinformaticians say they spend almost 80% of their time
          worrying about data grooming and file reformatting and only 20% on
          actual data analysis. Now, that’s insane, isn’t it? On Genestack you
          don’t have to worry about formats at all — our platform takes care of all the
          routine tasks so that you can focus on your work.


Reproducibility
---------------

Every file on Genestack "remembers" how it was made: all applications,
tool versions, parameters, and source files are recorded in the metadata
of each file.

This way, you can select any file in Genestack and visualize its "provenance",
a graph representing exactly what steps were taken to produce this file, which source
files and which applications were used to produce it.

Genestack hosts multiple tool versions at any given time in case you want to reproduce past results.

.. _ArrayExpress: https://www.ebi.ac.uk/arrayexpress/
.. _GEO: https://www.ncbi.nlm.nih.gov/geo/
.. _SRA: https://www.ncbi.nlm.nih.gov/sra/
.. _ENA: http://www.ebi.ac.uk/ena
.. _Ensembl: http://www.ensembl.org/index.html
.. _UCSC: https://genome.ucsc.edu/