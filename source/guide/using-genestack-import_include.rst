Importing data
--------------

Supported file types
~~~~~~~~~~~~~~~~~~~~

Here is a list of file types that can be imported into Genestack.
Note that gzippped (.gz) and zipped (.zip) files are also supported.

- **Microarray Data** — raw microarray data obtained from a microarray
  experiment (you can import Affymetrix (CEL), Agilent (TXT) or GenePix microarray data (GPR));
- **Infinium Microarray Data** — raw intensity data files for Illumina Infinium Microarrays (IDAT);
- **Infinium Methylation Beta Values** — methylation data matrices contained Beta-values
  (methylation ratios) for Illumina Infinium Microarrays (.tsv, .txt);
- **Methylation Array Annotation** — methylation chip annotation containing information about
  association of microarray probes to known genes (.tsv);
- **Raw Reads** — raw sequencing data (FASTQ, SRA or FASTA+QUAL);
- **Microarray Annotation** — annotation file containing information about
  association of microarray probes to biological entities like genes,
  transcripts and proteins (.txt, .csv);
- **Continuous Genomic Data** — contains information on continuous genome
  statistics, e.g. GC% content (WIGGLE, WIG);
- **Discrete Genomic Data** — information on discrete regions of the genome
  with an exact start and end position (BED);
- **Mapped Reads** — reads aligned to a specific reference genome (BAM or CRAM);
- **Ontology Files** — OWL, OBO or CSV files used to annotate metainfo;
- **Reference Genomes** — reference genome sequence for a specific organism
  with annotation (FASTA and GFF/GTF/GFF3);
- **Variation Files** — Genetic variations files, storing gene sequence
  variations (VCF);
- **Gene List** — the file includes the list of genes for a specific organism (.gpr, .txt, .tsv, .xls);
- **Gene Expression Signature** — the file includes the list of genes with expression pattern
  specific to an organism phenotype according to statistical significance (filtering based on p-value) (.gpr, .txt, .tsv, .xls);
- **Gene Signature Database** — a list of annotated gene sets, that can be used in enrichment analysis (.gmt).


.. Two new file types:
.. Gene List: store a list of genes with possibly additional annotation
.. Gene Expression Signature: store a list of genes and expression pattern (Log FC) with possibly additional annotation

.. Two formats are accepted:
.. .grp file of genes in separate lines. This is imported as Gene List.
.. .tsv file. If the file contains both gene names and log fold changes, it is imported
.. as Gene Expression Signature. If the file only contains gene names, it is imported
.. as Gene List. The importer will look at the headers of the .tsv file to try to detect
.. which columns may correspond to gene names or log fold changes (common variations are
.. supported such as ‘gene’/ ‘symbol’ for gene names, and ‘logFC’/’log fold change’ for log
.. fold changes). If it fails to detect them, the user will be asked to manually choose the
.. file type and specify the file headers corresponding to gene names or log fold changes.
.. Gene symbols and Ensembl/Entrez gene IDs are currently supported for gene names.

When you import files that are detected as raw sequencing or microarray data,
Genestack automatically creates a **dataset**, a special type of folder, and adds the assays to it.
Additional documents in any format (e.g. PDF, Word, text, etc.)
can be imported as **attachments** to a dataset. We will discuss the use of attachments below.
Some types of files, namely Reference Genome, Gene List, Gene Expression Signature,
Gene Signature Database, Genetic Variations, Ontology Files, Dictionary, Microarray Annotation,
Methylation Array Annotation, Infinium Beta Values, are not wrapped in
datasets on import because they are rarely uploaded and processed as batches.

When you perform any analysis on Genestack, other data types, which cannot be imported, can be created such as:

- **Affymetrix/Agilent/GenePix Microarrays Normalisation** — file with
  normalized Affymetrix/Agilent/GenePix microarrays data;
- **Differential Expression Statistics** — expression statistics for
  change in expression of individual genes or other genomic features between groups of samples,
  such as fold-changes, p-values, FDR, etc.;
- **Genome Annotations** — a technical file used for matching GO terms and
  gene symbols to gene coordinates;
- **Mapped Read Counts** — file is produced from Mapped Reads and contains the number of reads mapped to each feature of a reference
  sequence.

.. verify

Data import
~~~~~~~~~~~

There are several ways you can access the **Import** application:

- using the **Import data** link on the Dashboard;

.. image:: images/WP_import.png
   :scale: 90 %
   :align: center

- clicking the **Import** button in the File Manager;

.. image:: images/FM_import.png
   :scale: 90 %
   :align: center

- using an **import template**. We will describe what import template is and how to
  use it later in the guide.

.. image:: images/IT_import.png
   :scale: 90 %
   :align: center

Import data consists of 3 steps.
The first step results in temporary **Upload** files with your data created in the platform.
On the second step the biological data type is assigned to your imported data.
On the third step you can fill required metadata or import it from a text file.

Step 1: Getting data into the platform
+++++++++++++++++++++++

There are two ways to have your data imported into the platform:

1. **Upload data from your computer** — select or drag-and-drop files.

.. image:: images/import_start.png
   :scale: 80 %
   :align: center

2. **Import from URLs (FTP or HTTP/HTTPS)** — specify URLs for separate files or
   directories.

.. image:: images/URL_import.png
   :scale: 80 %
   :align: center

If you made a mistake during one of the following data import steps, you can
always reuse your previous **Upload** files: just select previously uploaded
data with **Use previous uploads** option and add more data if needed.
This way you will not have to upload your data for the second time
because everything recently uploaded is available as a temporary **Upload** file.

.. note:: **What is an Upload file?**

          *Describe the Upload file is a temporary file that stores your imported data
          and may be removed by the platform when it is no longer used. It's only purpose
          is to store your data during the Import Data steps and give you ability to ensure
          your Genestack Files are created and initialized correctly. After Genestack Files are
          created and initialized successfully, these temporary **Upload**
          files can be safely removed from the platform and no data will be lost.*

Uploading from you computer is done in multiple streams to increase upload speed.
Import from URLs is done in the background, which means that even while these files
are being uploaded, you can edit their metadata and use them in pipelines.

.. image:: images/uploading_step.png
   :scale: 80 %
   :align: center

If during uploading you lose your Internet connection, you will be able to
resume unfinished uploads later.

.. image:: images/resumed_uploads.png
   :scale: 85 %
   :align: center

Click the **Import files** button to proceed.

Step 2: Format recognition
++++++++++++++++++++++++++

After your data is uploaded, Genestack automatically recognizes file formats
and transforms them into biological data types: raw reads, mapped reads,
reference genomes, etc. Format conversions will be handled internally by
Genestack. You will not have to worry about formats at all.

.. image:: images/file_recognition.png
   :scale: 80 %
   :align: center

If files are unrecognized or recognized incorrectly, you can manually allocate
them to a specific data type: drag the **Upload** file and move it to the green
"Choose type" box at the top of the page.

.. image:: images/unrecognized_uploads.png
   :scale: 80 %
   :align: center

Choose the data type you find suitable:

.. image:: images/file_types_box.png
   :scale: 80 %
   :align: center

Click the **Create files** button to proceed.

Step 3: Editing metainfo
++++++++++++++++++++++++

In this step, you can describe uploaded data using an Excel-like spreadsheet.
Importantly, during this step, the import has already completed — you will
notice a message at the top of the page with the name of the folder where
the imported files are located (named "Imported on <date> <time>"), and
offering to share the data:

.. image:: images/import_edit_metainfo.png

By default, you see all metainfo fields available for files, you can fill them
or create new custom columns. Click the **Add column** button, name new metainfo
field and choose its type (Text, Integer, etc.):

.. image:: images/add_metainfo_field.png

You can also choose to apply a naming scheme. This allows you to generate
file names automatically based on other metainfo attributes.

.. image:: images/naming_scheme.png

Metainfo fields can be associated with specific dictionaries and
ontologies. We pre-uploaded some public dictionaries such as
the `NCBI Taxonomy`_ database for the "Organism" field, the Cellosaurus_ (a resource on cell lines),
the ChEBI_ for chemical compounds, and the `Cell Ontology`_ (cell types in animals).

.. _NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy
.. _Cellosaurus: http://web.expasy.org/cellosaurus/description.html
.. _ChEBI: https://www.ebi.ac.uk/chebi
.. _Cell Ontology: http://www.obofoundry.org/ontology/cl.html

We also created our own controlled vocabularies to cover Sex, Method and Platform fields.
You can find out more about ontologies in the :ref:`public-experiment-label` section.

Import with templates
~~~~~~~~~~~~~~~~~~~~~

You can create your own custom dictionary by importing it into the
platform as OWL, OBO or CSV file and attach it to the import template.

.. note:: **What is an import template?**

          Import templates allow you to select what metainfo attributes of your imported
          files will be tightly controlled (so you don’t lose any information in the
          process). Import templates allow you to set default fields for file metadata
          based on file type (e.g. Datasets, Discrete Genomic Data, Genetic
          Variations, etc.). Of course, if you’re only importing mapped reads, you don’t
          need to specify metainfo attributes for other data types.

You can select which import template to use in two ways: from the Dashboard,
or during the 3rd step of the import process by right-clicking on the
import template name ("Default template" is for the public one). You can create a copy of existing
import templates with **Make a copy** option in the context menu.

.. image:: images/copy-import-template.png
   :align: center
   :scale: 85 %

Genestack will attempt to fill metainfo fields automatically, but you can always
edit the contents manually during the import process. By using metainfo
templates you can make sure that all of your files will be adequately and
consistently described so you will not lose any valuable information. For
example, here is the list of metainfo attributes used by default to describe
Reference Genome data:

.. image:: images/default_import_template.png

**Import template editor** application allows to modify existing import templates and create
new ones with proper metainfo fields, requirements and controlled vocabularies. To access
the app right-click on a template's name and select the **Import template editor** from
the "Manage" submenu. To create new template on the basis of the default one you can also click
**Add import template** one the Dashboard.

.. image:: images/import_templates.png
   :scale: 45 %
   :align: center

Now let’s say you wish to create an import template to control
the metainfo attributes of raw reads (e.g. you always need to know the
tissue and sex of your samples). In order to do this, click on **Add import template**,
then look for the table related to Raw Reads and for the
fields "tissue" and "sex", change the required fields to *Yes*. As you can
see, the system controls what type of information can you put into your
metainfo fields. In this case, for tissue the system will map your entries to
the Uberon_ ontology (an integrative multi-species anatomy ontology) and
the metainfo type must be text.

.. _Uberon: http://uberon.github.io

.. image:: images/edit-template.png

If you want to **add other metainfo fields** that are not included in the table
already, you can do this at the bottom of the table where there are blank
spaces. For each entry, you must specify whether or not this field is
required and what is its metainfo type (e.g. text, yes/no, integer).

.. image:: images/metainfo_type_editor.png

If you are using a file kind that is not yet listed, you can add a new one by
clicking on the **Add file kind** button. Keep in mind that file kinds are
defined in Genestack — you will not be able to create a template entry for a
file kind that is not used on the platform.

When you are done, click on the blue **Import using this template** button.
This will take you to the **Import Data** app, where you can go through the three import
steps described above. You can find all the imported files in the "Imported" folder which can be accessed from the Dashboard and from the File
Manager.

Metadata import
~~~~~~~~~~~~~~~

Apart from editing metainformation manually, you can also import and validate the metainfo attached to the assays and
to the dataset on the platform.

.. image:: images/import_from_spreadsheet.png

Click **Import data from spreadsheet** button and select a local CSV or
Excel file containing metadata you would like to associate with the imported files.

.. image:: images/import_metainfo.png

Note that names in the first column in the file with metadata should exactly match names of the data
samples on the platform, based on the first "Name" column. For example, in our case metainfo
for the second sample does not match to any assays and is highlighted in red.

.. image:: images/import_metainfo_table_red.png

Use the **Select file** option to manually allocate the imported metadata to an appropriate
file.

.. image:: images/import_metainfo-select-file.png

Columns that are mapped to a metainfo field from the dataset's
template (by default data are imported with "Default" template) are highlighted in green.

.. image:: images/import_metainfo_table-green.png

On this step for each column you can specify whether it should be imported or not, and if it
should be mapped to some metainfo key from the import template, by clicking on the column header.

.. image:: images/metainfo-import-matching.png

Click **Import** when you finish editing the table. As a result, the table on the Metainfo Editor
page is filled in with metadata from the Excel-file.

.. image:: images/import_metainfo_complete.png

For instance, in this case we added new column
"Age" and filled "Organism", "Sex", "Tissue" and "Platform" columns that came from the default template.


Attachments
~~~~~~~~~~~

While importing a dataset into Genestack, you can also attach various files to it such as, for
example, a PDF file with the experiment plan or an R script, etc. When you open your newly-imported
datasets, all of the attachments will accompany it. They will be safely
stored on Genestack, so later you can download them from the platform, in case
they get lost on your computer.

**How to upload an attachment?**

Attachments should be uploaded together with the dataset. In the the Data Import application, choose
the attachments from your computer along with your dataset. The platform will
recognize the raw data, and all additional files that were unrecognised will
be added to the dataset as attachments.

.. image:: images/attachments.png

Besides, you can upload more attachments, or remove attachments in the Metainfo Editor.

.. image:: images/exp_attachments.png
