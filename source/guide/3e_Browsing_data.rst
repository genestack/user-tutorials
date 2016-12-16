Browsing Data
-------------

Efficient data search and browsing are at the core of Genestack. The
platform provides  rapid access to private, shared, and public data
analyses results; facilitates search for studies and assays across your
private, public, and shared data; and accepts queries using synonyms,
ontology expansions, and chemical similarity.

Experiment Browser
~~~~~~~~~~~~~~~~~~

Genestack Platform provides a rich collection of public experiments from SRA, ENA, GEO
and ArrayExpress, and Genestack synchronizes data from these databases regularly keeping
it up-to-date. There are currently more than 3 millions sequencing and microarray assays from over
100 000 experiments.
Experiment Browser app allows to browse these public datasets, as well as your private
experiments or the ones shared with you on Genestack. You can access the Experiment
Browser either from the Welcome Page or the Shortcuts Menu on the left-hand side.

You can search relevant data with **a free-text query**, and you can further
filter down experiments by **metadata attributes** using the checkboxes
on the left. These attributes are generated based on the metadata available for experiments.
For instance, you can set the filters 'Access', 'Method'
and 'Organism' to 'Public', 'RNA-Seq', 'Mus musculus', respectively,
to filter out publicly accessible data on mice obtained from  mouse RNA-Seq data.

|ExperimentBrowser|

Click **Save N matching assays** link to explore
the list of matching assays and save them into one folder.

|SaveMatchingAssays|

Moreover, Experiment Browser allows to find bioinformatic analyses results
associated with raw data. Indeed, if there are analysis performed on a given experiment,
and you have an access to these results (i.e. theses are yours or shared with you),
then under the experiment name you will see **"View N analysis results"** link.
Click it shows you the list of existing resulting files such as, for example, QC reports
or Genome Browser pages.

|AnalysisResults|

Clicking on the name of any of the found experiments will take you to
the **Metainfo Editor**, where you can learn more that experiment
and get some information on samples.

|EditMetainfo|

Besides that, on the Metainfo Editor page you can run the assays through a pipeline with
**"Star new data flow with application"**:
you can use existing data flows matching for the assays of choice or build pipeline step-by-step
 option.

|NewDF|

Click **"Use files in data flow"** button and
From the metainfo editor, you can also open the experiment in the **file manager** by clicking on
the experiment's name at the top of the page and selecting **Explore > File Manager**.

|fromMEtoFB|

File Manager
~~~~~~~~~~~~
**File Manager** is where you can easily access all of your private, public
and shared data. Here are other apps that help users better analyse their data and find
links between various results.

Clicking on the home icon will take you to the File Manager – a central
place on the platform, as it contains all of your files (you probably
got that already).

|FileManager|

The panel (tree view) on the left side is our file system navigator.
Here you can see many different folders. Let’s look at them in greater
detail:

**Created files** contains everything you have created on Genestack
Platform. Created a new import template? You’ll find it there. Processed
some of your files? You’ll find the results there. Created a new data
flow and want to share it? It will be in the Created files folder.

The files are organized by date, with oldest ones on top (however, you
can change this order to show the most recent ones - just click on the
header of the “Last Update” column). If you created a couple of files at
once using a data flow they will be located in one folder (called “Files
for XYZ data flow run <date>). In these folders you will find the very
result of your analysis (e.g. Genetic Variations file containing found
mutations), results of all intermediate analysis steps (e.g.
preprocessed reads, mapped reads etc. created by the apps participating
in your pipeline as you remember each contributing app creates a file),
as well as all original files (“Original Files for XYZ, a sub-folder in the “Dependencies” folder).

**Imported files** contains everything you have ever imported, organized by
date: all files imported at the same time (during one import action)
will be located in the same folder (until you move them around etc). 

Raw uploads contains all the files you’ve uploaded into Genestack -
fastq and bam files, pdf documents, excel tables etc.


.. note:: **What’s the difference between raw uploads and imported files?**

          When you have just started importing your files (in various formats like
          FASTQ, BAM etc), they all go to the specific storage area (“Raw uploads”
          folder). During import Genestack will recognize these uploaded files and
          allocate them to appropriate biological types (you can also do it
          manually), e.g. sequencing assays, mapped reads etc. These meaningful
          biological objects is what you work with on our platform and these are
          located in the “Imported files” folder.

**Exports** folder contains export files with download links. For example,
sets of exported microarrays. Get more information about Exporting Data from Genestack in
the Export Data section.

Below these four grouped folders, you will see two more: Shared with me
and Public Data.

**Shared with me** contains all files that other users have shared with
you or that you shared with other users. Our platform has collaboration
at its heart, but in order to keep things simple at this point, we’ll
talk about sharing at the very end of this guide (+ link to the guide
part about it).

**Public Data** contains all of the goodies we have preloaded the platform
with to make life a bit simpler for our users. This folder contains:

|PublicData|

#. **Codon tables**: currently 18 different tables such as yeast
   mitochondrial, vertebrate mitochondrial, blepharisma macronuclear
   etc;
#. **Dictionaries**: used for metainfo editing and curation, e.g. sex,
   sequencing platform, NCBI taxonomy. Read more about dictionaries in
   "Data and Metainfo Management" section;
#. **Example results**: so you can play around with our platform and see
   what types of visualizations are available;
#. **External databases**: sets of sequences with associated annotation;
   e.g. greengenes for 16S rRNA;
#. **Genome annotations**: for a range of different organisms and platforms
    (for WES Analysis);
#. **Microarray annotations**: annotation lists to be used as the
   translation table to link probes and common public domain sequences;
#. **Public analyses**: all files created during re-analysis of previously
   published data sets;
#. **Reference genomes**: various reference genomes for the most commonly
   analysed organisms;
#. **Public data flows**: all data flows available to our users, including
   tutorial data flows and the ones found on the Welcome page;
#. **Public experiments**: this is a feature we’re particularly proud of. We
   have preloaded the platform with thousands and thousands of publicly
   available experiments, from public repositories such as GEO,
   ArrayExpress, SRA, and ENA. Currently we have about 100,000
   experiments in our database. If you want to know more about a specific experiment use the
   Experiment Viewer app;
#. **Tutorials**: the folder contains files we use as examples during
   various tutorials. To read more on particular analysis types, go to
   https://genestack.com/tutorials/.

   Currently, we can offer you the following tutorials:

  -  `Getting Started With Genestack Platform`_
  -  `Testing Differential Gene Expression`_
  -  `Whole Genome Bisulfite Sequencing Analysis`_
  -  `Whole Exome Sequencing Analysis`_
  -  `Whole Genome Sequencing Analysis`_

To access the **context menu** for a given file, you can either right or left click
on the respective entry in the file browser. The topmost entry is the
app that was used to generate this file, or the app that should be used
to view it. The next four entries are submenus for each of the four different
types of apps that can be used on the file. Further down are options for
viewing and re-using the pipeline used to generate the file. The final
section allows you to manage file locations and names. For folders,
left-clicking opens the folder, while right-clicking opens the menu. You
can open file menus whenever you see a file name in link colors when
using the platform. The **Add to** option allows you to copy files while the
**Move to** option removes the original.

**Show all parent containers** gives you the option to quickly find all
copies of a file that are available to you. The **file accession** is a
unique identifier which allows you to find a file even when the file
name has changed.

|ParentContainers|

Above the file system navigator you can find the **Import button**. Clicking
it takes you to the Import app page, where you can upload your files,
import them into the platform and edit their metainfo. 

|import|

Next to the Import button, you can see a **New Folder button**. Using it
you will be able to create a new folder wherever you want. Another option
- **New folder with selection** - appears when you have selected files and
want to put all of them in a separate folder.

|NewFolder|

The **Preprocess, Analyse, Explore and Manage** menus at the top of the page
correspond to the four main actions you can undertake with your data.
These menus will become available when you select a file. 

|MatchingApps|

These apps are “clever” –  when you choose a file, the system will
suggest  apps which can work with the specific file type (e.g.
sequencing assay). However, you still need to think about the nature of
the data. For instance, if you want to align a raw WGBS sequencing assay
Genestack will suggest several mappers, but only the Bisulfite
Sequencing Mapping app will be suitable in this case. To figure out what
apps are recommended to process WGBS, WES, RNA-seq or other sequencing
data, go to the “Pipelines and applications” section of this guide.

**File search** in the top right corner allows you to search for files using
their metadata (names, organism, method). To limit the search by file
type or whether or not the file is shared with you, click on the little
triangle inside the search box.

|FileSearch|

Below the search box is a button to access your **briefcase**. Your
briefcase is a place where you can temporarily store files from various
folders. **To add files** to your briefcase hover over each
individual file and use the special “briefcase” button that appears or
select several files, right click on them and choose “Add to
briefcase...”. **To delete an item** from your briefcase hover over it and
click on the “x” button. **To clear all** items from the briefcase, select
“Clear all” option.

|BriefCase|

If you select a file, **three additional buttons** will show up, allowing
you to **share**, **delete** the file or **view metainfo** (an “eye”-icon) for the
file.

|3buttons1|

|3buttons2|

Use the **Share button** to share your
data with colleagues (the share button will not be available if you
are using a guest account).
Read more about sharing on Genestack in the "Data and Metainfo Management" part of the Guide.

|share|

The **Delete button** allows you to remove your files from the
system.

|delete|

**View metainfo** gives you more information about the file: technical (file
type, its accession and owner, when the file was created and modified,
etc.), biological (e.g. cell line, cell type, organism, etc.), and file
permissions.

|eye|

.. _Getting Started With Genestack Platform: https://genestack.com/tutorial/getting-started-with-genestack-platform/
.. _Testing Differential Gene Expression: https://genestack.com/tutorial/testing-differential-gene-expression-on-genestack-platform/
.. _Whole Genome Bisulfite Sequencing Analysis: https://genestack.com/tutorial/whole-genome-bisulfite-sequencing-analysis/
.. _Whole Exome Sequencing Analysis: https://genestack.com/tutorial/whole-exome-sequencing-data-analysis-on-genestack-platform/
.. _Whole Genome Sequencing Analysis: https://genestack.com/tutorial/wgs-analysis-on-genestack/
.. |SaveMatchingAssays| image:: images/save-matching-assays.png
.. |FileManager| image:: images/file-manager.png
.. |PublicData| image:: images/public-data.png
.. |ParentContainers| image:: images/parent-containers.png
.. |import| image:: images/import_start.png
.. |MatchingApps| image:: images/matching-apps.png
.. |FileSearch| image:: images/file-search.png
.. |BriefCase| image:: images/brief-case.png
.. |3buttons1| image:: images/3buttons-1.png
.. |3buttons2| image:: images/3buttons-2.png
.. |share| image:: images/share.png
.. |delete| image:: images/delete.png
.. |eye| image:: images/eye.png
.. |ExperimentBrowser| image:: images/experiment-browser.png
.. |NewFolder| image:: images/new-folder.png
.. |AnalysisResults| image:: images/analysis-results.png
.. |EditMetainfo| image:: images/DB-to-EditMetainfo.png
.. |fromMEtoFB| image:: images/From-ME-to-FB.png
.. |NewDF| image:: images/new-df.png