Building pipelines
------------------

Bioinformatic data analysis includes several steps, these steps vary depending on the type of
the data and your goals. For instance, WGS data analysis includes the following steps: check
the initial quality of raw reads, preprocessing of the data to improve the quality,
if it is needed, and alignment the reads onto a reference genome followed by identification
and annotation of genetic variants.

With Genestack you can either use one of the dataflows or build a pipeline
manually selecting customizable applications supported by the system.

Use the **Data Browser** to find a dataset you would like to analyse, click on it.
Then, on the **Metainfo Editor** page click on the button marked **Analyse** to
start creating a pipeline. If you want to analyse not the entire dataset but
some part of it, select the assays you wish to analyse and **Make a subset**.

So, select the first application you wish to see in your pipeline. For each
individual file the system suggests only applications that can be used to
analyse your data, considering its type and metadata.

Applications on the platform are divided in several categories:

- *Preprocess* to prepare the data for actual analysis;
- *Analyse* perform various kinds of analysis;
- *Explore* to visualise QC check or analysis results;
- *Manage* to operate with your files.

.. image:: images/pipeline_building.png
   :align: center
   :scale: 80 %

This will take you to the application page where you can:

- learn more about the application;
- view and edit application parameters;
- explore your results;
- add further steps to the file data flow (the pipeline).

.. image:: images/cla_page.png
   :align: center
   :scale: 80 %

To proceed click on **Add step** button that will show you the list of all the
matching applications.

.. image:: images/cla-add-step.png
   :align: center
   :scale: 80 %

Continue adding steps until you have finished building your pipeline. When
you add each of the steps, you create new files which end up in the **Created files**
and **My datasets** folders. However, these files are not yet ready to use —
they need to be initialized first.
