Importing data into the Genestack Platform
******************************************

.. raw:: html

    <iframe width="640" height="360" src="https://www.youtube.com/embed/eOl1uabctzI" frameborder="0" allowfullscreen="1">&nbsp;</iframe>

We have talked about the core concepts of Genestack and the geography of
the Platform. Now let's discuss importing data into the platform. On the `dashboard`_ you
can find an **Import data** option, and an **Import** button can be found in the File Manager.

|import file manager|

Once you click it, this will take you to the `Import Data`_ app
page. There are various options of importing your data. You can drag and
drop or select files from your computer, import data from URL or use
previous uploads.

|Import_step1|

After data is uploaded and imported, the platform automatically recognizes file
formats and transforms them into biological data types such as raw reads,
mapped reads, reference genomes and so on. This means you will not have to
worry about formats at all and this will most likely save you a lot of
time. If files are unrecognized, you can manually allocate them to a
specific data type by drag and drop.

|Import_step2|

On the next “Edit metainfo” step, you can describe uploaded data. Using an Excel-like spreadsheet you can
edit the file metainfo and add new attributes, for example cell type or
age.

|Import_step3|

Once this step is completed,
you can go to **Show files in File Manager** at the bottom of the page.
Take a look at a “kind” column ― there are no file formats, just
biological data types.

|File in file manager|

Additional option of importing your data is using import templates. On
the Dashboard you can find an `Add import template`_
option. Import templates allow you to specify required and optional
metainfo attributes for different file kinds. When you scroll down to
the bottom of the page, you will see an **Add import template** button.

|import welcome page|

.. |import file manager| image:: images/import-file-manager1.png
.. |Import_step1| image:: images/Import_step1.png
.. |Import_step2| image:: images/Import_step2.png
.. |Import_step3| image:: images/Import_step3.png
.. |File in file manager| image:: images/files_in_FM.png
.. |import welcome page| image:: images/import-welcome-page1.png
.. _dashboard: https://platform.genestack.org/endpoint/application/run/genestack/welcome
.. _Import Data: https://platform.genestack.org/endpoint/application/run/genestack/uploader
.. _Add import template: https://platform.genestack.org/endpoint/application/run/genestack/metainfotemplateeditorapp?action=openInBrowser
