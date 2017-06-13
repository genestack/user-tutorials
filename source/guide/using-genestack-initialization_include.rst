Initialising files
------------------

You can initialize files in different ways:

1. Using **Start initialization** option in the context menu.

Click on the name of your last created file at the top of the application page
and select "start initialization".

.. image:: images/start_initialization.png

2. Clicking **Start initialization now** in Data Flow Runner application.

If you want to save the pipeline and specific parameters you used here
to re-use again on other files you can create a new data flow. You need
to do this, before you start initialization. Click on the name of the
last created file, go to Manage and **Create new Data Flow**.

.. image:: images/create_new_data_flow.png

This will take you to the Data Flow Editor where you can introduce any last
changes to your pipeline. Click on **Run dataflow** button once you are done.

.. image:: images/data_flow_editor.png

This will take you to **Data Flow Runner** page where you can initialize the
computations (by clicking "Run Data Flow" in the last cell).

.. image:: images/run_data_flow.png

Choose **Start initialization now** option if you would like to run the
computations immediately or **Delay initialization till later**.

.. image:: images/start_initialization_now.png
   :scale: 65 %

This data flow, along with all your results (after computations are finished)
will be stored in the "Created files" folder.

3. Using **File Initializer** application.

Select the created uninitialized files (from data flow or File Manager), right
click on them, go to "Manage" and choose the "File Initializer" application.

.. image:: images/file_initializer_df.png

File Initializer reports the status of the files and allows you to initialize
those that need to be by clicking on their respective "Go!" buttons, or
"Initialize all" to do them all at once. Files do not need to be produced by
the same applications to be initialized together.

.. image:: images/file_initializer.png

4. Using **Start initialization** button in File Provenance.

Alternatively, you can click on the name of the last created file, go to Manage
and choose File Provenance application. The application displays the pipeline
and also allows you to run the computation using "Start initialization" button.
You’ll learn more about this application later in this section.

.. image:: images/file_provenance_init.png

You can track the progress of this task in Task Manager. Doing this will begin
initialization of all the files (including intermediate files) you have
created whilst building this pipeline.
