Task manager
------------

In the top-right corner of any page on Genestack, you can see a link called
**Tasks**. It will take you to the Task Manager, an application which allows you to
track the progress of your computations.

.. image:: images/task-manager.png

Besides, tasks can be *sorted* and *filtered* by application name, file name, accession, status,
a user who started tasks, last update and elapsed time.

Statuses in the Task Manager help you keep track of your tasks. Let’s look what
each status means:

-  *Starting* — the computation process has started to run;
-  *Done* — the task has finished successfully;
-  *Failed* — the computation has failed. To find out why click on **View logs**;
-  *Queued* — the task is pending for execution, for example when it is waiting for
   dependencies to complete initialization;
-  *Running* — your task is in progress;
-  *Blocked by dependency failure* — the computation cannot be completed
   because a task on which this one depends has failed;
-  *Killed* — the task has been canceled by the user.

You can also view *output and error logs* produced for each task.
Error logs tell you why a task has failed. Output logs contain information about
the exact details of what Genestack does with your files during the computation process, 
what specific tools and parameters are used, and so on. If the computations finished
successfully, error logs will be empty, but the logs can provide you with some
basic information about the output data.

.. image:: images/task-log.png

If you change your mind about a computation after it has started, remember that
you can *kill tasks* whenever you want by clicking the **Cancel** button, next to
the task status. To rerun an analysis click file name and select **Restart initialization**.



