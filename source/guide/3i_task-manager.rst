Task Manager
------------

In the upper right corner you can see a link called Tasks. It will take
you to the Task Manager, an application which allows you to track the
progress of your computations.

<screenshot “Task Manager”>|image26|

All your tasks can be sorted and filtered by file name, accession,
status, owner, last updated and elapsed time columns. Also you can ‘view
logs’ for each computation: the error log and the output log. Error logs
tell you why your task has failed. Output logs contain information about
the exact details of what Genestack does with your files during the
computation process, what specific tools and parameters are used, and so
on.

If the computations finished successfully, error logs will be empty, but
the logs can provide you with some basic statistics about the output
data, e.g. mapping statistics from the Unspliced Mapping with Bowtie2
app.

<screenshot “Output log in TM”>|image27|

If you change your mind about a computation after it has started
 remember that you can kill tasks whenever you want by clicking the
“Cancel” button.

Statuses in Task Manager help you keep track of your tasks. Let’s look
what each status means:

-  Created: a request for the computation has been  created and the task
   will be started soon;
-  Starting: the computation process has  started to run;
-  Done:  the task has finished successfully ;
-  Failed: the computation has  failed. To find out  why , click on
   “View logs”;
-  Queued: the task is waiting for dependencies to complete
   initialization or for computing resources to become available;
-  Running: your task is in progress;
-  Queueing: to fill
-  Blocked by dependency failure: the computation cannot be completed
   because a dependency has encountered an error
-  Killed: the task has been canceled .