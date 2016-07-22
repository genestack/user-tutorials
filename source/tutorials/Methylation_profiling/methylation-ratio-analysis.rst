Methylation ratio analysis
**************************

When sequencing reads for both murine phenotypes are mapped onto
reference genome with high enough quality, we can move on to the very
last step in our pipeline - determining * * DNA * * methylation
status of every single cytosine position on the both strands. In order
to do that, go back to the Data Flow Runner page. Click on the
"Methylation Ratio Analysis" to go to the app page where you can see
source files and command line options that could be easily changed. |DF
methylation analysis| Then return to the data flow click on “action”,
“add files”, chose the remaining merged mapped reads, and start
initialization. During this step we apply several options to remove
technical biases in WGBS data:

#. ****Trim N end-repairing fill-in bases set to “3”. **** This option
   allows to trim 3 bases from the read end to remove the DNA overhangs
   created during read end-repair in library preparation. It is
   important because this end repair procedure may introduce artefacts
   if the repaired bases contain methylated cytosines.
#. Report **only unique mappings**
#. ****Discard duplicated reads option** to remove duplicated reads
   which have identical sequences and could be the result of library
   preparation. These reads could be mapped to the same position and
   distort results of downstream analysis.**

The folder `Methylation Ratios for Rodriguez et al.,
2014 <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF968759&action=viewFile>`__
contains all the resulting files of methylation ratios estimation.

.. |DF methylation analysis| image:: https://genestack.com/wp-content/uploads/2015/12/DF-methylation-analysis.png
