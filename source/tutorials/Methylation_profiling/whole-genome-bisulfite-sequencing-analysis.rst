This post has videos accompanying it. You can watch them here: |Watch
all videos here | Or read full text and watch the videos separately:
Bisulfite sequencing approaches are currently considered a “gold
standard” for detecting DNA methylation. One of them, whole-genome
bisulfite sequencing (WGBS), provides single-base resolution of
methylated cytosines in genomic DNA. |DNA_methylation| Investigating
the methylation profile of DNA is extremely valuable because this type
of epigenetic modification controls gene expression, and is involved in
such processes as embryonic development, genomic imprinting,
X-chromosome inactivation and cell differentiation. Since methylation
takes part in so many cellular processes, it can be expected
that aberrant methylation can be associated with various diseases for
example different types of cancer. In this tutorial we will show you how
to bioinformatically analyse the whole genome DNA methylation with
Genestack applications.  The entire pipeline includes the following
steps:

#. `Setting up a WGBS experiment <#setting>`__
#. `Quality control of bisulfite sequencing
   reads <https://genestack.com/tutorial/quality-control-preprocessing-raw-reads/>`__
#. `Preprocessing of the raw reads: trimming adaptors, contaminants and
   low quality
   bases <https://genestack.com/tutorial/quality-control-preprocessing-raw-reads/>`__
#. `Bisulfite sequencing mapping of the preprocessed reads onto a
   reference
   genome <https://genestack.com/tutorial/mapping-sequencing-reads-merging-techinical-replicates/>`__
#. `Merging the mapped
   reads <https://genestack.com/tutorial/mapping-sequencing-reads-merging-techinical-replicates/>`__
#. `Quality control of the mapped
   reads <https://genestack.com/tutorial/quality-control-mapped-reads/>`__
#. `Methylation ratio
   analysis <https://genestack.com/tutorial/methylation-ratio-analysis/>`__
#. `Exploring the genome methylation levels in Genome
   Browser <https://genestack.com/tutorial/exploring-methylation-levels-genome-browser/>`__

`The data
flow <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF969172&action=viewFile&page=1>`__ used
in this tutorial, has been previously prepared by our team and put into
 `tutorial
folder <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF970554&action=viewFile&page=1>`__
for your convenience.

**** ****1. Setting up a WGBS experiment**** To go through all these
steps, we need to choose one of the WGBS experiments from Genestack
`Public
Experiments <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF070886&action=viewFile&page=1>`__ which
are the part of our `Public
data <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=public&action=viewFile&page=1>`__.
Feel free to reproduce the workflow on your own data uploaded with the
`Data
Importer <https://platform.genestack.org/endpoint/application/run/genestack/uploader>`__.
 What's more, you could open the folder
`"Tutorials" <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF000810&action=viewFile&page=1>`__  in
the Public data and select the "Whole-genome Bisulfite Sequencing Data
Analysis on Genestack Platform". In the folder you will find the WGBS
experiment used in this tutorial, processed data and all the other
needed files.  **Experiment by  `Rodriguez et al.,
2014 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49714>`__
 can serve as a clear example of the applications of WGBS to
genome-wide DNA 5-methylcytosine landscape profiling.** |public
experiments| To learn more just open the experiment in `Experiment
Viewer <https://platform.genestack.org/endpoint/application/run/genestack/experiment-viewer?a=GSF088374&action=viewFile>`__ :
|Experiment_Viewer| Briefly, the authors performed WGBS on DNA obtained
from mouse hematopoietic stem cells (HSCs) to investigate the mechanisms
that could promote changes in DNA methylation and contribute
to malignant transformation. They discovered extended DNA regions of low
methylation **—** “Canyons”**—** that are distinct from CpG islands and
shores and cover conserved domains frequently containing transcription
factors. Then, as DNA methyltransferase 3a (Dnmt3a) encoding gene is
often mutated in human leukemias, the authors also compared DNA
methylation patterns in purified wild type and Dnmt3a conditional
knockout mouse HSCs. And it was revealed that the loss of DNA Dnmt3a can
influence the Canyon size. Now let's start reproducing these results
with data flows pre-prepared by Genestack.

.. |Watch all videos here | image:: https://genestack.com/wp-content/uploads/2015/10/Zrzut-ekranu-2015-10-21-o-16.01.36-1024x109.png
   :class: aligncenter wp-image-3563 size-large
   :width: 604px
   :height: 64px
   :target: https://www.youtube.com/playlist?list=PLqGSwEO9VFw3ZfhBit9j2sTwTRiLvkJ6T
.. |DNA_methylation| image:: https://genestack.com/wp-content/uploads/2015/09/DNA_methylation-300x225.jpg
   :class: alignright wp-image-3052 size-medium
   :width: 300px
   :height: 225px
   :target: https://genestack.com/wp-content/uploads/2015/09/DNA_methylation.jpg
.. |public experiments| image:: https://genestack.com/wp-content/uploads/2015/12/public-experiments.png
   :class: alignnone wp-image-4182
   :width: 600px
   :height: 167px
   :target: https://genestack.com/wp-content/uploads/2015/12/public-experiments.png
.. |Experiment_Viewer| image:: https://genestack.com/wp-content/uploads/2015/08/Experiment_Viewer.png
   :class: aligncenter wp-image-2971 size-full
   :width: 701px
   :height: 685px
