**Quality control of raw reads** Raw sequencing reads can include PCR
primers, adaptors, low quality bases, duplicates and other contaminants
coming from the experimental protocols. That's why we recommend you to
check the quality of your raw data looking at such aspects as GC
percentage, per base sequence quality scores, and other quality
statistics. The easiest way to do this is to run \ `the Raw Reads QC
data
flow <https://platform.genestack.org/endpoint/application/run/genestack/dataflowrunner?a=GSF969011&action=createFromSources>`__:

After generating reports, you'll be able to review various statistics
and plots for each sample. Here is some of them for \ `Twin-N-iPSC
Rep3 <https://platform.genestack.org/endpoint/application/run/genestack/fastqc-report?a=GSF968986&action=viewFile#>`__:
**Per sequence GC content.** In a random library you expect a roughly
normal GC content distribution. An unusually shaped or shifted
distribution could indicate a contamination. |DGE\_per\_seq\_GC\_cont|
**Per base sequence quality** plot depicts the range of quality scores
for each position in the reads. A good sample will have qualities all
above 28: |DGE\_per\_base\_s| **Per sequence quality scores**
plot\ ** **\ shows the frequencies of quality scores in a sample. If the
reads are of good quality, the peak on the plot should be shifted to the
right as far as possible. In our case, the majority of Illumina reads
have a good quality - in the range from 35 to 40 (the best score is 41):
|DGE\_per\_seq\_quality\_sco| **Per base sequence content** shows the
four nucleotides’ proportions for each position. In a random library you
expect no nucleotide bias and the lines should be almost parallel with
each other: |DGE\_per\_base\_seq\_cont| There is a bias at the beginning
of the reads, which is common for RNA-Seq data. This occurs during
RNA-seq library preparation, when “random” primers are annealed to the
start of sequences. These primers are not truly random, and it leads to
a variation at the  beginning of the reads. **Sequence duplication
levels** plots represent the proportion of the library made up of
sequences with different duplication levels. Sequences with 1, 2, 3, 4,
etc duplicates are grouped to give the overall duplication level.
|DGE\_seq\_dupl\_lev| Looking at these plots, you may notice 15% of
sequences duplicated more than 10 times, 6% of sequences are repeated
more than 100 times, etc. The overall rate of  duplication is about 40
%. Nevertheless, while analysing transcriptome sequencing data, we
should not remove these duplicates because we don't know whether they
represent PCR duplicates or high gene expression of our samples. We have
run QC on all the data in the experiment and collected reports in `Raw
reads QC reports for Hibaoui et al
(2013) <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF000383&action=viewFile>`__ folder.

**3. Preprocessing of raw reads**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the quality of raw data has been checked, we can go back to the
main `Differential Gene Expression
Analysis <https://platform.genestack.org/endpoint/application/run/genestack/dataflowrunner?a=GSF968176&action=createFromSources>`__ data
flow and choose sources: |DGE\_data\_flow\_first\_step|

You can upload your samples directly into the data flow or select them
from the available datasets. Let's choose sources - 7 sequencing assays
from the experiment and a human reference genome, and create resulting
files in a specific folder.

QC reports can not only provide you with the information on the data
quality but can also help you to decide how to preprocess the data in
order to improve its quality and get more reliable results in further
analysis. There are various Genestack applications that allow you to do
preprocessing:\ |DGE\_preprocess\_apps| We'll run Trim Adaptors and
Contaminants app: |DGE\_trim\_ad| All resulting files are collected
in \ `Trimmed raw reads for Hibaoui et al
(2013) <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF967714&action=viewFile>`__ folder.

.. |DGE\_per\_seq\_GC\_cont| image:: https://genestack.com/wp-content/uploads/2015/08/DGE_per_seq_GC_cont.png
   :class: aligncenter size-full wp-image-2945
   :width: 539px
   :height: 349px
   :target: https://genestack.com/wp-content/uploads/2015/08/DGE_per_seq_GC_cont.png
.. |DGE\_per\_base\_s| image:: https://genestack.com/wp-content/uploads/2015/08/DGE_per_base_s.png
   :class: aligncenter size-full wp-image-2946
   :width: 726px
   :height: 277px
   :target: https://genestack.com/wp-content/uploads/2015/08/DGE_per_base_s.png
.. |DGE\_per\_seq\_quality\_sco| image:: https://genestack.com/wp-content/uploads/2015/08/DGE_per_seq_quality_sco.png
   :class: aligncenter size-full wp-image-2947
   :width: 545px
   :height: 355px
   :target: https://genestack.com/wp-content/uploads/2015/08/DGE_per_seq_quality_sco.png
.. |DGE\_per\_base\_seq\_cont| image:: https://genestack.com/wp-content/uploads/2015/08/DGE_per_base_seq_cont.png
   :class: aligncenter size-full wp-image-2948
   :width: 727px
   :height: 279px
   :target: https://genestack.com/wp-content/uploads/2015/08/DGE_per_base_seq_cont.png
.. |DGE\_seq\_dupl\_lev| image:: https://genestack.com/wp-content/uploads/2015/08/DGE_seq_dupl_lev.png
   :class: aligncenter size-full wp-image-2949
   :width: 736px
   :height: 413px
   :target: https://genestack.com/wp-content/uploads/2015/08/DGE_seq_dupl_lev.png
.. |DGE\_data\_flow\_first\_step| image:: https://genestack.com/wp-content/uploads/2015/07/DGE_data_flow_first_step.png
   :class: aligncenter size-full wp-image-2828
   :width: 405px
   :height: 575px
   :target: https://genestack.com/wp-content/uploads/2015/07/DGE_data_flow_first_step.png
.. |DGE\_preprocess\_apps| image:: https://genestack.com/wp-content/uploads/2015/07/DGE_preprocess_apps.png
   :class: aligncenter size-full wp-image-2860
   :width: 749px
   :height: 376px
.. |DGE\_trim\_ad| image:: https://genestack.com/wp-content/uploads/2015/07/DGE_trim_ad.png
   :class: aligncenter size-full wp-image-2895
   :width: 410px
   :height: 616px
