Mapping RNA-seq reads onto a reference genome
*********************************************

When all files were created, you can run the whole analysis here,
choosing Expression Navigator for genes. But first, let’s align RNA-seq
reads to the reference genome across splicing junctions and then compare
mappings in Genome Browser.\ |DGE\_spl\_mapping| Find Spliced Mapping
step, click on “7 files”. In “Explore” section choose “Genome Browser”
and start initialization there.

We run Spliced Mapping app with default parameters. To change them go to
the app page and choose "Edit parameters" button. If you want to learn
more about the app and its options, click on the app name and then on
"About application". |DGE\_spl\_map| Completed Mapped Reads files can be
found in `Mapped reads files for Hibaoui et al
(2013) <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF967837&action=viewFile>`__ folder. Let's
open some of them in Genome Browser to analyse reads \ `coverage on
chromosome
21 <https://platform.genestack.org/endpoint/application/run/genestack/genomeBrowser?a=GSF968535&action=viewFile&expired>`__\ on
the region chr21:30007376-40007694 (10 Mb): |DGE\_coverage\_21| Here
is a combined track for all trisomic and control samples:
|DGE\_GB\_combined\_track| As you see, the majority of chr21 genes are
indeed more expressed in the trisomic samples than in the euploid ones,
which is consistent with the overall up-regulation of chr21 genes in
individuals with Down syndrome.

**Quality control of mapped reads**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optional step is to check how mapping went using Mapped Reads QC
Report app. You can "generate reports" for each mapping separately or
just run `Mapped Reads Quality
Control <https://platform.genestack.org/endpoint/application/run/genestack/dataflowrunner?a=GSF968216&action=createFromSources>`__\ data
flow for multiple samples:

Output report includes mapping statistics such as:

#. **Mapped reads**: total reads which mapped to the reference genome;
#. **Unmapped reads**: total reads which failed to map to the reference
   genome;
#. **Mapped reads with mapped mate**: total paired reads where both
   mates were mapped;
#. **Mapped reads with partially mapped mate**: total paired reads where
   only one mate was mapped;
#. **Mapped reads with “properly” mapped mate**: total paired reads
   where both mates were mapped with the expected orientation;
#. **Mapped reads with “improperly” mapped mate**: total paired reads
   where one of the mates was mapped with unexpected orientation.

The **Coverage by chromosome** plot shows the percentage of bases
covered (y-axis) by at least N (x-axis)
reads.\ |Coverage\_by\_chromosome| For paired reads, you can look
at insert size statistics, such as median and mean insert sizes, median
absolute deviation and standard deviation of insert size. The **Insert
size distribution** plot is generated: |Insert\_size\_distribution| We
already prepared all QC reports for mapped reads and put them in `Mapped
reads QC reports for Hibaoui et al
(2013) <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF967840&action=viewFile>`__ folder.
You can open all of them in `Multiple QC Report
app <https://platform.genestack.org/endpoint/application/run/genestack/multiple-qc-plotter?a=GSF968715&action=viewFile>`__ to
view mapping statistics interactively: |DGE\_multiple\_qc\_plotter|
Overall, more than 80 % of reads are mapped. It includes properly and
partially mate pairs. Less than 11 % of reads are unmapped among the
samples. Additionally, you can sort your samples by QC statistics or
metainfo values. Read more what the app does in our blog post about
`i <https://genestack.com/blog/2014/12/10/interactive-sequencing-quality-control-reports/>`__\ `nteractive
sequencing quality control
reports <https://genestack.com/blog/2014/12/10/interactive-sequencing-quality-control-reports/>`__\ .

.. |DGE\_spl\_mapping| image:: https://genestack.com/wp-content/uploads/2015/07/DGE_spl_mapping.png
   :class: aligncenter size-full wp-image-2897
   :width: 401px
   :height: 613px
.. |DGE\_spl\_map| image:: https://genestack.com/wp-content/uploads/2015/08/DGE_spl_map-e1445441938143.png
   :class: aligncenter wp-image-2958 size-full
   :width: 600px
   :height: 729px
   :target: https://genestack.com/wp-content/uploads/2015/08/DGE_spl_map.png
.. |DGE\_coverage\_21| image:: https://genestack.com/wp-content/uploads/2015/07/DGE_coverage_21-e1445441975435.png
   :class: aligncenter wp-image-2899 size-full
   :width: 600px
   :height: 380px
   :target: https://genestack.com/wp-content/uploads/2015/07/DGE_coverage_21.png
.. |DGE\_GB\_combined\_track| image:: https://genestack.com/wp-content/uploads/2015/07/DGE_GB_combined_track-e1445442051712.png
   :class: aligncenter wp-image-2903 size-full
   :width: 600px
   :height: 397px
   :target: https://genestack.com/wp-content/uploads/2015/07/DGE_GB_combined_track.png
.. |Coverage\_by\_chromosome| image:: https://genestack.com/wp-content/uploads/2015/07/Coverage_by_chromosome-e1445442085712.png
   :class: aligncenter wp-image-2764 size-full
   :width: 600px
   :height: 400px
   :target: https://genestack.com/wp-content/uploads/2015/07/Coverage_by_chromosome.png
.. |Insert\_size\_distribution| image:: https://genestack.com/wp-content/uploads/2015/07/Insert_size_distribution-e1445442123895.png
   :class: aligncenter wp-image-2763 size-full
   :width: 600px
   :height: 398px
   :target: https://genestack.com/wp-content/uploads/2015/07/Insert_size_distribution.png
.. |DGE\_multiple\_qc\_plotter| image:: https://genestack.com/wp-content/uploads/2015/09/DGE_multiple_qc_plotter-e1445442157923.png
   :class: aligncenter wp-image-3083 size-full
   :width: 600px
   :height: 377px
   :target: https://genestack.com/wp-content/uploads/2015/09/DGE_multiple_qc_plotter.png
