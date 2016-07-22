Exploring the genome methylation levels in Genome Browser
*********************************************************

We can explore the distribution of genome methylation levels counted for
both murine phenotypes in  `Genome
Browser <https://platform.genestack.org/endpoint/application/run/genestack/genomeBrowser?a=GSF969175&action=viewFile>`__ .

As it was mentioned before, "Canyons" are the large unmethylated DNA
regions inside of a highly methylated locus that often harbour genes
playing a role in hematopoiesis and being deregulated in human leukemia.
Some Canyons can be exceptionally large, for example one associated with
the  *Pax6 * homeobox gene encoding a homeobox-containing protein
regulating  transcription is extended over 25 kb: |Genome Browser (Pax6;
only WTs)| Let’s compare our methylation ratios distribution in these
region with author’s results:

|Pax6-paper|

As DNMT3a is often mutated in
human leukemias, authors also examined the impact of loss of Dnmt3a on
the Canyon size. For this they compared low-methylated regions in HSCs
conditionally inactivated for Dnmt3a to WT HSCs: |Genome Browser (WT vs
KO)| This investigation revealed the fact that methylation loss in
Dnmt3a KO HSCs leads to the formation of new Canyons. Lack of Dnmt3a
does not affect regions inside Canyons but it results in changes of
Canyon boundaries: Canyon size can be decreased due to hypermethylation
or increased due to of hypomethylation. Moreover, at DNA regions
containing cluster of Canyons in WT HSCs,  larger Canyons (“Grand
Canyons”)  can be formed. We can see it on the example of  *HoxB*
regions in which Canyons are interrupted by short stretches of higher
methylation. All these findings suggest that Dnmt3a can be crucial for
maintaining methylation in the Canyon boundaries. **Now, let's take a
look at the original track for the same Canyon cluster to compare the
results:**

|Hox-paper2|

This experiment is a part of the large research
of changes in DNA methylation profile including different methodologies
such as, for example, whole genome bisulfite sequencing and CMS-seq to
reveal genome-wide distribution of mCs and hmCs, RNA-Seq to analyse
expression of Canyon-associated genes. This incredible work was turned
into a
research  `paper <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3920905/>`__,  and
the data sets can be found in our `Public
Experiments <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF070886&action=viewFile&page=1>`__ !
That's it for the tutorial, we hope you will enjoy working on your data
with Genestack! Later you can return back to the tutorial if necessary.
If you have any questions, suggestions etc, please leave them in
comments below or `email us <mailto:info@genestack.com>`__.

.. |Genome Browser (Pax6; only WTs)| image:: https://genestack.com/wp-content/uploads/2015/08/GB-Pax6-only-WTs.png
   :class: aligncenter wp-image-2989
   :width: 650px
   :height: 428px
   :target: https://genestack.com/wp-content/uploads/2015/08/GB-Pax6-only-WTs.png
.. |Pax6-paper| image:: https://genestack.com/wp-content/uploads/2015/08/Pax6-paper.png
   :class: size-full wp-image-2992 aligncenter
   :width: 517px
   :height: 288px
.. |Genome Browser (WT vs KO)| image:: https://genestack.com/wp-content/uploads/2015/08/GB-WT-vs-KO.png
   :class: aligncenter wp-image-2991
   :width: 650px
   :height: 344px
   :target: https://genestack.com/wp-content/uploads/2015/08/GB-WT-vs-KO.png
.. |Hox-paper2| image:: https://genestack.com/wp-content/uploads/2015/08/Hox-paper2.png
   :class: size-full wp-image-2990 aligncenter
   :width: 652px
   :height: 341px
