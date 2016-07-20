**GO-based enrichment analysis**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To further characterise the biological processes that might be affected
in trisomic samples, we performed `the downstream gene ontology (GO)
analysis <http://geneontology.org/>`__ of the DGE genes. For this, we'll
use GO Enrichment Analysis app, which performs the classic `Fisher's
exact test <https://en.wikipedia.org/wiki/Fisher%27s_exact_test>`__
based on gene counts, against GO annotations.  Open the app on one of
the completed Differential Expression Statistics files:

Changing the group and thresholds criteria in the Filter Options, you
can set what DE genes can be further used for enrichment analysis. Let's
run GO app twice, analysing down- and up-regulated genes separately.
|DGE\_GO\_filters| Analysing down-regulated genes we
observed significant enrichment for genes involved in multiple
developmental processes specifically in organ development and
morphogenesis, embryonic development and morphogenesis and system
development. Also, we found genes associated with nervous system-related
terms and specifically with nervous system development, brain
development, neurogenesis, generation of neurons, neuron differentiation
and axonogenesis. Moreover, there are some terms linked to cellular
adhesion (i.e. biological adhesion, cell adhesion, cell-cell adhesion)
and to the cadherin signalling pathway. Here is the first 8 of 1017 GO
terms related to the biological processes:
|DGE\_down\_DGE\_genes\_GO\_terms| Additionally, you can also
review molecular functions or cellular components that could be affected
in Twin-DS-iPSCs. Just change the "Top GO Terms" from Biological Process
category to the corresponding one. By comparison, this is the table from
the paper which listed the first 20 biological processes that might be
affected due to trisomy 21: |DGE\_down\_paper| Our GO analysis of the
genes up-regulated in trisomic samples revealed enrichment for functions
related to different metabolic and biological processes, regulation of
transcription and DNA-dependent transcription. Here is the list of the
first 8 of 339 biological process GO terms:
|DGE\_up\_regulated\_genes\_GO\_terms| Look at the GO terms associated
with up-regulated genes and reported in the paper: |DGE\_up\_paper| All
these biological processes can be found in our results. The difference
is in GO counts. But we expected it, because the ontologies are not
complete, they are being expanded constantly during the association of
gene products from the collaborating databases. If you'd like to check
it out, open differential expression statistics files stored in
folder \ `GO enrichment analysis for Hibaoui et al
(2013) <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF967843&action=viewFile>`__\ .
This is all for the tutorial. Why don’t you try repeating these steps
with your own data or using our
`public\ ** **\ experiments <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF070886&action=viewFile>`__\ ?
You can try it right now! Just open `the tutorial data
flow <https://platform.genestack.org/endpoint/application/run/genestack/dataflowrunner?a=GSF968015&action=createFromSources>`__ or
create your own one by adding new steps, changing sources and default
options. If you have any questions and comments, please submit them
below or email us at support@genestack.com\ .  Follow us on
Twitter: \ `@genestack <https://twitter.com/genestack>`__.

.. |DGE\_GO\_filters| image:: https://genestack.com/wp-content/uploads/2015/09/DGE_GO_filters-e1445441681370.png
   :class: aligncenter wp-image-3232 size-full
   :width: 600px
   :height: 434px
   :target: https://genestack.com/wp-content/uploads/2015/09/DGE_GO_filters.png
.. |DGE\_down\_DGE\_genes\_GO\_terms| image:: https://genestack.com/wp-content/uploads/2015/09/DGE_down_DGE_genes_GO_terms-e1445441710305.png
   :class: aligncenter wp-image-3233 size-full
   :width: 600px
   :height: 421px
   :target: https://genestack.com/wp-content/uploads/2015/09/DGE_down_DGE_genes_GO_terms.png
.. |DGE\_down\_paper| image:: https://genestack.com/wp-content/uploads/2015/08/DGE_down_paper.png
   :class: aligncenter size-full wp-image-2953
   :width: 613px
   :height: 403px
   :target: https://genestack.com/wp-content/uploads/2015/08/DGE_down_paper.png
.. |DGE\_up\_regulated\_genes\_GO\_terms| image:: https://genestack.com/wp-content/uploads/2015/09/DGE_up_regulated_genes_GO_terms-e1445441756258.png
   :class: aligncenter wp-image-3234 size-full
   :width: 600px
   :height: 418px
   :target: https://genestack.com/wp-content/uploads/2015/09/DGE_up_regulated_genes_GO_terms.png
.. |DGE\_up\_paper| image:: https://genestack.com/wp-content/uploads/2015/08/DGE_up_paper.png
   :class: aligncenter size-full wp-image-2955
   :width: 616px
   :height: 198px
   :target: https://genestack.com/wp-content/uploads/2015/08/DGE_up_paper.png
