Calculate read coverage for genes
*********************************

After mapping, we can count reads mapped to annotated features (genes,
exons, etc) running ** ** Quantify Raw Coverage in Genes app:
|DGE_quantify_genes| To run the app, click on "7 files" and then
"Start initialization". For our analysis we counted reads mapping within
exons, grouping them by gene_id and assigning reads to all exons they
overlap with. We calculated read coverage in all samples and collected
resulting files in `Raw gene counts for Hibaoui et al
(2013) <https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF967836&action=viewFile>`__ folder.

.. |DGE_quantify_genes| image:: images/DGE_quantify_genes.png
