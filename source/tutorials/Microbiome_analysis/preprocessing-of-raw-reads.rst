Preprocessing of raw reads
**************************

FastQC reports help you understand whether it is necessary to improve the
quality of your data by applying trimming, filtering, adaptor clipping and
other preprocessing steps. Here is the list of Genestack preprocess
applications available for raw reads: 

|Microbiome_preprocess_apps|

Once we have checked the quality of the raw reads, let's start building the
Microbiome Data Analysis pipeline:

.. Video - Building Microbiome Analysis pipeline
.. raw:: html

    <iframe width="640" height="360" src="" frameborder="0" allowfullscreen="1">&nbsp;</iframe>

Our preprocessing procedure will include two steps - adaptor trimming and
filtering out low quality reads. For this, choose all 8 samples in the
experiment and select "Trim Adaptors and Contaminants" application in
Preprocess section:

|Microbiome_trim_adaptors_and_contaminants|

This brings you to the appplication page

Trimmed reads are stored in `Trimmed raw reads for Alfano et al (2015)`_
folder.



Variant calling
***************

After duplicate removal, the next step is to identify different genomic
variants including SNVs, indels, MNVs, etc. For this, we'll use Variant
Calling application based on samtools mpileup:

.. raw:: html

    <iframe width="640" height="360" src="" frameborder="0" allowfullscreen="1">&nbsp;</iframe>

The app automatically scans every position along the genome, computes all the
possible genotypes from the aligned reads, and calculates the probability
that each of these genotypes is truly present in your sample. Then genotype
likelihoods are used to call the SNVs and indels.

|WES_variant_calling|

We run Variant Calling with default parameters, identifying multi-allelic
SNPs and indels, excluding non-variant sites and not considering anomalous
read pairs. Maximum read depth per position was set as 250 and minimum number
of gapped reads for an indel candidate is 1. Base alignment quality (BAQ)
recalculation is turned on by default. It helps to rule out false positive
SNP calls due to alignment artefacts near small indels. For more information
about the app and its options, click on the app name and then on 'About
application'.

When files will be complete, you can analyse `variants in Genome Browser`_:

|WES_variants_GB|

Genome Browser application allows you investigate the variants interactively:
how many mutations are in this particular gene or region, review some
information about detected variants such as average mapping quality and raw
read depth and compare variants enrichment between samples. Analysing variants
in Genome Browser, you can notice a large amount of both exome WES–specific and
WGS-specific SNVs. We identified variants for each sample separately and put
them in `Variants for Clark et al (2011)`_ folder.

.. |Microbiome_preprocess_apps| image:: images/Microbiome_preprocess_apps.png
.. |Microbiome_trim_adaptors_and_contaminants| image:: images/Microbiome_trim_adaptors_and_contaminants.png

.. |WES_variant_calling| image:: images/WES_variant_calling.png
.. |WES_variants_GB| image:: images/WES_variants_GB.png
.. _Filtered mapped reads for Clark et al (2011): https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF999208&action=viewFile&page=1
.. _variants in Genome Browser: https://platform.genestack.org/endpoint/application/run/genestack/genomeBrowser?a=GSF999281&action=viewFile
.. _Variants for Clark et al (2011): https://platform.genestack.org/endpoint/application/run/genestack/filebrowser?a=GSF999229&action=viewFile&page=1.. _
