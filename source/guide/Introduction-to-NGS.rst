Introduction to Next Generation Sequencing Data Analysis
========================================================

NGS technologies, such as WGS, RNA-Seq, WES, WGBS, ChIP-Seq etc,  generate lots of output data.
Before we start talking about various apps available on Genestack and how to choose appropriate
ones for your analysis, let’s take a moment to go through the basics of sequencing analysis. To
help you better understand the processes involved, we will use the example of genetic variant
analysis for Whole Exome Sequencing data. A typical WES data analysis pipeline includes: raw
reads quality control, preprocessing, mapping, post-alignment processing, variant calling,
annotation, and prioritization  (`Bao et al., 2010`_).

The first thing you need to do with sequencing data is quality control of raw reads.  Here you
should check metrics such as GC content of sequences, sequence duplication levels,
overrepresentation of sequences, and so on. After examining the QC report, you need to decide
whether any preprocessing is needed and if so, what type? For instance, if your sequencing data
is contaminated due to the sequencing process, you may choose to trim adaptors and contaminants
from your data. Or if your data shows high sequence duplication levels, you may want to discard
duplicated reads. Quality control and preprocessing are essential steps, because if you don’t
make sure your data is of good quality to begin with, you cannot fully rely on analysis  results.


After you’ve checked the quality of your data and if necessary, preprocessed it, the next step
is mapping, also called aligning, of your reads to a reference genome or reference transcriptome.
It allows to determine the nucleotide sequence of data being studied with no need of de novo
assembly because obtained reads are compared with a reference already existed in a data base.
For example, in our case, aligning WES reads allows you to discover nucleotides that vary
between reference and tested sequence and the accuracy of variant identification depends
on the mapping accuracy (`The 1000 Genomes Project Consortium, 2010`_).

After you have mapped your reads, it’s a good idea to check the mapping quality, as
some of the biases in the data only show up after the mapping step. Similarly to what you’ve
done before with raw sequencing reads, if you are unsatisfied with the mapping quality, you can
process the mapped reads and, for instance, remove duplicated mapped reads (which could be PCR
artifacts). Post-alignment processing is very important, as it can greatly  improve the accuracy
and quality of further variant analysis.
Once the sequence is aligned to a reference genome, the data needs to be analyzed in an experiment-
specific fashion. Here we will use the WES reads mapped against the reference genome to perform
variant analysis, including variant calling and predicting the effects  found variants produce
on known genes (e.g. amino acid changes or  frame shifts). In this step you compare your sequence
with the reference sequence, look at all the differences and try to establish how big of an
influence do these changes have on the gene. For instance, if it’s a synonymous variant, it will
probably have low influence on the gene as such a change causes a codon that produces the same
amino acid. However, if it’s a large deletion, you can assume that it will have a large effect
on the gene function.

When it comes to visualising your data: the standard tool for viewing of mapped reads and
identified variants is a Genome Browser. Since visualization is one of the concepts at the core
of our platform, on Genestack you will find a range of other useful tools that will help you
better understand your data considering their nature. For example for WES or WGS data we suggest
to use Variant Explorer which can be used to sieve through thousands of variants and  allow users
to focus on their most important findings.

In summary, the basic steps of bioinformatics analysis are: data preparation, mapping, difference
analysis, and result visualisation.

On Genestack we have divided the available applications into four categories based on what they do:
Preprocess, Analyse, Explore, and Manage. Apps in the preprocess section are used to prepare raw or
mapped data for future computations (by improving their quality) or to uncover any poor quality data
 or contaminants. Apps in the Analyse section are computational apps used to perform actual data
 analysis. Apps in the Explore section are interactive, visual applications allowing you to view
 the results of your analyses, identify the most interesting findings and interpret your results.
 The manage section apps allow you to keep track of data, create a template for your pipeline and
 explore and edit the metadata of your files.


.. _`Bao et al., 2010`: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4179624/&sa=D&ust=1480960531812000&usg=AFQjCNEanbRs0Pes_OFgveaUiLQ59pPQfQ
.. _`The 1000 Genomes Project Consortium, 2010`: http://www.nature.com/nature/journal/v467/n7319/full/nature09534.html&sa=D&ust=1480960531815000&usg=AFQjCNFuLcgFRjG08dkh7vuZnnQPJPmbsQ