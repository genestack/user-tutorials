Reference genomes
-----------------

One way or another, many bioinformatics analysis pipelines rely on the use of a reference genome.
For instance,  we use reference genomes in `DNA methylation analysis`_, in `differential gene
expression analysis`_, and in the `analysis of transcriptomic heterogeneity
within populations of cells`_. The choice of a reference genome can increase
the quality and accuracy of the downstream analysis or it can have a harmful
effect on it. For example, it has been shown that the choice of a gene
annotation has a big impact on RNA-seq data analysis, but also on `variant
effect prediction`_.

.. _DNA methylation analysis: http://genestack-user-tutorials.readthedocs.io/tutorials/Methylation_profiling/index.html
.. _differential gene expression analysis: http://genestack-user-tutorials.readthedocs.io/tutorials/DGE_analysis/index.html
.. _analysis of transcriptomic heterogeneity within populations of cells: https://genestack.com/blog/2014/09/24/single-cell-rna-seq-analysis-tutorial/
.. _variant effect prediction: http://genestack-user-tutorials.readthedocs.io/tutorials/WGS_data_analysis/index.html

On Genestack, you can find `several reference genomes`_ for some of the most
common model organisms. We are adding more and more reference genomes of model
organisms to this list regularly.

.. _several reference genomes: https://platform.genestack.org/endpoint/application/run/genestack/signin?original_url=%2Fendpoint%2Fapplication%2Frun%2Fgenestack%2Ffilebrowser%3Fa%3DGSF000018%26action%3DviewFile%26page%3D1
 
For some organisms we provide several genomes, e.g.  there are a couple of
reference genomes for *Homo sapiens*.

.. image:: images/public_reference_genomes.png

What are the differences between these reference
genomes? And how do you chose the correct one?  The answer is not so
straightforward and depends on several factors – let’s discuss each of them:

1. **Reference genome assembly and release version**

For instance: "Homo sapiens / GRCh37 release 75" vs "Homo sapiens / GRCh38
release 86".

The numbers correspond to versions (or “builds”) of the reference genome – the
higher the number, the more recent the version. We generally recommend you use
the latest version possible. One thing to remember is that for the newest
genome builds, it is likely that resources such as genome annotations and
functional information will be limited, as it takes time for Ensembl/ UCSC to
integrate additional genomic data with the new build. You can read more about
it a `blog post`_ from Genome Spot blog and in `this article`_ from Bio-IT.

.. _blog post: http://genomespot.blogspot.ru/2015/06/mapping-ngs-data-which-genome-version.html
.. _this article: http://www.bio-itworld.com/2014/1/27/getting-know-new-reference-genome-assembly.html

2. **One organism – many strains**

K12 and O103 are two different strains of *E.coli*. K12_ is an innocuous strain
commonly used in various labs around the world. O103_ is a pathogenic strain,
commonly isolated from human cases in Europe. Depending on your experiment, you
should choose a matching reference genome.

.. _K12: https://www.genome.wisc.edu/resources/strains.htm
.. _O103: http://aem.asm.org/content/79/23/7502.full

3. **Toplevel sequence or primary assembly**

- **Toplevel** reference genomes contain all chromosomes, sequence regions not
  assembled into chromosomes and padded haplotype/patch regions.

- **Primary assembly** genomes contain all toplevel sequence region excluding
  haplotypes and patches.

We are strongly recommend to use primary assembly reference genomes, since they
are best for performing sequence similarity searches while patches and
haplotypes would confuse analysis.

4. **DNA or cDNA**

- **DNA** - reference genome contains sequence of genomic DNA;
- **cDNA** reference genome consists of all transcripts sequences for actual and
  possible genes, including pseudogenes.

5. **Masked, soft-masked and unmasked genomes**

There are three types of Ensembl reference genomes: unmasked, soft-masked and
masked.

Masking is used to detect and conceal interspersed repeats and low complexity
DNA regions so that they could be processed properly by alignment tools.
Masking can be performed by special tools, like RepeatMasker_. The tool goes
through DNA sequence looking for repeats and low-complexity regions.

.. _RepeatMasker: http://www.repeatmasker.org/

There are two types of masked reference genomes: masked and soft-masked.

- **Masked** reference genomes are also known as hard-masked DNA sequences.
  Repetitive and low complexity DNA regions are detected and replaced with
  ‘N’s. The use of masked genome may adversely affect the analysis
  results, leading to wrong read mapping and incorrect variant calls.


.. note:: **When should you use a masked genome?**

          We generally do not recommend using masked genome, as it relates to the
          loss of information (after mapping, some "unique" sequences may not be
          truly unique) and does not guarantee 100% accuracy and sensitivity (e.g.
          masking cannot be absolutely perfect). Moreover, it can lead to the
          increase in number of falsely mapped reads.


- In **soft-masked** reference genomes, repeats and low complexity regions are
  also detected but in this case they are masked by converting to a lowercase
  variants of the base (e.g. acgt).


.. note:: **When should you use a soft-masked genome?**

          The soft-masked sequence does contain repeats indicated by lowercase
          letters, so the use of soft-masked reference could improve the quality
          of the mapping without detriment to sensitivity. But it should be noted
          that most of the alignment tools do not take into account soft-masked
          regions, for example BWA, tophat, bowtie2 tools always use all bases in
          alignment weather they are in lowercase nucleotides or not. That is why,
          there is no actual benefit from the use of soft masked genome in
          comparison with unmasked one.


- We recommend you use **unmasked** genomes when you do not want to lose any
  information. If you want to perform some sort of filtering, it is better to do
  so  after the mapping step.

Usually, reference genome name includes information about all these factors:
organism, genome assembly, release, primary assembly/toplevel, masking
procedure and molecule.

*Example*:

To perform Whole Exome Sequencing analysis, we recommend you use an unmasked
reference genome of the latest releases and assemblies (e.g. Homo sapiens /
GRCh38 release 85 (primary assembly, unmasked, DNA) for human samples).

The bioinformatics community is divided on the topic of the use of reference
genomes. It is our personal opinion that it is best to always use unmasked
genome and perform filtering after the mapping step. However, if you would like
to read more on the topic, we suggest taking a look at the following papers:

#. McCarthy DJ, Humburg P, Kanapin A, Rivas MA, Gaulton K, Cazier JB, Donnelly P.
   Choice of transcripts and software has a large effect on variant annotation.
   Genome Med. 2014;6(3):26. DOI: 10.1186/gm543;
#. Frankish A, Uszczynska B, Ritchie GR, Gonzalez JM, Pervouchine D, Petryszak R,
   et al. Comparison of GENCODE and RefSeq gene annotation and the impact of
   reference geneset on variant effect prediction. BMC Genomics. 2015;16 (Suppl
   8):S2. DOI: 10.1186/1471-2164-16-S8-S2.
