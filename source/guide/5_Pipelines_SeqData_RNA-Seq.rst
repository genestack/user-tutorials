RNA-seq
~~~~~~~

Mapping (also called alignment) of sequencing reads allows us to detect
variation in samples by comparing your data to the reference genome. By
doing this you can begin to analyse the relationship between variations
in genotype and phenotype in the population. Next generation sequencing
produces single-end or paired-end reads. For single-end sequence reads,
the sequencer reads the fragment only from one end and usually stops
before reaching the other. For paired-end reads, sequencing starts off
at one end, reads a specified numbers of base pairs, and then starts
another round of the reading from the opposite end of the fragment.
Paired-end sequencing improves the ability to detect genetic
rearrangements (e.g. deletions). This is due to the additional data
carried by pairing reads - they can only be a certain maximum distance
away from each other which limits the regions of the genome to which
they can be mapped. This is particularly useful for regions which are
repeated throughout the genome.

To compare your data to the reference genome, you need to find a
corresponding part of that sequence for each of the reads in our data
– this is the essence of sequence mapping. Following mapping, you will
be able to look at specific variations (SNPs, InDels etc).


Spliced Mapping with TopHat2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This application is used to map Raw Reads with transcriptomic data like
RNA-seq to a Reference Genome, taking into account splice junctions.

Let’s take a look at the app page and talk about various parameters:

|spliced mapping|

Details on various settings:

1)If you are using strand-specific RNA-seq data, the option
“Strand-specificity protocol” will let you choose between the “dUTP” and
“ligation” method. If you are not sure whether your RNA-seq data is
strand-specific or not, you can try using Subsample reads to make a
small subsample, map it with Spliced Mapping and check the coverage in
Genome Browser for genes on the two
strands.

2)By default, the application uses annotated transcripts from the
Reference Genome to distinguish between novel and known junctions. Using
the option “Rule for mapping over known annotation” you can restrict
mappings only across known junctions or infer splice junctions without
any reference annotation.

3)With default settings, the application will report the single best
mapping for each read, even if there are multiple valid mapping
positions. The option “Number of “best” mappings to report” lets you
increase the number of reported mappings. This can be used together with
“Rule for filtering mappings” to choose whether to keep reads mapping to
uniquely or to multiple positions, e.g. report up to 5 possible
mappings, and only for multi-hit reads. If you want to be stricter, you
can set the number of allowed mismatches from 2 to 1 or 0.

4)For paired reads, using the option “Disallow unique mappings of one
mate” you can discard pairs of reads where one mate maps uniquely and
the other to multiple positions. Selecting “Disallow discordant
mappings” will discard all mappings where the two mates map uniquely but
with unexpected orientation, or where the distance between two mapped
mates differs from and internally estimated fragment length, including
mates mapping to different chromosomes.

We used this app in the Testing Differential Gene Expression tutorial
that can be found
`here <https://www.google.com/url?q=https://genestack.com/tutorial/mapping-rna-seq-reads-onto-a-reference-genome/&sa=D&ust=1480960531934000&usg=AFQjCNFMSiaZdYZX9Sp1-nzMlTdCUM_5DA>`__ (link)

Spliced Mapping with STAR
^^^^^^^^^^^^^^^^^^^^^^^^^

Gene Quantification with RSEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

RSEM Report application uses STAR mapper [2] to align reads against
reference transcripts and applies `Expectation-Maximization
algorithm <https://www.google.com/url?q=https://en.wikipedia.org/wiki/Expectation%25E2%2580%2593maximization_algorithm&sa=D&ust=1480960531936000&usg=AFQjCNEBalCXSFMsrtTXF_2Onc0ebj2Djw>`__ [3]
to estimate gene and isoform expression levels from RNA-Seq data.

However, it’s important to know the fragment length distribution to
estimate expression levels from for single-end data accurately. In this
case, you need to specify the “Estimated average fragment length (for
single-end reads only)”. Typical Illumina libraries produce fragment
lengths ranging from 180–200 bp. By default it’s equal to 190. For
paired-end reads, the average fragment length can be directly estimated
from the reads.

Also, you can set the “Estimated standard deviation of fragment length
(for single-end reads only)” (the default value is 20). For paired-end
reads this value will be estimated from the input data.

“The RNA-Seq protocol used to generate the reads is strand specific”? If
yes, check it. By default, the app considers the reads as
non-strand-specific.

When the app did its job, click View report in Explore section to get
gene and isoform level expression estimates. The output report
represents a table with the following main columns:

-  transcript\_id - name of the transcript;
-  gene\_id - name of the gene which the transcript belongs to. If no
   gene information is provided, gene\_id and transcript\_id are the
   same;
-  length - transcript's sequence length (poly(A) tail is not counted);
-  effective\_length - counts only the positions that can generate a
   valid fragment. If no poly(A) tail is added, effective length is
   equal to transcript length - mean fragment length + 1. If one
   transcript's effective length is less than 1, this transcript's both
   effective length and abundance estimates are set to 0;
-  expected\_count - the sum of the posterior probability of each read
   comes from this transcript over all reads;
-  TPM - transcripts per million normalized by total transcript count in
   addition to average transcript length;
-  FPKM - fragments per kilobase of exon per million fragments mapped;
-  IsoPct - the percentage of the transcript's abundance over its parent
   gene's abundance. If the parent gene has only one isoform or the gene
   information is not provided, this field will be set to 100.

The application is based on the
`RSEM <https://www.google.com/url?q=http://deweylab.github.io/RSEM/&sa=D&ust=1480960531941000&usg=AFQjCNFjee9tNzlLAIYCCE5MdnaUL6BK6g>`__ program
and
`STAR <https://www.google.com/url?q=https://github.com/alexdobin/STAR&sa=D&ust=1480960531942000&usg=AFQjCNFAhrRRKB4LfbFwdzh7e47fEJTb0Q>`__ mapper,
which are distributed under the GPLv3 license.

References:

#. Li B and Dewey C N. “RSEM: accurate transcript quantification from
   RNA-Seq data with or without a reference genome.” BMC Bioinformatics
   2011 12:323, doi: 10.1186/1471-2105-12-323
#. Dobin A, Davis C A, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut
   P, Chaisson M and Gingeras T R. “STAR: ultrafast universal RNA-seq
   aligner.” Bioinformatics 2012 29(1): 15-21.
#. Do C B and Batzoglou S. “What is the expectation maximization
   algorithm?” Nature biotechnology, 2008 26(8): 897-899.

Gene Quantification with HTSeq-count
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Differential Gene Expression Analysis

Input: Mapped Reads and Reference Genome

Output: Mapped Read Counts (containing information about number of reads
overlapping each gene specified in the reference annotation)

Further apps to use: Test Differential Gene Expression

Depending on your tasks, you should specify the feature type for which
overlaps choosing from “exon”, “CDS” (coding DNA sequence), “3’UTR” (the
3’ untranslated region) or “5’UTR” (the 5’ untranslated region). For
example, you may consider each exon as a feature in order to check for
alternative splicing.

By default, the “gene-id” will be used as a feature identifier. If some
features will have the same feature identifier the application will
consider all these features as relating to the same feature.

You also need to choose a rule for overlaps that dictates how mapped
reads that overlap genomic features will be treated. There are three
overlap resolution modes: union, strict-intersection, and non-empty
intersection.

The first one - “union” - is the most recommended. It combines all cases
when the read (or read pair) at least partly overlaps the feature. The
“strict-intersection” mode is about strict intersection between the
feature and the read overlapping this feature. But if you are interested
in counting reads that are fully or partly intersected with the feature,
you should use the last mode. It’s important that the read will be
counted for feature if it overlaps precisely only one feature. If the
read overlaps with more than one feature, it will not be counted. 

|image47|

An additional useful option is “Strand-specific reads”. The application
takes into account the direction of the read and the reference, so that
a read from the wrong direction, even if it is mapped to the right
place, will not be counted. This option can be useful if your data is
strand-specific and you are interested in counting of reads overlapping
with feature regarding to whether these reads are mapped to the same or
the opposite strand as the feature. Choose “yes”, if the reads were
mapped to the same strand as the feature and “reverse” - if the reads
were mapped on the opposite strand as the feature. Specify “no”, if you
don’t consider strand-specificity.

Spliced Mapping and quantification with Kallisto
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Kallisto Report application quantifies abundances of transcripts from
RNA-Seq data without the need for alignment. It uses
`Expectation-Maximization
algorithm <https://www.google.com/url?q=http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html&sa=D&ust=1480960531950000&usg=AFQjCNFybXjTk5KPRxNbmjxlhLL-3DHqww>`__ [2]
on “pseudoalignments” to find a set of potential transcripts a read
could have come from.

Use “Strand-specificity protocol” parameter to specify how to process
the pseudoalignments. By default, the app doesn’t consider
strand-specificity (“none” value). To run the app in strand specific
mode, change this value to “forward” if you are interested only in
fragments where the first read in the pair pseudoaligns to the forward
strand of a transcript. If a fragment pseudoaligns to multiple
transcripts, only the transcripts that are consistent with the first
read are kept. The “reverse” is the same as “forward” but the first read
will be pseudomapped to the reverse strand of the transcript.

To correct the transcript abundances according to the model of sequences
specific bias, check “Enable sequence based bias correction” option.

In the case of single-end reads, the “Estimated average fragment length
(for single-end reads only)” option must be used to specify the average
fragment length. Typical Illumina libraries produce fragment lengths
ranging from 180–200 bp. By default it’s equal to 190. For paired-end
reads, the average fragment length can be directly estimated from the
reads.

Also, you can set the “Estimated standard deviation of fragment length
(for single-end reads only)” (the default value is 20). For paired-end
reads this value will be estimated from the input data.

Choose View report app in Explore section to review the Kallisto output
report. It represents the table with the following main columns:

-  target\_id - feature name, e.g. for transcript, gene;
-  length - feature length;
-  eff\_length - effective feature length, i.e. a scaling of feature
   length by the fragment length distribution;
-  est\_counts - estimated feature counts;
-  tpm - transcripts per million normalized by total transcript count in
   addition to average transcript length.

The application is based on the
`Kallisto <https://www.google.com/url?q=https://pachterlab.github.io/kallisto/&sa=D&ust=1480960531955000&usg=AFQjCNE1bXBQ3qF6jc_wDAJCNw22jHdrZg>`__`  <https://www.google.com/url?q=https://pachterlab.github.io/kallisto/&sa=D&ust=1480960531955000&usg=AFQjCNE1bXBQ3qF6jc_wDAJCNw22jHdrZg>`__`program <https://www.google.com/url?q=https://pachterlab.github.io/kallisto/&sa=D&ust=1480960531956000&usg=AFQjCNGGXcC1THxvDTJyvb4No9GcizX0NA>`__ which
is distributed under a non-commercial license you can find
`here <https://www.google.com/url?q=https://pachterlab.github.io/kallisto/download&sa=D&ust=1480960531957000&usg=AFQjCNGSAn4zmuZk6doICwwo5UJRomORqA>`__.

References:

#. Bray N L, Pimentel H, Melsted P and Pachter L. “Near-optimal
   probabilistic RNA-seq quantification.” Nature Biotechnology 2016
   34:525–527, doi:10.1038/nbt.3519
#. Do C B and Batzoglou S. “What is the expectation maximization
   algorithm?” Nature biotechnology, 2008 26(8): 897-899.

Links:
`https://pachterlab.github.io/kallisto/ <https://www.google.com/url?q=https://pachterlab.github.io/kallisto/&sa=D&ust=1480960531958000&usg=AFQjCNFwEnUvLs3G_D06Nm-namQCuBh5lA>`__ (kallisto
tool),
`https://pachterlab.github.io/kallisto/download <https://www.google.com/url?q=https://pachterlab.github.io/kallisto/download&sa=D&ust=1480960531959000&usg=AFQjCNHK84Jj-rL_8ko43QfPhd4F-BSIiw>`__ (licence),
`http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html <https://www.google.com/url?q=http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html&sa=D&ust=1480960531959000&usg=AFQjCNF-2DCoxnVl6lfofDiGhLfeJbq-1w>`__ (EM
algorithm)

Isoforms quantification with cuffQuant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Differential Isoform Expression Analysis

Input: Mapped Reads (corresponding to isoform alignment) and Reference
Genome

Output: Multiple output files corresponding to samples with different
biological conditions and isoforms, can be further processed together
for Differential Isoforms Expression analysis.

Action: The app is used to quantify isoform expression.

Further applications: Test Differential Isoform Expression

Specific genes can produce a range of different transcripts encoding
various isoforms, i.e. proteins of varying lengths containing different
segments of the basic gene sequence. Such isoforms can be generated, for
example, in the process of alternative splicing.

We use this application to calculate expression levels of these
isoforms. It takes the input Mapped Reads (corresponding to isoform
alignment) and Reference Genome files. Multiple output files
corresponding to samples with different biological conditions and
isoforms, can be further processed together for Differential Isoforms
Expression analysis using Test Differential Isoform Expression
application.

Before running the application, you can choose strand-specificity
protocol used for generating your reads. By default, the application
takes “none” strand-specific data, but this value can be changed to
“dUTP” or “RNA-ligation”.

Switch the “No correction by effective length” option if you’d like to
not apply effective length normalization to transcript FPKM (fragments
per kilo bases of exons for per million mapped reads).

The application always makes an initial estimation procedure to more
accurately weight reads mapping to multiple places in the genome.

This application is based on cuffQuant which is a part of
`Cufflinks <https://www.google.com/url?q=http://cufflinks.cbcb.umd.edu/&sa=D&ust=1480960531964000&usg=AFQjCNHU3aK5lX71_5lPCL820JdJ4BeLtw>`__.

Test Differential Gene Expression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Differential Gene Expression Analysis

Input: Mapped Read Counts (from Quantify Raw Coverage in Genestack app)

Action: The app performs differential gene expression analysis between
groups of samples. You can create these groups manually or apply auto
grouping when the application helps you to group your samples according
to experimental factor such as disease, tissue, sex, cell type, cell
line, treatment, organism, etc.

Further apps to use: Expression Navigator app

The application supports two statistical R packages - DESeq2 and edgeR
to perform normalization across libraries, fit negative binomial
distribution and likelihood ratio test (LRT) using generalized linear
model (GLM). With edgeR, one of the following types of dispersion
estimate is used, in order of priority and depending on the availability
of biological replicates: Tagwise, Trended, or Common. Also, edgeR is
much faster than DESeq2 for fitting GLM model, but it takes slightly
longer to estimate the dispersion. It’s important that edgeR gives
moderated fold changes for the extremely lowly DE genes which DESeq2
discards, showing that the likelihood of a gene being significantly
differentially expressed is related to how strongly it's expressed. So,
choose one of the packages according to your desires and run the
analysis.

For each group, a GLM LRT is carried out to find Differentially
Expressed (DE) genes in this group compared to the average of the other
groups. In the case of 2 groups, this reduces to the standard analysis
of finding genes that are differentially expressed between 2 groups.
Thus, for N groups, the application produces N tables of Top DE genes.
Each table shows the corresponding Log2(Fold Change), Log2(Counts per
Million), P-Value, and False Discovery Rate for each gene. Look at all
result tables and plots in Expression Navigator for Genes application.

1) Log2(Fold Change). Let’s assume, that we have two groups - with tumor
and with control samples. Then, for each gene in sample we know read
counts (output of Quantify Raw Coverage in Genes application). If we
divide read counts value for gene X (in the tumor sample) by the read
counts value for gene X (in the control sample) we’ll get Fold Change
value:

Fold Change = tumor/control

And if we apply Log2 transform for this value we’ll get Log2(Fold
Change).

2) Log2(Counts per Million). Dividing each read count by millions yields
counts per million (cpm), a simple measure of read abundance that can be
compared across libraries of different sizes. And if we apply
Log2 transform for this value we’ll get Log2(Counts per Million).

3) p-value. The application also counts p-value for each gene. A low
p–value is seen as evidence that the null hypothesis may not be true
(i.e., our gene is differentially expressed).

4) False discovery rate. FDR is the expected proportion of Type I errors
among the rejected hypotheses.

Expression Navigator for RNA-seq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Differential Gene expression analysis, Differential Isoform
expression analysis.

|image48|Used to: filter and view the results of differential gene
expression analyses, including isoform expression.

The Expression navigator page contains 4 sections.

The topmost section, “Groups Information”, is a summary of the groups
available for comparison. Size refers to the number of samples used to
generate each group. The drop-down selection menu lets you choose which
groups to compare.

The leftmost section allows you to filter and choose genes for
comparison. You can filter by maximum acceptable false discovery rate
(FDR), up or down regulation, minimum log fold change (LogFC), and
minimum log counts per million (LogCPM).

1) Log2(Fold Change). Let’s assume, that we have two groups - with tumor
and with control samples. Then, for each gene in a sample we know read
counts (output of Quantify Raw Coverage in Genes application). If we
divide read counts value for gene X (in the tumor sample) by the read
counts value for gene X (in the control sample) we’ll get the Fold
Change value:

Fold Change = tumor/control

And if we apply a Log2 transform for this value we’ll get Log2(Fold
Change). Genes with positive Log FC are considered to be up-regulated in
the selected group, ones with negative Log FC are down-regulated.

2) Log2(Counts per Million). Dividing each read count by millions yields
counts per million (cpm), a simple measure of read abundance that can be
compared across libraries of different sizes. And if we apply Log2
transform for this value we’ll get Log2(Counts per Million).

Counts per Million =  reads(gene)\^106/reads(all genes)

3) p-value. The application also counts p-value for each gene. A low
p–value is seen as evidence that the null hypothesis may not be true
(i.e., our gene is differentially expressed).

4) False discovery rate. FDR is the expected proportion of Type I errors
among the rejected null hypotheses. In other words, it’s the fraction of
genes for which a significant variation was identified incorrectly. You
can read more about it
`here <https://www.google.com/url?q=http://www.cbil.upenn.edu/PaGE/fdr.html&sa=D&ust=1480960531980000&usg=AFQjCNGB8RddgrwvSzTLjucTxGejSMgqEA>`__.

The buttons at the bottom of the section allow you to refresh the list
based on your filtering criteria or clear your selection.

The top right
section contains a box
plots of expression levels. Genes are listed on the x axis with one bar
present for each  selected group. Log normalized expression levels are
plotted on the y axis.

The bottom right section contains a search box for genes of interest.
You can search for one gene at a time with auto-complete functionality.
These genes do not need to be on the filtered list.

This application is based on two R packages -
`DESeq2 <https://www.google.com/url?q=http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html&sa=D&ust=1480960531982000&usg=AFQjCNHUCbU9X0qoNDWWc4dqZwMwGZOhNw>`__ and
`edgeR <https://www.google.com/url?q=http://www.bioconductor.org/packages/2.13/bioc/html/edgeR.html&sa=D&ust=1480960531983000&usg=AFQjCNGga-RmFYLwPva_sWSoZyblphb2ig>`__.

You can read more about this app in the following
`tutorial <https://www.google.com/url?q=https://genestack.com/tutorial/counting-reads-mapped-to-annotated-features/&sa=D&ust=1480960531983000&usg=AFQjCNGnBFzRjgLDGgRvjTDe0umXpihQ1w>`__.

References:

#. Love MI, Huber W and Anders S. “Moderated estimation of fold change
   and dispersion for RNA-seq data with DESeq2.” Genome Biology. 2014;
   15(12):550.
#. Robinson MD, McCarthy DJ and Smyth GK. “edgeR: a Bioconductor package
   for differential expression analysis of digital gene expression
   data.” Bioinformatics. 2010; 26(1):139-140.

Links:
`http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html <https://www.google.com/url?q=http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html&sa=D&ust=1480960531985000&usg=AFQjCNHmLvVJ1fGf6hrnWo7wzwmUMFX5aA>`__ (DESeq2)
and
`http://www.bioconductor.org/packages/2.13/bioc/html/edgeR.html <https://www.google.com/url?q=http://www.bioconductor.org/packages/2.13/bioc/html/edgeR.html&sa=D&ust=1480960531985000&usg=AFQjCNGhfjFBlDDWaXv1rE6MuQhQpFzWrA>`__ (edgeR);
`http://www.cbil.upenn.edu/PaGE/fdr.html <https://www.google.com/url?q=http://www.cbil.upenn.edu/PaGE/fdr.html&sa=D&ust=1480960531986000&usg=AFQjCNEL6qdg6agmqnsXX1OwJqfeK99CNQ>`__ (FDR)

Expression Navigator for splice isoforms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Testing Differential Isoform Expression Analysis

Input: Multiple Data files with FPKM isoform counts (produced by
Quantify FPKM Coverage in Isoforms app)

Action: The app performs differential expression (DE) analysis between
two groups of samples corresponding to different conditions. You can
create these groups manually or apply auto grouping when the application
helps you to group your samples according to experimental factor such as
disease, tissue, cell type, cell line, treatment etc. It’s important
that to run DE analysis you need to create only two groups.

Results exploration:  Expression Navigator for Isoforms

In “Program Options” section you can apply two options to run the
analysis. The first one - “Apply fragment bias correction” - if checked,
the application runs the bias detection and correction algorithm which
can significantly improve accuracy of transcript abundance estimates.
Use the second option “Apply multiple reads correction” if you’d like to
apply the multiple reads correction.

The application finds isoforms that are differentially expressed between
2 groups and produces 2 tables of Top DE transcripts. Each table shows
the corresponding Log2(Fold Change), Log2(Counts per Million), P-Value,
and False Discovery Rate for each isoform. To visualize your results run
Expression Navigator for Isoforms application.

1) Log2(Fold Change). Let’s assume, that we have two groups - with tumor
and with control samples. Then, for each transcript in sample we know
read counts (output of Quantify Raw Coverage in Genes
 application). If we divide read counts
value for transcript X (in the tumor sample) by the read counts value
for transcript X (in the control sample) we’ll get Fold Change value:

Fold Change = tumor/control

And if we apply Log2 transform for this value we’ll get Log2(Fold
Change).

2) Log2(Counts per Million). Dividing each read count by millions yields
counts per million (cpm), a simple measure of read abundance that can be
compared across libraries of different sizes. And if we apply
Log2 transform for this value we’ll get Log2(Counts per Million).

3) p-value. The application also counts p-value for each isoform. A low
p–value is seen as evidence that the null hypothesis may not be true
(i.e., our isoform is differentially expressed).

4) False discovery rate. FDR is the expected proportion of Type I errors
among the rejected hypotheses.

This application is based on cuffdiff which is a part of
`Cufflinks <https://www.google.com/url?q=http://cufflinks.cbcb.umd.edu/&sa=D&ust=1480960531994000&usg=AFQjCNH88na23xIz5CUAowl7LLWgSpC31A>`__.

Single cell RNA-seq analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Single-cell RNA-seq Analyser + Visualiser need to be merged

Good for: Single-cell RNA-seq Analysis

Input: single-cell RNA-seq data

Action: The app identifies heterogeneously-expressed (HE) genes across
cells, while accounting for technical noise.

Further apps to use : Single-cell RNA-seq Visualiser

The application supports two algorithms for HE analysis. The first uses
spike-in data (artificially introduced RNAs of known abundance) to
calibrate a noise model [1]. The second method is a non-parametric
algorithm based on smoothing splines and doesn’t require the presence of
spike-in data.

To identify highly variable genes you can try different options.
“Exclude samples with low coverage” option (switched by default) allows
you to exclude or include for analysis samples with low read
counts.

Set significance level for the p-value (-10log₁₀(p)). The application
will use the default of 1, which corresponds to selecting genes for
which P is smaller than 0.1.

The “Use spike-ins to calibrate noise” option determines whether or not
spike-in data should be taken into account. If you select only one
folder before running the app, you will use spike-free algorithm and
this option will be switched off by default. But if you select two
folders, one for biological and the other for spike-in data, you can use
the Brennecke algorithm [1]  which requires this option.  

The next three options will be available if spike-ins are included in
the experiment and “Use spike-ins to calibrate noise” option is
switched. You’ll be able to set “Expected
biological CV” which is the minimum threshold chosen for quantifying the
level of biological variability (CV - coefficient of variation) expected
in the null hypothesis of the model. The default value is 0.5.

The other two options - “Noise fit - proportion of genes with high CV²
to remove” and “Noise fit - proportion of genes with low mean expression
to remove” - enable us to exclude a fraction of spike-in genes to fit
the noise model, because extreme outliers tend to skew the fit. The
default values for these options are 0 and 0.85, consequently.

To look at the HE analysis results, open the created Single-cell RNA-seq
Analysis page in  Single-cell RNA-seq visualizer.

This application is based on such `R
packages <https://www.google.com/url?q=http://cran.r-project.org/&sa=D&ust=1480960532001000&usg=AFQjCNE1yhmhcF9OD882Ld6di-TrSBg14w>`__ as
DESeq, statmod, ape, flashClust and RSJONIO.

References:

#. Brennecke P, Anders S, Kim JK, Kolodziejczyk AA, Zhang X, Proserpio
   V, Baying B, Benes V, Teichmann SA, Marioni JC, Heisler MG.
   “Accounting for technical noise in single-cell RNA-seq experiments.”
   Nature Methods. 2013 Sep 22; 10(11):1093–1095.

Read more about single-cell RNA-seq analysis on Genestack here:
`https://genestack.com/blog/2016/02/22/visualisation-clustering-methods-single-cell-rna-seq-data/ <https://www.google.com/url?q=https://genestack.com/blog/2016/02/22/visualisation-clustering-methods-single-cell-rna-seq-data/&sa=D&ust=1480960532003000&usg=AFQjCNFAjkflTkJ-VOc9Pmyr7WT2N61K8Q>`__

Single-cell RNA-Seq Visualisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Single-cell RNA-seq Analysis

Used to: Explore cell-to-cell variability in gene expression in even
seemingly homogeneous cell populations based on scRNA-Seq datasets.

The application shows basic statistics such as the number of identified
highly variable genes across the analysed samples. It also provides
several quality control (QC) plots allowing to check the quality of raw
sequencing data, estimate and fit technical noise for the
Brennecke algorithm, and detect the genes with significantly high
variability in expression. Expression of the highly variable genes
across all cell samples is represented by an interactive clustered
heatmap. Finally, several plots in the Samples Visualisation section can
be used to detect cell subpopulations and identify novel cell
populations based on gene expression heterogeneity in the single-cell
transcriptomes.

QC plots are adopted from the original paper by Brennecke et al [1]. In
all the plots described below, gene expression levels are normalized
using the DESeq normalization procedure [2].

QC plots are adopted from the original paper by Brennecke et al [1].

The first plot describing the quality of raw data is the Scatter Plot of
Normalised Read Counts, which shows the cell-to-cell correlation of
normalized gene expression levels. Each dot represents a gene, its
x-coordinate is the normalized gene count in the first cell, and its
y-coordinate is the normalized gene count in the second cell. If
spike-ins were used during the analysis, separate plots will be rendered
for spike-in genes and for sample genes.

The Technical Noise Fit and Highly Variable Genes plots provide a visual
summary of the gene expression noise profile in your dataset across all
cells. They graph the squared coefficient of variation (CV2) against the
average normalized read counts across samples.  The Gene Expression
Variability QC plot allows you to visualize the genes whose expression
significantly varies across cells. A gene is considered as highly
variable if its coefficient of biological variation is significantly
higher than 50% (CV2 > 0.25)  and the biological part of its coefficient
of variation is significantly higher than a user-defined threshold (its
default value is 50%, and can be modified in the Single-cell
Analyser). The coefficient of variation is defined as the standard
deviation divided by the mean. It is thus a standardized measure of
variance.

If spike-ins were used to calibrate technical noise, then the separate
Technical Noise Fit plot is displayed. On this plot, each dot
corresponds to a “ technical gene” (spike-in gene).It plots the mean
normalized count across all samples on the x-coordinate and the squared
coefficient of variation (CV2) of the normalized counts across all
samples on the y-coordinate. The coefficient of variation is defined as
the standard deviation divided by the mean. It is thus a standardized
measure of variance. The plot also represents the fitted noise model as
a solid red line (with 95% confidence intervals as dotted red lines). It
allows you to check whether the noise model fits the data reasonably
well. If it is not the case, you should change the noise fitting
parameters in the Single-cell Analysis application.

The interactive heatmap depicts the log normalised read count of each
significant highly variable gene (rows) in each cell sample (columns).
Hierarchical clustering of molecular profiles from cell samples is based
on the similarity in gene expression of highly expressed genes and
allows identification of  molecularly distinct cell populations. The
heatmap is clustered both by columns and by rows, to identify clusters
of samples with similar gene expression profiles, and clusters of
potentially co-expressed genes. The bi-clustered heatmap is provided by
an open source interactive Javascript library
`InCHlib <https://www.google.com/url?q=http://openscreen.cz/software/inchlib/home/&sa=D&ust=1480960532013000&usg=AFQjCNGnCwLQvBZYAwnvVft_NSwJUYeZrg>`__ (Interactive
Cluster Heatmap library) [3].

The Samples Visualisation section provides interactive plots used to
cluster cell samples based on expression of highly variable genes.
Currently, two alternative methods are supported for visualisation and
clustering of samples: the first one is based on the t-distributed
Stochastic Neighbour Embedding (t-SNE) algorithm [4] and the second one
uses Principal Component Analysis (PCA). For automatic cluster
identification, the k-means clustering algorithm can be used in
combination with either  t-SNE or PCA.
K-means clustering requires you to supply
a number of clusters to look for ("k"). You can either enter it manually
using the dropdown menu or use the suggested value estimated using the
"elbow" method (choosing a value of k such that increasing the number of
clusters does not significantly reduce the average "spread" within each
cluster).

The Interactive Principal Component Analysis (PCA) scatter plot is
rendered using the
`NVD3 <https://www.google.com/url?q=http://nvd3.org/&sa=D&ust=1480960532015000&usg=AFQjCNGqXKChcZFjmqBSR5lfGkPjYLtq_A>`__ Javascript
library. The PCA features and k-means algorithm results are computed
using R's built-in functions
`prcomp <https://www.google.com/url?q=https://stat.ethz.ch/R-manual/R-patched/library/stats/html/prcomp.html&sa=D&ust=1480960532015000&usg=AFQjCNG0r7sbyWopaE14KyEE4d1vgwm92A>`__ and
`knn <https://www.google.com/url?q=https://stat.ethz.ch/R-manual/R-devel/library/class/html/knn.html&sa=D&ust=1480960532016000&usg=AFQjCNEqyNo-UhfT52yacNJBHNwelCFISA>`__.
The t-SNE transformation is computed using the
`Rtsne <https://www.google.com/url?q=http://cran.r-project.org/web/packages/Rtsne/index.html&sa=D&ust=1480960532017000&usg=AFQjCNGbgjxYIH_Ao0k-ARQ5A9JAqJLUwQ>`__ package.

You can read more about the app
`here <https://www.google.com/url?q=https://genestack.com/blog/2016/02/22/visualisation-clustering-methods-single-cell-rna-seq-data/&sa=D&ust=1480960532018000&usg=AFQjCNGwmsnPH2lWurlcrYcwrekhm-9OkQ>`__.

References: 

#. Brennecke P, Anders S, Kim JK, Kołodziejczyk AA, Zhang X, Proserpio
   V, Baying B, Benes V, Teichmann SA, Marioni JC and Heisler MG.
   “Accounting for technical noise in single-cell RNA-seq experiments.”
   Nature Methods 2013 10(11), 1093-1095. PMID: 24056876
#. Anders S and Huber W. "Differential expression analysis for sequence
   count data". Genome Biology 2010 11:R106  
#. Škuta C, Bartůněk P and Svozil D. “InCHlib–interactive cluster
   heatmap for web applications.” Journal of Cheminformatics 2014 6(1),
   1-9.
#. van der Maaten LJP and Hinton GE. “Visualizing High-Dimensional Data
   Using t-SNE.” Journal of Machine Learning Research 2008 9(11),
   2579-2605

Read more about single-cell RNA-seq analysis on Genestack here:
`https://genestack.com/blog/2016/02/22/visualisation-clustering-methods-single-cell-rna-seq-data/ <https://www.google.com/url?q=https://genestack.com/blog/2016/02/22/visualisation-clustering-methods-single-cell-rna-seq-data/&sa=D&ust=1480960532022000&usg=AFQjCNFHLp_YAJtq-t55uRJlHo1K1NAPwg>`__

NOTE: Reference Genomes

One way or another, most bioinformatics analysis pipelines, regardless
of the data type analysed, require the use of a reference genome. For
instance,  we use reference genomes in `DNA methylation
analysis <https://www.google.com/url?q=https://genestack.com/tutorial/whole-genome-bisulfite-sequencing-analysis/&sa=D&ust=1480960532024000&usg=AFQjCNEON1E936WzebWy5w5hCqDobFfbyQ>`__,
in `differential gene expression
analysis <https://www.google.com/url?q=https://genestack.com/tutorial/testing-differential-gene-expression-on-genestack-platform/&sa=D&ust=1480960532025000&usg=AFQjCNF8iK-m3LAGKdEi3YCpFxG4BQO4jg>`__,
and analysis of the `transcriptomic heterogeneity within populations of
cells <https://www.google.com/url?q=https://genestack.com/blog/2014/09/24/single-cell-rna-seq-analysis-tutorial/&sa=D&ust=1480960532025000&usg=AFQjCNF8rzNCeKOex8EvDd2Y0DVNHe855A>`__.
The choice of a reference genome can increase the quality and accuracy
of the downstream analysis or it can have a harmful effect on it. For
instance, it has been shown that the choice of a gene annotation has a
big impact on RNA-seq data analysis, but also on `variant effect
prediction <https://www.google.com/url?q=http://www.intechopen.com/books/references/next-generation-sequencing-advances-applications-and-challenges/impact-of-gene-annotation-on-rna-seq-data-analysis%23B34&sa=D&ust=1480960532026000&usg=AFQjCNHaCNQIeKNrxp0ot4JjVhZTFfN3fA>`__[
1, 2].

On Genestack, you can find `several reference
genomes <https://www.google.com/url?q=https://platform.genestack.org/endpoint/application/run/genestack/signin?original_url%3D%252Fendpoint%252Fapplication%252Frun%252Fgenestack%252Ffilebrowser%253Fa%253DGSF000018%2526action%253DviewFile%2526page%253D1&sa=D&ust=1480960532027000&usg=AFQjCNHhglPt1_5NPlkDPhTBlVLosmZJ6A>`__ for
some of the most common model organisms. We are adding more and more
reference genomes of model organisms to this list regularly.

 

For some organisms we provide several genomes, e.g.  there are 3
reference genomes for H. sapiens. What are the differences between these
reference genomes? And how do you chose the correct one?  The answer is
not so straightforward and depends on several factors – let’s discuss
each of them:

 

1) Versions of the reference genome

For instance:  Homo sapiens GRCh37.75 (unmasked) vs GRCh38.80 (unmasked)

The numbers correspond to versions (or “builds”) of the reference genome
– the higher the number, the more recent the version. We generally
recommend you use the latest version possible. One thing to remember is
that for the newest genome builds, it’s likely that resources such as
genome annotations and functional information will be limited, as it
takes time for Ensembl/ UCSC to integrate additional genomic data with
the new build. You can read more about it a `blog
post <https://www.google.com/url?q=http://genomespot.blogspot.co.uk/2015/06/mapping-ngs-data-which-genome-version.html&sa=D&ust=1480960532030000&usg=AFQjCNFBJPxQvY3k5N-9Vf16-S9qYj1Sqg>`__ from
Genome Spot blog and in
`this <https://www.google.com/url?q=http://www.bio-itworld.com/2014/1/27/getting-know-new-reference-genome-assembly.html&sa=D&ust=1480960532030000&usg=AFQjCNHF02fY7GpNpuUrkVxx4steDpXYng>`__`article <https://www.google.com/url?q=http://www.bio-itworld.com/2014/1/27/getting-know-new-reference-genome-assembly.html&sa=D&ust=1480960532030000&usg=AFQjCNHF02fY7GpNpuUrkVxx4steDpXYng>`__ from
Bio-IT.

 

2) One organism – many strains

K12 and O103:H2 are two different strains of E.coli.
`K12 <https://www.google.com/url?q=https://www.genome.wisc.edu/resources/strains.htm&sa=D&ust=1480960532032000&usg=AFQjCNH9An3zJ5ilwpvmVlEVsxbLsRLFFA>`__ is
an innocuous strain commonly used in various labs around the world.
`O103:H2 <https://www.google.com/url?q=http://aem.asm.org/content/79/23/7502.full&sa=D&ust=1480960532033000&usg=AFQjCNEHcWj9cqdYfeXTHPadk8p4twNGrQ>`__ is
a pathogenic strain, commonly isolated from human cases in Europe.
Depending on your experiment, you should choose a matching reference
genome.

 

3) Masked, soft-masked and unmasked genomes

There are three types of Ensembl reference genomes: unmasked,
soft-masked and masked.

Generally speaking, it’s recommended to use unmasked reference genome
builds for alignment.

Masking is used to detect and conceal interspersed repeats and low
complexity DNA regions so that they could be processed properly by
alignment tools.

Masking can be performed by special tools, like
`RepeatMasker <https://www.google.com/url?q=http://www.repeatmasker.org/&sa=D&ust=1480960532035000&usg=AFQjCNHVO50QRN7dUkikjWjvgoLYePqSgg>`__.
 The tool goes through DNA sequence looking for repeats and
low-complexity regions.

There are two types of masked reference genomes: masked and soft-masked.

MASKED

Masked reference genomes are also known as hard-masked DNA sequences.
Repetitive and low complexity DNA regions are detected and replaced with
‘N’s. The use of masked genome may adversely affect the analysis
results, leading to wrong read mapping and incorrect variant calls.

 

When should you use a masked genome?

We generally don’t recommend using masked genome, as it relates to the
loss of information (after mapping, some “unique” sequences may not be
truly unique) and does not guarantee 100% accuracy and sensitivity (e.g.
masking cannot be absolutely perfect). Moreover, it can lead to the
increase in number of falsely mapped reads.

SOFT-MASKED

 

In soft-masked reference genomes, repeats and low complexity regions are
also detected but in this case they are masked by converting to a
lowercase variants of the base (e.g. acgt).

 

When should you use a soft-masked genome?

The soft-masked sequence does contain repeats indicated by lowercase
letters, so the use of soft-masked reference could improve the quality
of the mapping without detriment to sensitivity. But it should be noted
that most of the alignment tools do not take into account soft-masked
regions, for example BWA, tophat, bowtie2 tools always use all bases in
alignment weather they are in lowercase nucleotides or not. That is why,
there is no actual benefit from the use of soft masked genome in
comparison with unmasked one.

 

Therefore, we recommend you use unmasked genomes when you don’t want to
lose any information. If you want to perform some sort of filtering,
it’s better to do so  after the mapping step.

 

Example:

To perform WES analysis, we recommend you use an unmasked reference
genome of the latest releases and assemblies (e.g. Homo sapiens /
GRCh38.80 (unmasked) for human samples).

 

The bioinformatics community is divided on the topic of the use of
reference genomes. It is our personal opinion that it is best to always
use unmasked genome and perform filtering after the mapping step.
However, if you would like to read more on the topic, we suggest taking
a look at the following papers:

#. `McCarthy DJ, Humburg P, Kanapin A, Rivas MA, Gaulton K, Cazier JB,
   Donnelly P. Choice of transcripts and software has a large effect on
   variant annotation. Genome Med.
   2014;6(3):26; <https://www.google.com/url?q=https://genomemedicine.biomedcentral.com/articles/10.1186/gm543&sa=D&ust=1480960532039000&usg=AFQjCNFQKTLLLg3B69W8VzPfoavzieNoow>`__

2. `Frankish A, Uszczynska B, Ritchie GR, Gonzalez JM, Pervouchine D,
   Petryszak R, et al. Comparison of GENCODE and RefSeq gene annotation
   and the impact of reference geneset on variant effect prediction. BMC
   Genomics. 2015;16 (Suppl
   8):S2 <https://www.google.com/url?q=http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-16-S8-S2&sa=D&ust=1480960532040000&usg=AFQjCNEhK7CXAJi8svzmvtqxfNNceHfm2w>`__
