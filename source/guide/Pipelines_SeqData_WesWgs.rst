Genome/ exome sequencing
~~~~~~~~~~~~~~~~~~~~~~~~

Unspliced Mapping with BWA
^^^^^^^^^^^^^^^^^^^^^^^^^^

On Genestack, you will find two Unspliced Mapping applications. This one
is based on the
`BWA <https://www.google.com/url?q=http://bio-bwa.sourceforge.net/&sa=D&ust=1480960532041000&usg=AFQjCNFYstLRthqjGP-BtyzLwe9HS6KRLg>`__tool
and is used to map exome sequencing reads to a Reference Genome. It is
meant to be used further with our Variant Calling application which is
in turn based on `samtools
mpileup <https://www.google.com/url?q=http://samtools.sourceforge.net/mpileup.shtml&sa=D&ust=1480960532042000&usg=AFQjCNGpkrHDwz5QYy5CU1RuQFJoCWqgIQ>`__.

|unspliced mapping with BWA|

BWA’s MEM algorithm will be used to map paired or single-ends reads from
70 bp up to 1Mbp (“mem” option in command line). For reads up to 70 bp
the algorithm called BWA-backtrack will be applied. This algorithm is
implemented with the “aln” command, which produces the suffix array (SA)
coordinates of the input reads. Then the application converts these SA
coordinates to chromosome coordinates using the “samse” command (if your
reads are single-end) or “sampe” (for paired-end reads).

We used this app in the Whole Exome Sequencing Data Analysis tutorial
that can be found
`here <https://www.google.com/url?q=https://genestack.com/tutorial/whole-exome-sequencing-data-analysis-on-genestack-platform/&sa=D&ust=1480960532043000&usg=AFQjCNEgMlyhiYZgyATe8MnVYwl2hoL55Q>`__.

Unspliced Mapping with Bowtie2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This application is based on the
`Bowtie2 <https://www.google.com/url?q=http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml&sa=D&ust=1480960532044000&usg=AFQjCNFOzcbBeA6op29d_stzX10eJKYp4w>`__tool
and is used to map sequencing libraries to a Reference Genome. Suitable
sequencing methods include DNA-seq and ChIP-seq.

|unspliced mapping with bowtie2|

Let’s talk a bit about the various settings:

1)By default the application will report the best mapping for one
mappable read. If you are interested in reads mapping to multiple
positions, switch off this option and set N mappable positions for one
read in the text box “Limit the number of mappings to search”.

2)You can apply a rule for filtering mappings to choose whether to keep
reads mapping uniquely or to multiple positions. If you want to be
stricter, you can change the maximum number of allowed mismatches, e.g.
if you set it to 1, any mapping with 2 or more mismatches won’t be
reported.

3)For paired reads, using the option “Disallow unique mappings of one
mate” you can discard pairs of reads where one mate maps uniquely and
the other to multiple positions. Selecting “Disallow discordant
mappings” will discard all mappings where the two mates map uniquely but
with unexpected orientation or where the distance between two mapped
mates differs from and internally estimated fragment length, including
mates mapping to different chromosomes.

Read more about differences between BWA and Bowtie2 on our
`forum <https://www.google.com/url?q=http://forum.genestack.org/t/unspliced-mapping-with-bwa-app-vs-unspliced-mapping-with-bowtie2-app/36/2&sa=D&ust=1480960532046000&usg=AFQjCNHMGtJKMz1PN9VHw-BLzEMS4G5bYw>`__.

Variant calling with samtools and bcftools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Variant Calling, Whole Exome Sequencing Analysis, Whole Genome
Sequencing Analysis

Input: Mapped Reads

Action: identifying genomic variants from Mapped Reads files.

Further apps to use:Effect Prediction and Genome Browser or Variant
Explorer for exploring results

The app uses samtools mpileup which automatically scans every position
supported by an aligned read, computes all the possible genotypes
supported by these reads, and then calculates the probability that each
of these genotypes is truly present in your sample.

As an example, let’s consider the first 1000 bases in a Reference Genome
file. Suppose the position 35 (in reference G) will have 27 reads with a
G base and two reads with a T nucleotide. Total read depth will be 29.
In this case, the application concludes with high probability that the
sample has a genotype of G, and the T reads are likely due to sequencing
errors. In contrast, if the position 400 in reference genome is T, but
it is covered by 2 reads with a C base and 66 reads with a G (total read
depth equal to 68), it means that the sample more likely will have G
genotype.

Then the application executes bcftools call which uses the genotype
likelihoods generated from the previous step to call genetic variants
and outputs the all identified variants in the Genetic Variations file.

By default, the application call both SNPs and indels, but if you’d like
to report only SNPs change “Variants to report” option to “SNPs only”
value. Also, you can tell the application to call only multi-allelic
variants, switching the “Call multi-allelic variants” option. The
multiallelic calling is recommended for most tasks.

To skip anomalous read pairs in variant calling, use option “Discard
anomalous read pairs” checked.

In some cases, it’ll be interested to report only potential variant
sites and exclude monomorphic ones (sites without alternate alleles) in
output Genetic Variation file. For this purpose, switch the option “Only
report variant sites”.

The application allows you to set maximum read depth to consider per
position and minimum number of gapped reads for an INDEL candidate. The
default value for the first option is 250 reads at the position per
input Mapped Reads file. For the second one, value is not set by
default.

Moreover, base alignment quality (BAQ) recalculation is turned on by
default. It helps to rule out false positive SNP calls due to alignment
artefacts near small indels.

The application will always write DP (number of reads covering
position), DV (number of high-quality variant reads), DP4 (number of
forward reference, reverse reference, forward non-reference and reverse
non-reference alleles used in variant calling) and SP (phred-scaled
strand bias P-value) tags in the output file.

You are also able to select chromosomes for analysis, using “Chromosome
to analyse” option and merge samples with the same metainfo key
(specify “Key to merge samples)”. The last option can be useful for
merging technical replicates.

The result Genetic Variations file can be opened in Genome Browser as a
separate  variation track, further annotated using Effect Prediction
application, or viewed immediately using Variant Explorer application.

This application is based on
`SAMtools <https://www.google.com/url?q=http://www.htslib.org/doc/samtools-1.1.html&sa=D&ust=1480960532055000&usg=AFQjCNFwdKm7yBHfHi6jm4j8pH433nu17Q>`__ and
`BCFtools <https://www.google.com/url?q=http://www.htslib.org/doc/bcftools-1.1.html&sa=D&ust=1480960532055000&usg=AFQjCNFOwJEgoQz7drG9vyiBT7c2nzCelQ>`__.

Variant effect prediction with SnpEff
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Whole Exome Sequencing Analysis, Whole Genome Sequencing
Analysis

Input: Genetic Variations file

Action: The app annotates these variants based on their genomic
locations and calculates the effects they produce on known genes.

Further apps to use: Variant Explorer, View Report

Use Variant Explorer application to know what effect is generated by
each separate variant. If you’d like to see the whole list of effects
and annotations as well as to get some general statistics (for example,
to know number of variants by chromosome, find out how many variants are
corresponding to SNP or insertions, to know number of effects by type
and region and some other information), just open this output annotated
Genetic Variations file in View Report application.

This application is based on open source
`SnpEff <https://www.google.com/url?q=http://snpeff.sourceforge.net/SnpEff_manual.html&sa=D&ust=1480960532059000&usg=AFQjCNFeW4EzcYHgiT0J3ml4QfiSuTPRxg>`__ tool.

Variant Association Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Variant Calling, Whole Exome Sequencing Analysis, Whole Genome
Sequencing Analysis

Input: Genetic Variations files

Action: The application extends input Genetic Variations with P values
of allelic association with variants found in 'control' Genetic
Variations file.

Further apps to use: Variant Explorer

The association test is based on comparing allele frequencies in groups
of samples using two-tailed `Fisher Exact
Test <https://www.google.com/url?q=http://en.wikipedia.org/wiki/Fisher%2527s_exact_test&sa=D&ust=1480960532062000&usg=AFQjCNF-F7myaweRgSBSsp1oC316tPZ9Xw>`__,
which compares tables of alternative allele count and reference allele
count in called genotypes.

For example, if a variant is found to be common in a 'control' group,
but is very rare or enriched in your cohort, than that variant is
associated with your population’s phenotype (e.g. trait or disease) and
have a correspondingly small p-value.

If 'control' file has an additional info of dividing samples into the
smaller groups, the values will be also calculated for the groups.

For example, to reduce spurious allelic association due to population
stratification in `1000 Genomes
Project <https://www.google.com/url?q=http://www.1000genomes.org/&sa=D&ust=1480960532063000&usg=AFQjCNF1ZanVn015hOAQ9yYbMRVrwc4X2g>`__ data,
samples are also divided to the main ancestries groups (American,
European, East Asian, and African).

So if you run the analysis with this file, FET is computed not only
against the combined 1000G population, but also against each of its main
ancestries groups.

To control the false discovery rate due to multiple FET testing,
Benjamini-Hochberg P Value adjustment is applied.

With Variant Explorer application the file can be filtered and sorted by
this values.

Variant Explorer
^^^^^^^^^^^^^^^^

|image51|


Good for: Any analysis type dealing with genetic variants
Used to: Interactively explore genetic variations such as SNPs, MNPs,
and indels at specific genomic positions. The app not only displays the
information about variants but also allows you to sort and filter by
various fields, such as mutation type, quality, locus, etc.

Variant Explorer takes as input a  Genetic Variations file which can be
imported or generated with the Variant Calling app. If you open it in
the app, you’ll see default DP (Raw read depth) and MQ (Average mapping
quality) columns (“Other” tab in “Columns” section).

Variants can be annotated with the Effect Prediction app that analyzes
genomic position of the variants and reveals the effects they produce on
known genes (such as amino acid changes, synonymous and nonsynonymous
mutations, etc.). For such variants the following information will be
shown (find it in “Effect prediction” tab):

-  Effect - effect predicted by SnpEff tool;
-  Impact - impact predicted by SnpEff tool;
-  Functional class - functional class of a region, annotated by SnpEff
   tool.

Moreover, the app calculates “Additional metrics” such as genotype
frequencies for homozygous samples with reference and alteration alleles
(GF HOM REF and GF HOM ALT columns correspondingly), reads depth for
homozygous samples with alteration allele (DP HOM ALT) and reads depth
for heterozygous samples (DP HET).

To change the default columns or add more columns, choose them in the
corresponding tabs in “Columns” section and “Save” your changes. After
that all selected columns will be displayed in Table viewer.

You can read more about this app in the following
`tutorial <https://www.google.com/url?q=https://genestack.com/tutorial/wgs-exploring-variants/%23filtering&sa=D&ust=1480960532068000&usg=AFQjCNFKsWZyvjtKfnI1DPMwtD0YvIw4KA>`__.

Intersect Genomic Features
^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: various analysis types

Input: Mapped Reads file or Genetic Variations file

Action: The app performs an intersection between several feature files
such as Mapped Reads files or Genetic Variations files. 

Output: Depending on input files, you can get different outputs, either
Mapped Reads or Genetic Variations files.

Further apps to use: depends on the analysis type
\

With default settings, the application will report overlapping features
(see option “Rule for filtering”). For example, you could isolate single
nucleotide polymorphisms (SNPs) that overlap with SNPs from another
file. For this, intersect two Genetic Variations files. But there are
cases when you’d like to know which features don’t overlap with other
ones. To get such outputs, use “Report non-overlapping features” filter.

The application has also other possibilities. For example, by setting
minimum overlapping fraction equal to 10 (default value), you can check
whether a feature of interest has at least 10% of its length overlapping
another feature.

The “Rule for overlap strandedness” option allows you to ignore overlaps
on the same strand or on the other strand. By default, overlapping
features are reported without respect to the strandedness.

This application is based on
`BEDtools <https://www.google.com/url?q=http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html&sa=D&ust=1480960532075000&usg=AFQjCNGU8dqh1cQxlk22wUALFNLXZK0Llg>`__.

DbNSFP Annotation
^^^^^^^^^^^^^^^^^

Good for: Whole Exome Sequencing Analysis, Whole Genome Sequencing
Analysis

Input: Genomic Variants

Action: The app processes variants adding annotations from `dbNSFP
database <https://www.google.com/url?q=https://sites.google.com/site/jpopgen/dbNSFP&sa=D&ust=1480960532076000&usg=AFQjCNECi4z5Eln9J4cljP2qQym75ATpeQ>`__,
computes FET (Fisher’s Exact Test) for the corresponding` 1000 Genomes
Project <https://www.google.com/url?q=http://www.1000genomes.org/&sa=D&ust=1480960532077000&usg=AFQjCNFvsi5JusJ5sf2fwpNkB5VlRUdniA>`__ data
and tests for HWE (Hardy-Weinberg equilibrium).

Further apps to use: Annotated variants can be further interactively
analysed in Variant Explorer

The app uses
`VCFtools <https://www.google.com/url?q=http://vcftools.sourceforge.net/&sa=D&ust=1480960532078000&usg=AFQjCNFVEJHrn3K8lJPf6pKO7Hwblhfxuw>`__.
