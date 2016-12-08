Pipelines and applications
==========================

Applications available on Genestack fall into four categories:
Preprocess, Analyse, Explore and Manage.

Preprocess contains all applications used to process files pre- or
post-alignment in order to increase the data quality.

Analyse contains all mappers and all other apps required to analyse
sequencing data.

Explore contains all interactive graphical interface applications
allowing users to view the results of their
computations. Apps for visualizing QC reports, the Genome Browser, apps
 for exploring  genomic variants, and many more.

Manage contains apps used to manage your data: apps dealing with data
flows, file provenance and so on.

An extended version (including information on licensing and references)
of every application description found in this guidebook can be found in
“About application” section of each of the individual apps.

How to find  the About Application section?

On the app page click on the name of the app in the upper left corner
and select  “About application”.

<About Application>|image32|

Sequencing data
---------------

Raw Reads preprocessing & QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequencing Assay Viewer app
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use this app to look through the content of a Sequencing Assay or a Raw
Reads file and search specific nucleotide sequences which can be exact,
reverse, complement or reverse complement to the sequence of interest.

<Sequencing Assay Viewer> |image33|

How to access the app?

Select the assay you are interested in, right click on it and from the
“Explore” section select the name of the app.

FastQC report
^^^^^^^^^^^^^

Let’s begin the analysis pipeline. The usual first step of any NGS data
analysis is quality control of raw sequencing reads. According to the
“garbage in, garbage out” rule, if we begin our analysis with poor
quality reads, we shouldn’t expect great results at the end. Luckily,
there are a few procedures that can be used to improve the data quality
if that proves to be unsatisfactory.

How can you check the quality of the data?

The tool used for raw reads quality check is FastQC Report app, based on
`FastQC
tool <https://www.google.com/url?q=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/&sa=D&ust=1480960531831000&usg=AFQjCNH02ePdMGERL56j74uXPHwHmKTndg>`__ developed
by Simon Andrews at the Babraham Institute.

The quickest way to perform the quality assessment of your data in
Genestack is via the use of one of the public data flows– Raw Reads
Quality Control Data Flow.

In order to use it, select all of your raw reads, right click on them
and from the dropdown menu select “Run dataflow on selection” and choose
the appropriate data flow.  Remember you need to initialize the
computation! On the Data Flow Runner page click on Run Data Flow and
select Start initialization now

You will have to wait for the results (you can track the progress of
your tasks in Task Manager). Once they are completed, you can find your
files in Created Files folder.

Since these files were created using a data flow, they will be located
in one folder (see Platform architecture for more details). To open up
one of these reports, click on the report and from the dropdown menu
select FastQC Report app.

If you don’t want to use a data flow, you can select all of your raw
reads, right click on them, go to “Explore” and select the appropriate
app. In this case you’ll also need to wait until the computations are
done and you will be able to track the progress of your tasks using Task
Manager.

On the FastQC Report page you can view both the result and the
provenance of the report file. At the top of the page you will see the
file name and the version of the FASTQC app used. The View parameters
button will show you the command line options used to generate the
report.  Below that you will see the File Data flow, in this case it
should only contain two app entries. Finally, the results can be viewed
in the Reports section. Here you will find  various graphs that
visualize the quality of your data. We’ll go through all of them one by
one and tell you:

a)how they should  look for data of perfect quality; 

b)how they may look if there’s something wrong with your data; 

c)what you can do if the quality is unsatisfactory.

The metrics table gives you quick indicators as to the status of each of
the quality metrics calculated.  Warnings: Yellow Triangles

Failures: Red X-es

Metrics

1) Basic Statistics

[FastQC 1]|image34|

Info on type and number of reads, GC content, and total sequence length.

2) Sequence Length Distribution

[FastQC 2]|image35|

Reports lengths of all sequences.

Warning

This report will raise an issue if the lengths are not identical, but
this can usually be ignored, as it is expected for some sequencing
platforms.

3) Per sequence GC content

[FastQC 3]|image36|

For data of good quality, the graph will show a normal, bell-shaped
distribution.

Warning

It raises a warning when the sum of the deviations from the normal
distribution represents more than 15% of the reads.

Warnings are usually caused by a presence of contaminants. Sharp peaks
may represent a presence of a very specific contaminant (e.g. an
adaptor). Broader peaks may indicate contamination with a range of
contaminants.

Improving data quality

Best solution: Run the Trim Adaptors and Contaminants preprocessing app.

4) Per base sequence quality plot

[FastQC 4]|image37|

For data of good quality, the median quality score per
base (Phred) should not drop below 20.

Failure 

A failure will be raised if the lower quartile for quality at any
base position is less than 5 or if the median for any base is less than
20.

Improving data quality
Best solution: If the quality of the library falls to a low level over
the course of a read, the blueprint solution is to perform quality
trimming of low quality bases or omitting low quality reads. This  can
be performed using Trim Low Quality Bases or Filter By Quality Score
apps respectively.

5) Per sequence quality scores plot 

[FastQC 5]|image38|

Ideally, we’d want to see a sharp peak at the very end of the graph
(meaning most frequently observed mean quality scores are above 27)

Warning

A warning is raised when the peak is shifted to the left, which means
the most frequently observed mean quality is below 27. This equals to a
0.2% error rate.

Improving data quality

Best solution: Perform quality-based trimming or selection using Trim
Low Quality Bases or Filter By Quality Score apps respectively.

6) Per base sequence content

[FastQC 6]|image39|

 

Ideally, in a random library we would see four parallel lines
representing the relative base composition. Fluctuations at the
beginning of reads in the tested sample may be caused by adapter
sequences or other contaminations of the library.

A bias at the beginning of the reads is common for RNA-Seq data. This
occurs during RNA-seq library preparation, when “random” primers are
annealed to the start of sequences. These primers are not truly random,
and it leads to a variation at the  beginning of the reads.

Warning

A warning will be raised  if the difference between A and T, or G and C
is greater than 10% at any position.

Improving data quality

If there is instability at the start of the read the consensus is that
no QC is necessary. If variation appears over the course of a read the
Trim to Fixed Length app may be used. If there is persistent variation
throughout the read it may be best to discard it. Some datasets may
trigger a warning due to the nature of the sequence. For example,
bisulfite sequencing data will have almost no Cytosines. Some species
may be unusually GC rich or poor and therefore also trigger a
warning.

7) Sequence duplication levels plots

[FastQC 7]|image40|

Reports total number of reads, number of distinct reads and mean
duplication rates.

Warning

This module will issue a warning if non-unique sequences make up more
than 20% of the total.

There are two potential types of duplicates in a library: technical
duplicates arising from PCR artefacts or biological duplicated which are
natural collisions where different copies of exactly the same sequence
are randomly selected. From a sequence level there is no way to
distinguish between these two types and both will be reported as
duplicates here.

Improving data quality

If the observed duplications are due to primer/adaptor contamination,
they can be removed using the Trim Adaptors and Contaminants app. Filter
Duplicated Reads can also be used for DNA sequencing data but will
distort expression data.

8) Overrepresented Sequences

[FastQC 8]|image41|

Shows the highly overrepresented sequences (more than 0.1% of total
sequence) in the sample

Warning

A warning will be raised  if any sequence is found to represent more
than 0.1% of the total.

There are several possible sources of overrepresented sequences:

–technical biases (one region was sequenced several times; PCR
amplification biases)

–feature of library preparation (e.g. for targeted sequencing)

–natural reasons (RNA-Seq libraries can naturally present high
duplication rates)

Overrepresented sequences should only worry you if you think they are
present due to technical biases.

Improving data quality

Procedures and caveats for improving data quality are the same as for
sequence duplication level.

Multiple QC Report
^^^^^^^^^^^^^^^^^^

You can also view a couple of reports at once using our Multiple QC
Report App. Go to the Created Files folder, select all the FastQC
reports you wish to compare, right click and select Multiple QC report.

Select from a range of QC keys to display on the plot (Total Nucleotide
Count (mate 1 and 2), GC Content % (mate 1 and 2), Number of distinct
reads (mate 1 and 2), number of reads (mate 1 and 2).

You can highlight the interesting reports and put them in a separate
folder.

<screenshot>|image42|

What are the signs that something is wrong with our data?

GC content that is far from 50% and read counts that are low compared to
other files in the dataset are ways of identifying which files
should not be used for further analysis.

Subsample Reads
^^^^^^^^^^^^^^^

Action:  used to create a random subset of raw reads.

The number of reads in the subset can be changed (default: 50,000). It
is also possible to specify a fraction of the original number of
reads.Changing the random seed value will let you create different
subsets with the same number of reads. Using the same random seed and
the same number of reads will result in identical subsets.

This application is based on
`Seqtk <https://www.google.com/url?q=https://github.com/lh3/seqtk&sa=D&ust=1480960531871000&usg=AFQjCNFaavr1xxB-goj-qyxDMaqTgd5njw>`__.

Best used when:

When the quality of the raw reads is unsatisfactory, several
preprocessing apps are available on the platform that can increase the
quality of your raw reads. Here we will walk you through each one and
give you a checklist to use when deciding which to select. After each of
the preprocessing steps, you can use the FastQC Report app again to
compare the quality pre- and post-processing (remember that in order to
do this, you need to run a different computation, this time inputting
processed data source files into the data flow).

Filter Duplicated Reads
^^^^^^^^^^^^^^^^^^^^^^^

Action: discards duplicated sequenced fragments from raw reads data. If
the sequence of two paired reads or a single read occurs multiple times
in a library, the output will include only one copy of that sequence.

The phred quality scores are created by keeping the highest score across
all identical reads for each position.

This tool is based on
`Tally <https://www.google.com/url?q=http://www.ebi.ac.uk/~stijn/reaper/tally.html&sa=D&ust=1480960531874000&usg=AFQjCNGSCUslmJdaVfxMgfxRfu6XqZ7B5w>`__.

Best used
when:

If you suspect contamination with primers, or some  other repetitive
sequence. This should be evident from Sequence duplication levels and
Overrepresented Sequences of the FastQC report. Keep in mind this app
should not be used with RNA-seq data as it will remove observed
differences in expression level.

After completing preprocessing, it’s a good idea to run a FastQC report
on the preprocessed files to see if the quality has improved.

Filter By Quality Score
^^^^^^^^^^^^^^^^^^^^^^^

Action: discards reads in a raw reads sample based on Phred33 quality
scores. You can change the minimum quality score, which is set to 20 by
default. A score of 20 means that there is a 1/100 probability that a
base was called incorrectly. In comparison, a score of 30 corresponds to
a 1/1000 probability.

You can also discard reads specifying a minimum percentage of bases to
be above the minimum quality score.

This tool is based on fastq\_quality\_filter, which is part of the
`FASTX-Toolkit <https://www.google.com/url?q=http://hannonlab.cshl.edu/fastx_toolkit/&sa=D&ust=1480960531878000&usg=AFQjCNFdpUyemH0OOfSQC7BusQ6otEFjmQ>`__.

Best used when:

If you have some low quality reads, but others are high-quality. You
should be able to tell if this is the case from the shape of the Per
sequence quality scores plot from FastQC. It may also be worth trying
this app if the per base sequence quality is low.

Trim Adaptors and Contaminants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Action: finds and trims adaptors and known contaminating sequences from
raw reads data. It is possible to specify the minimum length of trimmed
reads. Trimmed reads below the minimum length are discarded.

The app uses an internal list of sequences that can be considered as
contaminants. This list is based on the possible sequencing technologies
and platform used. For instance, it contains widely used PCR primers and
adaptors for Illumina, ABI etc. You can view the full list
`here <https://www.google.com/url?q=https://s3.amazonaws.com/bio-test-data/Genestack_adapters.txt&sa=D&ust=1480960531881000&usg=AFQjCNFst2bVH0ONqjijIMuLGMl02gh88g>`__.

This tool is based
on `fastq-mcf <https://www.google.com/url?q=https://code.google.com/p/ea-utils/wiki/FastqMcf&sa=D&ust=1480960531882000&usg=AFQjCNFm6647jAO33m4WZpSGH3Zvv6nn7A>`__,
one of the
`EA-Utils <https://www.google.com/url?q=https://code.google.com/archive/p/ea-utils/&sa=D&ust=1480960531883000&usg=AFQjCNHE_8KWOdIrCPTQ_lhxTFPRl2emWQ>`__ utilities.

Best used when:

You have irregularities in GC content, in base content at the start of
reads, duplicated reads. Since this QC app relies on sequence matching
it should be run first if used in conjunction with other QC apps

Trim Low Quality Bases
^^^^^^^^^^^^^^^^^^^^^^

Action: removes bases with a low phred33 quality score in raw reads
data. Note that a quality value of 3 means that there is a 50% chance
the base is wrong, and lower values represent even higher probabilities
of error. That’s why it can be useful to remove such bases from your
data.

Imagine you have a sequence:

Sequence:             C     G    T       A       G       A     C     T

Phred score          10   20   30      40     40      30    20   10

Error probability    .1   .01  .001  .0001 .0001 .001 .01   .1

The app will find the fragment of the read where the sum of all
probability errors will not be more than 0.01 (in our case).

In this case, the best sequence will be "TAGA" (.001\^2 + .0001\^2 =
.0022). Other fragments will have the sum of error probabilities more
than the cutoff (0.01)

Best used when:

If your per-base quality declines over the course of your reads the Trim
Low Quality Bases will select the highest quality region for each read.

This tool is based on the
`Seqtk <https://www.google.com/url?q=https://github.com/lh3/seqtk&sa=D&ust=1480960531888000&usg=AFQjCNFUVpRUIwwFfj5NUsDAZn_9jI1Mcg>`__ tool,
which uses the Phred algorithm.

Trim Reads to Fixed Length
^^^^^^^^^^^^^^^^^^^^^^^^^^

Action: trims a specific amount of bases from the extremities of all
reads in a sample.

You should specify the first base and the last base that should be
kept. For example, if you set 5 as the first base to keep and 30 as the
last base to keep, it means that the application trims all nucleotides
to the 5 position and all nucleotides from the 30th base.

This tool is based on fastx\_trimmer, which is part of the
`FASTX-Toolkit <https://www.google.com/url?q=http://hannonlab.cshl.edu/fastx_toolkit/&sa=D&ust=1480960531891000&usg=AFQjCNF1hob9o2h0-j49uKNqqhYZqPwV2g>`__.

Best used when: Trim to fix length is helpful when you want to obtain
reads of specific length (regardless of the quality).

Mapped Reads Preprocessing and QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mapped Reads QC Report
^^^^^^^^^^^^^^^^^^^^^^

In order to perform the mapped reads QC we follow a similar procedure to
the one used to generate FastQC reports. After selecting all the mapped
reads we wish to check the quality of, we can use the Mapped Reads QC
public data flow, initialize the computations, and then explore the
results. You can read more about the Mapped Reads QC Report app in the
“Explore” section of this guide.

An individual Mapped Reads QC report contains a range of mapping
statistics including:

#. Mapped reads: total number of reads which mapped to the reference
   genome;
#. Unmapped reads: total reads which failed to map to the reference
   genome;
#. Mapped reads with mapped mate: total paired reads where both mates
   were mapped;
#. Mapped reads with partially mapped mate: total paired reads where
   only one mate in the pair was mapped;
#. Mapped reads with “properly” mapped mate: total paired reads where
   both mates were mapped with the expected orientation;
#. Mapped reads with “improperly” mapped mate: total paired reads where
   one of the mates was mapped with an unexpected orientation.

→ what should we be on a lookout for here?

Large numbers of reads that are not properly mapped.|image43|

As well as two graphs.

1)Coverage by chromosome plot |image44|

This plot shows the percentage of reads covered by at least x reads. The
amount of coverage you are expecting varies with the experimental
techniques you are using. Normally you want similar coverage patterns
across all chromosomes, but this may not be the case if e.g. you are
dealing with advanced stage cancer. .

What should it look like normally?

What does it look like when data is of poor quality ( + what can we do
about it)

let's just imagine that we have a plot which shows coverage only for one
chromosome --> 1 line. On the x-axis we have the number of reads (e.g
100), on y-axis - percentage of chromosome bases covered by this number
of reads (e.g. 10%). So, it looks like we have 100-reads coverage for
10% of chromosome.

2) The insert size distribution plot

|image45|

What should it look like normally?

What does it look like when data is of poor quality ( + what can we do
about it)

This plot shows the  distribution of insert sizes. Inserts are the
distance between reads in mate pairs. Insert sizes can show e.g. indel
mutations if our data is from a specific genomic region.

Targeted Sequencing QC Report
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good to use during: Whole Exome Sequencing Analysis

Besides general quality control of mapped reads, you might also want to
assess whether the target capture has been successful, i.e. if most of
the reads actually fell on the target, if the targeted bases reached
sufficient coverage, etc. To do that, you can use Targeted Sequencing QC
Report.

By default the application allows you to compute enrichment statistics
for reads mapped only on exome. If you go to the app page, change the
value to “Both exome and target file” and select the appropriate target
annotation file, you get both exome and/or target enrichment statistics.

The following enrichment statistics are computed:

-  Number and proportion of mapped reads on target
-  Mean coverage on target with at least 2X coverage
-  Target bases with at least 2, 10, 20, 30, 40, and 50 x coverage

You can generate reports directly by choosing all of the files, right
clicking on them and choosing an appropriate app or  one of our
dedicated public data flows (Targeted Sequencing Quality Control public
data flow).

You can analyse the output for multiple reports at once using the
Multiple QC Report app.

Watch the video here: https://youtu.be/\_jHrtq\_3ya8

This application is based
on BED`tools <https://www.google.com/url?q=https://code.google.com/p/bedtools/&sa=D&ust=1480960531903000&usg=AFQjCNHFYsSqknf5t--ej96MWqPvN1jMEA>`__,
`Picard <https://www.google.com/url?q=http://broadinstitute.github.io/picard/&sa=D&ust=1480960531903000&usg=AFQjCNE7Nx1DN1A6MJS58mdncbZw3paNKQ>`__ `tools <https://www.google.com/url?q=http://broadinstitute.github.io/picard/&sa=D&ust=1480960531904000&usg=AFQjCNHQu-By-46lV8YOZ9fOB5PWZPMzGA>`__,
and SAMtools.

Mark Duplicated Mapped Reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Best used when: Duplicated reads are reads of identical sequence
composition and length, mapped to the same genomic position. Marking
duplicated reads can help speed up processing for specific apps, e.g.
the Variant Calling application, where processing additional identical
reads would lead to early PCR amplification effects (jackpotting)
contributing noise to the signal.

You can read more about Duplicated Mapped Reads in this excellent
`SeqAnswers
thread <https://www.google.com/url?q=http://seqanswers.com/forums/showthread.php?t%3D6854&sa=D&ust=1480960531906000&usg=AFQjCNEf4S1SCRUDkW22TsOHHRxjWD6Bvg>`__.

Action:goes through all reads in a Mapped Reads file, marking as
“duplicates” for paired or single reads where the orientation and the 5’
mapping coordinate are the same.

3’ coordinates are not considered due to two reasons:

#. The quality of bases generated by sequencers tends to drop down
   toward the 3’ end of a read. Thus its alignment is less reliable
   compared to the 5’ bases.
#. If reads are trimmed at 3’ low-quality bases before alignment, they
   will have different read lengths resulting in different 3’ mapping
   coordinates.

In such cases, when the distance between two mapped mates differs from
the internally estimated fragment length, including mates mapping to
different chromosomes, the application will not identify or use them but
will not fail due to inability to find the mate pair for the reads.

Marking duplicated reads can help speed up processing for specific apps,
e.g. the Variant Calling application.

This tool is based on MarkDuplicates, part of
`Picard <https://www.google.com/url?q=http://broadinstitute.github.io/picard/&sa=D&ust=1480960531908000&usg=AFQjCNFUTguXnVG8T-pHtUqYKTBvJRxSqQ>`__`  <https://www.google.com/url?q=http://broadinstitute.github.io/picard/&sa=D&ust=1480960531909000&usg=AFQjCNGhATTNeU1Rut4z-myvh2ew4jquEw>`__`tools <https://www.google.com/url?q=http://broadinstitute.github.io/picard/&sa=D&ust=1480960531909000&usg=AFQjCNGhATTNeU1Rut4z-myvh2ew4jquEw>`__.

Remove Duplicated Mapped Reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Best used when:The point of removing duplicated mapped reads is to try
to limit the influence of early PCR selection (jackpotting). Whether or
not you should remove duplicate mapped reads depends on the type of data
you have. If you are dealing with whole-genome sequencing data where
expected coverage is low and sequences are expected to be present in
similar amounts, removing duplicated reads will reduce processing time
and have little deleterious effect on analysis. If however you are
processing RNA-seq data, where the fold-variation in expression can be
up to 10^7, reads are relatively short, and your main point of interest
is the variation in expression levels, this probably isn’t the tool for
you.

You can read more about Duplicated Mapped Reads in this excellent
`SeqAnswers
thread <https://www.google.com/url?q=http://seqanswers.com/forums/showthread.php?t%3D6854&sa=D&ust=1480960531910000&usg=AFQjCNFadUu7kTaUIPWmbsa6k4trTNpkHA>`__.

Action: goes through all reads in a Mapped Reads file, marking as
“duplicates” paired or single reads where the orientation and the 5’
mapping coordinate are the same and discarding all except the “best”
copy.

3’ coordinates are not considered due to two reasons:

#. The quality of bases generated by sequencers tends to drop down
   toward the 3’ end of a read. Thus its alignment is less reliable
   compared to the 5’ bases.
#. If reads are trimmed at 3’ low-quality bases before alignment, they
   will have different read lengths resulting in different 3’ mapping
   coordinates. 

The app also takes into account interchromosomal read pairs.

In such cases, when the distance between two mapped mates differs from
the internally estimated fragment length, including mates mapping to
different chromosomes, the application  app cannot identify them but
will not fail due to inability to find the mate pair for the reads.

This tool is based on MarkDuplicates, part of the `Picard
tools <https://www.google.com/url?q=http://broadinstitute.github.io/picard/&sa=D&ust=1480960531914000&usg=AFQjCNH7a8doEzmn-2YlGchG7q_J-PR-YA>`__.

Subsample Reads
^^^^^^^^^^^^^^^

Best used when: For example, if you want to take a look at what your
final experimental results will look like, but don’t want to spend time
processing all your data right away.

Action: used to create a random subset of mapped reads.

Use subsampling ratio option to set a fraction of mapped reads you’d
like to extract (default: 50%). Changing random seed value will let you
produce different subsets with the same number of mapped reads. Using
the same random seed and the same subsampling ratio will result in
identical subsets.

This application is based on
`SAMtools <https://www.google.com/url?q=http://samtools.sourceforge.net/&sa=D&ust=1480960531916000&usg=AFQjCNFB4gFPcb-Qn-otAuuvXdgQxS-qew>`__.

Merge Mapped Reads
^^^^^^^^^^^^^^^^^^

Best used when: For example, if you have multiple replicates of the same
experiment and want to combine them before producing your final result.

Action: used to merge multiple Mapped Reads files, producing one single
output Mapped Reads file.

This application is based on
`SAMtools <https://www.google.com/url?q=http://samtools.sourceforge.net/&sa=D&ust=1480960531918000&usg=AFQjCNExyI1vxeDPJ4fJDe3oEq6iaUomvA>`__.

Merge Variants 
^^^^^^^^^^^^^^^

Best used when: Merging Genomic Variations files can be useful, when you
have, for example, one Genetic Variations file for SNPs and another one
for Indels. After their merging, the result Genetic Variations file will
separately contain information about SNPs and about Indels.

Action: allows you to merge two or more Genetic Variations files into a
single file.

This application is based on
`BCFtools <https://www.google.com/url?q=http://samtools.github.io/bcftools/bcftools.html&sa=D&ust=1480960531922000&usg=AFQjCNENqYzPwnsR_l1c-R1nKiaEfyV6JA>`__.

Concatenate Variants 
^^^^^^^^^^^^^^^^^^^^^

Best used when: Concatenation would be appropriate if you, for example,
have separate Genetic Variations files for each chromosome, and simply
wanted to join them 'end-to-end' into a single Genetic Variations file.

Action: allows you to join two or more Genetic Variations files by
concatenating them into a larger, single file.

The application always allows overlaps so that the first position at the
start of the second input will be allowed to come before the last
position of the first input. There is an option to remove duplicated
variants to make sure that there are no redundant results.

This application is based on
`BCFtools <https://www.google.com/url?q=http://samtools.github.io/bcftools/bcftools.html&sa=D&ust=1480960531926000&usg=AFQjCNFoChUsLd1NE-xsBd1GInhmlBtuHw>`__.
