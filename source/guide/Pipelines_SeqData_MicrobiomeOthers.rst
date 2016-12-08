Microbiome Analysis
~~~~~~~~~~~~~~~~~~~

Microbiome Analysis with QIIME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Microbiome Analysis

Input: Targeted Microbiome Sequencing Data

Action: Application report on the recorded microbial species and the
percentage composition of the sample.

Output: Clinical and Research reports

Further apps to use: None

The application can be used to create a  clinical or research
microbiology report with abundance plots and microbiological diversity
metrics.

Metrics include:

– counts for every taxonomic unit (how many reads match to a given
group)

– alpha diversity (within each sample, how rich the sample is e.g.
number of taxa identified)

– beta diversity (difference between a pair of samples) (heterogeneity
of samples)

Microbiome analysis is performed using QIIME (open source tool), using
“Greengenes 13.8” (for bacteria) and UNITE (for fungi) reference
databases to estimate the taxonomic composition of the microbial
communities.

The OTU picking step is performed using an open-reference procedure with
uclust. Taxonomy assignment is done using the blast algorithm. Any reads
that failed mapping to reference sequence are excluded. Tools used:
qiime 1.9.1

Others
~~~~~~

Genome Browser
^^^^^^^^^^^^^^

Good for: Variant Calling, Methylation Profiling, Whole Exome Sequencing
Analysis, Whole Genome Sequencing Analysis and many more.
Used to: View and explore different types of genomic data: mapped reads,
genetic variants, methylation ratios and others.

There are several tracks that can be visualized in Genome Browser:

-  Reference genome: displays annotated genes, transcripts, and their
   coordinates;

-  Coverage: represents the sequencing reads coverage for Mapped Reads

[Genome browser whole exome sequencing tutorial]];|image54|

-  Variation: shows genetic variants (SNPs, insertions etc.), their
   exact position on the genome, average mapping quality and raw read
   depth;

[genome browser methylation tutorial];|image55|

-  DNA methylation ratio: reflects the proportion of methylated and
   unmethylated cytosine residues.

Also you can manage tracks: add new ones, hide or delete them. When
manipulating with multiple tracks you can use the tracks mentioned above
to create Combined track or Formula track. On the combined track several
tracks are imposed and shown together, thereby comparing coverage for
different samples. Or you can apply some basic mathematical operations
and create formulas based on your genomic data, for example, quantify
average value between values corresponding to different samples. The
results of the computations will be shown on the formula track.

Moreover, each track can be personalised by changing its properties
(track color, normalized values, show only SNPs, maximum and minimum
values to be shown on a track, etc.). Use “Edit” button to change
properties for multiple tracks at once.

Genome Browser allows you to browse either a specific genomic position
(search by coordinates) or a specific feature (search by feature name).
You can navigate through the data to find a feature of interest or
explore regions surrounding the feature, and zoom in to nucleotide
resolution. The found feature can be marked with sticky notes (Shift +
click on the position on the track). When you share the Genome Browser
page with your collaborators, sticky notes will  help to focus their
attention on your findings.

You can see the Genome browser in action in this blog
`post <https://www.google.com/url?q=https://genestack.com/blog/2015/05/28/navigation-in-genestack-genome-browser/&sa=D&ust=1480960532101000&usg=AFQjCNE3r6NoPVzIZm7LOxsU0h9eopDvDQ>`__.
