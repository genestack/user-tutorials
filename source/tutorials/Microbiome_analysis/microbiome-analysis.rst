Lets now identify microbial species and their abundances in the microbiome
samples.

.. Video - Microbiome Analysis step
.. raw:: html

    <iframe width="640" height="360" src="" frameborder="0" allowfullscreen="1">&nbsp;</iframe>↵



Sample count: 8
Reads per sample: 377,823
Clustered reads per sample: 313,487


"Group samples by" option lets you change 


Display absolute OTU counts
Hide unidentified OTUs
Hide partially identified OTUs



The plot displays the relative abundance of OTUs at a highest taxonomic
resolution: genus (L6 tab) and species (L7 tab). You can change resolution to
the L2 level to see what phyla are the most abundant across the samples.

For example. our results shows that, at low taxonomic resolution (L2 tab), the
composition of microbial communities is similar between samples.
Bacteroidetes (8,30–86.73%), Firmicutes (1.46–650,49%) and Proteobacteria
(1,38–64,96%) are the most abundant phyla across most of the samples, followed
by Actinobacteria (0,02-15,14%) and Fusobacteria (0-10,33%).

.. image:: images/Microbiome_analysis_L2_level_plot.png

You can see these results in the table as well. Click on "Total" header in the
table to get the most abundant phyla across the samples:

.. image:: images/Microbiome_analysis_L2_level_table.png

Our findings are consistent with the paper results:

.. image:: images/Microbiome_analysis_L2_level_table_paper.png

At a higher taxonomic resolution (L7 tab), more than 820 microorganisms were
identified. The koala eye and the koala gastrointestinal tract are characterized
by distinct microbial communities. To measure the similarity between the bacterial
communities we used principal component analysis (PCA) based on Pearson
correlation coefficients:

.. image:: images/Microbiome_PCA.png

You may change the PCA type in the upper-left corner of the plot and try other
statistics to quantify the compositional dissimilarity between samples:
bray_curtis, abund_jaccard, euclidean, binary_pearson, binary_jaccard.

However, in comparison to the paper, authors used principal coordinate analysis
(PCoA) to show the similarity between the koala microbial communities:

.. image:: images/Microbiome_PCoA.png



What is the difference between the PCA and PCoA? Both
methods are used to visualize the data, but different mathematical approaches
are applied to the data. The purpose of PCA is to represent as much of the
variation as possible in the first few axes. For this, first, the variables are
centred to have a mean of zero and then the axes are rotated (and re-scaled).
In the end, we have two axes: the first one contains as much variation as
possible, the second one contains as much of the remaining variation as
possible, etc.



Congratulations! You've just gone through the entire tutorial!
