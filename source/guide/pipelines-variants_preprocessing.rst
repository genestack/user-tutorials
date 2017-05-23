Variants preprocessing
~~~~~~~~~~~~~~~~~~~~~~

While analysing variants, you also can preprocess them. Just select Genetic
Variations file and click on "Preprocess" section to see what applications
are available for you.

Merge variants
++++++++++++++

Merging variants can be useful, when you have, for example, one Genetic
Variations file for SNPs and another one for Indels. After their merging, the
result Genetic Variations file will separately contain information about SNPs
and about Indels.

**Action**: to merge two or more Genetic Variations files into a single file.

.. image:: images/merge_variants.png

This application is based on the `BCFtools`_.

.. _BCFtools: http://samtools.github.io/bcftools/bcftools.html

Concatenate variants
++++++++++++++++++++

Concatenation would be appropriate if you, for example, have separate Genetic
Variations files for each chromosome, and simply wanted to join them
"end-to-end" into a single Genetic Variations file.

**Action**: to join two or more Genetic Variations files by concatenating them
into a larger, single file.

.. image:: images/concatenate_variants.png

The application always allows overlaps so that the first position at the start
of the second input will be allowed to come before the last position of the
first input.

1. The **Remove duplicated variants** option checks for the duplicated variants and
   makes sure that there are no redundant results. (default: unchecked)

The application is based on the `BCFtools`_.

.. _BCFtools: http://samtools.github.io/bcftools/bcftools.html