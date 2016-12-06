Methylation Analysis
~~~~~~~~~~~~~~~~~~~~~

Bisulfite sequencing mapping with BSMAP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The application is based on the
`BSMAP <https://www.google.com/url?q=https://code.google.com/archive/p/bsmap/&sa=D&ust=1480960532079000&usg=AFQjCNFq0kN0aK1f-Wy2i7s1c83XjQg8IA>`__ tool
and is used to map high-throughput bisulfite reads at the level of the
whole genome.|bisulfite sequencing mapping|

Let’s talk a bit about various settings:

1)The option “Number of mismatches” lets you set the maximum number of
allowed mismatches per read. Changing this number you can affect
application runtime and percentage of mapped reads. There is an increase
in the percentage of mapped reads and in the application runtime when
increasing this value. For example, by default the read could be mapped
to the genome with no more than 5 mismatches.

2)By default, the application only reports unique hits for one mappable
read. But if your reads are mapped to multiple positions in the genome,
than you can change rule for multiple mappings to report one random
“best” mapping. This stops duplicated genome regions from being omitted
altogether.

3)Depending on the BS data generation protocol that was used to
construct the bisulfite converted library, BS reads need to  be analysed
in different ways. If  the “Lister” protocol was used, your reads will
be mapped to two forward strands. You can read more about this protocol
in Lister et al. [1]. If you Choose the “Cokus” protocol the application
will align your reads to all four strands. You can find more details
about this protocol in the original study by Cokus et al. [2].

We used this app in the Methylation Profiling Using Genestack Platform
tutorial that can be accessed
`here <https://www.google.com/url?q=https://genestack.com/tutorial/mapping-sequencing-reads-merging-techinical-replicates/&sa=D&ust=1480960532082000&usg=AFQjCNEzDwcTF01UsBP5l0UyOKnKYYJKIA>`__.

Reduced representation bisulfite sequencing mapping with BSMAP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The application is based on the
`BSMAP <https://www.google.com/url?q=https://code.google.com/archive/p/bsmap/&sa=D&ust=1480960532083000&usg=AFQjCNGrxXhzcONteprQELjc63McEx1vhg>`__ tool
and is used for mapping reduced representation bisulfite sequencing
reads to the specific digestion sites on the genome. |reduced
representation bisulfite sequencing mapping|

Let’s talk a bit about various settings:

1)You should set the “Enzyme sequence” which was recognized by by the
 restriction enzyme used to digest genomic DNA in the process of library
preparation. By default, the application uses the ‘C-CGG’ sequence which
is recognised in MspI restriction.

2)The option “Number of mismatches” lets you set the maximum number of
allowed mismatches per read. Decreasing this number you can reduce
application runtime and percentage of mapped reads. By default the
application aligns reads to the reference genome with no more than 5
mismatches.

3)By default the application only reports unique hits for one mappable
read. You can change the rule for multiple mappings to report one random
“best” mapping, if your reads are mapped to multiple positions in the
genome.

4) Choose the BS data generation protocol that was used to construct the
bisulfite converted library. If it is the Lister protocol [1], than your
reads will be mapped to two forward strands.  Reads generated using the
Cokus experimental protocol [2] will be aligned to all four strands.

 

References:

#. `Lister R, Pelizzola M, Dowen RH, Hawkins RD, Hon G, Tonti-Filippini
   J, Nery JR, Lee L, Ye Z, Ngo Q-M, Edsall L, Antosiewicz-Bourget J,
   Stewart R, Ruotti V, Millar AH, Thomson JA, Ren B, Ecker JR. “Human
   DNA methylomes at base resolution show widespread epigenomic
   differences.” <https://www.google.com/url?q=http://europepmc.org/abstract/MED/19829295&sa=D&ust=1480960532085000&usg=AFQjCNG66MkWxikJT0StWhOxW1ei40wiWQ>`__`Nature. <https://www.google.com/url?q=http://europepmc.org/abstract/MED/19829295&sa=D&ust=1480960532086000&usg=AFQjCNGpU4pwPyy6XnfY0z4BvUuolapZgw>`__`2009
   462(7271):315-22. <https://www.google.com/url?q=http://europepmc.org/abstract/MED/19829295&sa=D&ust=1480960532086000&usg=AFQjCNGpU4pwPyy6XnfY0z4BvUuolapZgw>`__
#. `Cokus SJ, Feng S, Zhang X, Chen Z, Merriman B, Haudenschild CD,
   Pradhan S, Nelson SF, Pellegrini M, Jacobsen SE. “Shotgun bisulphite
   sequencing of the Arabidopsis genome reveals DNA methylation
   patterning.” <https://www.google.com/url?q=http://europepmc.org/abstract/MED/18278030&sa=D&ust=1480960532086000&usg=AFQjCNGTngx6W4ckwk5HLaZRD1DR6crp2A>`__`Nature <https://www.google.com/url?q=http://europepmc.org/abstract/MED/18278030&sa=D&ust=1480960532087000&usg=AFQjCNF4zsutJJDSCWNBASaorGJJMoBK6Q>`__`.
   2008
   452(7184):215–219. <https://www.google.com/url?q=http://europepmc.org/abstract/MED/18278030&sa=D&ust=1480960532087000&usg=AFQjCNF4zsutJJDSCWNBASaorGJJMoBK6Q>`__

Methylation Ratio Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Good for: Methylation Profiling

Input: Mapped Reads

Action: The app allows to determine the percent methylation at each ‘C’
base in the Mapped Reads file.

Further apps to use: Genome Browser

To get results filtered by depth of coverage use “Minimum coverage”
option. By default, this value is not set. But raising it to a higher
value (e.g. 5) requires that at least five reads will cover the
position.

For paired-end mappings, you can trim from 1 to 240 fill-in nucleotides
in the DNA fragment end-repairing. By default, this “Trim N
end-repairing fill-in bases” option is switched off. For RRBS mappings,
the number of fill-in bases could be determined by the distance between
cuttings sites on forward and reverse strands. If you analyse WGBS
mappings, it’s recommended to set this number between 0~3.

Switch “Report loci with zero methylation ratios” option to report
positions with zero methylation. The application doesn’t apply this
option by default.

To combine CpG methylation ratio from both strands, set “Combine ratios
on both strands” option switched. By default, it is unchecked. If you
want to process only unique mappings, check “Only unique mappings”
option.

For paired reads, using the option “Discard discordant mappings” you can
discard all mappings where the two mates map uniquely but with
unexpected orientation, or where the distance between two mapped mates
differs from and internally estimated fragment length, including mates
mapping to different chromosomes.

Sometimes you need to remove duplicates from your Mapped Reads files.
For this purpose, use “Discard duplicated reads” option.

To ignore positions where there is a possible C/T SNPs detected, choose
“skip” value for “C/T SNPs filtering” option. If you want to correct the
methylation ratio according to the C/T SNP information estimated by the
G/A counts on reverse strand, set “correct” value. By default, the
application doesn’t consider C/T SNPs (“no-action” value).

This application is based on methratio.py script.
