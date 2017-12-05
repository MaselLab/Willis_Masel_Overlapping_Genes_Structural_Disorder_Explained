In this directory are two different version of the same general program.

These programs find the average ISD of an overlapping region by only considering a subset of the sequence being analyzed. 

IUPred is an algorithm that goes through amino acid sequences and scores each amino acid based on its propensity to exhibit structural disorder. The average is then taken over all scores to determine the average level of ISD of a gene.

IUPred was run on each of our amino acid sequences in our database. To determine the average ISD of each overlapping region, only a subset of the ISD scores were used. Specifically, only those corresponding to the amino acids in the overlapping regions.

The purpose of the programs in each directory is to take the averages of the individual ISD scores corresponding to regions of overlap from input sequences. 

The programs in each directory are as follows:

______________________________________________

SampleProgram:

This program comes with a tab-delimited file containing sample sequences and regions of overlap. It prints the UID and average overlapping ISD score in tab-delimited format to STDOUT.


______________________________________________

OverlappingISDFinder_WithRelativeAges

This program is roughly the same as the SampleProgram, except instead of printing the individual mean ISD values for each input gene, it groups genes based on whether they've been classified as novel or ancestral and by relative reading frame and prints the means and standard errors for each group as a whole.

This program comes with overlapping genes from our work. Included are the genes whose relative ages (with respect to their overlapping partners) have been classified using phylostratigraphy, codon adaptation index, or both. 


