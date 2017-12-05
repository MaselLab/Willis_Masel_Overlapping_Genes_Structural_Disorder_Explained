The main program in this directory is FindOverlappingSequence.py 
Its purpose is to find the nucleotide sequence that is shared by two overlapping genes. 

There are two plain text documents called “Sequence1.fasta” and “Sequence2.fast”. The user can copy and paste the nucleotide sequences of the overlapping genes they want to analyze into these two files in FASTA format. Sequence1.fasta should contain one of the genes in the overlapping pair, while Sequence2.fasta should contain the other. 

The program will then read them when it is run and will print the longest shared nucleotide sequence to stdout along with the length of the overlap. The two overlapping genes are also printed so that the user may visually inspect the results. Colors have been defined to print to STDOUT so that the overlapping regions are printed in a different color. This allows the user to easily see where the region of overlap is. 

The file LongestCommonSubstring.py contains the function longest_common_substring and is used by FindOverlappingSequence.py to determine the overlapping region. 