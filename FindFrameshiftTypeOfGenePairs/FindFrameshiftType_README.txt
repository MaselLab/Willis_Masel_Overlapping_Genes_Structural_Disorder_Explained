The file FindFrameshiftType.py is used to determine the relative frameshift of two genes in an overlapping pair. In a same-strand overlap, one gene will always be in the +1 frame relative to its overlapping partner. This forces its overlapping partner into the +2 frame. 

This program takes in overlapping genes from a flat file in tab-delimited format with the following entries:

	[1] Sequence UID
	[2] Pair UID (one per pair)
	[3] Species
	[4] Gene Name
	[5] Overlapping Sequence (the nucleotide sequence shared by both 
	both members of the overlap
	[6] Coding Gene (Nucleotide sequence of the gene in question)
	
This file is read in twice so that the gene being read in initially can be matched up with its overlapping partner.

The program functions by finding the location of the longest shared sequence in order to determine the relative frameshift.

The main file (run by the user) is FindFrameshiftType.py. It will print the output in tab-delimited format to STDOUT, which the user may redirect to a flat file if desired.

The output from a previous run is already saved in this directory as
	Frameshifts.txt

The output entries are as follows:
	[1] Sequence UID
	[2] Pair UID
	[3] Relative Frameshift (either +1 or +2)

The original input file used to generate Frameshifts.txt is also present in this directory as
	OverlappingVirusGenes_92PairsWithSequenceData.txt