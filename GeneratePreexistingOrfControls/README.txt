The purpose of the python script in this directory is to take in coding sequences and to find all open reading frames in the +1 and +2 reading frames. If multiple ORFs are found to terminate at a common stop codon, only the longest is returned to the user. 

There are two data files that are included with this script. 

	[1] - VirusGenes_Nonoverlapping_Controls_AncestralHomologs.txt
	[2] - VirusGenes_Nonoverlapping_Controls_Complete.txt

File [1] contains the pre-overlapping ancestrally-homologous genes for generating controls and File [2] contains the non-overlapping viral genes for generating controls. 

The program prints the ORFs to the terminal in tab-delimited format for the user. The user can choose to redirect the terminal output to a text file using the > command. 