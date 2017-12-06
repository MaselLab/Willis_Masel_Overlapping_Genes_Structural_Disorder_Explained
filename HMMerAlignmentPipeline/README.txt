This directory contains the files and programs that are needed to preliminarily assign homology group IDs to all overlapping and non-overlapping genes in our dataset.

In order for these programs to run, the user needs to have HMMer installed on their machine: http://hmmer.org/

The main program that will be primarily what the user will interact with is the following: FullHMMerPipeline.py

It utilizes the submodules ProteinConversion.py, SetConsolidator.py and Fasta_File_Maker.py to function. 

The plain text file Exceptions.txt and pairwiseExceptions.txt are used to exclude UIDs from homology groups after an initial HMMer run. In this work, all genes were fed through HMMer to get preliminary homology group IDs. The predicted homology groups were then aligned using Geneious. Any gene that had less than 50% sequence similarity (using a Blosum62 matrix and 1 as the similarity threshold) with all other genes in its homology group was added to the list of exceptions. 

The programs are located inside the director Full Pipeline, while the results from the runs are stored in this directory for easy viewing. Three files are generated from each run. 
	
	1 - ControlsHomologyGroups.txt - These are all the non-overlapping genes and their controls, printed in tab-delimited format with the entries 
		[1] - UID
		[2] - HomologyGroupUID

	2 - OverlappingHomologyGroups.txt - These are in the same format as 1, except instead of non-overlapping genes and controls, these are the overlapping genes. 1 & 2 are stored as separate flat files because overlapping genes and controls were stored as separate MySQL data tables.

	3 - HomologyGroups.txt - This file is for user inspection. It contains each homology group as a bracketed group. The UIDs in that group are shown in quotations while the homology group UIDs are without quotations and are the last entry in each list. Control UIDs are printed with a c in front of them to distinguish them from overlapping gene UIDs. 