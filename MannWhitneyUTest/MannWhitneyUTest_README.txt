MANN WHITNEY U TEST README

The purpose of this collection of programs is to look at the shared sequence of two overlapping genes in the two separate reading frames and to run a Mann Whitney U-Test on the collection of relative adaptedness values to get a p-value. This will determine whether the codon usage in the two reading frames is significantly different. 

Two example files are given. These are to be used together:
	[1] - TurnipYellowsVirus_88
	[2] - TurnipYellowsVirus_88_Designations

[1] Gives the relative adaptedness value for each codon appearing in the reading frame of the two genes being compared. Both the ancestral and novel gene values are combined in this file. 

[2] Partitions file [1] into either ancestral or novel by designating a value with either a 0 or a 1 to distinguish them. 

These two files are read into the script MannWhitneyUTest.r. 

Other files are included in this working directory so that the user can reproduce all data files used in this test.

	