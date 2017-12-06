This folder contains the functional program RSCU_and_RelativeAdaptedness_Calculator.py. The purpose of this program is to take in the raw codon count values from various viral species found in the file RawCodonCount.txt, which is included in this folder, and calculates the RSCU values and Relative adaptedness values. 

Definitions for how these are calculated can be found in the paper based on this research, as well as in the source material: 

Graur, D., Sater, A. K., & Cooper, T. F. (2016). Molecular and genome evolution. Sunderland, MA: Sinauer Associates, Inc.

Running the program can be done easily from the command line. 


## Note: The user may update this program in the future to a better design using dictionaries instead of lists (this design is old from when she was initially learning to code) ##

=========================================================

The input should be in tab-delimited format with no header

The first three rows should contain the following information:

	[1] UID
	[2] Virus
	[3] Designation - Only raw codon counts should be input into this file and should have the designation "Number".
	
Following these three entries, the raw counts for amino acids should be in the following order:

	UUU = [4]
        UUC = [5]
        UUA = [6]
        UUG = [7]
        UCU = [8]
        UCC = [9]
        UCA = [10]
        UCG = [11]
        UAU = [12]
        UAC = [13]
        UAA = [14]
        UAG = [15]
        UGU = [16]
        UGC = [17]
        UGA = [18]
        UGG = [19]
        CUU = [20]
        CUC = [21]
        CUA = [22]
        CUG = [23]
        CCU = [24]
        CCC = [25]
        CCA = [26]
        CCG = [27]
        CAU = [28]
        CAC = [29]
        CAA = [30]
        CAG = [31]
        CGU = [32]
        CGC = [33]
        CGA = [34]
        CGG = [35]
        AUU = [36]
        AUC = [37]
        AUA = [38]
        AUG = [39]
        ACU = [40]
        ACC = [41]
        ACA = [42]
        ACG = [43]
        AAU = [44]
        AAC = [45]
        AAA = [46]
        AAG = [47]
        AGU = [48]
        AGC = [49]
        AGA = [50]
        AGG = [51]
        GUU = [52]
        GUC = [53]
        GUA = [54]
        GUG = [55]
        GCU = [56]
        GCC = [57]
        GCA = [58]
        GCG = [59]
        GAU = [60]
        GAC = [61]
        GAA = [62]
        GAG = [63]
        GGU = [64]
        GGC = [65]
        GGA = [66]
        GGG = [67]
