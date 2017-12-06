This directory contains four data tables that were used in our research. They are in comma-delimited format and are as follows:

	1 - OverlappingVirusGenes_92PairsWithSequenceData.csv
	2 - OverlappingVirusGenes_PValue_0.035_Length_200.csv
	3 - VirusGenes_Nonoverlapping_Controls_Complete.csv
	4 - VirusCodonBias.csv

1 - Contains 184 genes, translating to 92 overlapping gene pairs, from viral genomes. These are the 184 genes pulled from the literature for which we could find sequence data.

2 - Contains the 47 overlapping gene pairs (94 genes) that we were able to classify the relative ages of (either using phylostratigraphy, codon bias, or both). 

3 - Contains 150 non-overlapping genes from the same viral species as the 92 overlapping pairs. Also contains 450 controls generated from these 150 genes (3 per control). The controls are as follows:
	- +1 artificially-frameshifted control
	- +2 artificially-frameshifted control
	- Randomly-scrambled nucleotide control

4 - Contains the raw codon count, RSCU values and relative adaptedness values for each codon for the viral species used in this work. These values were used to determine the codon bias for genes by comparing the codon usage of the gene with the codon usage of the genome as a whole. 

====================================================

For data tables 1 and 2, the headings are as follows. 

	[1] - Sequence UID
	[2] - Pair UID - one UID per overlapping pair
	[3] - HomologyGroupID - homology group ID assigned by single-link clustering using both HMMer and Geneious (>= 50% sequence similarity from Geneious alignment with BLOSUM62 matrix and similarity threshold 1)
	[4] - HomologyGroupID_HMMerOnly - similar to above, except the homology group ID was determined using HMMer only, not Geneious
	[5] - VirusType - single or double stranded DNA or RNA
	[6] - OrderOrGroup - Viral order or group, if available
	[7] - Family
	[8] - Genus
	[9] - Species
	[10] - PValue - Statistical-significance of CAI value (determined by running a Mann-Whitney U Test on the relative adaptedness values)
	[11] - GenomeAccessionNumber
	[12] - ProteinAccessionNumber 
	[13] - GeneName
	[14] - OverlappingGene
	[15] - GeneDesignation - Either Ancestral or Novel
	[16] - Verification - A poorly-titled entry which indicates whether the gene's relative age was determined in the literature using phylogenetic analysis. If it was, it is listed as "verified", if the relative ages were classified using CAI values only, it is listed as Tentative. 
	[17] - OverlapType - Either Internal (one gene contained entirely within the other), or Terminal (the 3' end of the upstream gene overlapping with the 5' end of the downstream)
	[18] - LengthOfOverlap - In nucleotides
	[19] - OverlappingSequence - The shared nucleotide sequence of two overlapping genes
	[20] - SequenceUsedForCodonAnalysis - The overlapping sequence trimmed so it's in the same reading frame as the listed gene
	[21] - CodingGene 
	[22] - ProteinSequence
	[23] - NoCysProteinSequence
	[24] - FrameshiftFromAncestral - The relative reading frame of the gene with respect to the ancestral gene of the pair. If the gene is the ancestral, then it is 0.
	[25] - FrameshiftFromNovel - same as above, except novel 
	[26] - FrameshiftRelativeToOverlappingGene - The relative reading frame of the gene with respect to its overlapping partner. This is always either 1 or 2. 
	[27] - NoCysIUPredMeanISD - The mean of all of the IUPred raw scores
	[28] - OverlappingOnly_NoCysIUPredMeanISD - The mean of the IUPred raw scores from the overlapping regions only (see: ISDOfOverlappingRegions in repository)
	[29] - NoCysOrigDataIUPredAminoAcidISDScores - The raw IUPred scores
	[30] - NoCysProteinSeqLength

====================================================

Data table 3 has fewer entries. The ones that are there are similar to those in tables 1 and 2. The ones that differ are:

	OldSequenceUID - This is different from SequenceUID. There exists one unique SequenceUID per entry in the datable. In contrast, there is a unique OldSequenceUID per protein-coding genes (and is the SequenceUID), and is also used for the three nongenic controls created from that gene. 

	GeneDesignation - either CodingGene, ScrambledByNucleotideControl, ShiftedReadingFrameControlPlusOne, or ShiftedReadingFrameControlPlusTwo. 

	GCPercentContent

====================================================

Data table 4's entries should be self-explanatory. The Designation column tells the reader whether they are looking at the raw codon count, the relative adaptedness values or the RSCU values. The description of relative adaptedness values and RSCU values can be found in our paper, as well as the source we used when calculating them: 

Graur, D., Sater, A. K., & Cooper, T. F. (2016). Molecular and genome evolution. Sunderland, MA: Sinauer Associates, Inc.
