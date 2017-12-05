# The purpose of this script is to take the ISD values output by IUPred and
# to find the mean ISD values of the overlapping protein section of two overlapping
# Genes

import csv
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats


# Below, the user may specify a p-value cutoff and the name of the file which will be used

filename = 'OverlappingVirusGenes_PValue_0.035_Length_200.txt'
pValueCutoff = 0.035

# A series of empty lists are defined where the relevant individual ISD scores  will be stored
# based on whether the ISD scores are from a gene designated as ancestral or novel, and
# based on which reading frame that gene is in with respect to its overlapping partner

AncestralISDScores = []
NovelISDScores = []

AncestralPValuesPlusOne = []
AncestralPValuesPlusTwo = []

NovelPValuesPlusOne = []
NovelPValuesPlusTwo = []

AncestralPhyloPlusOne = []
AncestralPhyloPlusTwo = []

NovelPhyloPlusOne = []
NovelPhyloPlusTwo = []

NovelPlusOne = []
NovelPlusTwo = []

AncestralPlusOne = []
AncestralPlusTwo = []

UnverifiedISD = []
UnverifiedOverlappingISD = []


# First the file with the relevant information is read in in tab-delimited format

with open('%s' %filename, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:

        # The entries are assigned variable names
        UID = row[0]
        Virus = row[1]
        pValue = row[2]
        GeneName = row[3]
        GeneDesignation = row[4]
        Verification = row[5]

        ####### IMPORTANT:######
        
        # the overlapping sequence that needs to be used
        # is the overlapping sequence that has been trimmed so that it
        # is in the same reading frame as the gene which is being analyzed.
        # This is usually designated as SequenceUsedForCodonAnalysis.
        # If this sequence is not used, then when the overlapping section
        # is translated into an amino acid sequence, the sequence cannot be
        # found in the overall protein sequence of the gene being analyzed
        # and an error will result
        
        OverlappingSequence = row[6]
        CodingSequence = row[7]
        AminoAcidSequence = row[8]
        NoCysAminoAcidSequence = row[9]
        FrameshiftFromAncestral = row[10]
        FrameshiftFromNovel = row[11]
        ISDScores = row[12]

        if OverlappingSequence == '':
            pass
        else:

        
            # The nucleotide sequence which is shared by two overlapping
            # sequences is translated into an amino acid sequence.
            #The cysteines are then removed (because of the ambiguity of how
            #they will affect the ISD in IUPRed, so the true length of the protein
            #sequence as it was analyzed is the length of the cysteine-less sequence)
            #Note: overlapping sequence is the one used for codon analysis. Not exact
            #Overlapping sequence, but rather one in the right reading frame
            
            NoCysAminoAcidSequence = NoCysAminoAcidSequence.replace('*','')
            OverlappingProteinSequence = Seq(OverlappingSequence).translate()
            NoCysOverlappingProteinSeq = str(OverlappingProteinSequence).replace('C','').replace('*','')
        
            #The StartIndex finds where the overlapping protein sequence begins
            #in the whole protein sequence. The StopIndex indicates where the
            #overlapping protein sequence ends in the whole protein sequence
            
            StartIndex = NoCysAminoAcidSequence.index(NoCysOverlappingProteinSeq)
            StopIndex = len(NoCysOverlappingProteinSeq) + StartIndex


            #The ISDList is the list of raw data output by IUPred, and includes
            #each individual score for each amino acid in the full amino acid
            #sequence. The list ISD then is the ISDList spliced so that only the raw
            #data for the overlapping sections are included
            
            ISDList = ISDScores.split(',')
            ISD = ISDList[StartIndex : StopIndex]

            #The list ISD then has all of its values converted from the standard
            #string format into float format so that the mean can be found

            TotalISD = [float(i) for i in ISDList]
            OverlappingISD = [float(i) for i in ISD]
        
            #MeanISD is then the mean value of all of the raw scores for each
            #of the amino acids in the overlapping section

            MeanISD = np.mean(OverlappingISD)
            TotalISD = np.mean(TotalISD)

            #If a p-value (corresponding to whether or not codon adaptation data
            #was available for the gene or not. if no value, then it was not,
            #if it has a p-value, then it was analyzed for codon degeneracy) isn't
            #available for the entry, then it is either sorted as a gene whose
            #age has been classified, and sorts out whether it's a +1 or +2 frameshift
            #
            #If it hasn't been verified by phylogenetics, then it's unclassified, and
            #its ISD is sorted into a list for unclassified genes
            #
            #If it has a p-value, it determines whether the p-value is under the
            #cutoff of .02. If it is, then it's classified and is sorted accordingly
            #if it's p-value is over the cutoff, it's unverified and is again sorted
            #accordingly
            #
            #Finally, the last catagory performs the same function, seeing if the gene
            #was verified by phylogenetics. This was due to a later discovery of a paper
            #with information on gene age, and so idiosynchratic labeling.
  
            if pValue == '':
                if Verification == 'Verified':
                    if FrameshiftFromAncestral == '+1':
                        NovelPhyloPlusOne.append(MeanISD)
                        NovelPlusOne.append(MeanISD)
                    else:
                        if FrameshiftFromAncestral == '+2':
                            NovelPhyloPlusTwo.append(MeanISD)
                            NovelPlusTwo.append(MeanISD)
                    if FrameshiftFromNovel == '+1':
                        AncestralPhyloPlusOne.append(MeanISD)
                        AncestralPlusOne.append(MeanISD)
                    else:
                        if FrameshiftFromNovel == '+2':
                            AncestralPhyloPlusTwo.append(MeanISD)
                            AncestralPlusTwo.append(MeanISD)
                else:
                    UnverifiedOverlappingISD.append(MeanISD)
                    UnverifiedISD.append(TotalISD)
                
            else:
                pValue = float(pValue)
                if pValue <= pValueCutoff:
                    if FrameshiftFromAncestral == '+1':
                        NovelPValuesPlusOne.append(MeanISD)
                        NovelPlusOne.append(MeanISD)
                    else:
                        if FrameshiftFromAncestral == '+2':
                            NovelPValuesPlusTwo.append(MeanISD)
                            NovelPlusTwo.append(MeanISD)
                    if FrameshiftFromNovel == '+1':
                        AncestralPValuesPlusOne.append(MeanISD)
                        AncestralPlusOne.append(MeanISD)
                    else:
                        if FrameshiftFromNovel == '+2':
                            AncestralPValuesPlusTwo.append(MeanISD)
                            AncestralPlusTwo.append(MeanISD)
                else:
                    if Verification == 'Verified':
                        if FrameshiftFromAncestral == '+1':
                            NovelPhyloPlusOne.append(MeanISD)
                            NovelPlusOne.append(MeanISD)
                        else:
                            if FrameshiftFromAncestral == '+2':
                                NovelPhyloPlusTwo.append(MeanISD)
                                NovelPlusTwo.append(MeanISD)
                        if FrameshiftFromNovel == '+1':
                            AncestralPhyloPlusOne.append(MeanISD)
                            AncestralPlusOne.append(MeanISD)
                        else:
                            if FrameshiftFromNovel == '+2':
                                AncestralPhyloPlusTwo.append(MeanISD)
                                AncestralPlusTwo.append(MeanISD)
                    else:
                        
                        UnverifiedOverlappingISD.append(MeanISD)
                        UnverifiedISD.append(TotalISD)
                    
                    
                    
                    

NovelTotal = NovelPlusOne + NovelPlusTwo
AncestralTotal = AncestralPlusOne + AncestralPlusTwo

print('\n\nOverall Overlapping Mean ISDs\n')
print('________________________________\n')
print('Relative Age - Relative Reading Frame - (Total Count): Mean +- Standard Error\n')
print('________________________________\n')
print('Ancestral - +1 - (%s): %s +- %s' %(len(AncestralPlusOne),np.mean(AncestralPlusOne), stats.sem(AncestralPlusOne)))
print('Ancestral - +2 - (%s): %s +- %s' %(len(AncestralPlusTwo),np.mean(AncestralPlusTwo), stats.sem(AncestralPlusTwo)))
print('Novel - +1 - (%s): %s +- %s' %(len(NovelPlusOne),np.mean(NovelPlusOne), stats.sem(NovelPlusOne)))
print('Novel - +2 - (%s): %s +- %s\n\n' %(len(NovelPlusTwo),np.mean(NovelPlusTwo), stats.sem(NovelPlusTwo)))



print('Overlapping Novel Gene Regions (%s): %s +- %s' %(len(NovelTotal),np.mean(NovelTotal),stats.sem(NovelTotal)))

print('Overlapping Ancestral Gene Regions (%s): %s +- %s' %(len(AncestralTotal),np.mean(AncestralTotal),stats.sem(AncestralTotal)))


