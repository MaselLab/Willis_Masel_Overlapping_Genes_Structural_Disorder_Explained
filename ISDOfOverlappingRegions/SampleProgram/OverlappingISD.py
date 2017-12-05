# The purpose of this script is to take the ISD values output by IUPred and
# to find the mean ISD values of the overlapping protein section of two overlapping
# Genes

import csv
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats


# Below, the user may specify a p-value cutoff and the name of the file which will be used

filename = 'ArtificialSequences.txt'


# First the file with the relevant information is read in in tab-delimited format

with open('%s' %filename, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:

        # The entries are assigned variable names.
        UID = row[0]
        NoCysOverlappingSequence = row[1]
        NoCysAminoAcidSequence = row[2]
        ISDScores = row[3]


        #The StartIndex finds where the overlapping protein sequence begins
        #in the whole protein sequence. The StopIndex indicates where the
        #overlapping protein sequence ends in the whole protein sequence
        StartIndex = NoCysAminoAcidSequence.index(NoCysOverlappingSequence)
        StopIndex = len(NoCysOverlappingSequence) + StartIndex
        

        #The ISDList is the list of raw data output by IUPred, and includes
        #each individual score for each amino acid in the full amino acid
        #sequence. The list ISD then is the ISDList spliced so that only the raw
        #data for the overlapping sections are included
        ISDList = ISDScores.split(',')
        ISD = ISDList[StartIndex : StopIndex]
        
        #The list ISD then has all of its values converted from the standard
        #string format into float format so that the mean can be found
        OverlappingISD = [float(i) for i in ISD]
        
        #MeanISD is then the mean value of all of the raw scores for each
        #of the amino acids in the overlapping section
        MeanISD = np.mean(OverlappingISD)
        
        
        output= [ UID,MeanISD]
        print(*output, sep = '\t')
            


