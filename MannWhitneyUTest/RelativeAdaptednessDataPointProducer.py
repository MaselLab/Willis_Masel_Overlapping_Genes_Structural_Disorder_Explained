import csv
import os
from pathlib import Path
import time


###### This is an edit of the program CAIFinder.py #######


# The purpose of this file is to go through a file containing gene segments
# which are under investigation, to count the number of codons in each segment,
# and to then print out the relative adaptedness numbers corresponding with
# each codon in that virus species


# Here we open the file where our gene segments are stored for analysis and
# label the UID, Species, Gene name and Sequence so that we can call on them as variables
# later in the script. The gene sequence that is needed is the section of overlapping gene
# from each gene pair, trimmed so that it's in the reading frame of the gene
# it's associated with

#GeneRelativeAdaptedness = []

GenesFile = 'OverlappingVirusGenes.txt'
RelativeAdaptednessFile = 'VirusCodonBias.txt'
AncestralOutputFile = Path('./Ancestral.txt')
NovelOutputFile = Path('./Novel.txt')

if AncestralOutputFile.is_file():
    os.remove('%s'%AncestralOutputFile)
if NovelOutputFile.is_file():
    os.remove('%s'%NovelOutputFile)
time.sleep(0.5)


#########################################################################################################
#########################################################################################################
#########################################################################################################

with open('%s' %GenesFile, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        UID = row[0]
        PairUID = row[1]
        Species = row[2]
        GeneName = row[3]
        GeneDesignation = row[4]
        Sequence = row[5]

            
        # Dictionaries are created where we can call on the data which pertain to the codons
        # under analysis. Each dictionary is of the same form and contains a number which
        # is used as a reference as the first entry in every sublist. This number can
        # be used later to navigate between dictionaries to group data pertaining to a particular
        # codon.
        
        CodonIndex = [[0,'TTT','TTC'],[1,'TTA','TTG','CTT','CTC','CTA','CTG'],[2,'ATT','ATC','ATA'],[3,'ATG'],[4,'GTT','GTC','GTA','GTG'],[5,'TCT','TCC','TCA','TCG'],[6,'CCT','CCC','CCA','CCG'],[7,'ACT','ACC','ACA','ACG'],[8,'GCT','GCC','GCA','GCG'],[9,'TAT','TAC'],[10,'TAA','TAG','TGA'],[11,'CAT','CAC'],[12,'CAA','CAG'],[13,'AAT','AAC'],[14,'AAA','AAG'],[15,'GAT','GAC'],[16,'GAA','GAG'],[17,'TGT','TGC'],[18,'TGG'],[19,'CGT','CGC','CGA','CGG','AGA','AGG'],[20,'AGT','AGC'],[21,'GGT','GGC','GGA','GGG']]

        CodonCount = [[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]


        RelativeAdaptedness =[[0,0,0],[1,0,0,0,0,0,0],[2,0,0,0],[3,0],[4,0,0,0,0],[5,0,0,0,0],[6,0,0,0,0],[7,0,0,0,0],[8,0,0,0,0],[9,0,0],[10,0,0,0],[11,0,0],[12,0,0],[13,0,0],[14,0,0],[15,0,0],[16,0,0],[17,0,0],[18,0],[19,0,0,0,0,0,0],[20,0,0],[21,0,0,0,0]]


        Index = ['Phe','Leu','Ile','Met','Val','Ser','Pro','Thr','Ala','Tyr','Stop','His','Gln','Asn','Lys','Asp','Glu','Cys','Trp','Arg','Ser','Gly']
        
        # The goal of the program is to now go through the sequences which need analyzing
        # and to count the number of codons of each type. First a counting variable is
        # defined "TotalCodonCount",
        TotalCodonCount = 0
        Sequence1Length = int(len(Sequence))
        Sequence1CodonNumber = int(Sequence1Length)
        
        # A for loop is created which goes through the nucleotide sequence under investigation
        # in groups of 3. The codon is cross-referenced with the dictionary above and 1 is added
        # to the value in the dictionary CodonCount in the same location. In this manner, each
        # type of codon is counted.
        for n in range(0,Sequence1CodonNumber,3):
            TotalCodonCount += 1
            codon = Sequence[n:n+3]
            for entry in CodonIndex:
                if codon in entry:
                    if codon == 'ATG':
                        pass
                    else:
                        if codon == 'TGG':
                            pass
                        else:
                            external_index = CodonIndex.index(entry)
                            internal_index = entry.index(codon)
                            CodonCount[external_index][internal_index]=CodonCount[external_index][internal_index]+1
                            
        # Now the program opens the file which contains all the relative adaptedness
        # data points which were generated from the online database numbers. The loop goes through
        # and finds the relative adaptedness values which pertain to the virus species in question
        # and labels all of the values with the codon which they correspond to.
        
        with open('%s' %RelativeAdaptednessFile, 'r') as f:
            reader = csv.reader(f, delimiter = '\t')
            for row in reader:
                UID = row[0]
                VirusName = row[1]
                Designation = row[2]
                if VirusName != Species:
                    pass
                else:
                    if Designation != 'RelativeAdapted':
                        pass
                    else:
                        UUU = float(row[3])
                        UUC = float(row[4])
                        UUA = float(row[5])
                        UUG = float(row[6])
                        UCU = float(row[7])
                        UCC = float(row[8])
                        UCA = float(row[9])
                        UCG = float(row[10])
                        UAU = float(row[11])
                        UAC = float(row[12])
                        UAA = float(row[13])
                        UAG = float(row[14])
                        UGU = float(row[15])
                        UGC = float(row[16])
                        UGA = float(row[17])
                        UGG = float(row[18])
                        CUU = float(row[19])
                        CUC = float(row[20])
                        CUA = float(row[21])
                        CUG = float(row[22])
                        CCU = float(row[23])
                        CCC = float(row[24])
                        CCA = float(row[25])
                        CCG = float(row[26])
                        CAU = float(row[27])
                        CAC = float(row[28])
                        CAA = float(row[29])
                        CAG = float(row[30])
                        CGU = float(row[31])
                        CGC = float(row[32])
                        CGA = float(row[33])
                        CGG = float(row[34])
                        AUU = float(row[35])
                        AUC = float(row[36])
                        AUA = float(row[37])
                        AUG = float(row[38])
                        ACU = float(row[39])
                        ACC = float(row[40])
                        ACA = float(row[41])
                        ACG = float(row[42])
                        AAU = float(row[43])
                        AAC = float(row[44])
                        AAA = float(row[45])
                        AAG = float(row[46])
                        AGU = float(row[47])
                        AGC = float(row[48])
                        AGA = float(row[49])
                        AGG = float(row[50])
                        GUU = float(row[51])
                        GUC = float(row[52])
                        GUA = float(row[53])
                        GUG = float(row[54])
                        GCU = float(row[55])
                        GCC = float(row[56])
                        GCA = float(row[57])
                        GCG = float(row[58])
                        GAU = float(row[59])
                        GAC = float(row[60])
                        GAA = float(row[61])
                        GAG = float(row[62])
                        GGU = float(row[63])
                        GGC = float(row[64])
                        GGA = float(row[65])
                        GGG = float(row[66])
                        
                        # Once all the of the relative adaptedness values are found and
                        # labeled, they are placed into a dictionary in the same form as
                        # the dictionaries above.
                        RelativeAdaptedness = [[0,UUU,UUC],[1,UUA,UUG,CUU,CUC,CUA,CUG],[2,AUU,AUC,AUA],[3,AUG],[4,GUU,GUC,GUA,GUG],[5,UCU,UCC,UCA,UCG],[6,CCU,CCC,CCA,CCG],[7,ACU,ACC,ACA,ACG],[8,GCU,GCC,GCA,GCG],[9,UAU,UAC],[10,UAA,UAG,UGA],[11,CAU,CAC],[12,CAA,CAG],[13,AAU,AAC],[14,AAA,AAG],[15,GAU,GAC],[16,GAA,GAG],[17,UGU,UGC],[18,UGG],[19,CGU,CGC,CGA,CGG,AGA,AGG],[20,AGU,AGC],[21,GGU,GGC,GGA,GGG]]
                    
                        RelativeAdaptednessOverall = []
                        GeneRelativeAdaptedness = []

                        numberCount = 0
                        for entry in RelativeAdaptedness:
                            CountIndex = entry[0]
                            entry = entry[1:]
                            for item in entry:
                                numberOfCodons = CodonCount[CountIndex][entry.index(item)+1]
                                if item == 0.0:
                                    item = 0.5
                                if numberOfCodons == 0:
                                    pass
                                else:
                                    for integer in range(numberOfCodons):
                                        GeneRelativeAdaptedness.append(item)

                        output_preliminary = [UID,PairUID,Species,GeneName]
                        GeneRelativeAdaptedness = [str(x) for x in GeneRelativeAdaptedness]
                        output = output_preliminary + GeneRelativeAdaptedness

                        if GeneDesignation == 'Novel':
                            NovelFileOutput = open('%s'%NovelOutputFile, 'a')
                            NovelFileOutput.write('\t'.join(output[0:]) + '\n')
                            NovelFileOutput.close()
                        elif GeneDesignation == 'Ancestral':
                            AncestralFileOutput = open('%s'%AncestralOutputFile, 'a')
                            AncestralFileOutput.write('\t'.join(output[0:]) + '\n')
                            AncestralFileOutput.close()


