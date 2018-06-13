import csv
import sys

# The purpose of this file is to act as an accessory to the R file
# MannWhitneyUTest.r which needs two separate files to run.

# This program takes in a file in tab-delimited format with the following designations:
# From the first file associated with the gene that's been designated as novel
# in the overlapping gene pair
# Entry 1: UID associated with the virus
# Entry 2: UID which identifies overlapping genes as a pair
# Entry 2: Name of the viral species
# Entry 3: Name of gene
# Entry 4: List of relative adaptedness values from the gene designated as the
# novel gene in a gene pair
#
# From the second file associated with the gene that's been designated as ancestral
# in the overlapping gene pair
# The file is scanned to find the same viral species as that found in the previous file
# and the relative adaptedness values are collected as a list corresponding to the
# tentatively-designated ancestral gene in the pair

TentativeNovelGenesFile = 'Novel.txt'
TentativeAncestralGenesFile = 'Ancestral.txt'

with open('%s' %TentativeNovelGenesFile, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        UID = row[0]
        PairUID = row[1]
        Species = row[2]
        Gene = row[3]
        NovelValues = row[4:]
        with open('%s' %TentativeAncestralGenesFile, 'r') as g:
            reader2 = csv.reader(g, delimiter = '\t')
            for row in reader2:
                if row[1] == PairUID:
                    AncestralValues = row[4:]
                else:
                    pass
    
        # A blank list is then created and will be used to designate the relative
        # adaptedness values as ancestral with a zero, or as novel with a one.
        # The relative adaptedness values are then concatenated to create one list
        # so they can be analyzed as a single dataset
        
        Designations = []
        for item in AncestralValues:
            Designations.append(0)
        for item in NovelValues:
            Designations.append(1)
        Values = AncestralValues + NovelValues



        # The relative adaptedness values are then printed in a form that can be read
        # by the R script
        orig_stdout = sys.stdout
        f = open('%s_%s' %(Species,PairUID),'w')
        sys.stdout = f
        print(*Values, sep='\t')
        sys.stdout=orig_stdout
        f.close()

        # A second file is then printed out in the same format as the previous file
        # used to differentiate the files. Both are required to run the
        # Mann-Whitney U-Test
        orig_stdout = sys.stdout
        f = open('%s_%s_Designations' %(Species,PairUID),'w')
        sys.stdout = f
        print(*Designations, sep='\t')
        sys.stdout=orig_stdout
        f.close()

