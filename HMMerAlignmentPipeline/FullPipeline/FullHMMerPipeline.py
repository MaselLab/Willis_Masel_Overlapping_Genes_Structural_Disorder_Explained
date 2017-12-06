# The purpose of this file is to read in two tab-delimited files
# that have been downloaded from a MySQL database. The first
# is a file containing overlapping genes, the second is a file
# containing non-overlapping controls.

# Format of input files: tab-delimited, no headings


# The Overlapping genes file will have the following entries:
# [1] UID
# [2] Gene Name
# [3] Overlapping Nucleotide Sequence (in the same frame as the
# overall gene)

# The controls file will have the following entries:
# [1] UID
# [2] Gene Designation (either Coding Gene, Scrambled, or
# Artificially Frameshifted)
# [3] Gene Name
# [4] Protein Sequence

# We want to run the overlapping sections of the overlapping genes,
# the frameshifted controls and the coding genes through the program
# HMMer to see if there are any homologous genes. This program will utilize
# phmmer to assign homology group IDs. It will then output two files:
# [1] The UIDs and homology group IDs for the overlapping genes
# in tab-delimited format
# [2] The same as the above, but for the non-overlapping controls.

# First, declare the two files you wish to use:

OverlappingGenesFilename = 'OverlappingVirusGenes_93PairsWithSequenceData.txt'
ControlsFilename = 'VirusGenes_Nonoverlapping_Controls_Complete.txt'

# This program allows a user to specify the UIDs of genes that are only very
# weakly homologous to the other genes in a homology group. Because the purpose
# of this program is to determine which genes are not independent datapoints when
# analyzing ISD, if a gene is only very weakly homologous, then we may be able
# to consider it an indepenent datapoint. The exceptions that are included in
# the file included with this program were determined manually using Geneious.
# They were considered after a pHMMer run. It was decided that genes that have
# a less than 50% positive match with their homologous partner
# are included in the exceptions list here.

ExceptionsFilename = 'Exceptions.txt'
PairwiseExceptionsFilename = 'pairwiseExceptions.txt'

# If the user does not wish to run the HMMer program and only wants to look at the
# output of a previous run, then NewHMMerRun can be set to False. When it's set to
# True, HMMer will be called to run on the input sequences

NewHMMerRun = True

############################################################
############################################################
############################################################
# The necessary modules and submodules are imported
import os
from pathlib import Path
import sys
import csv
import os.path
from ProteinConversion import ProteinConversion
from Fasta_File_Maker import Fasta_File_Maker
from SetConsolidator import SetConsolidator
import time 
# Before the program runs, it checks whether
# or not the program has already run and saved the output files
# to the target location. If it has, then those files are removed.
# This is so no new results are printed to the ends of preexisting files
# while loops are used instead of if statements to allow the files to be
# removed before proceeding. If the program runs too quickly and if statements
# are used, it can lead to the files not being properly removed
ControlHomGpFile = Path('../ControlsHomologyGroups.txt')
OverlapHomGpFile = Path('../OverlappingHomologyGroups.txt')
HomologyGpFile = Path('../HomologyGroups.txt')
while ControlHomGpFile.is_file():
    os.remove('%s'%ControlHomGpFile)
while OverlapHomGpFile.is_file():
    os.remove('%s'%OverlapHomGpFile)
while HomologyGpFile.is_file():
    os.remove('%s'%HomologyGpFile)

time.sleep(0.5)

# The program will load the UIDs of the genes that will be placed in their own
# homology groups so that this reader will know to exclude them while looking at
# the pHMMer output

Exceptions = []
pairwiseExceptions = []

with open('%s'%ExceptionsFilename, 'r') as f:
    reader = csv.reader(f, delimiter = ',')
    for row in reader:
        for item in row:
            UID = item
            Exceptions.append(UID)
with open('%s'%PairwiseExceptionsFilename, 'r') as f:
    reader = csv.reader(f,delimiter=',')
    for row in reader:
        pairException = [row[0],row[1]]
        sortedPairException = sorted(pairException)
        pairwiseExceptions.append(sortedPairException)
        

# ProteinConversion converts the overlapping nucleotide
# sequences into protein sequences, and then combines the overlapping genes
# file and the non-overlapping controls file into one large tab-delimited
# file

ProteinConversion(OverlappingGenesFilename,ControlsFilename)

# Fasta_File_Maker converts the large tab-delimited file into Fasta format
# so that HMMer can be used. Two separate fasta files are made. They are
# identical. One is used as the query file, the other serves as the
# database.

Fasta_File_Maker()

# phmmer is then used to compare each sequence in the query file to the
# database file. The output is saved as pHMMerOutput.txt. This is only done
# if the user specifies that a new pHMMer run is desired. If not, a previous
# pHMMer output file is used (which, so long as the user has not renamed the
# file, will be saved with the same name. 

if NewHMMerRun == True:
    os.system('phmmer --tblout pHMMerOutput.txt -E 0.0001 --domE 0.0001 --incE 0.0001 --incdomE 0.0001 FastaFile1.txt FastaFile2.txt')

# Once HMMer is run, the fasta files are removed from the directory

os.remove('./FastaFile1.txt')
os.remove('./FastaFile2.txt')

# The process of assigning homology group IDs starts here

HMMerOutputFile = 'pHMMerOutput.txt'
FileForComparison = 'AllSequencesFile.txt'

# Because of the formatting of the tblout file from hmmer, we need to specify
# unwanted characters so we can get rid of them. In this case, the unwanted
# characters are spaces, of which there are many

# The list homologousPairs will be used to store genes that hit one another
# pairwise (as output from HMMer)

# TotalList will be used to keep track of the unique IDs of the genes
# that had a hit in HMMer. This way, we can later assign unique IDs
# to all the sequences that only hit themselves 

# We also use the previously-defined Exceptions list which contains
# the UIDs of the gene pairs that are only weakly homologous (as determined
# manually using a Geneious alignment). Genes that appear in the Exceptions list
# will not be put in a homology group with another gene, and will instead get their
# own HomologyGroupID

unwanted = ['']
homologousPairs = []
TotalList = []

# Using the csv module, the file that was generated by a pHMMer run is opened
with open('%s' %HMMerOutputFile, 'r') as f:
    reader = csv.reader(f, delimiter = ' ')
    for row in reader:
        
        # An empty list is declared so we can reformat the hmmer output file
        newRow = []
        
        # we look at each element line by line in the file and if it's not
        # a blank space, then we add it to the newly reformatted row
        for element in row:
            if element not in unwanted:
                newRow.append(element)
                
        # If a line begins with a hash, then it's either a column heading
        # or notation in the hmmer output and it's ignored
        if newRow[0][0] != '#':
            
            # The first element in the new row is the UID of the query gene
            # The second is the gene that it hit
            # We combine these into a list containing both genes
            gene1 = newRow[0]
            gene2 = newRow[2]
            pair = [gene1,gene2]

            # We sort the genes so that we can compare lists (since lists
            # are order-dependent)
            pair = sorted(pair)
            
            # If the genes are the same, we discard them (since every gene
            # hits itself)
            if gene1 != gene2:
                
                # We make sure that the genes that were hit are not genes
                # that have been determined to be only weakly homologous
                # to the others in their homology group.
                if gene1 not in Exceptions:
                    if gene2 not in Exceptions:
                        if pair not in pairwiseExceptions:
                        
                            # The TotalList of gene IDs that had a hit then has
                            # these genes appended to it, so long as they passed the
                            # previous filters
                            TotalList = TotalList + pair
                            
                            # Since if A hits B, then B hits A, if we count all
                            # non-identical entries that hit, we double-count, so
                            # before we add the pair to the list of homologous pairs
                            # we check to see whether they're already in there.
                            if pair not in homologousPairs:
                                homologousPairs.append(pair)
             
# The output from SetConsolidator will be a list of disjoint
# lists containing UIDs that hit eachother in the HMMer run
homologyGroups_NoDuplicates = SetConsolidator(homologousPairs,pairwiseExceptions)
print(homologyGroups_NoDuplicates)


# We define homologyGroupID to be 0. We will sequentially add 1
# to it for each sublist in homologyGroups_NoDuplicates
homologyGroupID = 0

# We will also keep track of which UIDs have a homologyGroupID
# assigned to it. Any UIDs that do not have a homologyGroupID at
# the end of this program will be assigned a unique one
GenesWithHomologyGroupIDs = []
Genes_NoHomologs = []

# The output files are declared here. These will be where the desired
# outputs will be printed. One for the non-overlapping controls, one
# for the overlapping genes. The files will be saved one directory up
# so they're easy to locate. The output files are the following:
# [1] ControlFileFullName = tab-delimited file with the UIDs of all
# control genes and their homology group IDs (for easy upload to MySQL)
# [2] OverlappingFileFullName = same as above, except the genes are the
# overlapping genes, not the controls
# [3] HomologyGroupsOutputFullFilename = a basic file for the user to inspect
# to see what homology groups (with more than one entry) exist. The output format
# is one homology group per line with the following entries: [UID1, UID2, ..., HomologyGroupID]
save_path = '../'
ControlsFileFullName = os.path.join(save_path, 'ControlsHomologyGroups.txt')
OverlappingFileFullName = os.path.join(save_path, 'OverlappingHomologyGroups.txt')
HomologyGroupsOutputFullFilename = os.path.join(save_path, 'HomologyGroups.txt')

# Each sublist in homologyGroups_NoDuplicates counts as a homology
# group. Each sublist will have a unique homologyGroupID assigned to it
# as an integer in the last position of each sublist
for group in homologyGroups_NoDuplicates:
    homologyGroupID+=1
    group.append(homologyGroupID)

# The list of homology groups is printed to a reference file for the user to inspect
# only homology groups with more than one gene in it are included in this file
HomologyGroupsOutput = open('%s' %HomologyGroupsOutputFullFilename,'w')
for element in homologyGroups_NoDuplicates:
    HomologyGroupsOutput.write('%s\n'%str(element))
HomologyGroupsOutput.close()

# Next, we print the UIDs and their homology groups in tab-delimited
# format to output files.
for group in homologyGroups_NoDuplicates:
    for element in group:

        # We make sure that we're only selecting UIDs and not HomologyGroupIDs
        # The HomologyGroupIDs are always stored in the last position of each
        # sublist, so these are excluded.
        if group.index(element) != len(group)-1:

            # Because we had to combine the UIDs into one large file
            # to run HMMer, we added a 'c' to each of the control's
            # UIDs to distinguish them from the overlapping genes'
            # UIDs. As we're printing our output file, we use this
            # c to distinguish which gene is an overlapping gene
            # and which is a control. The c is removed for the
            # output file so the homology UIDs can be easily
            # uploaded into the MySQL database. The append function is used here so that
            # we can open the file more than once throughout this script without
            # overwriting previous entries
            if element[0] != 'c':
                output = [str(element),str(group[len(group)-1])]
                OverlappingFile = open('%s' %OverlappingFileFullName, 'a')
                OverlappingFile.write('\t'.join(output[0:]) + '\n')
                OverlappingFile.close()

            # If the UID starts with a 'c', then it is a control and it is appended
            # to a separate file
            else:
                output = [str(element[1:]),str(group[len(group)-1])]
                ControlsFile = open('%s'%ControlsFileFullName, 'a')
                ControlsFile.write('\t'.join(output[0:]) + '\n')
                ControlsFile.close()

            # Every UID that has a HomologyUID is appended to this
            # list so that any that doesn't appear will have a
            # unique Homology ID given to it
            GenesWithHomologyGroupIDs.append(element)

# finally, the large file with all the sequences is opened and
# the UIDs are looked at to see which ones need a unique homology group ID. 

with open('%s'%FileForComparison,'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        UID = row[0]
        GeneName = row[1]
        Sequence = row[2]

        if UID not in GenesWithHomologyGroupIDs:
            homologyGroupID+=1

            # The same files above are written in in the same way to get the rest
            # of the UIDs
            if UID[0] != 'c':
                output = [str(UID),str(homologyGroupID)]
                OverlappingFile = open('%s' %OverlappingFileFullName, 'a')
                OverlappingFile.write('\t'.join(output[0:]) + '\n')
                OverlappingFile.close()
            else:
                output = [ str(UID[1:]),str(homologyGroupID)]
                ControlsFile = open('%s'%ControlsFileFullName, 'a')
                ControlsFile.write('\t'.join(output[0:]) + '\n')
                ControlsFile.close()

# Finally, to ensure that our directory stays moderately organized, the large tab-delimited file
# is removed.
os.remove('./AllSequencesFile.txt')
os.remove('./pHMMerOutput.txt')



