# The purpose of this program is to take in data from a tab-delimited file,
# read in the genetic sequences, remove the start/stop codons from the sequence
# scramble the nucleotides so they are in random order. It will then scan the
# sequence to check for stop codons. If it finds any, it will select two
# nucleotides within the stop codon at random and will swap them. It will do this
# until all stop codons within the body of the scrambled sequence have been
# disassembled.

# Here the user specifies the file that 
filename = 'VirusGenes_Nonoverlapping_Controls_Complete.txt'

# The file should have entries which match the following:
# 1 - UID
# 2- Family
# 3 - Species
# 4 - Gene Designation
# 5 - Gene Name
# 6 - Coding Sequence


##########################


import random
from Bio.Seq import Seq
from Remove_Start_Stop_Function import removeStartStop
from Sequence_Scrambler_Function import sequenceScrambler
from Stop_Codon_Scrambler_Function import stopCodonScrambler
import csv

with open('%s' %filename, 'r') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:

        # The various entries in the table will be defined here for ease of use, and so if
        # the program needs to be modified, the indices will only need
        # to be changed here

        UID = row[0]
        Family = row[1]
        Species = row[2]
        GeneDesignation = row[3]
        GeneName = row[4]
        NucleotideSequence = row[5]

                              
        if GeneDesignation != 'CodingGene':
            pass

        # These next conditional statements check to see whether or not there is an
        # 'N' in the nucleotide sequence (which specifies unknown nucleotides).
        # If an N exists, it will be excluded and the program will move onto the
        # next sequence. Similarly, the next conditional statement checks the length
        # of the nucleotide sequence in question. If it is not a multiple of three
        # then the sequence is excluded.
        
        else:
            if 'N' in str(NucleotideSequence):
                pass
            else:
                if len(str(NucleotideSequence))%3 != 0:
                    pass
                
                # The following lines can be unhashed if the program gets interrupted
                # due to an error or by the user intentionally disrupting the
                # program and the user wishes to start the program again at a specific line

                else:
                    Input_nucleotide_sequence = str(NucleotideSequence)
        
                    list_containing_trimmed_sequence_start_and_stop = removeStartStop(Input_nucleotide_sequence)
                    trimmed_nucleotide_sequence = list_containing_trimmed_sequence_start_and_stop[0]
                    start_codon = list_containing_trimmed_sequence_start_and_stop[1]
                    stop_codon = list_containing_trimmed_sequence_start_and_stop[2]

                    # The function sequenceScrambler() takes in a nucleotide sequence,
                    # and returns a new string with the same nucleotides but in a random
                    # order
        
                    scrambled_sequence = sequenceScrambler(trimmed_nucleotide_sequence)
                    
                    # The function stopCodonScrambler() works to "disassemble" any
                    # stop codons it finds in the body of the sequence
        
                    scrambled_sequence_with_stop_codons_disassembled = stopCodonScrambler(scrambled_sequence)
        
                    #To complete the sequence, we take the scrambled sequence and put the
                    #start and stop codons back on the ends.
        
                    complete_scrambled_sequence = start_codon + scrambled_sequence_with_stop_codons_disassembled + stop_codon


                    
                    # The reason this final step is implemented has to do with the
                    # removeStartStop() function, which, if the nucleotide sequence
                    # it is given does not have a start and/or stop codon, it will
                    # return the string '***' instead. So that this does not show
                    # up in the final result, the replace function can be used.
        
                    final_sequence = complete_scrambled_sequence.replace("*","")

                    final_sequence_protein = Seq(final_sequence).translate()
                    final_sequence_protein_nocys = str(final_sequence_protein).replace("C","")

                    output = [UID, Family, Species, 'ScrambledByNucleotideControl', GeneName, final_sequence, final_sequence_protein, final_sequence_protein_nocys]
                    
                    print(*output, sep='\t')





