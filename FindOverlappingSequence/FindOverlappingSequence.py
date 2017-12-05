from Bio import SeqIO
from Bio.Seq import Seq
from LongestCommonSubstring import longest_common_substring

# Here we define bcolors so that we can print to different colors
# in stdout. This is so we can highlight the regions of overlap
# when we print the results
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# The purpose of this file is to take in two overlapping
# nucleotide sequences in FASTA format and to find the
# nucleotide sequence that is shared by both.

# The two FASTA files (found in the directory where this script
# is stored) are read in using the SeqIO function from the
# Bio module. The sequences are stored as sequence1 and
# sequence2 respectively.
for seq_record in SeqIO.parse('Sequence1.fasta','fasta'):
    sequence1 = seq_record.seq
    SeqRecord1 = seq_record.id

for seq_record in SeqIO.parse('Sequence2.fasta', 'fasta'):
    sequence2 = seq_record.seq
    SeqRecord2 = seq_record.id

# We import the submodule LongestCommonSubstring and use the
# function longest_common_substring, which takes in two strings
# and returns the longest shared substring. In this case, that's
# the nucleotide sequence shared by both genes.
Overlap = longest_common_substring(sequence2, sequence1)


# Here we determine the start and stop indices of the overlapping
# region within the two genes involved in the overlap.
# When we print the results, we'll highlight the regions of overlap
# in both regions so that they can be visually inspected.
start1 = str(sequence1).index(str(Overlap))
stop1 = len(str(Overlap))+start1
start2 = str(sequence2).index(str(Overlap))
stop2 = len(str(Overlap))+start2

print('Results:\n')
print('========================================\n\n')
print('Sequence 1: %s' %SeqRecord1)
print(sequence1[:start1] +bcolors.OKBLUE + sequence1 + bcolors.ENDC + sequence1[stop1:])
print('\n\nSequence 2: %s' %SeqRecord2)
print(sequence2[:start2] +bcolors.OKBLUE + sequence2[start2:stop2] + bcolors.ENDC + sequence2[stop2:])
print('\n\nOverlapping Sequence:')
print(bcolors.OKBLUE + Overlap + bcolors.ENDC)
print('\n\nLength of Overlap: %s' %len(Overlap))


