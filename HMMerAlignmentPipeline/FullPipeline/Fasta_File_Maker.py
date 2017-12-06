from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import csv
import sys

# The purpose of this file is to take in data from a tab-delimited file
# and to print it in fasta format so that it can be used by HMMer.

def Fasta_File_Maker():

    orig_stdout = sys.stdout 

    filename = 'AllSequencesFile.txt'
    #filename = 'OverlappingVirusGenes_93PairsWithSequenceData.txt'
    with open('%s' %filename, 'r') as f:
        reader = csv.reader(f, delimiter = '\t')
        for row in reader:

            UID = row[0]
            GeneName = row[1]
            Sequence = row[2]
            record=SeqRecord(Seq(Sequence, generic_protein), id = '%s' %UID, name = '%s' %GeneName)
            sys.stdout=open('FastaFile1.txt','a')
            print(record.format('fasta'))
            sys.stdout = open('FastaFile2.txt','a')
            print(record.format('fasta'))
    sys.stdout.close()
    sys.stdout=orig_stdout

#Fasta_File_Maker()



            
