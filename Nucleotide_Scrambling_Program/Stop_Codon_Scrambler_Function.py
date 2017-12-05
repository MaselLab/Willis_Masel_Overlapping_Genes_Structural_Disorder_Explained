import random

# The purpose of this subroutine is to look for stop codons
# within the body of a nucleotide sequence. If one is found, then
# one of the nucleotides within the body of the stop codon and another nucleotide
# within the body of the sequence are selected.
# We then swap them. This process is repeated until all stop codons are eliminated




def stopCodonScrambler(scrambled_and_trimmed_nucleotide_sequence):
    true = 1
    false = 0
    stop_codons_in_body = true


    while stop_codons_in_body == true:
        
        # count is used to determine whether or not a stop codon has been found
        # If nothing is added to the count
        # then no stop codon was found and at the end of this loop
        # stop_codons_in_body is set to false. 
        count = 0

        # 'n' will be used to look at the indices of our sequence.
        for n in range(len(scrambled_and_trimmed_nucleotide_sequence)-1):

            
            # For there to be a stop codon, the codon must start with T, which means
            # there must be a T at an index which is a multiple of three. If this is true
            # we move on to further investigations of the codon. If not, then we continue
            # scanning the sequence
            if scrambled_and_trimmed_nucleotide_sequence[n] == 'T' and n%3 ==0:
                
                # The next if statements look to see what follows the T at the start of a codon. If any
                # of the three if statements below are true, a stop codon exists. If not, there is no
                # stop codon.
                
                if scrambled_and_trimmed_nucleotide_sequence[n+1] == 'A' and scrambled_and_trimmed_nucleotide_sequence[n+2] == 'A':
                    # TAA is a stop codon, so if this is true, one has been found.
                    # A random nucleotide within the codon is selected by its index
                    # A second nucleotide within the body of the sequence needs to be found (at random) to swap it with.
                    random_integer_to_select_nucleotide = random.randint(n,n+2)
                    second_random_integer_to_select_nucleotide = random.randint(0,len(scrambled_and_trimmed_nucleotide_sequence)-1)
                    # A while loop is used to make sure that the two randomly selected indices are not the same
                    while random_integer_to_select_nucleotide == second_random_integer_to_select_nucleotide:
                        second_random_integer_to_select_nucleotide = random.randint(0,len(scrambled_and_trimmed_nucleotide_sequence)-1)

                    # First the relative sizes of the indices are determined for ease of swapping items within
                    # a string environment
                    if random_integer_to_select_nucleotide < second_random_integer_to_select_nucleotide:
                        smaller_Index = random_integer_to_select_nucleotide
                        larger_Index = second_random_integer_to_select_nucleotide
                    else:
                        smaller_Index = second_random_integer_to_select_nucleotide
                        larger_Index = random_integer_to_select_nucleotide

                    # The nucleotides associated with the indices are determined
                    smaller_Index_nucleotide = scrambled_and_trimmed_nucleotide_sequence[smaller_Index]
                    larger_Index_nucleotide = scrambled_and_trimmed_nucleotide_sequence[larger_Index]

                    # and their orders are swapped within the body of the string
                    scrambled_and_trimmed_nucleotide_sequence = scrambled_and_trimmed_nucleotide_sequence[:smaller_Index]+larger_Index_nucleotide + scrambled_and_trimmed_nucleotide_sequence[smaller_Index+1 : larger_Index] + smaller_Index_nucleotide + scrambled_and_trimmed_nucleotide_sequence[larger_Index+1:]

                    # One is added to the count indicating that codon was found
                    count += 1

                # The same procedure is carried out for the other two possible stop codons
                # CODON TAG
                else:
                    if scrambled_and_trimmed_nucleotide_sequence[n+1] == 'A' and scrambled_and_trimmed_nucleotide_sequence[n+2] == 'G':
                        random_integer_to_select_nucleotide = random.randint(n,n+2)
                        second_random_integer_to_select_nucleotide = random.randint(0,len(scrambled_and_trimmed_nucleotide_sequence)-1)
                        while random_integer_to_select_nucleotide == second_random_integer_to_select_nucleotide:
                            second_random_integer_to_select_nucleotide = random.randint(0,len(scrambled_and_trimmed_nucleotide_sequence)-1)
                        if random_integer_to_select_nucleotide < second_random_integer_to_select_nucleotide:
                            smaller_Index = random_integer_to_select_nucleotide
                            larger_Index = second_random_integer_to_select_nucleotide
                        else:
                            smaller_Index = second_random_integer_to_select_nucleotide
                            larger_Index = random_integer_to_select_nucleotide
                        smaller_Index_nucleotide = scrambled_and_trimmed_nucleotide_sequence[smaller_Index]
                        larger_Index_nucleotide = scrambled_and_trimmed_nucleotide_sequence[larger_Index]
                        scrambled_and_trimmed_nucleotide_sequence = scrambled_and_trimmed_nucleotide_sequence[:smaller_Index]+larger_Index_nucleotide + scrambled_and_trimmed_nucleotide_sequence[smaller_Index+1 : larger_Index] + smaller_Index_nucleotide + scrambled_and_trimmed_nucleotide_sequence[larger_Index+1:]
                        count = count + 1
                
                    # CODON TGA
                    else:
                        if scrambled_and_trimmed_nucleotide_sequence[n+1] == 'G' and scrambled_and_trimmed_nucleotide_sequence[n+2] == 'A':
                            random_integer_to_select_nucleotide = random.randint(n,n+2)
                            second_random_integer_to_select_nucleotide = random.randint(0,len(scrambled_and_trimmed_nucleotide_sequence)-1)
                            while random_integer_to_select_nucleotide == second_random_integer_to_select_nucleotide:
                                second_random_integer_to_select_nucleotide = random.randint(0,len(scrambled_and_trimmed_nucleotide_sequence)-1)
                            if random_integer_to_select_nucleotide < second_random_integer_to_select_nucleotide:
                                smaller_Index = random_integer_to_select_nucleotide
                                larger_Index = second_random_integer_to_select_nucleotide
                            else:
                                smaller_Index = second_random_integer_to_select_nucleotide
                                larger_Index = random_integer_to_select_nucleotide
                            smaller_Index_nucleotide = scrambled_and_trimmed_nucleotide_sequence[smaller_Index]
                            larger_Index_nucleotide = scrambled_and_trimmed_nucleotide_sequence[larger_Index]
                            scrambled_and_trimmed_nucleotide_sequence = scrambled_and_trimmed_nucleotide_sequence[:smaller_Index]+larger_Index_nucleotide + scrambled_and_trimmed_nucleotide_sequence[smaller_Index+1 : larger_Index] + smaller_Index_nucleotide + scrambled_and_trimmed_nucleotide_sequence[larger_Index+1:]
                            count = count + 1
                        else:
                            pass

        #If the count is zero, then no stop codons were found
        if count == 0:
            stop_codons_in_body =false

    return scrambled_and_trimmed_nucleotide_sequence
