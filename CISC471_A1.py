"""
Assignment 1 - CISC 471
Written by Emily Gauvreau
20074874

I have included a comment for each function for it's purpose and have only included comments
for lines like list comprehension where their function might not be immediately seen

Additionally, I have made the assumption that the sequences will be provided in uppercase 
letters as per the examples in the assignment description. As well that all sequence letters will be one of ACGT
"""
import itertools

nucleotides = "ACGT"

# Question 1: Calculate the percentage of GC present in the sequence
def calculateGC(DNAsequence):
    size = len(DNAsequence)
    totalGC = 0
    for i in DNAsequence:
        if i == "G" or i == "C":
            totalGC += 1

    percent = round(totalGC/size * 100, 2)
    return percent

# Question 2: Return the reverse compliment of a strand of DNA
def reverseCompliment(DNAsequence):
    # Switched the letters to they corresponding pair (A&T, G&C)
    compliment = ""
    for i in DNAsequence:
        if i == "G":
            compliment += "C"
        elif i == "C":
            compliment += "G"
        elif i == "A":
            compliment += "T"
        else:
            compliment += "A"
        
    # Use string splicing to reverse the string starting at the back
    compliment = compliment[::-1]

    return compliment

# def rnaToProtein(DNAsequence):
    
#     codons = {
#         "UUC": "F", "UUC": "F", "UUA": "L", "UUG": "L",
#         "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
#         "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
#         "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
#         "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
#         "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
#         "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
#         "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
#         "UAU": "Y", "UAC": "Y", "CAU": "H", "CAC": "H",
#         "AAU": "N", "AAC": "N", "GAU": "D", "GAC": "D",
#         "CAA": "Q", "CAG": "Q", "AAA": "K", "AAG": "K",
#         "GAA": "E", "GAG": "E", "UGU": "C", "UGC": "C",
#         "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
#         "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
#         "GGU": "G", "GCG": "G", "GGA": "G", "GGG": "G",
#         "UGG": "W", "UAA": "Stop", "UAG": "Stop", "UGA": "Stop",
#     }
    
#     if len(DNAsequence) % 3 == 0:
#         rnaStrand = ""
#         for i in DNAsequence:
#             if i == "G":
#                 rnaStrand += "C"
#             elif i == "C":
#                 rnaStrand += "G"
#             elif i == "A":
#                 rnaStrand += "U"
#             else:
#                 rnaStrand += "A"

#         protein = ""
#         for i in range(0, len(rnaStrand), 3):
#             singleCodon = rnaStrand[i:i + 3]
#             protein += codons[singleCodon]
        
#         return (rnaStrand, protein)
    

# Question 4: Returns the Hamming Distance (# of mismatches) between two sequences
def hammingDistance(sequence1, sequence2):

    if len(sequence1) == len(sequence2):
        count = 0
        for i in range(len(sequence1)):
            if sequence1[i] != sequence2[i]:
                count += 1
        return count
    else:
        return "Lengths are not equal cannot calculate distance."


# Question 5: Return the number of times a pattern is found within a sequence
def count(DNAsequence, pattern):

    count = 0
    patternLen = len(pattern) 
    for i in range(0, len(DNAsequence)-(patternLen-1)):
        if DNAsequence[i:i+patternLen] == pattern: 
            count += 1
    return count

# Question 6: Finds all k-mers and the number of times they appear, returns the max occurence 
def mostFrequentKMer(DNAsequence, numK):
    
    kmers = {}
    for i in range(0, len(DNAsequence)-(numK-1)):
    # for each element starting a zero until k length before the end
    # ensures that the search does not go out of index range
        pattern = DNAsequence[i:i+numK]
        if pattern in kmers:
            count = kmers[pattern][0]
            indexes = kmers[pattern][1]
            indexes.append(i)
            kmers[pattern] = (count+1, indexes)

        else:
            kmers[pattern] = (1, [i])

    # determine the max value by parsing the first index of the dictionary values.
    maxValue = max(kmers.values(), key=lambda x : x[0])[0]

    # create a new dictionary {key: indices} if the key is a max occurrence
    patterns = {k: v[1] for k, v in kmers.items() if v[0] == maxValue}

    return patterns


# HELPER: given a length k the function will return all possible combinations of ACGT that are k length alphabetically
def generateAllKmers(numK):
    kmers = map(''.join, itertools.product(nucleotides, repeat=numK))
    return list(sorted(kmers))

# Question 7: Returns the list of neighbours that are at most d mismatches away
def neighbors(pattern, d):

    size = len(pattern)
    if d <= size:
        neighborList = []
        kmers = generateAllKmers(size)
        for var in kmers:
            if hammingDistance(pattern, var) <= d:
                neighborList.append(var)
        return neighborList
    else: 
        return "The value of d cannot be longer than the pattern itself."

# Question 8: Returns the list of all kmers that are only d mismatches away from the pattern with their index
def approxOccurrence(DNAsequence, pattern, d):
    
    kmers = []
    sizeSequence = len(DNAsequence)
    sizePattern = len(pattern)
    for i in range(0, sizeSequence-(sizePattern-1)):
        alter = DNAsequence[i:i+sizePattern]
        if hammingDistance(pattern, alter) <= d:
            kmers.append([alter, i])
    return kmers

# Question 9: Returns the list of k-mers that appear (L, t)-clump
def clumps(genome, k, L, t):

    uniqueKmers = []
    sizeGenome = len(genome)
    for i in range(0, sizeGenome-(L-1)):
        kmers = {}
        currentClump = genome[i:i+L]
        for j in range(0, L-(k-1)):
            pattern = currentClump[j:j+k]
            if pattern in kmers:
                kmers[pattern] += 1
            else:
                kmers[pattern] = 1
        patterns = [k for k, v in kmers.items() if v >= t]
        uniqueKmers.extend(patterns)
    
    return set(uniqueKmers)

# Question 10:
def sharedPairs(k, sequence1, sequence2):
    pass

# Question 12: Returns the frequency array for all k-mers in a sequence
# Uses functions previously defined to obtain values.
def freqTable(DNAsequence, k):

    table = [[],[],[]]
    kmers = generateAllKmers(k) # already sorted
    frequency = []
    for i in kmers:
        frequency.append(count(DNAsequence, i))

    table[0] = kmers
    table[1] = list(range(0, len(kmers)))
    table[2] = frequency

    return table


def main():
    # Test Functions with the data taken from either the assignment or the corresponding Rosalind Question

    # GCContent = calculateGC("ATGCTTAGGACT")
    # print("GC Content Percentage: ", GCContent)

    # compliment = reverseCompliment("ATGCTTAGGACT")
    # print("Reverse compliment: ", compliment)

    # rnaStrand, protein = rnaToProtein("ATGCTTAGGACT")
    # print("Transcribed RNA Strand: ", rnaStrand)
    # print("Protein sequence: ", protein)

    # hamD = hammingDistance("ATGCTTAGGACT", "ATGCTTACCACT")
    # print("Hamming Distance: ", hamD)

    # patternCount = count("ACAACTATGCATACTATCGGGAACTATCCT","ACTAT")
    # print("Total k-mer pattern: ", patternCount)

    # mostFreq = mostFrequentKMer("ACAACTATGCATACTATCGGGAACTATCCT", 3)
    # print("Most Freq K-Mer with indices: ", mostFreq)

    # numNeigh = neighbors("ACTAT", 4)
    # print("Total List of Neighbours: ", numNeigh)

    # occurrences = approxOccurrence("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC", "ATTCTGGA", 3)
    # print("Approximate Occurrence: ", occurrences)

    # clumpVal = clumps("CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC", 5, 75, 4)
    # print("K-Mer Clumps: ", clumpVal)

    # frequencyArray = freqTable("ACGCGGCTCTGAAA", 2)
    # print(frequencyArray)

if __name__ == "__main__":
    main()