import Parsing_Helper as parse

def validSubstitution(FASTA, mutation):
    enzyme = "Invalid"
    newFASTA = parse.FASTA(FASTA)[0]
    position = parse.FASTA(FASTA)[1]
    if (newFASTA[position] == "C" and mutation == "T"):
        enzyme = "CBE"
        newString = newFASTA
        strand = "+"
    elif (newFASTA[position] == "A" and mutation == "G"):
        enzyme = "ABE"
        newString = newFASTA
        strand = "+"
    elif (newFASTA[position] == "T" and mutation == "C"):
        enzyme = "ABE"
        newString = parse.reverser(FASTA, mutation)[0]
        mutation = parse.complementDict[mutation]
        position = parse.FASTA(newString)[1]
        newString = parse.FASTA(newString)[0]
        strand = "-"
    elif (newFASTA[position] == "G" and mutation == "A"):
        enzyme = "CBE"
        newString = parse.reverser(FASTA, mutation)[0]
        mutation = parse.complementDict[mutation]
        position = parse.FASTA(newString)[1]
        newString = parse.FASTA(newString)[0]
        strand = "-"
    if (enzyme == "Invalid"):
        return enzyme
    return newString, position, mutation, enzyme, strand


# FASTA = "ACCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTGTATCATCATCTTAGCGACGGT(G)TAGCATGCTCTATCATCATCTCATGCTCTATCATCATCTGCATACGCATGCTCTATCATCATCTGTTAAATATAT"
# mutation = "A"
# print (validSubstitution(FASTA, mutation))