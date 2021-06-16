# print ("Enter your FASTA Sequence: ")
# FASTA = input()
# print ("Enter the numeric position of mutation you want to induce: ")
# position = input()
# print (NG_Finder(FASTA, position))

def NG_Finder (FASTA, position):
    optimal_string = ""
    answer = False
    listByPos = []
    if len(FASTA) < 11:
        return print("ERROR: FASTA given is too small, please enter a FASTA that is 11 or more base pairs.")
    if ((position - 7) < 0) or (position + 4 > len(FASTA)):
        return print("ERROR: Selected mutation site is out of bounds for editing in this FASTA file. Add more base pairs to the ends for editing.")
    for x in range (position-7, position+4):
        optimal_string = optimal_string + FASTA[x]
    counter = position - 7 
    print (optimal_string)
    for x in optimal_string:
        if x == "G":
            if (counter != position-7):
                tuple = (counter, FASTA[counter-1]+"G")
                listByPos.append(tuple)
                answer = True
        counter += 1
    if answer == False:
        return print("No available NG PAM sites for given mutation.")
    return listByPos

FASTA = "ACTAGCTACGACTACGCATACGCATGCTCTATCATCATCTGTTAAATATAT"
print (NG_Finder(FASTA,12))

