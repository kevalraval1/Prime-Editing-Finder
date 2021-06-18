# print ("Enter your FASTA Sequence: ")
# FASTA = input()
# NG_Finder(FASTA)

def NG_Finder (FASTA):
    counter = 0
    position = 0
    new_string = ""
    for char in FASTA:
        if (char != "(") and (char != ")"):
            new_string = new_string + char
            counter += 1
        elif (char == ")"):
            continue
        else:
            position = counter + 1
            continue
    optimal_string = ""
    answer = False
    listByPos = []
    if len(new_string) < 11:
        return print("ERROR: FASTA given is too small, please enter a FASTA that is 11 or more base pairs.")
    if ((position - 7) < 0) or (position + 4 > len(new_string)):
        return print("ERROR: Selected mutation site is out of bounds for editing in this FASTA file. Add more base pairs to the ends for editing.")
    for x in range (position-7, position+4):
        optimal_string = optimal_string + new_string[x]
    counter = position - 7 
    print (optimal_string)
    for x in optimal_string:
        if x == "G":
            if (counter != position-7):
                tuple = (counter, new_string[counter-1]+"G")
                listByPos.append(tuple)
                answer = True
        counter += 1
    if answer == False:
        return print("No available NG PAM sites for given mutation.")
    if answer == True:
        return print(listByPos)

FASTA = "ACTAGCTACGA(C)TACGCATACGCATGCTCTATCATCATCTGTTAAATATAT"
FASTA1 = "TAAATATAT"
FASTA2 = "ACTAGCTACGACTACGCATACGCATGCTCTATCATCATCTGTTAAATA(T)AT"
FASTA3 = "ACTA(G)CTACGACTACGCATACGCATGCTCTATCATCATCTGTTAAATATAT"
FASTA4 = "AAAAAAAAAA(A)AAAAAAAAAAAAAAAAAAAAA"

NG_Finder(FASTA)
NG_Finder(FASTA1)
NG_Finder(FASTA2)
NG_Finder(FASTA3)
NG_Finder(FASTA4)




