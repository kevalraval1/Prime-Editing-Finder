# print ("Enter your FASTA Sequence: ")
# FASTA = input()
# NG_Finder(FASTA)

from termcolor import colored

def reverser(FASTA, mutation):
    global newFASTA
    global newMutation
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '(': ')', ')': '('} 
    newFASTA = ''.join(complement.get(base, base) for base in reversed(FASTA))
    newMutation = ''.join(complement.get(base, base) for base in reversed(mutation))

def parsedFASTA (FASTA):
    global newString
    global position
    newString = ""
    counter = 0
    position = 0
    for char in FASTA:
        if (char != "(") and (char != ")"):
            newString = newString + char
            counter += 1
        elif (char == ")"):
            continue
        else:
            position = counter + 1
            continue
    if len(newString) < 11:
        return print("ERROR: FASTA given is too small, please enter a FASTA that is 11 or more base pairs.")
    if ((position - 7) < 0) or (position + 4 > len(newString)):
        return print("ERROR: Selected mutation site is out of bounds for editing in this FASTA file. Add more base pairs to the ends for editing.")
    if ((position - 27) < 0) or ((position + 11) > len(newString)):
        return print("ERROR: Some Spacer or Extension sequences may not be made, PAM may be out of bounds of the FASTA file.")
    return print("Successfully parsed FASTA")

def NG_Finder (newString, position):
    global listByPos
    listByPos = []
    optimal_string = ""
    answer = False
    for x in range (position-7, position+4):
        optimal_string = optimal_string + newString[x]
    counter = position - 7 
    print ("Mutation at 7th position in string: " + optimal_string)
    for x in optimal_string:
        if x == "G":
            if (counter != position-7):
                tuple = (counter, newString[counter-1]+"G")
                listByPos.append(tuple)
                answer = True
        counter += 1
    if answer == False:
        return print("No available NG PAM sites for given mutation.")
    if answer == True:
        return print(listByPos)

def spacer(newString, listByPos):
    global listOfSpacers
    listOfSpacers = []
    for tup in listByPos:
        spacerSequence = ""
        if (newString[tup[0]-21]) != "G":
            spacerSequence = spacerSequence + "G"
        for bases in range ((tup[0]-21), (tup[0]-1)):
            spacerSequence = spacerSequence + newString[bases]
        listOfSpacers.append(spacerSequence)
    return print(listOfSpacers)

def extension(newString, listByPos, mutation):
    global listOfExtensions
    listOfExtensions = []
    for tup in listByPos:
        extensionSequence = ""
        for bases in range ((tup[0]-17), (tup[0]+9)):
            if bases == position-1:
                extensionSequence = extensionSequence + mutation.lower()
                continue
            extensionSequence = extensionSequence + newString[bases]
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
        reverse_complement = ''.join(complement.get(base, base) for base in reversed(extensionSequence))
        listOfExtensions.append(reverse_complement)
    return (print(listOfExtensions))

def analysisPrinter(listByPos, listOfSpacers, listOfExtensions):
    for count in range (0, len(listByPos)):
        print ("-------------------")
        if (listByPos[count][0] + 1 == position):
            print (colored(("** PAM DESTROYED **"), 'blue'))
            print (colored("PAM " + str(count + 1) + ": " + str(listByPos[count][1]), 'blue'))
            print (colored("Position: " + str(listByPos[count][0]), 'blue'))
            print (colored("Spacer sequence Top: " + "cacc" + listOfSpacers[count] + "gtttt", 'blue'))
            print (colored("Extension sequence Top: " + "gtgc" + listOfExtensions[count], 'blue'))
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
            reverse_complement_spacer = ''.join(complement.get(base, base) for base in reversed(listOfSpacers[count]))
            reverse_complement_extension = ''.join(complement.get(base, base) for base in reversed(listOfExtensions[count]))
            print (colored("Spacer sequence Bottom: " + "ctctaaaac" + reverse_complement_spacer, 'blue'))
            print (colored("Extension sequence Bottom: " + "aaaa" + reverse_complement_extension, 'blue'))
            continue
        print ("PAM " + str(count + 1) + ": " + str(listByPos[count][1]))
        print ("Position: " + str(listByPos[count][0]))
        print ("Spacer sequence Top: " + "cacc" + listOfSpacers[count] + "gtttt")
        print ("Extension sequence Top: " + "gtgc" + listOfExtensions[count])
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
        reverse_complement_spacer = ''.join(complement.get(base, base) for base in reversed(listOfSpacers[count]))
        reverse_complement_extension = ''.join(complement.get(base, base) for base in reversed(listOfExtensions[count]))
        print ("Spacer sequence Bottom: " + "ctctaaaac" + reverse_complement_spacer)
        print ("Extension sequence Bottom: " + "aaaa" + reverse_complement_extension)
    print ("----------------------------------")

def main(FASTA, mutation):
    parsedFASTA(FASTA)
    NG_Finder(newString, position)
    spacer(newString, listByPos)
    extension(newString, listByPos, mutation)
    print ("Successfully found spacer and extension sequences for all PAMs.")
    print (colored("\n*PLUS STRAND ANALYSIS*\n", 'red'))
    analysisPrinter(listByPos, listOfSpacers, listOfExtensions)
    print (colored("\n*MINUS STRAND ANALYSIS*\n", 'red'))
    reverser(FASTA, mutation)
    parsedFASTA(newFASTA)
    NG_Finder(newString, position)
    spacer(newString, listByPos)
    extension(newString, listByPos, newMutation)
    print ("Successfully found spacer and extension sequences for all PAMs.")
    analysisPrinter(listByPos, listOfSpacers, listOfExtensions)

FASTA = "ACCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTGTATCATCATCTTAGCGACGT(G)TAGCATGCTCTATCATCATCTCATGCTCTATCATCATCTGCATACGCATGCTCTATCATCATCTGTTAAATATAT"

mutation = "T"
main(FASTA, mutation)







