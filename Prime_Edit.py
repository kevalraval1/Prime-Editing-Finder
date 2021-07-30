from tkinter import *
import sys, os, regex as re

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
    position = 0 # position is the index of the mutation in the new string
    for char in FASTA:
        if (char != "(") and (char != ")"):
            newString += char
            counter += 1
        elif (char == ")"):
            continue
        else:
            position = counter # Imagine mutation at position 1
            continue
    if len(newString) < 11:
        return print("ERROR: FASTA given is too small, please enter a FASTA that is 11 or more base pairs.")
    if ((position - 7) < 0) or (position + 4 > len(newString)):
        return print("ERROR: Selected mutation site is out of bounds for editing in this FASTA file. Add more base pairs to the ends for editing.")
    if ((position - 27) < 0) or ((position + 11) > len(newString)):
        return print("ERROR: Some Spacer or Extension sequences may not be made, PAM may be out of bounds of the FASTA file.")
    return print("Successfully parsed FASTA")

def regexCompiler (input):
    compile = ""
    for x in input:
        if x == 'A':
            compile += 'A'
        elif x == 'C':
            compile += 'C'
        elif x == 'G':
            compile += 'G'
        elif (x == 'T') or (x == 'U'):
            compile += 'T'
        elif x == 'M':
            compile += '[A|C]'
        elif x == 'R':
            compile += '[A|G]'
        elif x == 'W':
            compile += '[A|T]'
        elif x == 'S':
            compile += '[C|G]'
        elif x == 'Y':
            compile += '[C|T]'
        elif x == 'K':
            compile += '[G|T]'
        elif x == 'V':
            compile += '[A|C|G]'
        elif x == 'H':
            compile += '[A|C|T]'
        elif x == 'D':
            compile += '[A|G|T]'
        elif x == 'B':
            compile += '[C|G|T]'
        elif x == 'N':
            compile += '[G|A|T|C]'
    print(compile)
    return re.compile(compile)

def pamDestroyed(inputPAM, mutation, listByPos):
    searchString = ""
    for bases in range (position - len(inputPAM) + 1, position + len(inputPAM)):
        searchString += newString[bases]
    set1 = set()
    for match in pattern.finditer(searchString, overlapped=True):
        set1.add(match.start() + position - len(inputPAM) + 1)
    searchString2 = ""
    for bases in range (position - len(inputPAM) + 1, position + len(inputPAM)):
        if bases == position:
            searchString2 += mutation.upper()
        else:
            searchString2 += newString[bases]
    set2 = set()
    for match in pattern.finditer(searchString2, overlapped=True):
        set2.add(match.start() + position - len(inputPAM) + 1)
    checkSet = set2.symmetric_difference(set1)
    for elements in checkSet:
        if elements not in set1: # PAM created
            tempTuple = (elements, searchString2[elements - (position - len(inputPAM) + 1):(len(inputPAM) + elements - (position - len(inputPAM) + 1))], 2)
            listByPos.append(tempTuple)
        if elements not in set2: # PAM destroyed
            counter = 0
            for x in listByPos:
                if x[0] == elements:
                    tempList = list(x)
                    tempList[2] = 1
                    x = tuple(tempList)
                    listByPos[counter] = x
                counter += 1

def sequenceFinder (newString, position, inputPAM):
    global listByPos
    global pattern
    listByPos = []
    optimal_string = ""
    for x in range (position - (4 + len(inputPAM)), position + (3 + len(inputPAM))): # position = index of mutation, 4 bases from end of first possible PAM
        optimal_string += newString[x]                                               # 2 bases from the start of last possible PAM, + 1 for range function (ends 1 before)
    print ("Mutation at optimal position in string: " + optimal_string)
    pattern = regexCompiler(inputPAM)
    for match in pattern.finditer(optimal_string, overlapped=True):
        tempTuple = (match.start() + (position - (4 + len(inputPAM))), match.group(), 0) # List by Pos = [(start index in newstring, PAM sequence, 0/1/2 where 0 = PAM untouched 1 = destroyed 2 = created)]
        listByPos.append(tempTuple)
    if len(listByPos) == 0:
        return print("No available PAM sites for given mutation.")
    else:
        return print(listByPos)

def spacer(newString, listByPos):
    global listOfSpacers
    listOfSpacers = []
    for tup in listByPos:
        spacerSequence = ""
        if (newString[tup[0] - 20]) != "G":
            spacerSequence += "G"
        for bases in range ((tup[0] - 20), (tup[0])):
            spacerSequence += newString[bases]
        listOfSpacers.append(spacerSequence) # List of Spacers = [spacer 1 for PAM 1, spacer 2 for PAM 2]
    return print(listOfSpacers)

def extension(newString, listByPos, mutation):
    global listOfExtensions
    listOfExtensions = []
    for tup in listByPos:
        extensionSequence = ""
        for bases in range ((tup[0] - (PBSlength + 3)), (tup[0] + 10)): # PBS Length always ends 3 BP before PAM, RTT always starts 3 BP before PAM and is 13 BP long
            if bases == position:
                extensionSequence += mutation.lower()
                continue
            extensionSequence += newString[bases]
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
        reverse_complement = ''.join(complement.get(base, base) for base in reversed(extensionSequence))
        listOfExtensions.append(reverse_complement) # List of Extensions = [Extension 1 for PAM 1, Extension 2 for PAM 2]
    return (print(listOfExtensions)) 

def ngRNA(newString, mutation):
    global listOfngRNA
    listOfngRNA = []
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '(': ')', ')': '('}
    tempList = list(newString)
    tempList[position] = mutation.lower()
    newString = "".join(tempList)
    reversedString = ''.join(complement.get(base, base) for base in reversed(newString))
    newPosition = len(newString) - position - 1
    optimalString = ""
    for base in range (newPosition + 1, newPosition + 8):
        optimalString += reversedString[base]
    match = re.search(pattern, optimalString)
    if (match != None):
        ngRNA1 = ""
        if reversedString[newPosition + 1 + match.start() - 20] != "G":
            ngRNA1 += "G"
        for bases in range (newPosition + 1 + match.start() - 20, newPosition + 1 + match.start()):
            ngRNA1 += reversedString[bases]
        ngRNA2 = ''.join(complement.get(base, base) for base in reversed(ngRNA1))
        tuple = (ngRNA1, ngRNA2)
        listOfngRNA.append(tuple) # List of ngRNA = [(top ngRNA, bot ngRNA), (top ngRNA, bot ngRNA)] plus strand then minus strand
    else:
        tuple = ("N/A", "N/A")
        listOfngRNA.append(tuple)

def analysisPrinter(listByPos, listOfSpacers, listOfExtensions, file1):
    if (len(listByPos) == 0):
        file1.write("NO PAM SITES AVAILABLE IN GIVEN INPUT")
        return
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
    printList = [f"ngRNA Top: cacc{listOfngRNA[0][0]}", 
                f"ngRNA Bottom: aaac{listOfngRNA[0][1]}\n", 
                ]
    for count in range(0, len(listByPos)):
        reverse_complement_spacer = ''.join(complement.get(base, base) for base in reversed(listOfSpacers[count]))
        reverse_complement_extension = ''.join(complement.get(base, base) for base in reversed(listOfExtensions[count]))
        tempList = ["-------------------", 
                    f"PAM {str(count + 1)}: {str(listByPos[count][1])}", 
                    f"Position: {str(listByPos[count][0] + 1)}",
                    f"Spacer sequence Top: cacc{listOfSpacers[count]}gtttt", 
                    f"Spacer sequence Bottom: ctctaaaac{reverse_complement_spacer}",
                    f"Extension sequence Top: gtgc{listOfExtensions[count]}",
                    f"Extension sequence Bottom: aaaa{reverse_complement_extension}\n"]
        if ((listByPos[count][2]) == 1):
            tempList[0] = "-------------------\n** PAM DESTROYED **"
            tempString = "\n".join(tempList)
            printList.insert(2, tempString)
            continue
        if ((listByPos[count][2] == 2)):
            tempList[0] = "-------------------\n** PAM CREATED **"
            tempString = "\n".join(tempList)
            printList.insert(3, tempString)
            continue
        tempString = "\n".join(tempList)
        printList.append(tempString)
    printList.append("----------------------------------")
    addString = "\n".join(printList)
    file1.write(addString)

def main(FASTA, mutation, filename, inputPAM, PBS):
    global PBSlength
    PBSlength = int(PBS)
    #For testing code:
    completename = os.path.join(os.path.dirname("Prime_Edit.py"), (filename + ".txt"))
    #For testing executable:
    # completename = os.path.join(os.path.dirname(sys.executable), (filename + ".txt"))
    parsedFASTA(FASTA)
    sequenceFinder(newString, position, inputPAM)
    pamDestroyed(inputPAM, mutation, listByPos)
    spacer(newString, listByPos)
    extension(newString, listByPos, mutation)
    ngRNA(newString, mutation)
    print ("Successfully found spacer and extension sequences for all PAMs.")
    file1 = open(completename, "w")
    file1.write("\n*PLUS STRAND ANALYSIS*\n\n")
    analysisPrinter(listByPos, listOfSpacers, listOfExtensions, file1)
    file1.write("\n\n*MINUS STRAND ANALYSIS*\n\n")
    reverser(FASTA, mutation)
    parsedFASTA(newFASTA)
    sequenceFinder(newString, position, inputPAM)
    pamDestroyed(inputPAM, newMutation, listByPos)
    spacer(newString, listByPos)
    extension(newString, listByPos, newMutation)
    ngRNA(newString, mutation)
    print ("Successfully found spacer and extension sequences for all PAMs.")
    analysisPrinter(listByPos, listOfSpacers, listOfExtensions, file1)
    print ("Exiting...")
    file1.close()
    sys.exit()

# FASTA = "ACCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTGTATCATCATCTTAGCGACGGT(G)TAGCATGCTCTATCATCATCTCATGCTCTATCATCATCTGCATACGCATGCTCTATCATCATCTGTTAAATATAT"
