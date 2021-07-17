''' 
TODO : Make PAMs variable
- try 2 BP PAMs with ambiguity codes
    - REGEX characters search
- try adjusting windows depending on PAM sizes
    - global variable of PAM size, adjust indexes on that size
'''
from tkinter import *
import sys, os, regex as re

window = Tk()
window.title("Prime Editing: spG Peg Design")

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
            for x in listByPos:
                if x[0] == elements:
                    tempList = list(x)
                    tempList[2] = 1
                    x = tuple(tempList)

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
        return print("No available NG PAM sites for given mutation.")
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
    for base in range (newPosition, newPosition + 7):
        optimalString += reversedString[base]
    match = re.search(pattern, optimalString)
    if (match != None):
        ngRNA1 = ""
        if reversedString[match.start() - 20] != "G":
            ngRNA1 += "G"
        for bases in range (match.start() - 20, match.start()):
            ngRNA1 += reversedString[bases]
        ngRNA2 = ''.join(complement.get(base, base) for base in reversed(ngRNA1))
        tuple = (ngRNA1, ngRNA2)
        listOfngRNA.append(tuple) # List of ngRNA = [(top ngRNA, bot ngRNA), (top ngRNA, bot ngRNA)] plus strand then minus strand
    else:
        tuple = ("N/A", "N/A")
        listOfngRNA.append(tuple)

def analysisPrinter(listByPos, listOfSpacers, listOfExtensions, file1):
    addString = ""
    addString += ("ngRNA Top: " + "cacc" + listOfngRNA[0][0] + "\n")
    addString += ("ngRNA Bottom: " + "aaac" + listOfngRNA[0][1] + "\n\n")
    for count in range (0, len(listByPos)):
        if (listByPos[count][2] == 1):
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
            reverse_complement_spacer = ''.join(complement.get(base, base) for base in reversed(listOfSpacers[count]))
            reverse_complement_extension = ''.join(complement.get(base, base) for base in reversed(listOfExtensions[count]))
            addString += ("-------------------\n")
            addString += ("** PAM DESTROYED **\n")
            addString += ("PAM " + str(count + 1) + ": " + str(listByPos[count][1]) + "\n")
            addString += ("Position: " + str(listByPos[count][0] + 1) + "\n")
            addString += ("Spacer sequence Top: " + "cacc" + listOfSpacers[count] + "gtttt" + "\n")
            addString += ("Spacer sequence Bottom: " + "ctctaaaac" + reverse_complement_spacer + "\n")
            addString += ("Extension sequence Top: " + "gtgc" + listOfExtensions[count] + "\n")
            addString += ("Extension sequence Bottom: " + "aaaa" + reverse_complement_extension + "\n")
            break
        if (listByPos[count][2] == 2):
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
            reverse_complement_spacer = ''.join(complement.get(base, base) for base in reversed(listOfSpacers[count]))
            reverse_complement_extension = ''.join(complement.get(base, base) for base in reversed(listOfExtensions[count]))
            addString += ("-------------------\n")
            addString += ("** PAM CREATED **\n")
            addString += ("PAM " + str(count + 1) + ": " + str(listByPos[count][1]) + "\n")
            addString += ("Position: " + str(listByPos[count][0] + 1) + "\n")
            addString += ("Spacer sequence Top: " + "cacc" + listOfSpacers[count] + "gtttt" + "\n")
            addString += ("Spacer sequence Bottom: " + "ctctaaaac" + reverse_complement_spacer + "\n")
            addString += ("Extension sequence Top: " + "gtgc" + listOfExtensions[count] + "\n")
            addString += ("Extension sequence Bottom: " + "aaaa" + reverse_complement_extension + "\n")
            break
    for count in range (0, len(listByPos)):
        if (listByPos[count][2] == 1) or (listByPos[count][2] == 2):
            continue
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
        reverse_complement_spacer = ''.join(complement.get(base, base) for base in reversed(listOfSpacers[count]))
        reverse_complement_extension = ''.join(complement.get(base, base) for base in reversed(listOfExtensions[count]))
        addString += ("-------------------\n")
        addString += ("PAM " + str(count + 1) + ": " + str(listByPos[count][1]) + "\n")
        addString += ("Position: " + str(listByPos[count][0] + 1) + "\n")
        addString += ("Spacer sequence Top: " + "cacc" + listOfSpacers[count] + "gtttt" + "\n")
        addString += ("Spacer sequence Bottom: " + "ctctaaaac" + reverse_complement_spacer + "\n")
        addString += ("Extension sequence Top: " + "gtgc" + listOfExtensions[count] + "\n")
        addString += ("Extension sequence Bottom: " + "aaaa" + reverse_complement_extension + "\n")
    addString += ("----------------------------------")
    file1.write(addString)

def main():
    FASTA = FASTAEntry.get()
    mutation = mutationEntry.get()
    filename = filenameEntry.get()
    inputPAM = pamEntry.get()
    global PBSlength
    PBSlength = int(PBSlengthEntry.get())
    #For testing code:
    completename = os.path.join(os.path.dirname("Prime_Edit.py"), (filename + ".txt"))
    #For testing executable:
    # completename = os.path.join(os.path.dirname(sys.executable), (filename + ".txt"))
    file1 = open(completename, "w")
    parsedFASTA(FASTA)
    sequenceFinder(newString, position, inputPAM)
    pamDestroyed(inputPAM, mutation, listByPos)
    spacer(newString, listByPos)
    extension(newString, listByPos, mutation)
    ngRNA(newString, mutation)
    print ("Successfully found spacer and extension sequences for all PAMs.")
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

# FASTA = "ACCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTGTATCATCATCTTAGCGACGT(G)TAGCATGCTCTATCATCATCTCATGCTCTATCATCATCTGCATACGCATGCTCTATCATCATCTGTTAAATATAT"

# mutation = "T"
# main(FASTA, mutation)

canvas = Canvas(window, height = 200, width = 600)
canvas.pack()

frame = Frame(window,relief = 'groove')
frame.place(relx = 0.1, rely = 0.1, relwidth = 0.8, relheight = 0.8)

welcome = Label(frame, text = "Welcome to the spG Peg Design Program", fg = "Black")
welcome.pack(side = "top")


FASTAEntry = Entry(frame, width = 50)
FASTAEntry.pack(side = "top")
FASTAEntry.insert(0, "Please enter the DNA sequence")

mutationEntry = Entry(frame, width = 50)
mutationEntry.pack(side = "top")
mutationEntry.insert(0, "Please enter the desired mutation")

pamEntry = Entry(frame, width = 50)
pamEntry.pack(side = "top")
pamEntry.insert(0, "Please enter the desired PAM sequence")

filenameEntry = Entry(frame, width = 50)
filenameEntry.pack(side = "top")
filenameEntry.insert(0, "Please enter the desired .txt filename")

PBSlengthEntry = Entry(frame, width = 50)
PBSlengthEntry.pack(side = "top")
PBSlengthEntry.insert(0, "Please enter PBS Length (from 7-17)")


enterButton = Button(window, text = "Start", padx = 10, pady = 5, fg = "Black", bg = "gray", command = main)
enterButton.pack()

window.mainloop()
