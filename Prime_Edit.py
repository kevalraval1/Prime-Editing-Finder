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

def sequenceFinder (newString, position, inputPAM):
    global listByPos
    listByPos = []
    optimal_string = ""
    for x in range (position - (4 + len(inputPAM)), position + (3 + len(inputPAM))): # position = index of mutation, 4 bases from end of first possible PAM
        optimal_string += newString[x]                                               # 2 bases from the start of last possible PAM, + 1 for range function (ends 1 before)
    print ("Mutation at optimal position in string: " + optimal_string)
    pattern = regexCompiler(inputPAM)
    for match in pattern.finditer(optimal_string, overlapped=True):
        tempTuple = (match.start() + (position - (4 + len(inputPAM))), match.group()) # List by Pos = [(index in newstring, PAM sequence)]
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
        for bases in range ((tup[0]-(PBSlength + 4)), (tup[0] + 9)):
            if bases == position:
                extensionSequence += mutation.lower()
                continue
            extensionSequence += newString[bases]
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
        reverse_complement = ''.join(complement.get(base, base) for base in reversed(extensionSequence))
        listOfExtensions.append(reverse_complement)
    return (print(listOfExtensions))

def ngRNA(newString, mutation):
    global listOfngRNA
    listOfngRNA = []
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '(': ')', ')': '('}
    tempList = list(newString)
    tempList[position-1] = mutation.lower()
    newString = "".join(tempList)
    reversedString = ''.join(complement.get(base, base) for base in reversed(newString))
    counter = -1
    newPosition = len(newString) - position + 1
    for base in range (newPosition, newPosition + 7):
        if reversedString[base] == "G":
            counter = base - 1
            break
    if (counter != -1):
        ngRNA1 = ""
        if reversedString[counter - 20] != "G":
            ngRNA1 = ngRNA1 + "G"
        for bases in range (counter - 20, counter):
            ngRNA1 = ngRNA1 + reversedString[bases]
        ngRNA2 = ''.join(complement.get(base, base) for base in reversed(ngRNA1))
        tuple = (ngRNA1, ngRNA2)
        listOfngRNA.append(tuple)
    if (counter == - 1):
        tuple = ("N/A", "N/A")
        listOfngRNA.append(tuple)

def analysisPrinter(listByPos, listOfSpacers, listOfExtensions, file1):
    addString = ""
    addString = addString + ("ngRNA Top: " + "cacc" + listOfngRNA[0][0] + "\n")
    addString = addString + ("ngRNA Bottom: " + "aaac" + listOfngRNA[0][1] + "\n\n")
    for count in range (0, len(listByPos)):
        if (listByPos[count][0] + 1 == position):
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
            reverse_complement_spacer = ''.join(complement.get(base, base) for base in reversed(listOfSpacers[count]))
            reverse_complement_extension = ''.join(complement.get(base, base) for base in reversed(listOfExtensions[count]))
            addString = addString + ("-------------------\n")
            addString = addString + ("** PAM DESTROYED **\n")
            addString = addString + ("PAM " + str(count + 1) + ": " + str(listByPos[count][1]) + "\n")
            addString = addString + ("Position: " + str(listByPos[count][0] + 1) + "\n")
            addString = addString + ("Spacer sequence Top: " + "cacc" + listOfSpacers[count] + "gtttt" + "\n")
            addString = addString + ("Spacer sequence Bottom: " + "ctctaaaac" + reverse_complement_spacer + "\n")
            addString = addString + ("Extension sequence Top: " + "gtgc" + listOfExtensions[count] + "\n")
            addString = addString + ("Extension sequence Bottom: " + "aaaa" + reverse_complement_extension + "\n")
            break
    for count in range (0, len(listByPos)):
        if (listByPos[count][0] + 1 == position):
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
