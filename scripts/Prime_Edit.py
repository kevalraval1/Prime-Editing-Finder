import sys, os, regex as re
import Parsing_Helper as parse

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

def pamDestroyed(position, newString, inputPAM, mutation, listByPos, pattern):
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

def sequenceFinder (pattern, newString, position, inputPAM):
    listByPos = []
    optimal_string = ""
    for x in range (position - (4 + len(inputPAM)), position + (3 + len(inputPAM))): # position = index of mutation, 4 bases from end of first possible PAM
        optimal_string += newString[x]                                               # 2 bases from the start of last possible PAM, + 1 for range function (ends 1 before)
    print ("Mutation at optimal position in string: " + optimal_string)
    for match in pattern.finditer(optimal_string, overlapped=True):
        tempTuple = (match.start() + (position - (4 + len(inputPAM))), match.group(), 0) # List by Pos = [(start index in newstring, PAM sequence, 0/1/2 where 0 = PAM untouched 1 = destroyed 2 = created)]
        print(match)
        listByPos.append(tempTuple)
    if len(listByPos) == 0:
        return print("No available PAM sites for given mutation.")
    else:
        print (listByPos)
        return listByPos

def spacer(newString, listByPos):
    listOfSpacers = []
    for tup in listByPos:
        spacerSequence = ""
        if (newString[tup[0] - 20]) != "G":
            spacerSequence += "G"
        for bases in range ((tup[0] - 20), (tup[0])):
            spacerSequence += newString[bases]
        listOfSpacers.append(spacerSequence) # List of Spacers = [spacer 1 for PAM 1, spacer 2 for PAM 2]
    print(listOfSpacers)
    return listOfSpacers

def extension(position, newString, listByPos, mutation, PBSlength):
    listOfExtensions = []
    for tup in listByPos:
        extensionSequence = ""
        for bases in range ((tup[0] - (PBSlength + 3)), (tup[0] + 10)): # PBS Length always ends 3 BP before PAM, RTT always starts 3 BP before PAM and is 13 BP long
            if bases == position:
                extensionSequence += mutation.lower()
                continue
            extensionSequence += newString[bases]
        reverse_complement = ''.join(parse.complementDict.get(base, base) for base in reversed(extensionSequence))
        listOfExtensions.append(reverse_complement) # List of Extensions = [Extension 1 for PAM 1, Extension 2 for PAM 2]
    print(listOfExtensions)
    return listOfExtensions

def ngRNA(position, newString, mutation, pattern):
    listOfngRNA = []
    tempList = list(newString)
    tempList[position] = mutation.lower()
    newString = "".join(tempList)
    reversedString = ''.join(parse.complementDict.get(base, base) for base in reversed(newString))
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
        ngRNA2 = ''.join(parse.complementDict.get(base, base) for base in reversed(ngRNA1))
        tuple = (ngRNA1, ngRNA2)
        listOfngRNA.append(tuple) # List of ngRNA = [(top ngRNA, bot ngRNA), (top ngRNA, bot ngRNA)] plus strand then minus strand
    else:
        tuple = ("N/A", "N/A")
        listOfngRNA.append(tuple)
    return listOfngRNA

def analysisPrinter(listByPos, listOfSpacers, listOfExtensions, listOfngRNA, file1):
    if (len(listByPos) == 0):
        file1.write("NO PAM SITES AVAILABLE IN GIVEN INPUT")
        return
    printList = [f"ngRNA Top: cacc{listOfngRNA[0][0]}", 
                f"ngRNA Bottom: aaac{listOfngRNA[0][1]}\n", 
                ]
    for count in range(0, len(listByPos)):
        reverse_complement_spacer = ''.join(parse.complementDict.get(base, base) for base in reversed(listOfSpacers[count]))
        reverse_complement_extension = ''.join(parse.complementDict.get(base, base) for base in reversed(listOfExtensions[count]))
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
    return

def main(FASTA, mutation, filename, inputPAM, PBS):
    PBSlength = int(PBS)
    newString = parse.FASTA(FASTA)[0]
    position = parse.FASTA(FASTA)[1]
    pattern = regexCompiler(inputPAM)
    listByPos = sequenceFinder(pattern, newString, position, inputPAM)
    pamDestroyed(position, newString, inputPAM, mutation, listByPos, pattern)
    listOfSpacers = spacer(newString, listByPos)
    listOfExtensions = extension(position, newString, listByPos, mutation, PBSlength)
    listOfngRNA = ngRNA(position, newString, mutation, pattern)
    print ("Successfully found spacer and extension sequences for all PAMs.")
    #For testing code:
    completename = os.path.join(os.path.dirname("Main_Interface.py"), (filename + ".txt"))
    #For testing executable:
    # completename = os.path.join(os.path.dirname(sys.executable), (filename + ".txt"))
    file1 = open(completename, "w")
    file1.write("\n*PLUS STRAND ANALYSIS*\n\n")
    analysisPrinter(listByPos, listOfSpacers, listOfExtensions, listOfngRNA, file1)
    FASTA = parse.reverser(FASTA, mutation)[0]
    mutation = parse.reverser(FASTA, mutation)[1]
    newString = parse.FASTA(FASTA)[0]
    position = parse.FASTA(FASTA)[1]
    pattern = regexCompiler(inputPAM)
    listByPos = sequenceFinder(pattern, newString, position, inputPAM)
    pamDestroyed(position, newString, inputPAM, mutation, listByPos, pattern)
    listOfSpacers = spacer(newString, listByPos)
    listOfExtensions = extension(position, newString, listByPos, mutation, PBSlength)
    listOfngRNA = ngRNA(position, newString, mutation, pattern)
    print ("Successfully found spacer and extension sequences for all PAMs.")
    file1.write("\n\n*MINUS STRAND ANALYSIS*\n\n")
    analysisPrinter(listByPos, listOfSpacers, listOfExtensions, listOfngRNA, file1)
    print ("Exiting...")
    file1.close()
    sys.exit()

# FASTA = "ACCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTGTATCATCATCTTAGCGACGGT(G)TAGCATGCTCTATCATCATCTCATGCTCTATCATCATCTGCATACGCATGCTCTATCATCATCTGTTAAATATAT"
