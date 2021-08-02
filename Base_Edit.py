import Parsing_Helper as parse
import os, sys

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

def findgRNA(workingStrand, position):
    listOfgRNA = []
    for times in range (5):
        gRNA = ""
        PAM = ""
        for bases in range ((position - 3 - times), (position + 17 - times)):
            if (bases == (position - 3 - times) and workingStrand[bases] != "G"):
                gRNA += "G"
            elif (bases == position):
                gRNA += workingStrand[bases].lower()
                continue
            gRNA += workingStrand[bases]
        for bases in range (position + 17 - times, position + 20 - times):
            PAM += workingStrand[bases]
        tempTuple = (gRNA, PAM)
        listOfgRNA.append(tempTuple)
    return listOfgRNA

def bystanderMutations(listOfgRNA, enzyme):
    for index, tempTuple in enumerate(listOfgRNA):
        count = 0
        for bases in range (3, 8):
            if enzyme == "CBE":
                if (tempTuple[0][bases]) == "C":
                    count += 1
            elif enzyme == "ABE":
                if (tempTuple[0][bases]) == "A":
                    count += 1
        tempList = list(tempTuple)
        tempList.append(count)
        newTuple = tuple(tempList)
        listOfgRNA[index] = newTuple
    return listOfgRNA

def analysisPrinter(listOfgRNA, enzyme, strand, file1):
    printList = [f"Enzyme to use: {enzyme}",
                f"Working Strand: {strand}"]
    for count, tuple in enumerate(listOfgRNA):
        tempList = ["------------------------------",
                    f"guideRNA Top {count + 1}: cacc{tuple[0]}",
                    f"guideRNA Bottom: {parse.reverser(tuple[0])[0]}",
                    f"PAM: {tuple[1]}",
                    f"Number of Bystander mutations: {tuple[2]}"]
        tempString = '\n'.join(tempList)
        printList.append(tempString)
    printList.append("------------------------------")
    addString = '\n'.join(printList)
    file1.write(addString)
    return

def main(FASTA, mutation, fileName):
    tempTuple = validSubstitution(FASTA, mutation)
    newString = tempTuple[0]
    position = tempTuple[1]
    newMutation = tempTuple[2]
    enzyme = tempTuple[3]
    strand = tempTuple[4]
    listOfgRNA = findgRNA(newString, position)
    listOfgRNA = bystanderMutations(listOfgRNA, enzyme)
    #For testing code:
    completename = os.path.join(os.path.dirname("Main_Interface.py"), (fileName + ".txt"))
    #For testing executable:
    # completename = os.path.join(os.path.dirname(sys.executable), (fileName + ".txt"))
    file1 = open(completename, "w")
    analysisPrinter(listOfgRNA, enzyme, strand, file1)
    file1.close()
    sys.exit()



# FASTA = "ACCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTCTATCATCATCTCATGCTGTATCATCATCTTAGCGACGGT(G)TAGCATGCTCTATCATCATCTCATGCTCTATCATCATCTGCATACGCATGCTCTATCATCATCTGTTAAATATAT"
# mutation = "A"
# print (validSubstitution(FASTA, mutation))