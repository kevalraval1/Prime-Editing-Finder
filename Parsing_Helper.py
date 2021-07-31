complementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '(': ')', ')': '('}

def reverser(FASTA, mutation):
    newFASTA = ''.join(complementDict.get(base, base) for base in reversed(FASTA))
    newMutation = ''.join(complementDict.get(base, base) for base in reversed(mutation))
    return newFASTA, newMutation

def FASTA (FASTA):
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
    print("Successfully parsed FASTA")
    return (newString, position)
