import argparse
import math
import sys
from copy import deepcopy

# Parse input parameters
parser = argparse.ArgumentParser(description='Analyze the norm of a given transcorrelated FCIDUMP file.')
parser.add_argument('-f', '--fcidmp', type=str, action='store', required=True, help='Input FCIDUMP file')
parser.add_argument('-o', '--output', type=str, action='store', required=True, help='Output FCIDUMP file')
parser.add_argument('-t', '--threshold', type=float, required=False, help='L1 truncation threshold for the Hamiltonian')
parser.add_argument('-l', '--l2threshold', type=float, required=False, help='L2 truncation threshold for the Hamiltonian')
parser.add_argument('--threebody_pruning', action='store_true', help='If present, returns also the FCIDUMP where only the three-body part is pruned')

args = parser.parse_args()
inputFileName = args.fcidmp
outputFileName = args.output
l1threshold = args.threshold
l2threshold = args.l2threshold
generate3BIntegrals = args.threebody_pruning

if l1threshold and l2threshold:
    sys.exit("Please just set only one between the L1 and L2 thresholds")

inFile = open(inputFileName, mode='r')

oneBodyLNorm = 0.
twoBodyLNorm = 0.
threeBodyLNorm = 0.

oneBodyTerms = 0
twoBodyTerms = 0
threeBodyTerms = 0

listOfTerms = []
listOfOneBodyTerms = []
listOfTwoBodyTerms = []
listOfThreeBodyTerms = []

# == Parsing of the integral file ==

for i in inFile:
    indices = [int(j) for j in i.split()[1:]]
    factor = float(i.split()[0])
    # One-body 
    if indices[2] == 0:
        oneBodyTerms += 1
        oneBodyLNorm += factor*factor
        if l1threshold:
            listOfTerms.append([indices, abs(factor), factor])
            listOfOneBodyTerms.append([indices, abs(factor), factor])
        if l2threshold:
            listOfTerms.append([indices, factor*factor, factor])
            listOfOneBodyTerms.append([indices, factor*factor, factor])
    # Two-body
    elif indices[4] == 0:
        twoBodyTerms += 1
        twoBodyLNorm += factor*factor
        if l1threshold:
            listOfTerms.append([indices, abs(factor), factor])
            listOfTwoBodyTerms.append([indices, abs(factor), factor])
        if l2threshold:
            listOfTerms.append([indices, factor*factor, factor])
            listOfTwoBodyTerms.append([indices, factor*factor, factor])
    else:
        threeBodyTerms += 1
        threeBodyLNorm += factor*factor
        if l1threshold:
            listOfTerms.append([indices, abs(factor), factor])
            listOfThreeBodyTerms.append([indices, abs(factor), factor])
        if l2threshold:
            listOfTerms.append([indices, factor*factor, factor])
            listOfThreeBodyTerms.append([indices, factor*factor, factor])

if l1threshold:
    print("L1 norm of the one-body terms :"+str(oneBodyLNorm))
    print("L1 norm of the two-body terms :"+str(twoBodyLNorm))
    print("L1 norm of the three-body terms :"+str(threeBodyLNorm))

if l2threshold:
    print("L2 norm of the one-body terms :"+str(math.sqrt(oneBodyLNorm)))
    print("L2 norm of the two-body terms :"+str(math.sqrt(twoBodyLNorm)))
    print("L2 norm of the three-body terms :"+str(math.sqrt(threeBodyLNorm)))

listOfPrunedTerms = deepcopy(listOfTerms)
listOfPrunedOneBodyTerms = deepcopy(listOfOneBodyTerms)
listOfPrunedTwoBodyTerms = deepcopy(listOfTwoBodyTerms)
listOfPrunedThreeBodyTerms = deepcopy(listOfThreeBodyTerms)

if l1threshold or l2threshold:
    listOfPrunedTerms.sort(key=lambda x: x[1])
    listOfPrunedOneBodyTerms.sort(key=lambda x: x[1])
    listOfPrunedTwoBodyTerms.sort(key=lambda x: x[1])
    listOfPrunedThreeBodyTerms.sort(key=lambda x: x[1])
    sumOfTerms = 0.
    sumOfOneBodyTerms = 0.
    sumOfTwoBodyTerms = 0.
    sumOfThreeBodyTerms = 0.
    # == L1 truncation ==
    if l1threshold:
        while sumOfOneBodyTerms < l1Threshold and len(listOfPrunedOneBodyTerms) != 0:
            sumOfOneBodyTerms += listOfPrunedOneBodyTerms.pop(0)[1]
        print("Number of one-body terms after L1 truncation :"+str(len(listOfPrunedOneBodyTerms)))
        while sumOfTwoBodyTerms < l1Threshold and len(listOfPrunedTwoBodyTerms) != 0:
            sumOfTwoBodyTerms += listOfPrunedTwoBodyTerms.pop(0)[1]
        print("Number of two-body terms after L1 truncation :"+str(len(listOfPrunedTwoBodyTerms)))
        while sumOfThreeBodyTerms < l1Threshold and len(listOfPrunedThreeBodyTerms) != 0:
            sumOfThreeBodyTerms += listOfPrunedThreeBodyTerms.pop(0)[1]
        print("Number of three-body terms after L1 truncation :"+str(len(listOfPrunedThreeBodyTerms)))
        while sumOfTerms < l1Threshold and len(listOfPrunedTerms) != 0:
            sumOfTerms += listOfPrunedTerms.pop(0)[1]
        print("Number of coupling terms after L1 truncation :"+str(len(listOfPrunedTerms)))
    # == L2 truncation ==
    if l2threshold:
        while math.sqrt(sumOfOneBodyTerms) < l2threshold and len(listOfPrunedOneBodyTerms) != 0:
            sumOfOneBodyTerms += listOfPrunedOneBodyTerms.pop(0)[1]
        print("Number of one-body terms after L2 truncation :"+str(len(listOfPrunedOneBodyTerms)))
        while math.sqrt(sumOfTwoBodyTerms) < l2threshold and len(listOfPrunedTwoBodyTerms) != 0:
            sumOfTwoBodyTerms += listOfPrunedTwoBodyTerms.pop(0)[1]
        print("Number of two-body terms after L2 truncation :"+str(len(listOfPrunedTwoBodyTerms)))
        while math.sqrt(sumOfThreeBodyTerms) < l2threshold and len(listOfPrunedThreeBodyTerms) != 0:
            sumOfThreeBodyTerms += listOfPrunedThreeBodyTerms.pop(0)[1]
        print("Number of three-body terms after L2 truncation :"+str(len(listOfPrunedThreeBodyTerms)))
        while math.sqrt(sumOfTerms) < l2threshold and len(listOfPrunedTerms) != 0:
            sumOfTerms += listOfPrunedTerms.pop(0)[1]
        print("Number of coupling terms after L2 truncation :"+str(len(listOfPrunedTerms)))
    # == In all cases, the FCIDUMP must be written with what's left ==
    # Note that the DMRG calculation is done by truncating the *whole* Hamiltonian (so, one-, two-, and three-body parts)
    if outputFileName:
        outputFile = open(outputFileName, mode="w")
        templateString = "  {value:19.12E}   {index[0]:2d}   {index[1]:2d}   {index[2]:2d}   {index[3]:2d}   {index[4]:2d}   {index[5]:2d}\n"
        for i in listOfPrunedTerms:
            outputFile.write(templateString.format(value=i[2], index=i[0]))
        outputFile.close()
    # If specified, generates also the FCIDUMP with only the three-body term pruned
    if generate3BIntegrals:
        outputFile = open(outputFileName+"_ThreeBodyPruning", mode="w")
        templateString = "  {value:19.12E}   {index[0]:2d}   {index[1]:2d}   {index[2]:2d}   {index[3]:2d}   {index[4]:2d}   {index[5]:2d}\n"
        for i in listOfOneBodyTerms:
            outputFile.write(templateString.format(value=i[2], index=i[0]))
        for i in listOfTwoBodyTerms:
            outputFile.write(templateString.format(value=i[2], index=i[0]))
        for i in listOfPrunedThreeBodyTerms:
            outputFile.write(templateString.format(value=i[2], index=i[0]))
        outputFile.close()

inFile.close()

