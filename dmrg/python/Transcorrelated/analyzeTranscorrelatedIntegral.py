import argparse
import math

# Parse input parameters
parser = argparse.ArgumentParser(description='Analyze the norm of a given transcorrelated FCIDUMP file.')
parser.add_argument('-f', '--fcidmp', type=str, action='store', required=True, help='Input FCIDMP file')
parser.add_argument('-t', '--threshold', type=float, required=False, default=0., help='Threshold for the Hamiltonian terms')
parser.add_argument('-l', '--l2threshold', type=float, required=False, help='L2 truncation threshold for the Hamiltonian')

args = parser.parse_args()
inputFileName = args.fcidmp
threshold = args.threshold
l2Threshold = args.l2threshold

inFile = open(inputFileName, mode='r')

# Analysis part
oneBodyNorm = 0.
twoBodyNorm = 0.
threeBodyNorm = 0.

oneBodyLSquaredNorm = 0.
twoBodyLSquaredNorm = 0.
threeBodyLSquaredNorm = 0.

oneBodyTerms = 0
twoBodyTerms = 0
threeBodyTerms = 0

listOfTerms = []
listOfTwoBodyTerms = []
listOfThreeBodyTerms = []

for i in inFile:
    indices = [int(j) for j in i.split()[1:]]
    factor = float(i.split()[0])
    # Conventional thresholding
    if abs(factor) > threshold:
        if indices[2] == 0:
            oneBodyTerms += 1
            oneBodyNorm += abs(factor)
            oneBodyLSquaredNorm += factor*factor
        elif indices[4] == 0:
            twoBodyTerms += 1
            twoBodyNorm += abs(factor)
            twoBodyLSquaredNorm += factor*factor
            if l2Threshold:
                listOfTerms.append([indices, factor*factor])
                listOfTwoBodyTerms.append([indices, factor*factor])
        else:
            threeBodyTerms += 1
            threeBodyNorm += abs(factor)
            threeBodyLSquaredNorm += factor*factor
            if l2Threshold:
                listOfTerms.append([indices, factor*factor])
                listOfThreeBodyTerms.append([indices, factor*factor])

print("Number of one-body terms :"+str(oneBodyTerms)+" with L1 norm "+str(oneBodyNorm))
print("Number of two-body terms :"+str(twoBodyTerms)+" with L1 norm "+str(twoBodyNorm))
print("Number of three-body terms :"+str(threeBodyTerms)+" with L1 norm "+str(threeBodyNorm))

print("L2 norm of the one-body terms :"+str(math.sqrt(oneBodyLSquaredNorm)))
print("L2 norm of the two-body terms :"+str(math.sqrt(twoBodyLSquaredNorm)))
print("L2 norm of the three-body terms :"+str(math.sqrt(threeBodyLSquaredNorm)))

if l2Threshold:
    listOfTerms.sort(key=lambda x: x[1])
    listOfTwoBodyTerms.sort(key=lambda x: x[1])
    listOfThreeBodyTerms.sort(key=lambda x: x[1])
    sumOfTerms = 0.
    sumOfTwoBodyTerms = 0.
    sumOfThreeBodyTerms = 0.
    while math.sqrt(sumOfTwoBodyTerms) < l2Threshold and len(listOfTwoBodyTerms) != 0:
        sumOfTwoBodyTerms += listOfTwoBodyTerms.pop(0)[1]
    print("Number of two-body terms after L2 truncation :"+str(len(listOfTwoBodyTerms)))
    while math.sqrt(sumOfThreeBodyTerms) < l2Threshold and len(listOfThreeBodyTerms) != 0:
        sumOfThreeBodyTerms += listOfThreeBodyTerms.pop(0)[1]
    print("Number of three-body terms after L2 truncation :"+str(len(listOfThreeBodyTerms)))
    while math.sqrt(sumOfTerms) < l2Threshold and len(listOfTerms) != 0:
        sumOfTerms += listOfTerms.pop(0)[1]
    print("Number of coupling terms after L2 truncation :"+str(len(listOfTerms)))

inFile.close()
