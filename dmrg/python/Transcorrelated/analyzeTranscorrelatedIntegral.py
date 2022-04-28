import argparse

# Parse input parameters
parser = argparse.ArgumentParser(description='Analyze the norm of a given transcorrelated FCIDUMP file.')
parser.add_argument('-f', '--fcidmp', type=str, action='store', required=True, help='Input FCIDMP file')
parser.add_argument('-t', '--threshold', type=float, required=False, default=0., help='Threshold for the Hamiltonian terms')

args = parser.parse_args()
inputFileName = args.fcidmp
threshold = args.threshold
inFile = open(inputFileName, mode='r')

# Analyziz part
oneBodyNorm = 0.
twoBodyNorm = 0.
threeBodyNorm = 0.
oneBodyTerms = 0
twoBodyTerms = 0
threeBodyTerms = 0

for i in inFile:
    indices = [int(j) for j in i.split()[1:]]
    factor = float(i.split()[0])
    if abs(factor) > threshold:
        if indices[2] == 0:
            oneBodyTerms += 1
            oneBodyNorm += abs(factor)
        elif indices[4] == 0:
            twoBodyTerms += 1
            twoBodyNorm += abs(factor)
        else:
            threeBodyTerms += 1
            threeBodyNorm += abs(factor)

print("Number of one-body terms :"+str(oneBodyTerms)+" with L1 norm "+str(oneBodyNorm))
print("Number of two-body terms :"+str(twoBodyTerms)+" with L1 norm "+str(twoBodyNorm))
print("Number of three-body terms :"+str(threeBodyTerms)+" with L1 norm "+str(threeBodyNorm))

inFile.close()
