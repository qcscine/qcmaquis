#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#/**
# * @file
# * @copyright This code is licensed under the 3-clause BSD license.
# *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
# *            See LICENSE.txt for details.
# */

import argparse
from VibrationalResultFile import ResultFileVibrationaleMeasurement
from typing import List

def pruneFcidumpFile(listOfListModals: List[ List[int] ], fileName: str, outputFileName: str):
  with open(fileName, mode="r") as fileObject:
    with open(outputFileName, mode="w") as outputFile:
      for iLine in fileObject:
        splittedLine = iLine.split()
        splittedLine.pop()
        accepted = True
        for iElement in splittedLine:
          iSplittedElement = [int(j) for j in iElement.split("-")]
          assert(len(iSplittedElement) == 2)
          if (not iSplittedElement[1] in listOfListModals[iSplittedElement[0]-1]):
            accepted = False
            break
        if accepted:
          outputFile.write(iLine)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--resultfile", type=str, help="""
                                                           Input result files.
                                                           If more than a single file is provided, the list of pruned modes is the
                                                           union of the pruned set obtained from each file.
                                                           """, required=True, nargs="+")
  parser.add_argument("-f", "--fcidump", type=str, action="store", help="If present, prunes directly the FCIDUMP file")
  parser.add_argument("-o", "--output", type=str, action="store", help="""
                                                                       Name of the output file.
                                                                       Required if the [fcidump] option is provided.
                                                                       """)
  parser.add_argument("-t", "--threshold", type=float, action="store", help="Pruning threshold for the one-modal entropy", required=True)
  args = parser.parse_args()
  fcidumpName = args.fcidump
  outputFile = args.output
  threshold = args.threshold
  # Creates the result file object and does the pruning
  resultFileList = [ResultFileVibrationaleMeasurement(iFile) for iFile in args.resultfile]
  assert(all([iRes.getLatticeSize() == resultFileList[0].getLatticeSize() for iRes in resultFileList]))
  lstOfAcceptedModes = [[] for i in range(resultFileList[0].getNumberOfModes())]
  for iRes in resultFileList:
    oneModalEntropies = iRes.getSeparatedOneModalEntropy()
    # Loop over the modals of each mode and, if the entropy is larger than the threshold
    # at least for one of the result files, include the mode.
    for iMode, entropyList in enumerate(oneModalEntropies):
      for idx, iModalEntropy in enumerate(entropyList):
        if abs(iModalEntropy) > threshold and idx not in lstOfAcceptedModes[iMode]:
          lstOfAcceptedModes[iMode].append(idx)
  #
  if (fcidumpName and outputFile):
    pruneFcidumpFile(lstOfAcceptedModes, fcidumpName, outputFile)
  else:
    print(" == LIST OF ACCEPTED MODES == ")
    for iMode, iList in enumerate(lstOfAcceptedModes):
      print(" Mode "+str(iMode))
      listOfModes = " ".join([str(i) for i in sorted(iList)])
      print(listOfModes)
