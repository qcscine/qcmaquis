#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2019- Laboratory for Physical Chemistry, ETH Zurich
#*               2019-2020 by Annina Lieberherr
#*               2022- by Alberto Baiardi <abaiardi@ethz.ch>
#*
#* This software is part of the ALPS Applications, published under the ALPS
#* Application License; you can use, redistribute it and/or modify it under
#* the terms of the license, either version 1 or (at your option) any later
#* version.
#*
#* You should have received a copy of the ALPS Application License along with
#* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
#* available from http://alps.comp-phys.org/.
#*
#* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
#* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
#* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
#* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#* DEALINGS IN THE SOFTWARE.
#*
#*****************************************************************************

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
        if len(splittedLine) > 2:
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
  parser.add_argument("-r", "--resultfile", type=str, help="Input result file", required=True, nargs="+")
  parser.add_argument("-f", "--fcidump", type=str, action="store", help="If present, prunes directly the FCIDUMP file")
  parser.add_argument("-o", "--output", type=str, action="store", help="""
                                                                       Name of the output file.
                                                                       Required if the [fcidump] option is provided.
                                                                       """)
  parser.add_argument("-t", "--threshold", type=float, action="store", help="Pruning threshold for the one-particle modals", required=True)
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
      listOfModes = " ".join([str(i) for i in iList])
      print(listOfModes)

  
