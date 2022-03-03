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
import math
import re
from sys import dont_write_bytecode
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import lines
from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
from VibrationalResultFile import ResultFileVibrationaleMeasurement

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--resultfile", type=str, help="Input result file", required=True)
  parser.add_argument("-f", "--fcidump", type=str, action="store", help="If present, prunes directly the FCIDUMP file")
  parser.add_argument("-t", "--threshold", type=float, action="store", help="Pruning threshold for the one-particle modals", required=True)
  args = parser.parse_args()
  inputFileName = args.resultfile
  fcidumpName = args.fcidump
  threshold = args.threshold
  # Creates the result file object and does the pruning
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  oneModalEntropies = resultFile.getSeparatedOneModalEntropy()
  lstOfAcceptedModes = []
  for entropyList in oneModalEntropies:
    acceptedModes = []
    for idx, iModalEntropy in enumerate(entropyList):
      if (abs(iModalEntropy) > threshold):
        acceptedModes.append(idx)
    lstOfAcceptedModes.append(acceptedModes)
  if (fcidumpName):
    pass
  else:
    print(" == LIST OF ACCEPTED MODES == ")
    for iMode, iList in enumerate(lstOfAcceptedModes):
      print(" Mode "+str(iMode))
      listOfModes = " ".join([str(i) for i in iList])
      print(listOfModes)

  
