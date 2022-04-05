#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2022- Laboratory for Physical Chemistry, ETH Zurich
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
import matplotlib.pyplot as plt
import numpy as np
from vibrationalEntanglementDiagram import plotMutualInformation
from VibrationalResultFile import ResultFileVibrationaleMeasurement

def extractOneBodyEnergyFromFCIDUMP(resultFile):
  listOfSingleParticleEnergies = []
  for i in fcidumpFile:
    lineAsList = i.split()
    if len(lineAsList) == 3:
      firstMode = lineAsList[0].split("-")[0]
      secondMode = lineAsList[1].split("-")[0]
      firstModal = lineAsList[0].split("-")[1]
      secondModal = lineAsList[1].split("-")[1]
      if firstMode == secondMode and firstModal == secondModal:
        listOfSingleParticleEnergies.append([firstMode, firstModal, float(lineAsList[2])])
  return listOfSingleParticleEnergies

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-f", "--fcidump", type=str, action="store", help="If present, prunes directly the FCIDUMP file")
  args = parser.parse_args()
  fcidumpFile = open(args.fcidump, mode="r")
  listOfSingleParticleEnergies = extractOneBodyEnergyFromFCIDUMP(fcidumpFile)
  listOfSingleParticleEnergiesSorted = listOfSingleParticleEnergies[:]
  listOfSingleParticleEnergiesSorted.sort(key=lambda x: x[2])
  sorting = ""
  for i in listOfSingleParticleEnergiesSorted:
    sorting += str(listOfSingleParticleEnergies.index(i))
    sorting += ","
  print("Sorting of the modals in increasing energy order")
  print(sorting[:-1])
