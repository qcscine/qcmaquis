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

def pruneIntegralFile(listOfThresholds, integralFile, outputFile: str):
  inFile = open(integralFile, mode="r")
  outFile = open(outputFile, mode="w")
  for i in inFile:
    tmpList = [[int(k) for k in j.split("-")] for j in i.split()[:-1]]
    if all([j[1] <= listOfThresholds[j[0]-1] for j in tmpList]):
      outFile.write(i)
  outFile.close()

def calculateMaximumModalIndex(inputList, threshold: int):
  numberOfModes = max([i[0] for i in inputList])
  listOfThresholds = []
  for i in range(numberOfModes):
    listOfThresholds.append(max([idx[1] for idx in inputList[:threshold] if idx[0] == i+1]))
  return listOfThresholds

def extractOneBodyEnergyFromFCIDUMP(resultFile):
  listOfSingleParticleEnergies = []
  for i in resultFile:
    lineAsList = i.split()
    if len(lineAsList) == 3:
      firstMode = int(lineAsList[0].split("-")[0])
      secondMode = int(lineAsList[1].split("-")[0])
      firstModal = int(lineAsList[0].split("-")[1])
      secondModal = int(lineAsList[1].split("-")[1])
      if firstMode == secondMode and firstModal == secondModal:
        listOfSingleParticleEnergies.append([firstMode, firstModal, float(lineAsList[2])])
  return listOfSingleParticleEnergies

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-f", "--fcidump", type=str, action="store", help="If present, prunes directly the FCIDUMP file")
  parser.add_argument("-t", "--threshold", type=int, action="store", required=False, help="Threshold on the number of modals")
  args = parser.parse_args()
  fcidumpFile = open(args.fcidump, mode="r")
  threshold = args.threshold
  # Prints header
  print("")
  print(" =========================================== ")
  print("  SORTING MODALS IN INCREASING ENERGY ORDER  ")
  print(" =========================================== \n")
  # Extracts relevant data
  listOfSingleParticleEnergies = extractOneBodyEnergyFromFCIDUMP(fcidumpFile)
  listOfSingleParticleEnergiesSorted = listOfSingleParticleEnergies[:]
  listOfSingleParticleEnergiesSorted.sort(key=lambda x: x[2])
  sorting = ""
  for i in listOfSingleParticleEnergiesSorted:
    sorting += str(listOfSingleParticleEnergies.index(i))
    sorting += ","
  print(" Output sorting for the full modal set: ")
  print(" "+sorting[:-1])
  # If a threshold is provided, cleans the input file.
  if threshold:
    print("")
    print(" +-- Pruning the basis such that "+str(threshold)+" modal basis functions are retained --+\n")
    listOfThresholds = calculateMaximumModalIndex(listOfSingleParticleEnergiesSorted, threshold)
    for idx, i in enumerate(listOfThresholds):
      print(" "+str(i+1)+" modal basis functions retained for mode "+str(idx))
    newIntegralFileName = args.fcidump+"_Pruned"
    pruneIntegralFile(listOfThresholds, args.fcidump, newIntegralFileName)
    print("")
    print(" New integral file generated with name "+newIntegralFileName)
    print("")
    # Calculated pruned sorting
    listOfSingleParticleEnergiesPruned = [i for i in listOfSingleParticleEnergies if i[1] <= listOfThresholds[i[0]-1]]
    listOfSingleParticleEnergiesPrunedSorted = listOfSingleParticleEnergiesPruned[:]
    listOfSingleParticleEnergiesPrunedSorted.sort(key=lambda x: x[2])
    sortingPruned = ""
    for i in listOfSingleParticleEnergiesPrunedSorted:
      sortingPruned += str(listOfSingleParticleEnergiesPruned.index(i))
      sortingPruned += ","
    print(" Output sorting for the pruned modal set: ")
    print(" "+sortingPruned[:-1])

