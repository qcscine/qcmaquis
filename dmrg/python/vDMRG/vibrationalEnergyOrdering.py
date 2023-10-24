#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#/**
# * @file
# * @copyright This code is licensed under the 3-clause BSD license.
# *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
# *            See LICENSE.txt for details.
# */

import argparse
import matplotlib.pyplot as plt
import numpy as np
from vibrationalEntanglementDiagram import plotMutualInformation
from VibrationalResultFile import ResultFileVibrationaleMeasurement

def pruneIntegralFile(listOfThresholds, integralFile, outputFile: str):
  inFile = open(integralFile, mode="r")
  outFile = open(outputFile, mode="w")
  localText = ""
  for i in inFile:
    if all([j[1] <= listOfThresholds[j[0]-1] for j in [[int(k) for k in j.split("-")] for j in i.split()[:-1]]]):
      localText += i
  outFile.write(localText)
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
  parser.add_argument("-f", "--fcidump", type=str, action="store", help="Name of the FCIDUMP file")
  parser.add_argument("-t", "--threshold", type=int, action="store", required=False, help="Threshold on the number of modals")
  parser.add_argument("-p", "--prunefcidump", action="store_true", help="If present, prunes the FCIDUMP file and generates a new one") 
  args = parser.parse_args()
  fcidumpFile = open(args.fcidump, mode="r")
  threshold = args.threshold
  writeNewFcidump = args.prunefcidump
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
    if writeNewFcidump:
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
