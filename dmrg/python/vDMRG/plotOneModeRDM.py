#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#/**
# * @file
# * @copyright This code is licensed under the 3-clause BSD license.
# *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
# *            See LICENSE.txt for details.
# */

import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import re
import sys
from typing import Dict, List, Set
from VibrationalResultFile import ResultFileVibrationaleMeasurement

def extractIntegrals(inputFileName: str, numModes: int, numBasis: List[int]):
  """
  Extracts the integrals from a file and puts them in a proper
  data structure.
  """
  integralFile = open(inputFileName, mode="r")
  # Data structure initialization
  dictOneBody = {}
  dictTwoBody = {}
  for i in range(numModes):
    dictOneBody[i] = np.zeros([numBasis[i], numBasis[i]])
    for j in range(i):
      dictTwoBody[(i, j)] = np.zeros([numBasis[i], numBasis[i], numBasis[j], numBasis[j]])
  # Data extraction
  for i in integralFile:
    splittedLine = i.split()
    if len(splittedLine) == 3:
      iMode  = int(splittedLine[0].split("-")[0])-1
      iModal = int(splittedLine[0].split("-")[1])
      jMode  = int(splittedLine[1].split("-")[0])-1
      jModal = int(splittedLine[1].split("-")[1])
      assert iMode == jMode
      dictOneBody[iMode][iModal, jModal] = float(splittedLine[2])
    elif len(splittedLine) == 5:
      iMode  = int(splittedLine[0].split("-")[0])-1
      iModal = int(splittedLine[0].split("-")[1])
      jMode  = int(splittedLine[1].split("-")[0])-1
      jModal = int(splittedLine[1].split("-")[1])
      kMode  = int(splittedLine[2].split("-")[0])-1
      kModal = int(splittedLine[2].split("-")[1])
      lMode  = int(splittedLine[3].split("-")[0])-1
      lModal = int(splittedLine[3].split("-")[1])
      assert iMode == jMode and kMode == lMode
      dictTwoBody[(max(iMode, kMode), min(iMode, kMode))][iModal, jModal, kModal, lModal] = float(splittedLine[4])
    else:
      print(splittedLine)
      sys.exit("Pruning for three-mode Hamiltonians NYI")
  integralFile.close()
  return dictOneBody, dictTwoBody


def plotOneModeRDM(inputFileName: str):
  """
  Plots the one-mode entropy as a 2D plot.
  """
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  numberOfModals = resultFile.getNumberOfModals()
  oneModeRDM = resultFile.getOneModeRDM()
  print("Diagonal elements of the one-modal RDM")
  for i in range(resultFile.getLatticeSize()):
    print("Modal "+str(i)+" RDM = "+str(oneModeRDM[i, i]))
  fix, ax = plt.subplots()
  ax.matshow(oneModeRDM, cmap=plt.cm.Blues)
  offset = 0
  for i in numberOfModals:
    rect = plt.Rectangle((offset-0.5, offset-0.5), i, i, fill=False)
    ax.add_patch(rect)
    offset += i
  plt.show()


def getNaturalOrbitals(resultFileName: str, threshold: float):
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  oneModeRDM = resultFile.getOneModeRDM()
  numberOfModals = resultFile.getNumberOfModals()
  eigenValuesDict = {}
  eigenVectorsDict = {}
  # Main loop
  overallCont = 0
  for iCount, i in enumerate(numberOfModals):
    modeSpecificRDM = np.copy(oneModeRDM[overallCont: overallCont+i, overallCont: overallCont+i])
    eigenValues, eigenVectors = np.linalg.eigh(modeSpecificRDM)
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    # Pruning
    found = False
    i = 0
    while not found: 
      if abs(eigenValues[i]) < threshold or i == eigenValues.shape[0]-1:
        found = True
        i += 1
      else:
        i += 1
    eigenValuesDict[iCount] = eigenValues[:i]
    eigenVectorsDict[iCount] = eigenVectors[:i,:]
    overallCont += i
  return eigenValuesDict, eigenVectorsDict


def rotateOneBody(resultFileName: str, oneModeDict: Dict[int, np.array], threshold: float):
  """
  Transform the integrals (taken as a dict) based on the MO --> NO transformation
  calculated form a result file of a DMRG calculation.
  """
  ret = {}
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  numberOfModes = resultFile.getNumberOfModes()
  eigenValuesDict, eigenVectorsDict = getNaturalOrbitals(resultFileName, threshold)
  for i in range(numberOfModes):
    ret[i] = np.einsum("ij,jk,lk->il", eigenVectorsDict[i], oneModeDict[i], eigenVectorsDict[i])
  return ret


def rotateTwoBody(resultFileName: str, twoModeDict: Dict[Set, np.array], threshold: float):
  """
  Transform the integrals (taken as a dict) based on the MO --> NO transformation
  calculated form a result file of a DMRG calculation.
  """
  ret = {}
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  numberOfModes = resultFile.getNumberOfModes()
  eigenValuesDict, eigenVectorsDict = getNaturalOrbitals(resultFileName, threshold)
  for i in range(numberOfModes):
    for j in range(i):
      intermediate1 = np.einsum("jkbc,ab->jkac", twoModeDict[(i, j)], eigenVectorsDict[j])
      intermediate2 = np.einsum("jkac,dc->jkad", intermediate1, eigenVectorsDict[j])
      intermediate3 = np.einsum("ij,jkad->ikad", eigenVectorsDict[i], intermediate2)
      ret[(i, j)] = np.einsum("lk,ikad->ilad", eigenVectorsDict[i], intermediate3)
      #rotatedOneModeRDM = np.einsum("ji,kl,jkbc,ba,cd->iald", eigenVectorsDict[i], eigenVectorsDict[i], twoModeDict[(i, j)],
      #                                                        eigenVectorsDict[j], eigenVectorsDict[j])
  return ret


def printNewFCIDUMP(oneModeDict: Dict[int, np.array], twoModeDict: Dict[Set, np.array],
                    newFileName: str, printThreshold: float):
  """
  Dump new integrals to file
  """
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  numberOfModals = resultFile.getNumberOfModals()
  newFile = open(newFileName, mode="w")
  # One-body part
  for iKey, iVal in oneModeDict.items():
    nModals = numberOfModals[iKey]
    for iRow in range(iVal.shape[0]):
      for iCol in range(iVal.shape[1]):
        if abs(iVal[iRow, iCol]) > printThreshold:
          newLine = str(iKey+1)+"-"+str(iRow)+"    "+str(iKey+1)+"-"+str(iCol)+"      "+str(iVal[iRow, iCol])+"\n"
          newFile.write(newLine)
  # Two-body part
  for iKey, iVal in twoModeDict.items():
    nModalsI = numberOfModals[iKey[0]]
    nModalsJ = numberOfModals[iKey[1]]
    for iRow in range(iVal.shape[0]):
      for iCol in range(iVal.shape[1]):
        for jRow in range(iVal.shape[2]):
          for jCol in range(iVal.shape[3]):
            if abs(iVal[iRow, iCol, jRow, jCol]) > printThreshold:
              newLine = str(iKey[0]+1)+"-"+str(iRow)+"    "+str(iKey[0]+1)+"-"+str(iCol)+"    "\
                       +str(iKey[1]+1)+"-"+str(jRow)+"    "+str(iKey[1]+1)+"-"+str(jCol)+"    "+str(iVal[iRow, iCol, jRow, jCol])+"\n"
              newFile.write(newLine)
  newFile.close()


def plotNOONs(inputFileName: str):
  """
  Plot the NOONs on logarithmic scale
  """
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  oneModeRDM = resultFile.getOneModeRDM()
  numberOfModals = resultFile.getNumberOfModals()
  overallCont = 0
  fix, ax = plt.subplots()
  for i in numberOfModals:
    modeSpecificRDM = np.copy(oneModeRDM[overallCont: overallCont+i, overallCont: overallCont+i])
    eigenValues, eigenVectors = np.linalg.eigh(modeSpecificRDM)
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    xVal = np.arange(0, i)
    labelSpecific="Mode "+str(i)
    ax.set_yscale('log')
    ax.plot(xVal, eigenValues, label=labelSpecific)
    overallCont += i
  plt.show()


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--resultfile", type=str, help="Input result file", required=True)
  parser.add_argument("-n", "--nummodes", type=int, help="Number of vibrational modes", required=True)
  parser.add_argument("-b", "--basis",  type=int, nargs="+", help="Number of basis functions per mode", required=True)
  parser.add_argument("-f", "--fcidump", type=str, help="Name of the FCIDUMP file")
  parser.add_argument("-t", "--threshold", type=float, help="Threshold for the truncation of the NO basis", default=0.)
  parser.add_argument("-p", "--print", action="store_true", help="If present, prints the one-particle RDM and the corresponding eigenvalues")
  parser.add_argument("--printThreshold", type=float, help="Threshold for the printin", default=1.0E-16)
  args = parser.parse_args()
  inputFileName = args.resultfile
  fcidumpFile = args.fcidump
  threshold = args.threshold
  numModes = args.nummodes
  numBasis = args.basis
  doPrint = args.print
  printThreshold = args.printThreshold
  printingThreshold = args.printThreshold
  #
  if doPrint:
    plotOneModeRDM(inputFileName)
    plotNOONs(inputFileName)
  # Pruning
  if fcidumpFile:
    print("Extracting integrals")
    oneBodyDict, twoBodyDict = extractIntegrals(fcidumpFile, numModes, numBasis)
    print("Rotating one-body")
    oneBodyRotated = rotateOneBody(inputFileName, oneBodyDict, threshold)
    print("Rotating two-body")
    twoBodyRotated = rotateTwoBody(inputFileName, twoBodyDict, threshold)
    print("Printing FCIDUMP")
    printNewFCIDUMP(oneBodyRotated, twoBodyRotated, fcidumpFile+"_NaturalOrbitals", printThreshold)
