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

def getMutualInformationLaplacian(iMatrix):
  """
  Gets the laplacian of the mutual information graph 
  """
  laplacian = np.zeros(iMatrix.shape)
  idx = 0
  while idx < laplacian.shape[0]:
    jdx = 0
    while jdx < laplacian.shape[1]:
      if abs(iMatrix[idx, jdx]) > 1.0E-10:
        laplacian[idx, idx] += abs(iMatrix[idx, jdx])
      jdx += 1
    idx += 1
  #
  for i in range(laplacian.shape[0]):
    for j in range(laplacian.shape[1]):
      if abs(iMatrix[i, j]) > 1.0E-10:
        laplacian[i, j] -= np.abs(iMatrix[i, j])
  return laplacian

def fiedler(laplacian):
  """
  Gets the Fiedler ordering
  """
  eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
  print(eigenvalues)
  idx = eigenvalues.argsort()
  sort_ev = eigenvalues[idx]
  sort_vec = eigenvectors[:, idx]
  iRow = 1
  for i in range(laplacian.shape[0]-1):
    if sort_ev[iRow+i] > 1.0E-10:
      iRow = i+1
      break
  return sort_vec[:, iRow]

def reorder(s1, iMatrix, ordering):
  """
  Reorders the single-orbital entropy and mutual information based on
  a new ordering.
  """
  news1 = s1[ordering]
  newMutualInformation = np.zeros(iMatrix.shape)
  idx = 0
  while idx < ordering.shape[0]:
    jdx = 0
    while jdx < ordering.shape[0]:
      newMutualInformation[idx, jdx] = iMatrix[ordering[idx], ordering[jdx]]
      jdx += 1
    idx += 1
  return news1, newMutualInformation

def penaltyFunction(mutualInfo):
  """Calculates the penalty function (that is minimized to get the Fiedler ordering"""
  rowBegin=0
  rowEnd = mutualInfo.shape[0]
  colBegin=0
  colEnd = mutualInfo.shape[1]
  cost = 0.
  idx = rowBegin
  while idx < rowEnd:
    jdx = colBegin
    while jdx < colEnd:
      cost = cost + mutualInfo[idx, jdx] * abs(idx - jdx)**2
      jdx += 1
    idx += 1
  return cost


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--resultfile", type=str, help="Input result file", required=True)
  args = parser.parse_args()
  resultFile = ResultFileVibrationaleMeasurement(args.resultfile)
  # Calculates the Fiedler ordering
  # plot mutual information without ordering
  oneModalEntropy = resultFile.getOneModalEntropy()
  mutualInformation = resultFile.getMutualInformation()
  laplacian = getMutualInformationLaplacian(mutualInformation)
  fiedlerVector = fiedler(laplacian)
  # Order
  order = fiedlerVector.argsort()
  ofv = fiedlerVector[order]
  # plotMutualInformation(mutualInformation, oneModalEntropy, resultFile.getNumberOfModals(), resultFile.getLatticeSize())
  plt.matshow(mutualInformation, cmap='GnBu')
  new_s1, new_mutinf = reorder(oneModalEntropy, mutualInformation, order)
  # plotMutualInformation(new_mutinf, new_s1, resultFile.getNumberOfModals(), resultFile.getLatticeSize())
  plt.matshow(new_mutinf, cmap='GnBu')
  cost_old = penaltyFunction(mutualInformation)
  cost_new = penaltyFunction(new_mutinf)
  print("Old cost ", cost_old, " new cost ", cost_new)
  originalOrder = [i for i in range(len(order))]
  new_order = [originalOrder[order[i]] for i in range(len(order))]
  outputString = ""
  for iString in new_order:
    outputString += str(iString)+","
  print("Fiedler ordering: ", outputString[:-1])
