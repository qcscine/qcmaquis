#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2022 Laboratory for Physical Chemistry, ETH Zurich
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
import numpy as np
import matplotlib.pyplot as plt
from VibrationalResultFile import ResultFileVibrationaleMeasurement

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

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--resultfile", type=str, help="Input result file", required=True)
  args = parser.parse_args() 
  inputFileName = args.resultfile
  plotOneModeRDM(inputFileName)


