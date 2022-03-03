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

def plotThreeDimensionalOneModal(inputFileName: str, view1=30, view2=240):
  """
  Plots the one-modal entropy as a 3D plot.
  This plot format is more conventient than the entanglement diagram
  because it separates the different modes.
  """
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.view_init(view1, view2)
  oneModalEntropies = resultFile.getSeparatedOneModalEntropy()
  #
  for idx, iModeEntropies in enumerate(oneModalEntropies):
    xAxis = range(len(iModeEntropies))
    labelName = "Mode "+str(idx)
    ax.plot(xAxis, iModeEntropies, zs=idx, zdir='y', ls='-', label=labelName)
  # ax.legend(fontsize=30, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.15))
  ax.set_xlabel('\n\n\nMode')
  ax.set_ylabel('\n\n\nModal index')
  ax.set_zlabel('\n\n\nSingle modal entropy')
  # Forces x and y axes to be integers
  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  ax.yaxis.set_major_locator(MaxNLocator(integer=True))
  # Set the boundaries for the axes
  ax.set_xlim((0, max([len(i) for i in oneModalEntropies])))
  ax.set_ylim((0, len(oneModalEntropies)))
  #ax.set_zlim((-5, 1))
  ax.tick_params(axis='both', which='major', pad=15)
  # Uses logarithmic scale for z axis
  # ax.set_zscale('log')
  plt.show()

def plotTwoDimensionalOneModal(inputFileName: str, view1=30, view2=240):
  """
  Plots the one-modal entropy as a logarithmic 2D plot.
  """
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  oneModalEntropies = resultFile.getSeparatedOneModalEntropy()
  # Matplotlib
  fig = plt.figure()
  ax = fig.add_subplot(111)
  #
  for idx, iModeEntropies in enumerate(oneModalEntropies):
    xAxis = range(len(iModeEntropies))
    labelName = "Mode "+str(idx)
    ax.plot(xAxis, iModeEntropies, ls='-', label=labelName, marker="o")
  # ax.legend(fontsize=30, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.15))
  ax.set_xlabel('\n\n\nModal index')
  ax.set_ylabel('\n\n\nSingle modal entropy')
  # Forces x and y axes to be integers
  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  # Uses logarithmic scale for z axis
  ax.set_yscale('log')
  # Plots the legend
  plt.legend()
  plt.show()


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--resultfile", type=str, help="Input result file", required=True)
  parser.add_argument("--threedimensional", action="store_true", help="If present, do a 3D plot of the vibrational entanglement")
  parser.add_argument("--twodimensional", action="store_true", help="If present, do a 2D logarithmic plot of the entanglement")
  args = parser.parse_args() 
  inputFileName = args.resultfile
  doDiagram = args.diagram
  doTwoDimensional = args.twodimensional
  doThreeDimensional = args.threedimensional
  #
  if doThreeDimensional:
    plotThreeDimensionalOneModal(inputFileName)
  if doTwoDimensional:
    plotTwoDimensionalOneModal(inputFileName)
