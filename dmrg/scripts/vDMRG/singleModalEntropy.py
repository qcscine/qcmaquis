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
  doTwoDimensional = args.twodimensional
  doThreeDimensional = args.threedimensional
  #
  if doThreeDimensional:
    plotThreeDimensionalOneModal(inputFileName)
  if doTwoDimensional:
    plotTwoDimensionalOneModal(inputFileName)
