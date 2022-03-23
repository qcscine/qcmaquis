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
from VibrationalResultFile import ResultFileVibrationaleMeasurement

def plotMutualInformation(mutualInformation, oneModalEntropy, numModals, L,
                          scaling_area=400, alpha_scale=0.5):
    """
    Generates the plot associated with the entanglement diagram.
    """
    # Sets general variables
    plt.figure()
    theta = np.zeros(L)
    r = np.zeros(L)
    labels = np.zeros(L)
    area = np.zeros(L)
    slice_ = 2*math.pi/L
    nGridForCircles = 2000
    rInnerCircle = 0.9
    rOuterCircle = 1.1
    # Prepares the data for the single-orbital entropies
    for i in range(L):
        theta[i] = i * slice_ + slice_/2
        r[i] = 1.0
        labels[i] = i
        area[i] = oneModalEntropy[i]*scaling_area
    for iCont in range(len(area)):
      if abs(area[iCont]) < 1.0E-10:
        area[iCont] = 0
    ax = plt.subplot(111, polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi/2.0)
    # Grid that generates the Outer and Inner circles
    r_outer = np.zeros(nGridForCircles) + rInnerCircle
    r_inner = np.zeros(nGridForCircles) + rOuterCircle
    theta_lines = np.linspace(0, 2*math.pi, nGridForCircles)
    theta_marks = []
    theta_width = []
    alphas = []
    modesCont = 0
    iCont = 0
    # Data for the boundaries between the modals for different modes
    for i in numModals:
        theta_marks.append(2*modesCont*math.pi/L)
        theta_width.append(2.*i*math.pi/L)
        modesCont += i
        iCont += 1
        alphas.append(iCont*alpha_scale/len(numModals))
    for i in range(len(theta_marks)):
        plt.bar(theta_marks[i], 0.2, width=theta_width[i], bottom=0.9, color=(0, 0.45, 0.), alpha=alphas[i], zorder=0, align='edge')
    # black lines to separate the green areas
    plt.polar(theta_lines, r_outer, c="Black", linewidth=2, zorder=0)
    plt.polar(theta_lines, r_inner,c="Black", linewidth=2, zorder=0)
    for i in theta_marks:
        ax.plot((i, i), (rInnerCircle, rOuterCircle), c="Black", linewidth = 1, zorder =0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(b=False)
    plt.scatter(theta, r, c=(0.9, 0.17, 0.0044), s=area, zorder=1)
    legendlines = {}
    # # Put number for each mode
    for i in range(L):
        # We avoid printing the number so that the plot is not overcrowded
        # plt.text(theta[i], (r[i]+0.18), int(labels[i]), size='large', ha='center', va='center')
        for j in range(i, L):
            x = [theta[i], theta[j]]
            y = [1,1]
            Iij = mutualInformation[i, j]
            if Iij >= 0.1:
                line = lines.Line2D(x, y, linewidth=Iij, color='black', linestyle='-', alpha=1, label='0.1')
                legendlines['0.1'] = line
                ax.add_line(line)
            elif Iij >= 0.01:
                line = lines.Line2D(x, y, linewidth=10*Iij, color='gray', linestyle='-', alpha=1, label='0.01')
                legendlines['0.01'] = line
                ax.add_line(line)
            elif Iij >= 0.001:
                line = lines.Line2D(x, y, linewidth=0.1, color='grey', linestyle=':', alpha=1, label='0.001')
                legendlines['0.001'] = line
                ax.add_line(line)
    plt.tight_layout(h_pad=0.5)
    plt.subplots_adjust(bottom=0.2)
    ax.legend(legendlines.values(),[l.get_label() for l in legendlines.values()], bbox_to_anchor=(0.00,1.0),
              fancybox=True, shadow=True)
    plt.show()

def plotEntanglementDiagram(inputFileName: str):
  resultFile = ResultFileVibrationaleMeasurement(inputFileName)
  oneModalEntropy = resultFile.getOneModalEntropy()
  mutualInformation = resultFile.getMutualInformation()
  numberOfModals = resultFile.getNumberOfModals()
  L = resultFile.getLatticeSize()
  plotMutualInformation(mutualInformation, oneModalEntropy, numberOfModals, L)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-r", "--resultfile", type=str, help="Input result file", required=True)
  args = parser.parse_args() 
  inputFileName = args.resultfile
  plotEntanglementDiagram(inputFileName)
