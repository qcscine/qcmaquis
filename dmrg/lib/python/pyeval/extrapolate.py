#!/usr/bin/env python
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2015 Laboratory for Physical Chemistry, ETH Zurich
#*               2012-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
#*
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

# extrapolate the truncated weight to zero
# Usage: extrapolate.py [list of hdf5 result files]

import sys
import os
import subprocess
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

import plotutil
import pyalps.plot
import extrapolate_base


#rcParams['text.usetex'] = True
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}
matplotlib.rc('font', **font)
#rcParams.update({'figure.autolayout': True})
 
def plot(flist):

    extr = extrapolate_base.energy2y(flist)
    xdata = extr.xdata()[1:]
    ydata = extr.ydata()[1:]

    fit = extr.linfit()
    xf = np.linspace(1e-16,max(xdata),100)
    yf = fit(xf)

    print "extrapolated to", fit(0) 

    fig = plt.figure()
    ax = fig.add_subplot(111)
    dots = plt.plot(xdata, ydata, 'o', xf, yf, '-')

    #plt.ticklabel_format(style='sci', axis='x', scilimits=(3,4), useOffset=False)
    xfmt = plt.ScalarFormatter(useOffset=False, useMathText=True)
    xfmt.set_scientific(True)
    ax.yaxis.set_major_formatter(xfmt)
    
    plt.xlabel('truncation error $\\varepsilon$')
    plt.ylabel('Energy [Hartree]')
    fig.subplots_adjust(left=0.2)
    plt.legend()

    def autolabel(xd,yd,labels):
        # attach some text labels
        shift = (max(yd) - min(yd)) * 0.04
        for x,y,l in zip(xd, yd, labels):
            item = ax.text(x, y + shift, 'm=%d'%l,
                        ha='center', va='bottom')
            item.set_fontsize(12)

    autolabel(xdata, ydata, extr.m())
    plt.savefig('e.svg')



if __name__=='__main__':
    files = sys.argv[1:]
    plot(files)
