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
from argparse import ArgumentParser

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
#from matplotlib import rcParams

import extrapolate_base


#rcParams['text.usetex'] = True
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 15}
#matplotlib.rc('font', **font)
#rcParams.update({'figure.autolayout': True})

rc('axes', unicode_minus=False) # supposed to eliminate Glyph warning, but doesn't work for TeX
rc('font', **{'family':'serif','serif':['ComputerModernRoman'],'size':20})
rc('text', usetex=True)
 
def plot(flist, sweepnr=None):

    extr = extrapolate_base.energy2y(flist, sweepnr)
    xdata = extr.xdata()
    ydata = extr.ydata()

    fit = extr.linfit()
    xf = np.linspace(1e-16,max(xdata),100)
    yf = fit(xf)

    print "{:30s} {:13.9f}   {:12.1f}    {:5s}  {:3s}".format("extrapolation", fit(0), 0, "infty", "")
    print "{:30s} {:14.9f}".format("extrapolation error", fit(0)-ydata[0])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    dots = plt.plot(xdata, ydata, 'o', xf, yf, '-')

    plt.ticklabel_format(style='sci', axis='x', scilimits=(3,4), useOffset=False)
    xfmt = plt.ScalarFormatter(useOffset=False, useMathText=True)
    xfmt.set_scientific(True)
    ax.yaxis.set_major_formatter(xfmt)
    
    plt.xlabel(r'\textbf{truncation error} $\varepsilon$')
    plt.ylabel(r'\textbf{Energy [Hartree]}')
    #plt.xlabel('truncation error $\\varepsilon$')
    #plt.ylabel('Energy [Hartree]')
    fig.subplots_adjust(left=0.2)

    def autolabel(xd,yd,labels):
        # attach some text labels
        yshift = (max(yd) - min(yd)) * 0.05
        xshift = (max(xd) - min(xd)) * 0.2
        for x,y,l in zip(xd, yd, labels):
            #item = ax.text(x - xshift, y, '$m=%d$'%l,
            #            ha='right', va='center')
            item = ax.text(x, y, '$m=%d$'%l,
                        ha='center', va='bottom')
            item.set_fontsize(18)

    autolabel(xdata, ydata, extr.m())

    # add the extrapolated energy to the plot"
    ext_note = ax.text(max(xf)/2.0 + (max(xf)-min(xf)) * 0.05, fit(0),
                       'extrapolation:\n$%.6f$'%fit(0), ha='left', va='center')
    # position below is yaxis 
    #ext_note = ax.text(0, fit(0), '$%.6f$'%fit(0), ha='right', va='top')

    plt.savefig('e.svg')



if __name__=='__main__':
    parser = ArgumentParser(description='extrapolate DMRG results to zero truncation error')
    parser.add_argument('files', nargs='+', help="files to be loaded")
    parser.add_argument('-s', '--sweep', default=None, type=int, help="analyze at specific sweep")
    args=parser.parse_args()

    plot(args.files, args.sweep)
