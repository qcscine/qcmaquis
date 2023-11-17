#!/usr/bin/env python2
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


# usage: sweeps.py h5-result-file
# small interactive script to look at energies in sweeps


import pyalps
import pydmrg
import sys 
import numpy as np
import matplotlib.pyplot as plt 
import pyalps.plot
import plotutil

def plot(fnames):

    ret = pydmrg.LoadDMRGSweeps(fnames,['Energy'])

    plot_data = []
    for i,f in enumerate(ret):
        sweeps = []
        for sw in f:
            sweeps += list(sw[0].y)

        print "total number of sweep values", len(sweeps)
        props = f[0][0].props
        L = props['L']
        print "L", L
        swl = 2.0*(L-1)
        xdata = np.arange(len(sweeps))/swl
        ydata = np.array(sweeps)
      
        #print ydata
        print "Minimum energy:", np.min(ydata), "\n"

        pdata = pyalps.DataSet()
        pdata.x = xdata
        pdata.y = ydata
        pdata.props['label'] = "Energy " + props["resultfile"]
        plot_data.append(pdata)

    fig = plt.figure()
    plt.xlabel('Sweeps')
    plt.ylabel('Energy')
    pyalps.plot.plot(plot_data)

    #locs,labels = plt.xticks()
    #plt.xticks(locs, map(lambda x: "%g" % x, locs))
    #locs,labels = plt.yticks()
    #plt.yticks(locs, map(lambda x: "%.4f" % x, locs))

    #plt.gca().set_yticklabels(plt.gca().get_yticks())
    #plt.gca().ticklabel_format(useOffset=False)
    plt.legend()

    for ax in fig.axes:
        ax.callbacks.connect('xlim_changed', plotutil.on_xlim_changed)

    plt.show()
    

if __name__=='__main__':
    rfile = sys.argv[1:]
    plot(rfile)




