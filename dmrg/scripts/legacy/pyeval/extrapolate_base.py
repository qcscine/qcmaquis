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

# base class for extrapolations

import pyalps
import pyalps.fit_wrapper as fw
import numpy as np
import math
import os
import glob
import pydmrg

import warnings
warnings.filterwarnings("ignore")

def getData(flist, what = ['Energy'], sweepnr = None):

    # [file][sweep][observable]
    truncs = pydmrg.LoadDMRGSweeps(flist, ['TruncatedWeight'])
    energs = pydmrg.LoadDMRGSweeps(flist, what)
    props = truncs[0][0][0].props

    if sweepnr is None: sweepnr = min([len(f) - 1 for f in truncs])

    xdata = []
    ydata = []
    mdata = []
    print "{:30s} {:15s} {:15s}   m     sweep".format("result file", "energy", "truncation error")
    print "-----------------------------------------------------------------------------"
    for tr,e in zip(truncs, energs):
        props = tr[0][0].props
        sweeps = []
        bonds = []
        try:
            # check for reverse m data
            sweep_bond_dims = props['sweep_bond_dimensions']
            bond_dims = map(int, sweep_bond_dims.split(','))
            for i,b in enumerate(bond_dims[1:]+[0]):
                if bond_dims[i] > b and i < len(tr):
                    sweeps.append(i)
                    bonds.append(bond_dims[i])
                
        except KeyError:
            if len(tr[sweepnr][0].y) != 2 * (props['L'] - 1):
                print "\nWARNING: data in sweep", sweepnr, "incomplete\n"
            sweeps.append(len(tr)-1)
            bonds.append(int( e[sweepnr][0].props['max_bond_dimension'] + 0.5))

        for s,b in zip(sweeps, bonds):
            xdata.append(max(tr[s][0].y))
            ydata.append(min( e[s][0].y))
            mdata.append(b)
            print "{:30s} {:13.9f}   {:.6e}    {:5d}  {:3d}".format(props['resultfile'], ydata[-1], xdata[-1], mdata[-1], s)

    print ""

    return (xdata, ydata, mdata)

class extrapolator(object):

    def __init__(self, flist, what = ['Energy'], sweepnr = None):
        self.extflag = True
        self.flist = flist

        self._xdata, self._ydata, self._m = getData(flist, what, sweepnr)

        self._xdata.sort()
        self._ydata.sort()
        self._m.sort(reverse=True)

    def m(self):
        return self._m

    def linfit(self, num_points=None):
        if num_points is None:
            num_points = len(self._ydata)
        ff = lambda self, x, pars: pars[0]()*x + pars[1]()
        pars = [fw.Parameter(1.), fw.Parameter(-1.)]
        fw.fit(None, ff, pars, np.array(self._ydata[0:num_points]), np.array(self._xdata[0:num_points]))
        
        fit_func = lambda x: pars[0].get()*x + pars[1].get()
        return fit_func

    def expfit(self, num_points=None):
        if num_points is None:
            num_points = len(self._ydata)

        lf = self.linfit()

        ff = lambda self, x, pars: pars[0]()*(x**pars[1]()) + pars[2]()
        #pars = [fw.Parameter(1.0), fw.Parameter(1.01), fw.Parameter(self._ydata[0])]
        pars = [fw.Parameter(lf(1)-lf(0)), fw.Parameter(1.01), fw.Parameter(lf(0))]
        fw.fit(None, ff, pars, np.array(self._ydata[0:num_points]), np.array(self._xdata[0:num_points]))
        
        fit_func = lambda x: pars[0].get()*(x**pars[1].get()) + pars[2].get()
        return fit_func

    data = []
    _xdata, _ydata = [], []

class energy2y(extrapolator):
    def __init__(self, flist, sweepnr = None):
        super(energy2y, self).__init__(flist, ['Energy'], sweepnr)
    
    def get(self, fdata):
        return fdata[0].y[0]

    def ydata(self):
        return self._ydata

    def xdata(self):
        return self._xdata


