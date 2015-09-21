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

# base class for extrapolations

import pyalps
import pyalps.fit_wrapper as fw
import numpy as np
import math
import os
import glob
import pydmrg

class extrapolator(object):

    def __init__(self, flist, what = ['Energy'], sweepnr = None):
        self.extflag = True
        self.flist = flist

        print "files:", self.flist

        truncs = pydmrg.LoadDMRGSweeps(self.flist, ['TruncatedWeight'])

        if sweepnr is None:
            sweepnr = min([len(f) - 1 for f in truncs])

        print "extrapolating data at sweep", sweepnr

        props = truncs[0][0][0].props

        # TODO : find out why _xdata survives even if derived class is deleted
        # load truncated weights
        self._xdata = []
        for d in truncs:
            if len(d[sweepnr][0].y) != 2 * (props['L'] - 1):
                print "\nWARNING: data in sweep", sweepnr, "incomplete\n"

            self._xdata.append( max(d[sweepnr][0].y) )

        # load energies and bond dimensions
        self._ydata = []
        self._m = []
        energs = pydmrg.LoadDMRGSweeps(flist, what)
        for f in energs:
            # get min energy from last sweep
            self._ydata.append(min(f[sweepnr][0].y))
            self._m.append(int(f[sweepnr][0].props['max_bond_dimension'] + 0.5))

        self._xdata.sort()
        self._ydata.sort()
        self._m.sort(reverse=True)

        print "truncs", self._xdata
        print "energies", self._ydata

    def m(self):
        return self._m

    def ext(self, x=None):
        return self.extrapolate_to_zero(self._xdata, self._ydata)

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

    def extrapolate_to_zero(self, xdata, ydata):
        x0 = xdata[0]
        y0 = ydata[0]
        dx = xdata[1] - x0
        dy = ydata[1] - y0
        b = y0 - dy/dx*x0
        print "extrapol diff", abs(b-y0)
        return b

    data = []
    _xdata, _ydata = [], []

class energy2y(extrapolator):
    def __init__(self, flist, sweepnr = None):
        super(energy2y, self).__init__(flist, ['Energy'], sweepnr)
    
    def get(self, fdata):
        return fdata[0].y[0]

    def ydata(self):
        return [self.ext()] + self._ydata

    def xdata(self):
        return [0.0] + self._xdata


