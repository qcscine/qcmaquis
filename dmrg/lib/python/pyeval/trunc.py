#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
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

# small interactive script to examin truncated weights in sweeps
# Usage: trunc.py h5-result-file

import pyalps
import pydmrg
import numpy as np
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot
import plotutil


def plot(fname):

    ret = pydmrg.LoadDMRGSweeps([fname],['TruncatedWeight'])

    sweeps = []
    for sw in ret[0]:
        sweeps += list(sw[0].y)

    print "total number of values in sweep", len(sweeps)

    # Get length of 1 sweep
    props = ret[0][0][0].props
    L = props['L']
    swl = 2*(L-1)

    xdata = np.arange(len(sweeps), dtype=np.float64)/swl
    ydata = np.array(sweeps)

    pdata = pyalps.DataSet()
    pdata.x = xdata
    pdata.y = ydata
    pdata.props['label'] = "Truncated weight"

    fig = plt.figure()
    plt.ylabel('Truncated weight')
    plt.xlabel('sweep')
    pyalps.plot.plot(pdata)
    plt.legend()

    for ax in fig.axes:
        ax.callbacks.connect('xlim_changed', plotutil.on_xlim_changed)

    plt.show()

if __name__=='__main__':
    rfile = sys.argv[1]
    plot(rfile)
