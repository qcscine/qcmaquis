#!/usr/bin/env python2
#-*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
#*               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
#*               2014-2015 by Yingjin Ma <yma@ethz.ch>
#*               2018 by Leon Freitag <lefreita@ethz.ch>
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


## This script imports RDM derivatives from ALPS/QCMaquis results.h5 file
## and saves them to files with name (one|two)rdmderivative.[MPS parameter index]
## The output of this script is currently designed to be compatible with
## Yingjin's MPS linear response program, but this may change in the future.

import sys
import pyalps
import numpy as np

# dependency on pytools

# dependency on rdmsave_su2_onlyone.py


def save_1rdmderiv(rdm,filename_prefix):
    fmt = '%14.14e'

    L = int(rdm.props["L"])

    # dictionary of files. the file list for each derivative is maintained
    # in a dictionary

    f = {}

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[3]
        j = lab[4]

        filename="%s.%i-%i-%i" % (filename_prefix, lab[0], lab[1], lab[2])

        # open the file if f[filename] doesn't exist yet
        # and write L into it

        if (f.get(filename) == None):
            f[filename] = f.get(filename,open(filename,'w'))
            f[filename].write(str(L)+'\n')

        f[filename].write(str(i)+' '+str(j)+' '+str(fmt%val)+'\n')
        if (i!=j):
            f[filename].write(str(j)+' '+str(i)+' '+str(fmt%val)+'\n')

    # close all files
    for fi in f.itervalues():
        fi.close()

    return L

def load_1rdmderiv(inputfile,meas):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what=meas)[0][0]
    return rdm

def load_2rdmderiv(inputfile,meas):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what=meas)[0][0]
    rdm.y[0] = 0.5 * rdm.y[0]
    return rdm

def save_2rdmderiv(rdm,L,filename_prefix):
    fmt = '%14.14e'

    # dictionary of files. the file list for each derivative is maintained
    # in a dictionary

    f = {}

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[3]
        j = lab[4]
        k = lab[5]
        l = lab[6]

        filename="%s.%i-%i-%i" % (filename_prefix, lab[0], lab[1], lab[2])

        if (f.get(filename) == None):
            f[filename] = f.get(filename,open(filename,'w'))
            f[filename].write(str(L)+'\n')
        f[filename].write(str(i)+' '+str(j)+' '+str(k)+' '+str(l)+' '+str(fmt%val)+'\n')

    # close all files
    for fi in f.itervalues():
        fi.close()

if __name__ == '__main__':
    inputfile = sys.argv[1]

    for c in ['L','R']:
      rdm1=load_1rdmderiv(inputfile,'onerdmderiv%s'%c)
      L=save_1rdmderiv(rdm1,'onerdmderivative%s'%c)

      rdm2 = load_2rdmderiv(inputfile,'twordmderiv%s'%c)
      save_2rdmderiv(rdm2,L,'twordmderivative%s'%c)

