#!/usr/bin/env python
#-*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
#*               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

import sys
import pyalps

import numpy as np

#import numpy as np
def load_2rdm(inputfile):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what='twoptdm')[0][0]
    rdm.y[0] = 0.5 * rdm.y[0]
    # uncomment for CASPT2 comparison
    # rdm.y[0] = rdm.y[0]
    return rdm

def load_2rdm_matrix(rdm):
    L = rdm.props['L']
    odm = np.zeros([L,L,L,L])

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[0]
        j = lab[1]
        k = lab[2]
        l = lab[3]

        odm[i,j,k,l] = val

        if l!=k:
            odm[j,i,l,k] = val

        if not min(i,j) == min(l,k):
            odm[k,l,i,j] = val
            if k!=l:
                odm[l,k,j,i] = val

    return odm

def print_2rdm(rdm):
    fmt = '% -016.10E'
    #fmt = '%e'

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[0]
        j = lab[1]
        k = lab[2]
        l = lab[3]

        if abs(val) == 0: continue
 
        print i+1,j+1,k+1,l+1, fmt%val

        # print duplicates
        if l!=k:
            print j+1,i+1,l+1,k+1, fmt%val

        if not min(i,j) == min(l,k):
            print k+1,l+1,i+1,j+1, fmt%val
            if k!=l:
                print l+1,k+1,j+1,i+1, fmt%val

if __name__ == '__main__':
    inputfile = sys.argv[1]

    rdm = load_2rdm(inputfile)
    print_2rdm(rdm)
