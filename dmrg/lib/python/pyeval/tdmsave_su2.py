#!/usr/bin/env python2
#-*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
#*               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
#*               2014-2015 by Yingjin Ma <yma@ethz.ch>
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

# dependency on pytools

# dependency on rdmsave_su2_onlyone.py

def save_1rdm(rdm):
    fmt = '%14.14e'

    L = int(rdm.props["L"])
    mat = np.zeros((L,L))

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[0]
        j = lab[1]
        mat[i,j] = val;
        #mat[j,i] = val;

    f=open('oneparticle.tdm','w')
    f.write(str(L)+'\n')

    for i in range(L):
        for j in range(L):
            f.write(str(i)+' '+str(j)+' '+str(fmt%mat[i,j])+'\n')
    f.close()

    return L

def load_1rdm(inputfile):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what='transition_oneptdm')[0][0]
    return rdm

def load_2rdm(inputfile):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what='transition_twoptdm')[0][0]
    rdm.y[0] = 0.5 * rdm.y[0]
    return rdm

def save_2rdm(rdm,L):
    fmt = '%14.14e'

    f=open('twoparticle.tdm','w')
    f.write(str(L)+'\n')
    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[0]
        j = lab[1]
        k = lab[2]
        l = lab[3]
        # work around an indexing bug in SU2U1 evaluation
        if (((i == k) and (i != j) and (i != l)) != ((j == l) and (i != j) and (j != k))):
            f.write(str(k)+' '+str(l)+' '+str(i)+' '+str(j)+' '+str(fmt%val)+'\n')
            if (k != l):
                f.write(str(l)+' '+str(k)+' '+str(j)+' '+str(i)+' '+str(fmt%val)+'\n')
        else:
            f.write(str(i)+' '+str(j)+' '+str(k)+' '+str(l)+' '+str(fmt%val)+'\n')
            if (l != k):
                f.write(str(j)+' '+str(i)+' '+str(l)+' '+str(k)+' '+str(fmt%val)+'\n')


    f.close()

if __name__ == '__main__':
    inputfile = sys.argv[1]

    rdm1 =load_1rdm(inputfile)
    L=save_1rdm(rdm1)

    rdmT = load_2rdm(inputfile)
    save_2rdm(rdmT,L)

