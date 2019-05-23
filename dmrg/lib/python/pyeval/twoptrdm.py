#!/usr/bin/env python2
#-*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
#*               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
#*               2018-2019 by Stefan Knecht    <stknecht@ethz.ch>
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

import input as DmrgInput
import numpy as np

def load_2rdm(inputfile):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what='twoptdm')[0][0]
    rdm.y[0] = 0.5 * rdm.y[0]
    # uncomment for CASPT2 comparison
    # rdm.y[0] = rdm.y[0]
    return rdm

def load_2rdm_dbg(inputfile):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what='twoptdm')[0][0]
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

def print_2rdm_dbg(rdm,tag):
    #fmt = '% -016.10E'
    #fmt = '%e'
    tag1 = tag
    tag2 = tag

    f=open('twoparticle.rdm.%s.%s' % (tag1,tag2),'w')
    b=open('extDMRG_%s_%s.rdm2' % (tag1,tag2),'w')

    for lab, val in zip(rdm.x, rdm.y[0]):
        m = lab[0]
        n = lab[1]
        o = lab[2]
        p = lab[3]

        #print "raw data --> ", m,n,o,p, "\t", (val.real, val.imag)

        if abs(val.real) == 0 and abs(val.imag) == 0: continue

        # dump element (and permutations) for DIRAC
        dump_element(f, val,m,n,o,p)
        dump_element(f, val,o,p,m,n)
        dump_element(f,-val,m,p,o,n)
        dump_element(f,-val,o,n,m,p)

        # dump element for BAGEL as i+ j+ k l --> ikjl
        # permuattions are taken care of internally in BAGEL
        dump_element(b, val,m+1,o+1,n+1,p+1)

    f.close()
    b.close()

def dump_element(f,value,i,j,k,l):
    
    #print i,j,k,l, "\t", (value.real, value.imag)
    fmt  = '% -020.14E'
    f.write(str(i)+' '+str(j)+' '+str(k)+' '+str(l)+'  '+str(fmt%value.real)+'  '+str(fmt%value.imag)+'\n')

if __name__ == '__main__':
    inputfile = sys.argv[1]

    props = DmrgInput.loadProperties(inputfile)

    if props["symmetry"] == "u1dg":
        tag = sys.argv[2]
        rdm = load_2rdm_dbg(inputfile)
        print_2rdm_dbg(rdm,tag)
    else:
        rdm = load_2rdm(inputfile)
        print_2rdm(rdm)


