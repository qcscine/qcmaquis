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

#import numpy as np
def load_2rdm(inputfile):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what='twoptdm')[0][0]
    return rdm

def print_2rdm(rdm,tag):
    #fmt = '% -016.10E'
    #fmt = '%e'
    tag1 = tag
    tag2 = tag

    f=open('twoparticle.rdm.%s.%s' % (tag1,tag2),'w')

    for lab, val in zip(rdm.x, rdm.y[0]):
        m = lab[0]
        n = lab[1]
        o = lab[2]
        p = lab[3]

        if abs(val.real) and abs(val.imag) == 0: continue
 
        dump_element(f, val,m,n,o,p)
        dump_element(f, val,o,p,m,n)
        dump_element(f,-val,m,p,o,n)
        dump_element(f,-val,o,n,m,p)

    f.close()

def dump_element(f,value,i,j,k,l):
    
    print i,j,k,l, "\t", (value.real, value.imag)
    fmt  = '% -020.14E'
    f.write(str(i)+' '+str(j)+' '+str(k)+' '+str(l)+'  '+str(fmt%value.real)+'  '+str(fmt%value.imag)+'\n')

if __name__ == '__main__':
    inputfile = sys.argv[1]
    tag       = sys.argv[2]

    rdm = load_2rdm(inputfile)
    print_2rdm(rdm,tag)
