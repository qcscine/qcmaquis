#!/usr/bin/env python
#-*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2015 Laboratory for Physical Chemistry, ETH Zurich
#*               2015-2015 by Stefan Knecht <stefan.knecht@phys.chem.ethz.ch>
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
def load_3rdm(inputfile):
    # load data from the HDF5 result file
    rdm =  pyalps.loadEigenstateMeasurements([inputfile], what='transition_threeptdm')[0][0]
    rdm.y[0] = rdm.y[0]
    return rdm

def print_3rdm(rdm,tag1,tag2):
    fmt = '% -020.14E'

    # fix correct naming with tag1
    f=open('threeparticle.tdm.%s.%s' % (tag1,tag2),'w')
    print 'start of 3-TDM for states %s --> %s' %(tag1,tag2)

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[0]+1
        j = lab[1]+1
        k = lab[2]+1
        l = lab[3]+1
        m = lab[4]+1
        n = lab[5]+1

        if abs(val) == 0: continue
        if i == j and i == k: continue 
        if l == m and l == n: continue

        print i,j,k,l,m,n, fmt%val

        # 6 permutations (dealing with the transition 3-RDM)
        dump_element(f,val,i,j,k,l,m,n)
        dump_element(f,val,i,k,j,l,n,m)
        dump_element(f,val,j,i,k,m,l,n)
        dump_element(f,val,j,k,i,m,n,l)
        dump_element(f,val,k,i,j,n,l,m)
        dump_element(f,val,k,j,i,n,m,l)

    f.close()
    print 'end of 3-TDM for states %s --> %s' %(tag1,tag2)

def dump_element(f,value,i,j,k,l,m,n):
    
    fmt  = '% -020.14E'
    f.write(str(i)+' '+str(j)+' '+str(k)+' '+str(l)+' '+str(m)+' '+str(n)+' '+str(fmt%value)+'\n')

if __name__ == '__main__':

    rdm = load_3rdm(sys.argv[1])
    print_3rdm(rdm,sys.argv[2],sys.argv[3])
