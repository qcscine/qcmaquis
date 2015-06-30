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
def load_4rdm(inputfile):
    # load data from the HDF5 result file
    return pyalps.loadEigenstateMeasurements([inputfile], what='fourptdm')[0][0]

def print_4rdm(rdm):
    fmt = '% -016.10E'
    #fmt = '%e'

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[0]
        j = lab[1]
        k = lab[2]
        l = lab[3]
        m = lab[4]
        n = lab[5]
        o = lab[6]
        p = lab[7]

        if abs(val) == 0: continue

        print i+1,j+1,k+1,l+1,m+1,n+1,o+1,p+1, fmt%val

        valm05 = val*(-0.5)
        valm2  = val*(-2.0)

        # print permutations
        if i==j and k==l: # case 1: i==j k==l
           if m==n and o==p:
               print i+1,j+1,k+1,l+1,o+1,o+1,m+1,m+1, fmt%val

               print i+1,j+1,k+1,l+1,m+1,o+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,o+1,m+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,m+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,o+1,m+1, fmt%valm05
           elif m==n and o!=p:
               # regular permutations
               print i+1,j+1,k+1,l+1,m+1,m+1,p+1,o+1, fmt%val
               print i+1,j+1,k+1,l+1,o+1,p+1,m+1,m+1, fmt%val
               print i+1,j+1,k+1,l+1,p+1,o+1,m+1,m+1, fmt%val
               # remaining permutations
               print i+1,j+1,k+1,l+1,o+1,m+1,m+1,p+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,p+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,p+1,m+1,m+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,p+1,m+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,o+1,p+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,p+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,p+1,m+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,o+1,m+1,p+1, fmt%valm05
           elif m!=n and o==p:
               # regular permutations
               print i+1,j+1,k+1,l+1,n+1,m+1,o+1,o+1, fmt%val
               print i+1,j+1,k+1,l+1,o+1,o+1,m+1,n+1, fmt%val
               print i+1,j+1,k+1,l+1,o+1,o+1,n+1,m+1, fmt%val
               # remaining permutations
               print i+1,j+1,k+1,l+1,m+1,o+1,o+1,n+1, fmt%valm05
               print i+1,j+1,k+1,l+1,n+1,o+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,o+1,n+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,n+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,n+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,n+1,m+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,o+1,n+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,n+1,o+1,m+1,o+1, fmt%valm05
           elif m!=n and o!=p:
               if n == o:
                    print i+1,j+1,k+1,l+1,n+1,n+1,m+1,p+1, fmt%valm2
                    print i+1,j+1,k+1,l+1,n+1,n+1,p+1,m+1, fmt%valm2
                    print i+1,j+1,k+1,l+1,m+1,p+1,n+1,n+1, fmt%valm2
                    print i+1,j+1,k+1,l+1,p+1,m+1,n+1,n+1, fmt%valm2

                    print i+1,j+1,k+1,l+1,p+1,n+1,n+1,m+1, fmt%val
                    print i+1,j+1,k+1,l+1,p+1,n+1,m+1,n+1, fmt%val
                    print i+1,j+1,k+1,l+1,m+1,n+1,p+1,n+1, fmt%val
                    print i+1,j+1,k+1,l+1,n+1,p+1,m+1,n+1, fmt%val
                    print i+1,j+1,k+1,l+1,n+1,m+1,p+1,n+1, fmt%val
                    print i+1,j+1,k+1,l+1,n+1,m+1,n+1,p+1, fmt%val
                    print i+1,j+1,k+1,l+1,n+1,p+1,n+1,m+1, fmt%val
               else:
                    if o > p:
                        print i+1,j+1,k+1,l+1,m+1,n+1,p+1,o+1, fmt%val
                        print i+1,j+1,k+1,l+1,n+1,m+1,p+1,o+1, fmt%val
                        print i+1,j+1,k+1,l+1,n+1,m+1,o+1,p+1, fmt%val
                        print i+1,j+1,k+1,l+1,o+1,p+1,m+1,n+1, fmt%val
                        print i+1,j+1,k+1,l+1,o+1,p+1,n+1,m+1, fmt%val
                        print i+1,j+1,k+1,l+1,p+1,o+1,n+1,m+1, fmt%val
                        print i+1,j+1,k+1,l+1,p+1,o+1,m+1,n+1, fmt%val

                        print i+1,j+1,k+1,l+1,m+1,p+1,n+1,o+1, fmt%valm2
                        print i+1,j+1,k+1,l+1,m+1,p+1,o+1,n+1, fmt%valm2
                        print i+1,j+1,k+1,l+1,p+1,m+1,o+1,n+1, fmt%valm2
                        print i+1,j+1,k+1,l+1,p+1,m+1,n+1,o+1, fmt%valm2
                        print i+1,j+1,k+1,l+1,n+1,o+1,p+1,m+1, fmt%valm2
                        print i+1,j+1,k+1,l+1,n+1,o+1,m+1,p+1, fmt%valm2
                        print i+1,j+1,k+1,l+1,o+1,n+1,m+1,p+1, fmt%valm2
                        print i+1,j+1,k+1,l+1,o+1,n+1,p+1,m+1, fmt%valm2
                    elif p > o:
                        print i+1,j+1,k+1,l+1,m+1,n+1,p+1,o+1, fmt%val
                        print i+1,j+1,k+1,l+1,n+1,m+1,p+1,o+1, fmt%val
                        print i+1,j+1,k+1,l+1,n+1,m+1,o+1,p+1, fmt%val
                        print i+1,j+1,k+1,l+1,o+1,p+1,m+1,n+1, fmt%val
                        print i+1,j+1,k+1,l+1,o+1,p+1,n+1,m+1, fmt%val
                        print i+1,j+1,k+1,l+1,p+1,o+1,n+1,m+1, fmt%val
                        print i+1,j+1,k+1,l+1,p+1,o+1,m+1,n+1, fmt%val

        elif i==j and k!=l: # case 2: i==j k!=l
           if m==n and o==p: # 2x2 equal
               print i+1,j+1,k+1,l+1,o+1,o+1,m+1,m+1, fmt%val

               print i+1,j+1,k+1,l+1,m+1,o+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,m+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,o+1,m+1,o+1, fmt%valm05
           if m==n and o!=p: # 1 equal (first 2)
               print i+1,j+1,k+1,l+1,m+1,o+1,m+1,p+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,p+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,m+1,p+1, fmt%valm05
               print i+1,j+1,k+1,l+1,p+1,m+1,o+1,m+1, fmt%valm05
           if m!=n and o==p: # 1 equal (latter 2)
               print i+1,j+1,k+1,l+1,n+1,m+1,o+1,o+1, fmt%val
           if m!=n and n!=o and o!=p and m!=p: # all different
               print i+1,j+1,k+1,l+1,n+1,m+1,o+1,p+1, fmt%val

        elif i!=j and j==k and k!=l: # case 3: i!=j j==k k!=l
           if m==n and o==p: # 2x2 equal
               print i+1,j+1,k+1,l+1,o+1,o+1,m+1,m+1, fmt%val
               print i+1,j+1,k+1,l+1,o+1,m+1,o+1,m+1, fmt%val
               print i+1,j+1,k+1,l+1,m+1,o+1,m+1,o+1, fmt%val

               print i+1,j+1,k+1,l+1,m+1,o+1,o+1,m+1, fmt%valm2
               print i+1,j+1,k+1,l+1,o+1,m+1,m+1,o+1, fmt%valm2
           if m!=n and o==p: # 1 equal (latter 2)
               print i+1,j+1,k+1,l+1,o+1,o+1,m+1,n+1, fmt%val
               print i+1,j+1,k+1,l+1,o+1,m+1,o+1,n+1, fmt%val
               print i+1,j+1,k+1,l+1,m+1,o+1,n+1,o+1, fmt%val

               print i+1,j+1,k+1,l+1,m+1,o+1,o+1,n+1, fmt%valm2
           if m!=n and o!=p and m==p: # 1 equal (first and last)
               print i+1,j+1,k+1,l+1,m+1,o+1,n+1,m+1, fmt%val
           if m!=n and n!=o and o!=p and m!=p: # all different
               print i+1,j+1,k+1,l+1,m+1,o+1,n+1,p+1, fmt%val

        elif i!=j and k==l: # case 4: i!=j k==l
           if m==n and o==p: # 2x2 equal
               print i+1,j+1,k+1,l+1,o+1,o+1,m+1,m+1, fmt%val

               print i+1,j+1,k+1,l+1,m+1,o+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,m+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,m+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,o+1,m+1,o+1, fmt%valm05
           if m==n and o!=p: # 1 equal (first 2)
               print i+1,j+1,k+1,l+1,m+1,m+1,p+1,o+1, fmt%val
           if m!=n and o==p: # 1 equal (latter 2)
               print i+1,j+1,k+1,l+1,m+1,o+1,n+1,o+1, fmt%valm05
               print i+1,j+1,k+1,l+1,m+1,o+1,o+1,n+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,n+1,o+1,m+1, fmt%valm05
               print i+1,j+1,k+1,l+1,o+1,n+1,m+1,o+1, fmt%valm05
           if m!=n and n!=o and o!=p and m!=p: # all different
               print i+1,j+1,k+1,l+1,m+1,n+1,p+1,o+1, fmt%val

        elif i!=j and k!=l and j!=k: # case 5: i!=j j!=k k!=l
           if (m==o and n==p) or (m==p and n==o):
               print i+1,j+1,k+1,l+1,n+1,m+1,p+1,o+1, fmt%val
           if m==n and o==p:
               print i+1,j+1,k+1,l+1,o+1,o+1,m+1,m+1, fmt%val

def print_4rdm_compressed(rdm,tag1):
    fmt = '% -020.14E'
    #fmt = '%e'

    # fix correct naming with tag1
    f=open('fourparticle.rdm.%s' % (tag1),'w')

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[0]+1
        j = lab[1]+1
        k = lab[2]+1
        l = lab[3]+1
        m = lab[4]+1
        n = lab[5]+1
        o = lab[6]+1
        p = lab[7]+1

        if abs(val) == 0: continue

        dump_element(f,val,i,j,k,l,m,n,o,p)

        valm05 = val*(-0.5)
        valm2  = val*(-2.0)

        # print permutations
        if i==j and k==l: # case 1: i==j k==l
           if m==n and o==p:
               dump_element(f,val,i,j,k,l,o,o,m,m)

               dump_element(f,valm05,i,j,k,l,m,o,o,m)
               dump_element(f,valm05,i,j,k,l,m,o,m,o)
               dump_element(f,valm05,i,j,k,l,o,m,m,o)
               dump_element(f,valm05,i,j,k,l,o,m,o,m)
           elif m==n and o!=p:
               # regular permutations
               dump_element(f,val,i,j,k,l,m,m,p,o)
               dump_element(f,val,i,j,k,l,o,p,m,m)
               dump_element(f,val,i,j,k,l,p,o,m,m)
               # remaining permutations
               dump_element(f,valm05,i,j,k,l,o,m,m,p)
               dump_element(f,valm05,i,j,k,l,o,m,p,m)
               dump_element(f,valm05,i,j,k,l,p,m,m,o)
               dump_element(f,valm05,i,j,k,l,p,m,o,m)
               dump_element(f,valm05,i,j,k,l,m,o,p,m)
               dump_element(f,valm05,i,j,k,l,m,p,o,m)
               dump_element(f,valm05,i,j,k,l,m,p,m,o)
               dump_element(f,valm05,i,j,k,l,m,o,m,p)
           elif m!=n and o==p:
               # regular permutations
               dump_element(f,val,i,j,k,l,n,m,o,o)
               dump_element(f,val,i,j,k,l,o,o,m,n)
               dump_element(f,val,i,j,k,l,o,o,n,m)
               # remaining permutations
               dump_element(f,valm05,i,j,k,l,m,o,o,n)
               dump_element(f,valm05,i,j,k,l,n,o,o,m)
               dump_element(f,valm05,i,j,k,l,o,m,o,n)
               dump_element(f,valm05,i,j,k,l,o,n,o,m)
               dump_element(f,valm05,i,j,k,l,o,m,n,o)
               dump_element(f,valm05,i,j,k,l,o,n,m,o)
               dump_element(f,valm05,i,j,k,l,m,o,n,o)
               dump_element(f,valm05,i,j,k,l,n,o,m,o)
           elif m!=n and o!=p:
               if n == o:
                    dump_element(f,valm2,i,j,k,l,n,n,m,p)
                    dump_element(f,valm2,i,j,k,l,n,n,p,m)
                    dump_element(f,valm2,i,j,k,l,m,p,n,n)
                    dump_element(f,valm2,i,j,k,l,p,m,n,n)

                    dump_element(f,val,i,j,k,l,p,n,n,m)
                    dump_element(f,val,i,j,k,l,p,n,m,n)
                    dump_element(f,val,i,j,k,l,m,n,p,n)
                    dump_element(f,val,i,j,k,l,n,p,m,n)
                    dump_element(f,val,i,j,k,l,n,m,p,n)
                    dump_element(f,val,i,j,k,l,n,m,n,p)
                    dump_element(f,val,i,j,k,l,n,p,n,m)
               else:
                    if o > p:
                        dump_element(f,val,i,j,k,l,m,n,p,o)
                        dump_element(f,val,i,j,k,l,n,m,p,o)
                        dump_element(f,val,i,j,k,l,n,m,o,p)
                        dump_element(f,val,i,j,k,l,o,p,m,n)
                        dump_element(f,val,i,j,k,l,o,p,n,m)
                        dump_element(f,val,i,j,k,l,p,o,n,m)
                        dump_element(f,val,i,j,k,l,p,o,m,n)

                        dump_element(f,valm2,i,j,k,l,m,p,n,o)
                        dump_element(f,valm2,i,j,k,l,m,p,o,n)
                        dump_element(f,valm2,i,j,k,l,p,m,o,n)
                        dump_element(f,valm2,i,j,k,l,p,m,n,o)
                        dump_element(f,valm2,i,j,k,l,n,o,p,m)
                        dump_element(f,valm2,i,j,k,l,n,o,m,p)
                        dump_element(f,valm2,i,j,k,l,o,n,m,p)
                        dump_element(f,valm2,i,j,k,l,o,n,p,m)
                    elif p > o:
                        dump_element(f,val,i,j,k,l,m,n,p,o)
                        dump_element(f,val,i,j,k,l,n,m,p,o)
                        dump_element(f,val,i,j,k,l,n,m,o,p)
                        dump_element(f,val,i,j,k,l,o,p,m,n)
                        dump_element(f,val,i,j,k,l,o,p,n,m)
                        dump_element(f,val,i,j,k,l,p,o,n,m)
                        dump_element(f,val,i,j,k,l,p,o,m,n)

        elif i==j and k!=l: # case 2: i==j k!=l
           if m==n and o==p: # 2x2 equal
               dump_element(f,val,i,j,k,l,o,o,m,m)

               dump_element(f,valm05,i,j,k,l,m,o,o,m)
               dump_element(f,valm05,i,j,k,l,o,m,m,o)
               dump_element(f,valm05,i,j,k,l,o,m,o,m)
               dump_element(f,valm05,i,j,k,l,m,o,m,o)
           if m==n and o!=p: # 1 equal (first 2)
               dump_element(f,valm05,i,j,k,l,m,o,m,p)
               dump_element(f,valm05,i,j,k,l,m,p,o,m)
               dump_element(f,valm05,i,j,k,l,o,m,m,p)
               dump_element(f,valm05,i,j,k,l,p,m,o,m)
           if m!=n and o==p: # 1 equal (latter 2)
               dump_element(f,val,i,j,k,l,n,m,o,o)
           if m!=n and n!=o and o!=p and m!=p: # all different
               dump_element(f,val,i,j,k,l,n,m,o,p)

        elif i!=j and j==k and k!=l: # case 3: i!=j j==k k!=l
           if m==n and o==p: # 2x2 equal
               dump_element(f,val,i,j,k,l,o,o,m,m)
               dump_element(f,val,i,j,k,l,o,m,o,m)
               dump_element(f,val,i,j,k,l,m,o,m,o)

               dump_element(f,valm2,i,j,k,l,m,o,o,m)
               dump_element(f,valm2,i,j,k,l,o,m,m,o)
           if m!=n and o==p: # 1 equal (latter 2)
               dump_element(f,val,i,j,k,l,o,o,m,n)
               dump_element(f,val,i,j,k,l,o,m,o,n)
               dump_element(f,val,i,j,k,l,m,o,n,o)

               dump_element(f,valm2,i,j,k,l,m,o,o,n)
           if m!=n and o!=p and m==p: # 1 equal (first and last)
               dump_element(f,val,i,j,k,l,m,o,n,m)
           if m!=n and n!=o and o!=p and m!=p: # all different
               dump_element(f,val,i,j,k,l,m,o,n,p)

        elif i!=j and k==l: # case 4: i!=j k==l
           if m==n and o==p: # 2x2 equal
               dump_element(f,val,i,j,k,l,o,o,m,m)

               dump_element(f,valm05,i,j,k,l,m,o,o,m)
               dump_element(f,valm05,i,j,k,l,o,m,m,o)
               dump_element(f,valm05,i,j,k,l,o,m,o,m)
               dump_element(f,valm05,i,j,k,l,m,o,m,o)
           if m==n and o!=p: # 1 equal (first 2)
               dump_element(f,val,i,j,k,l,m,m,p,o)
           if m!=n and o==p: # 1 equal (latter 2)
               dump_element(f,valm05,i,j,k,l,m,o,n,o)
               dump_element(f,valm05,i,j,k,l,m,o,o,n)
               dump_element(f,valm05,i,j,k,l,o,n,o,m)
               dump_element(f,valm05,i,j,k,l,o,n,m,o)
           if m!=n and n!=o and o!=p and m!=p: # all different
               dump_element(f,val,i,j,k,l,m,n,p,o)

        elif i!=j and k!=l and j!=k: # case 5: i!=j j!=k k!=l
           if (m==o and n==p) or (m==p and n==o):
               dump_element(f,val,i,j,k,l,n,m,p,o)
           if m==n and o==p:
               dump_element(f,val,i,j,k,l,o,o,m,m)

    f.close()

def get_indx_contr(i,j,k,l):

    return (((i-1)*(i)*(i+1)*(i+2))/24
               +((j-1)*(j)*(j+1))/6
               +((k-1)*(k))/2
               +  l
               );

def dump_element(f,value,i,j,k,l,m,n,o,p):
    
    fmt  = '% -020.14E'
    ijkl = get_indx_contr(i,j,k,l)
    f.write(str(ijkl)+' '+str(m)+' '+str(n)+' '+str(o)+' '+str(p)+' '+str(fmt%value)+'\n')

if __name__ == '__main__':

    rdm       = load_4rdm(sys.argv[1])

    if sys.argc > 2: 
        print_4rdm_compressed(rdm,sys.argv[2])
    else:
        print_4rdm(rdm)

