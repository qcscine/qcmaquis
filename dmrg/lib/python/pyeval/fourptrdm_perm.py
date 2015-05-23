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

if __name__ == '__main__':
    inputfile = sys.argv[1]

    rdm = load_4rdm(inputfile)
    print_4rdm(rdm)
