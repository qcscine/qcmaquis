#!/usr/bin/env python
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2013 Laboratory for Physical Chemistry, ETH Zurich
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

import numpy as np

def assemble_halfcorr(diag, triang):
    """From the diagonal and upper triangle, construct a symmetric matrix
       diag: diagonal
       triang: upper triangle, sequential reversed rows"""

    L = len(diag)
    assert((L-1)*L/2 == len(triang.y[0]))

    ret = np.zeros((L,L))

    # set the diagonal
    for i in range(L):
        ret[i,i] = diag[i]

    for lab, val in zip(triang.x, triang.y[0]):
        i = lab[0]
        j = lab[1]
        ret[i,j] = val
        ret[j,i] = val

    # set upper triangle
    #offset = 0
    #for i in range(L):
    #    for j,s in zip(range(i+1,L), reversed(range(offset, offset + L-1-i))):
    #        ret[i,j] = triang[s]
    #        ret[j,i] = triang[s]

    #    offset += L-1-i

    return ret

def pretty_print(mat):
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            print "{0: .5f}".format(mat[i,j]),

        print ""
