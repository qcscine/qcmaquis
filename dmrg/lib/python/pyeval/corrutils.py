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

    ret = np.zeros((L,L))

    # set the diagonal
    for i in range(L):
        ret[i,i] = diag[i]

    for lab, val in zip(triang.x, triang.y[0]):
        i = lab[0]
        j = lab[1]
        ret[i,j] = val
        ret[j,i] = val

    return ret

def assemble_halfcorr_complex(diag, triang):
    """From the diagonal and upper triangle, construct a symmetric complex matrix
       diag: diagonal (real)
       triang: upper triangle, sequential reversed rows, complex"""

    L = len(diag)

    ret_real = np.zeros((L,L))
    ret_cplx = np.zeros((L,L))

    # set the diagonal
    for i in range(L):
        ret_real[i,i] = diag[i]

    for lab, val in zip(triang.x, triang.y[0]):
        #print val
        i = lab[0]
        j = lab[1]
        ret_real[i,j] = val.real
        ret_real[j,i] = val.real
        ret_cplx[i,j] = val.imag
        ret_cplx[j,i] = val.imag

    return (ret_real,ret_cplx)

def assemble_vector(dataset):
    """From the diagonal and upper triangle, construct a symmetric matrix
       diag: diagonal
       triang: upper triangle, sequential reversed rows"""

    L = len(dataset.y[0])

    ret = np.zeros(L)

    for lab, val in zip(dataset.x, dataset.y[0]):
        i = lab
        ret[i] = val

    return ret

def pretty_print(mat):
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            print "{0: .5f}".format(mat[i,j]),

        print ""
