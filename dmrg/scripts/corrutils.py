#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#/**
# * @file
# * @copyright This code is licensed under the 3-clause BSD license.
# *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
# *            See LICENSE.txt for details.
# */

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

    if len(triang.x) == 1:
        i = 0
        j = 1
        ret[i,j] = triang.y
        ret[j,i] = triang.y
    else:
        for lab, val in zip(triang.x, triang.y):
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

    ret = np.zeros((L,L),dtype=np.complex_)

    # set the diagonal
    for i in range(L):
        ret[i,i] = diag[i]

    for lab, val in zip(triang.x, triang.y):
        #print(val)
        i = lab[0]
        j = lab[1]
        ret[i,j] =  val[0]+val[1]*1j # complex assignment, converting from numpy array to a complex value
        # dv(j,i) = dv(i,j)*
        ret[j,i] =  np.conjugate(ret[i,j])

    return ret

def assemble_vector(dataset):
    """From the diagonal and upper triangle, construct a symmetric matrix
       diag: diagonal
       triang: upper triangle, sequential reversed rows"""

    L = len(dataset.y)

    ret = np.zeros(L)

    for lab, val in zip(dataset.x, dataset.y):
        i = lab
        ret[i] = val

    return ret

def merge_transpose(diag, obs1, obs2):
    L = len(diag)

    ret = np.zeros((L,L))

    if len(obs1.x) == 1:
        i = 0
        j = 1
        ret[i,j] = obs1.y
        ret[j,i] = obs1.y
    else:
        for lab, val in zip(obs1.x, obs1.y):
            i = lab[0]
            j = lab[1]
            ret[i,j] = val

    if len(obs2.x) == 1:
        i = 0
        j = 1
        ret[i,j] = obs2.y
        ret[j,i] = obs2.y
    else:
        for lab, val in zip(obs2.x, obs2.y):
            i = lab[0]
            j = lab[1]
            ret[j,i] = val

    return ret

def pretty_print(mat):
    for i in range(len(mat)):
        for j in range(len(mat[i])):
            if (abs(mat[i,j]) > 1e-5):
                print("{0: .5f}".format(mat[i,j]),end="")
            else:
                print("  .     ",end="")

        print("")
