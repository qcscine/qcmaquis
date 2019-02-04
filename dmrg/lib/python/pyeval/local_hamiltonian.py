#!/usr/bin/env python2
#-*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
#*               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
#*               2014-2015 by Yingjin Ma <yma@ethz.ch>
#*               2018      by Leon Freitag <lefreita@ethz.ch>
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

# This scripts extracts the local Hamiltonian from QCMaquis HDF5 file and
# saves it into an ASCII file named local_hamiltonian.txt

import sys
import pyalps
import numpy as np


def load_meas(inputfile,meas):
    # load data from the HDF5 result file
    m = pyalps.loadEigenstateMeasurements([inputfile], what=meas)[0][0]
    return m


def save_hamiltonian(m,diag,sigmavec=False):
    #fmt = '%14.14e'

    # get indices

    # convert the compound MPS indices (i.e. 0-1-2) to a single running index

    # first, the indices are converted into an intermediate form
    # itm_index = a*(m+1)^2+b*(m+1)+c where m is the max. index present

    order = max([ max([ m.x[i][j] for i in range(len(m.x)) ]) for j in range(3) ])

    def itm_index(index_list):
        return sum([index_list[j]*(order+1)**(2-j) for j in range(3)])

    # create an array of intermediate indices via a list comprehension
    # from the first three indices
    intermediate_index = [ itm_index(i) for i in [ [ m.x[k][j] for j in range(3) ] for k in range(len(m.x)) ] ]

    # now convert intermediate index to a (final) consecutive index with a dictionary

    idxdict = {}

    idx = 0
    for i in intermediate_index:
        if (idxdict.get(i) == None): # since intermediate indices contain duplicates, we consider only one unique index at a time
            idxdict[i] = idx
            idx += 1

    # now idxdict contains the mapping of the intermediate indices to the final ones

    # determine the matrix size from final indices
    n = max(idxdict.values())+1

    # create a numpy array

    mat = np.zeros((n,n)) if (not diag) else np.zeros(n)

    thresh = 1.0e-19

    # and fill it with the values
    # the zip expression returns a tuple, with the first element as the list of two original index lists and the second element as the values
    if (sigmavec):
        for lab, val in zip([ [ m.x[k][j] for j in range(3) ] for k in range(len(m.x)) ], m.y[0]):
                mat[idxdict[itm_index(lab)]] = val if abs(val) > thresh else 0.0
    else:
        if (not diag):
            for lab, val in zip([ [ [ m.x[k][j] for j in range(3) ], [ m.x[k][j] for j in range(3,6) ] ] for k in range(len(m.x)) ], m.y[0]):
                mat[idxdict[itm_index(lab[0])],idxdict[itm_index(lab[1])]] = val if abs(val) > thresh else 0.0
        else:
            for lab, val in zip([ [ [ m.x[k][j] for j in range(3) ], [ m.x[k][j] for j in range(3,6) ] ] for k in range(len(m.x)) ], m.y[0]):
                mat[idxdict[itm_index(lab[0])]] = val if abs(val) > thresh else 0.0

    # save the array to file

    filename = "local_hamiltonian.txt"
    np.savetxt(filename, mat,fmt='%.19e')


if __name__ == '__main__':
    # diag indicates if we measured only the diagonal
    diag = False
    sigmavec = False
    if (len(sys.argv) > 2):
        if (sys.argv[2] == "-d"): # local Hamiltonian diagonal
            diag = True
        if (sys.argv[2] == "-s"): # sigma vector
            diag = True
            sigmavec = True
    inputfile = sys.argv[1]
    if diag:
      if sigmavec:
        save_hamiltonian(load_meas(inputfile,"sigma_vector"),diag,sigmavec=True)
      else:
        save_hamiltonian(load_meas(inputfile,"local_hamiltonian_diag"),diag)
    else:
        save_hamiltonian(load_meas(inputfile,"local_hamiltonian"),diag)


