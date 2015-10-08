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

import numpy as np

from corrutils import pretty_print,assemble_halfcorr

def load_1spd(inputfile):
    """From the diagonal and upper triangle, construct a symmetric matrix
       diag: diagonal
       triang: upper triangle, sequential reversed rows"""

    diagup     =  pyalps.loadEigenstateMeasurements([inputfile], what='Nup')[0][0]
    diagdown   =  pyalps.loadEigenstateMeasurements([inputfile], what='Ndown')[0][0]
    triangup   =  pyalps.loadEigenstateMeasurements([inputfile], what='dm_up')[0][0]
    triangdown =  pyalps.loadEigenstateMeasurements([inputfile], what='dm_down')[0][0]

    # Create the full matrix from the diagonal (nup.y[0]) and upper triangle (dmup)
    dmu = assemble_halfcorr(diagup.y[0], triangup)
    dmd = assemble_halfcorr(diagdown.y[0], triangdown)

    # this is the spin-density matrix
    ds = dmu-dmd

    return ds

def print_1spdm(rdm):
    #fmt = '% -016.10E'
    fmt = '%e'

    L = int(rdm.props["L"])
    mat = np.zeros((L,L))

    for lab, val in zip(rdm.x, rdm.y[0]):
        i = lab[0]
        j = lab[1]

        mat[i,j] = val;
        mat[j,i] = val;

    pretty_print(mat)

if __name__ == '__main__':
    inputfile = sys.argv[1]

    spdm = load_1spdm(inputfile)
    print_1spdm(spdm)
