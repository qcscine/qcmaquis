#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2013 Laboratory for Physical Chemistry, ETH Zurich
#*               2012-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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
import scipy.linalg as sl

from copy import deepcopy

from corrutils import assemble_halfcorr as assy_hc
from corrutils import assemble_vector as assy_vec
from corrutils import pretty_print

import input as DmrgInput

#from corrutils import assemble_halfcorr_complex

def read_irreps(orbfile):
    of = open(orbfile, 'r')
    orbstring = of.readlines(1000)[1]

    substring = orbstring.split('=')[1]
    ret = map(int, [ s for s in substring.split(',') if s.strip()])
    return ret

def diag_dm(matrix):

    evals = sl.eigh(matrix)[0][::-1]
    for e in evals:
            print "{0: .5f}".format(e)

class oneptdm:
    def __init__(self, inputfile):
        # load data from the HDF5 result file
        self.nup    = assy_vec(pyalps.loadEigenstateMeasurements([inputfile], what='Nup')[0][0])
        self.ndown  = assy_vec(pyalps.loadEigenstateMeasurements([inputfile], what='Ndown')[0][0])
        self.dmup   = assy_hc(self.nup, pyalps.loadEigenstateMeasurements([inputfile], what='dm_up')[0][0])
        self.dmdown = assy_hc(self.ndown, pyalps.loadEigenstateMeasurements([inputfile], what='dm_down')[0][0])

    def rdm(self):
        return self.rdm_a() + self.rdm_b()

    def rdm_a(self):
        return deepcopy(self.dmup)

    def rdm_b(self):
        return deepcopy(self.dmdown)

def assemble_complex_dm(triang):
    """From the upper triangle, construct a symmetric complex matrix
       triang: upper triangle, sequential reversed rows, complex"""

    L = int(triang.props["L"])

    ret_real = np.zeros((L,L))
    ret_cplx = np.zeros((L,L))

    for lab, val in zip(triang.x, triang.y[0]):
        #print val
        i = lab[0]
        j = lab[1]
        ret_real[i,j] =  val.real
        ret_cplx[i,j] =  val.imag
        # dv(j,i) = dv(i,j)*
        ret_real[j,i] =  val.real
        ret_cplx[j,i] = -val.imag

    return (ret_real,ret_cplx)


def print_rdm1(inputfile,tag):

    tag1 = tag
    tag2 = tag

    f=open('oneparticle.rdm.%s.%s' % (tag1,tag2),'w')

    # load data from the HDF5 result file
    dm = pyalps.loadEigenstateMeasurements([inputfile], what='oneptdm')[0][0]

    # old way
    #n  = pyalps.loadEigenstateMeasurements([inputfile], what='N')[0][0]
    #dm = pyalps.loadEigenstateMeasurements([inputfile], what='dm')[0][0]

    # Create the full matrix from the upper triangle (dm)
    (dm_real, dm_imag) = assemble_complex_dm(dm)
    # old way
    #(dm_real, dm_imag) = assemble_halfcorr_complex(n.y[0], dm)

    spinors = int(dm.props["L"])
    for j in range(spinors):
        for i in range (spinors):
            dump_element(f,dm_real[i,j],dm_imag[i,j],i,j)

    f.close()

def dump_element(f,val_real,val_imag,i,j):
    
    print (val_real, val_imag), "\t", i, j
    fmt  = '% -020.14E'
    f.write(str(i)+' '+str(j)+'  '+str(fmt%val_real)+'  '+str(fmt%val_imag)+'\n')

if __name__ == '__main__':

    inputfile = sys.argv[1]

    props = DmrgInput.loadProperties(inputfile)

    if props["symmetry"] == "u1dg":
        print_rdm1(inputfile,sys.argv[2])
    else:
        dm_ = oneptdm(inputfile)
        dm = dm_.rdm()
        pretty_print(dm)

        blocks = [len(dm)]
        if len(sys.argv) == 3:
            orbfile = sys.argv[2]
            irrlist = read_irreps(orbfile)
            blocks = np.bincount(irrlist)[1:]

        print "Point group blocks", blocks

        bstart = 0
        for l in blocks:
            subblock = dm[bstart:bstart+l, bstart:bstart+l]
            diag_dm(subblock)
            bstart += l
