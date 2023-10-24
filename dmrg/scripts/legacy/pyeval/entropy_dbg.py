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

from corrutils import assemble_halfcorr_complex as assy_hc
from corrutils import pretty_print

class MaquisMeasurement:

    def __init__(self, inputfile):
        self.loc_n = pyalps.loadEigenstateMeasurements([inputfile], what='N')[0][0].y[0]
        self.norb = len(self.loc_n)
        DMRG_Parms = pyalps.getParameters([inputfile])
        orbital_order = map(int, DMRG_Parms[0]['orbital_order'].split(','))
        inv_order = []
        for i in range(self.norb):
			inv_order.append(orbital_order.index(i+1))

        self.orb_order = inv_order
        empty_diag = np.zeros(self.norb)

        self.corr_cdag_c   = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='dm')[0][0])
        self.corr_docc     = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='doccdocc')[0][0])

    def two_orb_rdm(self, p, q):
        pq_dm_matrix = np.zeros((4,4),dtype=np.complex_)

        pq_dm_matrix[0,0] = 1 - self.loc_n[p] - self.loc_n[q] + self.corr_docc[p,q]
        pq_dm_matrix[1,1] = self.loc_n[q] - self.corr_docc[p,q]
        pq_dm_matrix[2,2] = self.loc_n[p] - self.corr_docc[p,q]
        pq_dm_matrix[3,3] = self.corr_docc[p,q]
        pq_dm_matrix[1,2] = self.corr_cdag_c[p,q]
        pq_dm_matrix[2,1] = np.conjugate(pq_dm_matrix[1,2])
        
        #print "pqmatrix[3,3]{0},{1}".format(p,q), "and", "pqmatrix[12,12]{0},{1}".format(p,q), pq_dm_matrix[3,3], pq_dm_matrix[12,12]
        return pq_dm_matrix

    def s1(self):
        n = self.loc_n
        
        ret = np.zeros(len(n))
        for i in range(len(n)):
            m11 = 1 - n[i]
            m22 = n[i]

            # make sure that in the unlikely case that one eigenvalue is 0, log(x) gives something meaningful...
	        # meaning nothing since log(1) = 0.
       	    if m11 == 0: m11 = 1
    	    if m22 == 0: m22 = 1

            ret[i] = -sum(map(lambda x: x*np.log(x), [m11, m22]))

        return ret[self.orb_order]


    def s2(self):
        ret = np.zeros((self.norb,self.norb))

        for p in range(self.norb):
            for q in range(p+1, self.norb):
                rdm = self.two_orb_rdm(p,q)
                evals = [ x for x in sl.eigh(rdm)[0] if x > 1e-12]
                ret[p,q] = -sum(map(lambda x: x*np.log(x), evals))
                ret[q,p] = ret[p,q]

                #rdm = self.two_orb_rdm(q,p)
                #evals = [ x for x in sl.eigh(rdm)[0] if x > 1e-12]
                #ret[q,p] = -sum(map(lambda x: x*np.log(x), evals))

        tmp = ret[:,self.orb_order]
        tmp2 = tmp[self.orb_order,:]
        return tmp2

    def I(self):
        ret = np.zeros((self.norb,self.norb))
        s1 = self.s1()
        s2 = self.s2()
        for p in range(self.norb):
            for q in range(p+1, self.norb):
                ret[p,q] = 0.5 * (s1[p] + s1[q]- s2[p,q])
                ret[q,p] = ret[p,q]

        return ret


    def dump_raw(self):
        for k,v in zip(self.__dict__.keys(), self.__dict__.values()):
            try:
                print k
                pretty_print(v)
            except:
                print v

            print ""
                

if __name__ == '__main__':
    inputfile = sys.argv[1]

    guinea_pig = MaquisMeasurement(inputfile)

    print "s1 matrix"
    print guinea_pig.s1()

    print "s2 matrix"
    s2m = guinea_pig.s2()
    pretty_print(s2m)

    print "I (mutual information)"
    pretty_print(guinea_pig.I())

