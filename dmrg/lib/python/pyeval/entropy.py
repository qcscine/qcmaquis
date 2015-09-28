#!/usr/bin/env python
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
from corrutils import pretty_print

class MaquisMeasurement:

    def __init__(self, inputfile):
        self.loc_nup = pyalps.loadEigenstateMeasurements([inputfile], what='Nup')[0][0].y[0]
        self.loc_ndown = pyalps.loadEigenstateMeasurements([inputfile], what='Ndown')[0][0].y[0]
        self.loc_docc = pyalps.loadEigenstateMeasurements([inputfile], what='Nupdown')[0][0].y[0]

        self.norb = len(self.loc_nup)
        empty_diag = np.zeros(self.norb)

        self.corr_cdag_up_c_up = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='dm_up')[0][0])
        self.corr_cdag_down_c_down = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='dm_down')[0][0])

        self.corr_nupnup =     assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='nupnup')[0][0])
        self.corr_nupndown =   assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='nupndown')[0][0])
        self.corr_ndownnup =   assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='ndownnup')[0][0])
        self.corr_ndownndown = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='ndownndown')[0][0])

        self.corr_docc =       assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='doccdocc')[0][0])
        self.corr_trans_up =   assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='transfer_up_while_down')[0][0])
        self.corr_trans_down = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='transfer_down_while_up')[0][0])
        self.corr_trans_up_down2 = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='transfer_up_while_down_at_2')[0][0])
        self.corr_trans_up_down1 = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='transfer_up_while_down_at_1')[0][0])
        self.corr_trans_down_up2 = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='transfer_down_while_up_at_2')[0][0])
        self.corr_trans_down_up1 = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='transfer_down_while_up_at_1')[0][0])
        self.corr_trans_pair = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='transfer_pair')[0][0])

        self.corr_spinflip = assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='spinflip')[0][0])

        self.corr_nupdocc =    assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='nupdocc')[0][0])
        self.corr_ndowndocc =  assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='ndowndocc')[0][0])
        self.corr_doccnup =    assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='doccnup')[0][0])
        self.corr_doccndown =  assy_hc(empty_diag, pyalps.loadEigenstateMeasurements([inputfile], what='doccndown')[0][0])

    def one_pt_dm(self):
        
        dmu = deepcopy(self.corr_cdag_up_c_up)
        for i in range(self.norb):
            dmu[i,i] = self.loc_nup[i]

        dmd = deepcopy(self.corr_cdag_down_c_down)
        for i in range(self.norb):
            dmd[i,i] = self.loc_ndown[i]
            
        return dmu+dmd

    def two_orb_rdm(self, p, q):
        pq_dm_matrix = np.zeros((16,16))

        pq_dm_matrix[ 0, 0] = 1 + self.loc_docc[p] + self.loc_docc[q] + self.corr_docc[p,q] - self.loc_ndown[p] \
                              - self.corr_ndowndocc[p,q] - self.loc_ndown[q] - self.corr_doccndown[p,q] \
                              + self.corr_ndownndown[p,q] - self.loc_nup[p] - self.corr_nupdocc[p,q] + self.corr_nupndown[p,q] \
                              - self.loc_nup[q] - self.corr_doccnup[p,q] + self.corr_ndownnup[p,q] + self.corr_nupnup[p,q]
        pq_dm_matrix[ 1, 1] = -self.loc_docc[p] - self.corr_docc[p,q] + self.loc_ndown[p] + self.corr_ndowndocc[p,q] \
                            +  self.corr_doccndown[p,q] - self.corr_ndownndown[p,q] + self.corr_doccnup[p,q] - self.corr_ndownnup[p,q]

        pq_dm_matrix[ 2, 2] = -self.loc_docc[p] - self.corr_docc[p,q] + self.corr_doccndown[p,q] + self.loc_nup[p] \
                            +  self.corr_nupdocc[p,q] - self.corr_nupndown[p,q] + self.corr_doccnup[p,q] - self.corr_nupnup[p,q]

        pq_dm_matrix[ 3, 3] =  self.loc_docc[p] - self.corr_doccndown[p,q] - self.corr_doccnup[p,q] + self.corr_docc[p,q]
        pq_dm_matrix[ 4, 4] = -self.loc_docc[q] - self.corr_docc[p,q] + self.corr_ndowndocc[p,q] + self.loc_ndown[q] \
                            +  self.corr_doccndown[p,q] - self.corr_ndownndown[p,q] + self.corr_nupdocc[p,q] - self.corr_nupndown[p,q]

        pq_dm_matrix[ 5, 5] = self.corr_ndownndown[p,q] - self.corr_ndowndocc[p,q] - self.corr_doccndown[p,q] + self.corr_docc[p,q]
        pq_dm_matrix[ 6, 6] = self.corr_nupndown[p,q] - self.corr_doccndown[p,q] - self.corr_nupdocc[p,q] + self.corr_docc[p,q]
        pq_dm_matrix[ 7, 7] = self.corr_doccndown[p,q] - self.corr_docc[p,q]
        pq_dm_matrix[ 8, 8] = -self.loc_docc[q] - self.corr_docc[p,q] + self.corr_ndowndocc[p,q] + self.corr_nupdocc[p,q]
\
                             + self.loc_nup[q] + self.corr_doccnup[p,q] - self.corr_ndownnup[p,q] - self.corr_nupnup[p,q]

        pq_dm_matrix[ 9, 9] = self.corr_ndownnup[p,q] - self.corr_ndowndocc[p,q] - self.corr_doccnup[p,q] + self.corr_docc[p,q]
        pq_dm_matrix[10,10] = self.corr_nupnup[p,q] - self.corr_nupdocc[p,q] - self.corr_doccnup[p,q] + self.corr_docc[p,q]
        pq_dm_matrix[11,11] = self.corr_doccnup[p,q] - self.corr_docc[p,q]
        pq_dm_matrix[12,12] = self.loc_docc[q] - self.corr_nupdocc[p,q] - self.corr_ndowndocc[p,q] + self.corr_docc[p,q]
        pq_dm_matrix[13,13] = self.corr_ndowndocc[p,q] - self.corr_docc[p,q]
        pq_dm_matrix[14,14] = self.corr_nupdocc[p,q] - self.corr_docc[p,q]
        pq_dm_matrix[15,15] = self.corr_docc[p,q]

        pq_dm_matrix[1,4] = self.corr_cdag_down_c_down[p,q] - self.corr_trans_down_up1[p,q] - self.corr_trans_down_up2[p,q] \
                          + self.corr_trans_down[p,q]
        pq_dm_matrix[4,1] = pq_dm_matrix[1,4]

        pq_dm_matrix[2,8] = self.corr_cdag_up_c_up[p,q] - self.corr_trans_up_down1[p,q] - self.corr_trans_up_down2[p,q] \
                          + self.corr_trans_up[p,q]
        pq_dm_matrix[8,2] = pq_dm_matrix[2,8]

        pq_dm_matrix[3,6] = self.corr_trans_down_up1[p,q] - self.corr_trans_down[p,q]
        pq_dm_matrix[6,3] = pq_dm_matrix[3,6]

        pq_dm_matrix[3,9] = -self.corr_trans_up_down1[p,q] + self.corr_trans_up[p,q]
        pq_dm_matrix[9,3] = pq_dm_matrix[3,9]

        pq_dm_matrix[6,9] = self.corr_spinflip[p,q]                
        pq_dm_matrix[9,6] = pq_dm_matrix[6,9]

        pq_dm_matrix[3,12] = self.corr_trans_pair[p,q]
        pq_dm_matrix[12,3] = pq_dm_matrix[3,12]

        pq_dm_matrix[6,12] = -self.corr_trans_up_down2[p,q] + self.corr_trans_up[p,q]
        pq_dm_matrix[12,6] = pq_dm_matrix[6,12]

        pq_dm_matrix[9,12] = self.corr_trans_down_up2[p,q] - self.corr_trans_down[p,q]
        pq_dm_matrix[12,9] = pq_dm_matrix[9,12]

        pq_dm_matrix[7,13] = self.corr_trans_up[p,q]
        pq_dm_matrix[13,7] = pq_dm_matrix[7,13]

        pq_dm_matrix[11,14] = self.corr_trans_down[p,q]
        pq_dm_matrix[14,11] = pq_dm_matrix[11,14]
        
        #print "pqmatrix[3,3]{0},{1}".format(p,q), "and", "pqmatrix[12,12]{0},{1}".format(p,q), pq_dm_matrix[3,3], pq_dm_matrix[12,12]
        return pq_dm_matrix

    def s1(self):
        nu = self.loc_nup
        nd = self.loc_ndown
        nud = self.loc_docc

        assert(len(nu)==len(nd)==len(nud))
        ret = np.zeros(len(nu))
        for i in range(len(nu)):
            m11 = nu[i] - nud[i]
            m22 = nd[i] - nud[i]
            m33 = 1 - nu[i] - nd[i] + nud[i]
            m44 = nud[i]

            ret[i] = -sum(map(lambda x: x*np.log(x), [m11, m22, m33, m44]))

        return ret        


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

        return ret

    def I(self):
        ret = np.zeros((self.norb,self.norb))
        s1 = self.s1()
        s2 = self.s2()
        for p in range(self.norb):
            for q in range(p+1, self.norb):
                ret[p,q] = 0.5 * (s1[p] + s1[q]- s2[p,q])
                ret[q,p] = ret[p,q]

        return ret


    def corr_func(self, f1,f2):
        n = self.norb
        cf = np.zeros((n, n))

        for i in range(n):
            for j in range(i+1,n):
                cf[i,j] = self.two_orb_rdm(i,j)[f1,f2]

        return cf

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

    #ozone.dump_raw()
