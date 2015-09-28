#!/usr/bin/env python
#-*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2013 Laboratory for Physical Chemistry, ETH Zurich
#*               2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

class spinMeasurement:

    def __init__(self, inputfile):
    self.loc_nup          = pyalps.loadEigenstateMeasurements([inputfile], what='Nup')[0][0].y[0]
    self.loc_ndown        = pyalps.loadEigenstateMeasurements([inputfile], what='Ndown')[0][0].y[0]
    self.loc_nupdown      = pyalps.loadEigenstateMeasurements([inputfile], what='Nupdown')[0][0].y[0]

    self.loc_nup_nup      = pyalps.loadEigenstateMeasurements([inputfile], what='nupnup')[0][0].y[0]
    self.loc_ndown_nup    = pyalps.loadEigenstateMeasurements([inputfile], what='nupndown')[0][0].y[0]
    self.loc_nup_ndown    = pyalps.loadEigenstateMeasurements([inputfile], what='ndownnup')[0][0].y[0]
    self.loc_ndown_ndown  = pyalps.loadEigenstateMeasurements([inputfile], what='ndownndown')[0][0].y[0]
    self.loc_splus_sminus = pyalps.loadEigenstateMeasurements([inputfile], what='splus_sminus')[0][0].y[0]

    def sminusplus(self):
        nupdown   = self.loc_nupdown
        ndown     = self.loc_ndown
        splusminus = self.loc_splus_sminus

        return (sum(ndown) - sum(nupdown) + 2*sum(splus_sminus))

    def sz2(self):
        nup         = self.loc_nup
        ndown       = self.loc_ndown
        nupdown     = self.loc_nupdown
        nup_nup     = self.loc_nup_nup
        ndown_ndown = self.loc_ndown_ndown
        nup_ndown   = self.loc_nup_ndown
        ndown_nup   = self.loc_ndown_nup

        return ( 0.25*sum(nup) - 0.5*sum(nupdown)  + 0.25*sum(ndown) \
          + 2*0.25*(sum(nup_nup) - sum(nup_ndown) - sum(ndown_nup) + sum(ndown_ndown)))

    def sz(self):
        nup   = self.loc_nup
        ndown = self.loc_ndown

        return (0.5*(sum(nup) - sum(ndown)))

    S2 = sminusplus + s_z2 + s_z


if __name__ == '__main__':
    inputfile = sys.argv[1]

    guinea_pig = spinMeasurement(inputfile)

    S2 = guinea_pig.sminusplus() + guinea_pig.sz2() + guinea_pig.sz()

    print "  <psi|S_z|psi>:", guinea_pig.sz()
    print "<psi|S_z^2|psi>:", guinea_pig.sz2()
    print "  <psi|S^2|psi>:", S2
    print "       S       :", (-1 + (1+4*S2)**0.5)/2

    #print "sum(splus_sminus)", sum(splus_sminus)
    #print "sum(ndown) - sum(nupdown)", sum(ndown) - sum(nupdown)
    #nup_c_nup = sum(np.array([[x*y for y in nup] for x in nup]).flatten())
