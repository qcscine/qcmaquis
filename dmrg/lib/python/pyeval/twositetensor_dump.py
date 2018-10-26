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

# This scripts extracts the two-site tensor dump from a QCMaquis HDF5 file and
# saves it into an ASCII file named twositetensordump.X 
# X is supplied via the command line

import sys
import pyalps


def load_meas(inputfile,meas):
    # load data from the HDF5 result file
    m = pyalps.loadEigenstateMeasurements([inputfile], what=meas)[0][0]
    return m


def save_tst(m,stateno):
    with open('twositetensordump.%s' % stateno, 'w') as f:
        for lab, val in zip(m.x, m.y[0]):
            f.write('%020.14E %i %i %i\n' % (val,lab[1],lab[2],lab[3])) # the site information is ignored for now

if __name__ == '__main__':
    if (len(sys.argv) == 3):
        inputfile = sys.argv[1]
        save_tst(load_meas(inputfile,"tstdump"), sys.argv[2])
    else:
        print "Usage: ", sys.argv[0], "h5file state-number"
    

