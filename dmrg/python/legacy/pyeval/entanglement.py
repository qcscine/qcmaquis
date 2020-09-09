#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2013 Laboratory for Physical Chemistry, ETH Zurich
#*               2012-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
#*               2014-2015 by   Yingjin Ma   <yingjin.ma@phys.ethz.ch>
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

# All the calculating staff in entropy
from   entropy import MaquisMeasurement
from      copy import deepcopy

if __name__ == "__main__":

    filename = sys.argv[1]	
    simulation = MaquisMeasurement(filename)

    norb = len(simulation.s1())

    s1 = simulation.s1()	
    f1=open('Entropy_s1.out','w')
    f1.write("orbital   1'-Entropy"+'\n')
    for p in range(norb):
	if p<9:    
          f1.write("  "+str(p+1)+"       ")
	elif p<99:
          f1.write(" " +str(p+1)+"       ")
	elif p<999:
          f1.write(""  +str(p+1)+"       ")
        f1.write(str(s1[p])+'\n')
 	
    s2 = simulation.s2()	
    f2=open('Entropy_s2.out','w')
    f2.write("orb  orb     2'-Entropy"+'\n')
    for p in range(norb):
        for q in range(norb):
            if p<9:    
               if q<9:    
                  f2.write("  "+str(p+1)+"    "+str(q+1)+"     ")
  	       elif q<99:
                  f2.write("  "+str(p+1)+"   " +str(q+1)+"     ")
	       elif q<999:
                  f2.write("  "+str(p+1)+"  "  +str(q+1)+"     ")
	    elif p<99:
               if q<9:    
                  f2.write(" " +str(p+1)+"    "+str(q+1)+"     ")
  	       elif q<99:
                  f2.write(" " +str(p+1)+"   " +str(q+1)+"     ")
	       elif q<999:
                  f2.write(" " +str(p+1)+"  "  +str(q+1)+"     ")
	    elif p<999:
               if q<9:    
                  f2.write(""  +str(p+1)+"    "+str(q+1)+"     ")
  	       elif q<99:
                  f2.write(""  +str(p+1)+"   " +str(q+1)+"     ")
	       elif q<999:
                  f2.write(""  +str(p+1)+"  "  +str(q+1)+"     ")
            f2.write(str(s2[p][q])+'\n')

    mu = simulation.I()	
    f3=open('Entropy_mu.out','w')
    f3.write("orb  orb     Mutual inf."+'\n')
    for p in range(norb):
        for q in range(norb):
            if p<9:    
               if q<9:    
                  f3.write("  "+str(p+1)+"    "+str(q+1)+"     ")
  	       elif q<99:
                  f3.write("  "+str(p+1)+"   " +str(q+1)+"     ")
	       elif q<999:
                  f3.write("  "+str(p+1)+"  "  +str(q+1)+"     ")
	    elif p<99:
               if q<9:    
                  f3.write(" " +str(p+1)+"    "+str(q+1)+"     ")
  	       elif q<99:
                  f3.write(" " +str(p+1)+"   " +str(q+1)+"     ")
	       elif q<999:
                  f3.write(" " +str(p+1)+"  "  +str(q+1)+"     ")
	    elif p<999:
               if q<9:    
                  f3.write(""  +str(p+1)+"    "+str(q+1)+"     ")
  	       elif q<99:
                  f3.write(""  +str(p+1)+"   " +str(q+1)+"     ")
	       elif q<999:
                  f3.write(""  +str(p+1)+"  "  +str(q+1)+"     ")
            f3.write(str(mu[p][q])+'\n')



