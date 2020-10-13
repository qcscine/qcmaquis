#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#*****************************************************************************
#*
#* ALPS MPS DMRG Project
#*
#* Copyright (C) 2013 Laboratory for Physical Chemistry, ETH Zurich
#*               2014-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
#*               2015-2015 by Christopher Stein <steinc@phys.chem.ethz.ch>
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


#test implementation of the fiedler vector ordering


import sys
import numpy as np
#from pylab import *
from scipy import *
from numpy import *

import input as DmrgInput
import entropy

#calculate graph laplacian
def get_laplacian(mat_I):
    L = np.zeros(mat_I.shape)
    i = 0
    while i<L.shape[0]:
        j = 0
        while j<L.shape[0]:
            L[i,i] += mat_I[i,j]
            j += 1
        i += 1
    L = L-mat_I
    return L

#extract Fiedler vector
def fiedler(L):
    eigenvalues, eigenvectors = np.linalg.eig(L)
    idx = eigenvalues.argsort()
    sort_ev = eigenvalues[idx]
    sort_vec = eigenvectors[:,idx]
    fiedler_vec = sort_vec[:,1]
    return fiedler_vec


#reorder mutinf according to new ordering
def reorder(s1, mat_I, ordering):
    s1 = s1[ordering]
    new_mutinf = np.zeros(mat_I.shape)
    i = j = 0
    while i < ordering.shape[0]:
        j = 0
        while j < ordering.shape[0]:
            new_mutinf[i,j] = mat_I[ordering[i],ordering[j]]
            j += 1
        i += 1
    return s1, new_mutinf


#calculate cost measure
def cost_meas(mutinf, row_begin = 0, row_end = None, col_begin = 0, col_end = None):
    if row_end is None: row_end = mutinf.shape[0]
    if col_end is None: col_end = mutinf.shape[0]
    eta = 2.
    cost = 0.
    i = row_begin
    while i < row_end:
        j = col_begin
        while j < col_end:
            cost = cost + mutinf[i,j]*abs(i-j)**eta
            j += 1
        i += 1
    return cost


def condense(mutinf,nsym,occ_vec):
    "sum blocks of mutual information matrix and return condensed mutual information matrix of dimension nsym x nsym"
    cond_mutinf = np.zeros((nsym,nsym))
    for i in range(nsym):
        min_row = occ_vec[i]
        max_row = occ_vec[i+1]
        for j in range(nsym):
            sum_elements = 0.0
            min_col = occ_vec[j]
            max_col = occ_vec[j+1]
            k = min_row
            while k < max_row:
                l = min_col
                while l < max_col:
                    cond_mutinf[i,j] = cond_mutinf[i,j]+mutinf[k,l]
                    l = l+1
                k = k+1

    return cond_mutinf

def block_order(dim,occ,order):
    occ_new = np.zeros((occ.shape[0]-1)*2)
    orbvec = np.zeros(dim,dtype=int)
    for i in range((occ.shape[0]-1)):
        occ_new[2*i] =  occ[order[i]]
        occ_new[2*i+1] = occ[order[i]+1]
    i = 0
    for n in range(0,occ_new.shape[0],2):
        j = occ_new[n]
        while j < occ_new[n+1]:
            orbvec[i] = j
            j += 1
            i += 1
    return orbvec, occ_new

def bfo_gfv(occ_new, global_fo):
    """block-fiedler ordering graph-fiedler ordering"""
    orbvec = np.zeros(global_fo.shape[0],dtype=int)
    i = 0
    for j in range(0,occ_new.shape[0],2):
        k = 0
        for k in range(global_fo.shape[0]):
            if global_fo[k] > (occ_new[j]-1) and global_fo[k] < occ_new[j+1]:
                orbvec[i] = global_fo[k]
                i += 1
    return orbvec

def useful_printout_order(order):
    order_string = "orbital_order = \""
    for el in order:
        order_string = order_string + str(el) + ","
    order_string = order_string[:-1]
    order_string = order_string + "\""
    return order_string

def plot_all(mat_I,vec_s1,t1,t2,t3,original_order,fiedler_order,blocked_f_order):
    import matplotlib.pyplot as plt
    import mutinf
    mutinf.plot_mutinf(mat_I, vec_s1, original_order,  title=t1)
    mutinf.plot_mutinf(mat_I, vec_s1, fiedler_order,   title=t2)
    mutinf.plot_mutinf(mat_I, vec_s1, blocked_f_order, title=t3)
    #plotting
    plt.show()

#main program
if __name__ == '__main__':
    inputfile = sys.argv[1]

    entropies = entropy.MaquisMeasurement(inputfile)
    props = DmrgInput.loadProperties(inputfile)


    #control output
    full_output = False
    if (len(sys.argv) > 2 and sys.argv[2] == '-f'):
        full_output = True


    #mutual information matrix:
    mat_I = entropies.I()
    #single entropy vector:
    vec_s1 = entropies.s1()

    original_order = map(int, props["orbital_order"].split(','))

    #plot mutual information without ordering
    t1 = 'Mutual information plot for previous ordering'

    #primitive Fiedler ordering
    L = get_laplacian(mat_I)
    fiedler_vec = fiedler(L)

    #order
    order = fiedler_vec.argsort()
    ofv = fiedler_vec[order]

    fiedler_s1, fiedler_mutinf = reorder(vec_s1, mat_I, order)
    standard_cost = cost_meas(mat_I)
    fiedler_cost = cost_meas(fiedler_mutinf)

    fiedler_order = [original_order[order[i]] for i in range(len(order)) ]
    for el in fiedler_order:
       el += 1


    #plot mutual information with primitive Fiedler ordering
    t2 = 'Mutual information plot for standard Fiedler ordering'


    #Fiedler ordering with symmetry blocks: variant 1 -> in-block ordering according to global Fiedler vector
    #Leon: added the function to filter out empty elements in the list, otherwise it would crash
    site_types = [int(x) for x in filter(lambda l: l != '',props["site_types"].split(','))]
    site_types = [site_types[original_order[i]-1] for i in range(len(site_types))]

    blocked_cost = fiedler_cost
    blocked_order = fiedler_order
    blocked_mutinf = fiedler_mutinf
    blocked_s1 = fiedler_s1

    if max(site_types) > 1:
        occ = np.array([0])
        for i in range(len(site_types)-1):
            if site_types[i] != site_types[i+1]:
                occ = np.append(occ, np.array([i+1]))
        occ = np.append(occ, np.array(len(site_types)))


        #condense mutual information for each symmtery block
        cond_mutinf = condense(mat_I, len(occ) - 1, occ)
        L_cond = get_laplacian(cond_mutinf)
        cond_fv = fiedler(L_cond)
        order1 = cond_fv.argsort()
        #reorder the condensed mutual information
        c_new_s1, c_new_mutinf = reorder(cond_mutinf[:,1], cond_mutinf, order1)

        #reorder blocks in full mutual information
        block_reord, occ_new = block_order(mat_I.shape[0], occ, order1)
        block_s1, block_mutinf = reorder(vec_s1, mat_I, block_reord)

        #reorder according to global Fiedler vector
        #this is stupid -> better order within each block!
        blocked_f_order = bfo_gfv(occ_new, order)
        blocked_s1, blocked_mutinf = reorder(vec_s1, mat_I, blocked_f_order)

	#print '\nnumber of irreducible representations in symmetry-respecting Fiedler ordering:'
        #print [site_types[blocked_f_order[i]] for i in range(len(site_types))]
        blocked_f_order += 1
        t3 = 'Mutual information plot for symmetry-respecting Fiedler ordering'

        blocked_cost = cost_meas(blocked_mutinf)

    if full_output:
        plot_all(mat_I,vec_s1,t1,t2,t3,original_order,fiedler_order,blocked_f_order)
        print "\ncost measures:"
        print "standard ordering: ", standard_cost
        print "Fiedler ordering: ", fiedler_cost
        print "Block Fiedler ordering: ", blocked_cost

        print "\nFiedler ordering: "
        print fiedler_order

        print "\nBlock-Fiedler ordering: "
        print blocked_f_order

        print "\nrecommended ordering: "


    if max(site_types) > 1:
        print useful_printout_order(blocked_f_order)
    else:
        print useful_printout_order(fiedler_order)
