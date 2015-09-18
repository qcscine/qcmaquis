#!/usr/bin/env python
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
import s2
import numpy as np


if __name__ == '__main__':
    inputfile = sys.argv[1]

    guinea_pig = s2.MaquisMeasurement(inputfile)


import matplotlib.pyplot as plt
from matplotlib import lines
from pylab import *

#define functions

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
    print len(L)
    print len(L[0])
    eigenvalues, eigenvectors = np.linalg.eig(L)
    idx = eigenvalues.argsort()
    sort_ev = eigenvalues[idx]
    sort_vec = eigenvectors[:,idx]
    fiedler_vec = sort_vec[:,1]
    return fiedler_vec


#reorder mutinf according to new ordering
def reorder(s1, mat_I, order):
    s1 = s1[order]
    new_mutinf = np.zeros(mat_I.shape)
    i = j = 0
    while i < order.shape[0]:
        j = 0
        while j < order.shape[0]:
            new_mutinf[i,j] = mat_I[order[i],order[j]]
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
            cost = cost + mutinf[i,j]*abs(i-j)**2
            j += 1
        i += 1
    return cost


def condense(mutinf,nsym,occ_vec):
    cond_mutinf = np.zeros((nsym,nsym))
    for i in range(nsym):
        for j in range(nsym):
            cond_mutinf[i,j] = cost_meas(mutinf,occ_vec[i],occ_vec[i+1],occ_vec[j],occ_vec[j+1])
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

#plot of new mutual information
def plot_mutinf(new_mutinf,s1,ord,nr,title,cost):
    plt.figure(nr)
    N = len(new_mutinf)
    theta = np.zeros(N)
    r=np.zeros(N)
    labels=np.zeros(N)
    area=np.zeros(N)
    for i in range(N):
        theta[i]=-2*np.pi/N*i+np.pi/2
        r[i]=1.0
        labels[i]=ord[i]+1
        area[i]=s1[i]*250

    ax = plt.subplot(111,polar=True)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(b=False)
    c = plt.scatter(theta,r,c="Red",s=area)

    plt.title(title+'\nResults file: '+sys.argv[1]+'\n entanglement distance: '+str(cost))

#this is dummy:
    c1 = plt.scatter(theta,(r+0.1),c="red",s=0)

    legendlines = {}

    for i in range(N):
#  plt.annotate(int(labels[i]),xy=(theta[i],(r[i]+0.2)),size='xx-large',)
        plt.text(theta[i],(r[i]+0.18),int(labels[i]),size='medium',ha='center',va='center')
  
    
        for j in range(i,N):
             x =[theta[i],theta[j]]
             y =[1,1]
             if new_mutinf[i,j] >= 0.1:
                line = lines.Line2D(x, y, linewidth=2*5*new_mutinf[i,j], color='black',linestyle='-', alpha=1,label='0.1')
                legendlines['0.1'] = line
                ax.add_line(line)
             elif new_mutinf[i,j] >=0.01:
                 line = lines.Line2D(x, y, linewidth=2*15*new_mutinf[i,j], color='gray',linestyle='--', alpha=1,label='0.01')
                 legendlines['0.01'] = line
                 ax.add_line(line)
             elif new_mutinf[i,j] >=0.001:
                 line = lines.Line2D(x, y, linewidth=1.5, color='lime',linestyle=':', alpha=1,label='0.001')
                 legendlines['0.001'] = line
                 ax.add_line(line)


    ax.legend(legendlines.values(),[l.get_label() for l in legendlines.values()],bbox_to_anchor=(0.00,1.0),fancybox=True,shadow=True)
    #plt.savefig(str(nr)+'.pdf')

#main program

#Input from result file

#mutual information matrix:
mat_I=guinea_pig.I()
#single entropy vector:
vec_s1=guinea_pig.s1()

#occ = np.array([0])

#primitive Fiedler ordering
L = get_laplacian(mat_I)
fiedler_vec = fiedler(L)

#order
order = fiedler_vec.argsort()
ofv = fiedler_vec[order]


new_s1, new_mutinf = reorder(vec_s1,mat_I,order)
cost_old = cost_meas(mat_I)
cost_new = cost_meas(new_mutinf)


#Fiedler ordering with symmetry blocks: variant 1 -> in-block ordering according to global Fiedler vector

#var1_mutinf = mat_I
#final_order1 = range(order.shape[0])
#final_mutinf1 = mat_I
#final_s1_1 = vec_s1
#final_cost1 = cost_old
#j = 0
#while j < 5:
#    cond_mutinf = condense(var1_mutinf,occ.shape[0]-1,occ)
#    L_cond = get_laplacian(cond_mutinf)
#    cond_fv = fiedler(L_cond)
#    order1 = cond_fv.argsort()
#    c_new_s1, c_new_mutinf = reorder(cond_mutinf[:,1],cond_mutinf,order1)
#    block_reord, occ_new = block_order(mat_I.shape[0],occ,order1)
#    block_s1, block_mutinf = reorder(vec_s1,mat_I,block_reord)
#    blocked_f_order = bfo_gfv(occ_new, order)
#    var1_s1, var1_mutinf = reorder(vec_s1,mat_I,blocked_f_order)
#    
#    if cost_meas(var1_mutinf) < final_cost1:
#        final_cost1 = cost_meas(var1_mutinf)
#        final_mutinf1 = var1_mutinf
#        final_order1 = blocked_f_order
#        final_s1_1 = var1_s1
#    else: j += 1
    
    
    

#print cost_old, cost_new, final_cost1
print cost_old, cost_new
print order+1

#print final_order1+1




#plotting
t1 = 'mutual information plot from HF ordering'
plot_mutinf(mat_I,vec_s1,range(order.shape[0]),1,t1,cost_old)
t2 = 'mutual information plot from primitive Fiedler ordering'
plot_mutinf(new_mutinf,new_s1,order,2,t2,cost_new)
#t3 = 'mutual information plot from Fiedler block ordering'
#plot_mutinf(final_mutinf1,final_s1_1,final_order1,3,t3,final_cost1)
plt.show()




