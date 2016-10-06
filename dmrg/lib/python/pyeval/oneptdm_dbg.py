#!/usr/bin/env python

import sys
import pyalps

import numpy as np

from corrutils import assemble_halfcorr_complex

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

    print_rdm1(sys.argv[1],sys.argv[2])
