#!/usr/bin/env python

import sys
import pyalps

import numpy as np
import scipy.linalg as sl

from corrutils import assemble_halfcorr
from corrutils import assemble_halfcorr_complex
from corrutils import pretty_print

def print_rdm1(inputfile,tag):

    tag1 = tag
    tag2 = tag

    f=open('oneparticle.rdm.%s.%s' % (tag1,tag2),'w')

    # load data from the HDF5 result file
    n  = pyalps.loadEigenstateMeasurements([inputfile], what='N')[0][0]
    dm = pyalps.loadEigenstateMeasurements([inputfile], what='dm')[0][0]

    # Create the full matrix from the diagonal (nup.y[0]) and upper triangle (dmup)
    (dm_real, dm_imag) = assemble_halfcorr_complex(n.y[0], dm)

    spinors = len(dm_real)
    for j in range(spinors):
        for i in range (spinors):
            dump_element(f,dm_real[i,j],dm_imag[i,j],i,j)

    f.close()

def dump_element(f,val_real,val_imag,i,j):
    
    print i, j, "\t", (val_real, val_imag)
    fmt  = '% -020.14E'
    f.write(str(i)+' '+str(j)+'  '+str(fmt%val_real)+'  '+str(fmt%val_imag)+'\n')

if __name__ == '__main__':

    inputfile = sys.argv[1]
    tag       = sys.argv[2]

    print_rdm1(inputfile,tag)
