#!/usr/bin/env python

import sys
import pyalps

import numpy as np
import scipy.linalg as sl

from corrutils import assemble_halfcorr
from corrutils import assemble_halfcorr_complex
from corrutils import pretty_print

def print_rdm1(inputfile):

    # load data from the HDF5 result file
    n  = pyalps.loadEigenstateMeasurements([inputfile], what='N')[0][0]
    dm = pyalps.loadEigenstateMeasurements([inputfile], what='dm')[0][0]

    # Create the full matrix from the diagonal (nup.y[0]) and upper triangle (dmup)
    (dm_real, dm_imag) = assemble_halfcorr_complex(n.y[0], dm)

    spinors = len(dm_real)
    for j in range(spinors):
        for i in range (spinors):
            print i+1, j+1, "\t", (dm_real[i,j], dm_imag[i,j])

    return spinors

if __name__ == '__main__':

    inputfile = sys.argv[1]
    spinors   = print_rdm1(inputfile)

