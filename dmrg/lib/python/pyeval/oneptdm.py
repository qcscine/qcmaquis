#!/usr/bin/env python

import sys
import pyalps

import numpy as np
import scipy.linalg as sl

from corrutils import assemble_halfcorr
from corrutils import pretty_print

def read_irreps(orbfile):
    of = open(orbfile, 'r')
    orbstring = of.readlines(1000)[1]

    substring = orbstring.split('=')[1]
    ret = map(int, [ s for s in substring.split(',') if s.strip()])
    return ret

def diag_dm(matrix):

    evals = sl.eigh(matrix)[0][::-1]
    for e in evals:
            print "{0: .5f}".format(e)


if __name__ == '__main__':
    inputfile = sys.argv[1]

    # load data from the HDF5 result file
    nup = pyalps.loadEigenstateMeasurements([inputfile], what='Nup')[0][0]
    ndown = pyalps.loadEigenstateMeasurements([inputfile], what='Ndown')[0][0]
    dmup = pyalps.loadEigenstateMeasurements([inputfile], what='dm_up')[0][0]
    dmdown = pyalps.loadEigenstateMeasurements([inputfile], what='dm_down')[0][0]

    # Create the full matrix from the diagonal (nup.y[0]) and upper triangle (dmup)
    dmu = assemble_halfcorr(nup.y[0], dmup)
    dmd = assemble_halfcorr(ndown.y[0], dmdown)

    # this is the density matrix
    dm = dmu+dmd
    pretty_print(dm)

    blocks = [len(dm)]
    if len(sys.argv) == 3:
        orbfile = sys.argv[2]
        irrlist = read_irreps(orbfile)
        blocks = np.bincount(irrlist)[1:]

    print "Point group blocks", blocks

    bstart = 0
    for l in blocks:
        subblock = dm[bstart:bstart+l, bstart:bstart+l]
        diag_dm(subblock)
        bstart += l
