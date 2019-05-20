from numpy import *
import pyalps
from pyalps.ngs import archive

infile  = 'ref/nevol_pure_L3_J1.0_U1.0_dt0.02.dat'
outfile = 'ref/pure_bosons.h5'

observable = 'Local density'
meas_at = [1., 2., 3., 4., 5.]
parms = {
            't'     : 1.,
            'U'     : 1.,
        }

def create_archive(infile, outfile, observable, meas_at, parms):
    data = loadtxt(infile)
    sel = array( [where( abs(data[:,0]-ti)<1e-6 )[0][0] for ti in meas_at] )
    measurements = data[sel, 1:]

    ar = archive(outfile, 'w')
    ar['/parameters'] = parms
    for i, m in enumerate(measurements):
        path = '/simulation/sweep%s/results/%s/mean/value' % (i, observable)
        ar[path] = [m]


create_archive(   infile     = 'ref/nevol_pure_L3_J1.0_U1.0_dt0.02.dat'
                , outfile    = 'ref/pure_bosons.diag.h5'
                , observable = 'Local density'
                , meas_at    = [1., 2., 3., 4., 5.]
                , parms      = {
                                't'     : 1.,
                                'U'     : 1.,
                                }
              )

create_archive(   infile     = 'ref/nevol_pure_L3_J1.0_U1.0_Gamma1a0.2_dt0.02.dat'
                , outfile    = 'ref/pure_bosons_diss.diag.h5'
                , observable = 'Local density'
                , meas_at    = [1., 2., 3., 4., 5.]
                , parms      = {
                                't'         : 1.,
                                'U'         : 1.,
                                'Gamma1a'   : 0.2,
                                }
              )

create_archive(   infile     = 'ref/nevol_coherent_L3_J1.0_U1.0_Gamma1a0.2_dt0.02.dat'
                , outfile    = 'ref/coherent_bosons_diss.diag.h5'
                , observable = 'Local density'
                , meas_at    = [1., 2., 3., 4., 5.]
                , parms      = {
                                't'         : 1.,
                                'U'         : 1.,
                                'Gamma1a'   : 0.2,
                                }
              )
