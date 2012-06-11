import sys,os
from os.path import expanduser

import pyalps
from pyalps.apptest import createTest
from numpy import *
import argparse
from subprocess import check_call

# Settings
testname    = 'bh'
prog_dmrg   = '/Users/michele/PhD/hp2c/src_stable/build-merge/applications/dmrg/dmrg'
prog_meas   = '/Users/michele/PhD/hp2c/src_stable/build-merge/applications/dmrg/dmrg_meas'

# Input parameters
parms = {
			'nsweeps'					: 5,
			'nmainsweeps'				: 1,
			'ngrowsweeps'				: 1,

			'max_bond_dimension'		: 200,

			'truncation_final'			: 1e-10,

			'resultfile'                : testname+'.out.h5',
            'chkpfile'                  : testname+'.out.ckp.h5',
			
			'symmetry'                  : 'u1',
	    }
model = {
        	'LATTICE'                   : 'open chain lattice',
            'L'                         : 10,

            'MODEL'                     : 'boson Hubbard',
            'Nmax'                      : 2,
            't'                         : 1,
            'U'                         : 8,

            'CONSERVED_QUANTUMNUMBERS'  : 'N',
            'N_total'                   : 5,
            
            'MEASURE_LOCAL[Local density]' : 'n',
	    }


def run_dmrg(prog):
    pyalps.writeParameterFile(testname+'.parms', parms)
    pyalps.writeParameterFile(testname+'.model', model)
    
    if os.path.exists(testname+'.out.h5'):
        os.remove(testname+'.out.h5')
    if os.path.exists(testname+'.out.log'):
        os.remove(testname+'.out.log')
    
    cmd = [expanduser(prog), testname+'.parms', testname+'.model']
    fp = open(testname+'.out.log', 'w')
    return check_call(cmd, stdout=fp, stderr=fp)    
def run_meas(prog):
    pyalps.writeParameterFile(testname+'.parms', parms)
    pyalps.writeParameterFile(testname+'.model', model)
    
    if os.path.exists(testname+'.out.h5'):
        os.remove(testname+'.out.h5')
    if os.path.exists(testname+'.out.meas.log'):
        os.remove(testname+'.out.meas.log')
    
    cmd = [expanduser(prog), testname+'.parms', testname+'.model']
    fp = open(testname+'.out.meas.log', 'w')
    return check_call(cmd, stdout=fp, stderr=fp)    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test script for MAQUIS DMRG')
    parser.add_argument('--application', metavar='A', default='dmrg',
                       help='path to \'dmrg\' application')
    parser.add_argument('--measure', metavar='M', default='dmrg_meas',
                   help='path to \'dmrg_meas\' application')
    parser.add_argument('--create', dest='action', action='store_const',
                      const='create', default='run',
                      help='create test reference results (default: run dmrg)')

    
    args = parser.parse_args()
    
    if args.action == 'create':
        createTest(sys.argv[0], outputs=[testname+'.out.h5'])
    else:
        run_dmrg(prog_dmrg)
        run_meas(prog_meas)
