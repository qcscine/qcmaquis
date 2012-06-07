import sys,os
from os.path import expanduser

import pyalps
from pyalps.apptest import createTest
from numpy import *
import argparse
from subprocess import check_call

# Settings
testname    = 'fermiladder'
prog_dmrg   = 'dmrg'
prog_meas   = 'dmrg_meas'

# Input parameters
parms = {
			'nsweeps'					: 1,
			'nmainsweeps'				: 1,
			'ngrowsweeps'				: 1,

			'max_bond_dimension'		: 400,

			'truncation_initial'		: 1e-16,
			'truncation_final'			: 1e-16,

			'alpha_initial'				: 0.001,
			'alpha_main'				: 1e-6,
			'alpha_final'				: 0,
			
			'resultfile'                : testname+'.out.h5',
            'chkpfile'                  : testname+'.out.h5',
            'donotsave'                 : 1,
			
			'symmetry'                  : '2u1',
	    }
model = {
        	'LATTICE'                   : 'open ladder',
            'L'                         : 10,

            'MODEL'                     : 'fermion Hubbard',
            't'                         : 1,
            "t'"                        : 1,
            'U'                         : 8,

            'CONSERVED_QUANTUMNUMBERS'  : 'Nup,Ndown',
            'Nup_total'                 : 10,
            'Ndown_total'               : 10,
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
