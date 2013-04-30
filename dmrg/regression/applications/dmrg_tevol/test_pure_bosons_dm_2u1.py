#!/usr/bin/env alpspython

import sys

from pydmrg import apptest
import sys, os

testname       = os.path.splitext( os.path.basename(sys.argv[0]) )[0]
reference_dir  = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'ref/' )
reference_file = os.path.join(reference_dir, testname+'.h5')

parms = {
            'nsweeps_img'                : 0,
            'nsweeps'                    : 50,
            
            'max_bond_dimension'         : 100,
            
            'truncation_final'           : 1e-10,
            
            'dt'                         : 0.1,
            'te_type'                    : 'nn',
            'te_optim'                   : 1,
            'expm_method'                : 'heev',
            
            'measure_each'               : 10,
            'chkp_each'                  : 50,
            
            'resultfile'                 : testname+'.out.h5',
            'chkpfile'                   : testname+'.out.ckp.h5',
            
            'init_state'                 : 'basis_state_generic',
            'init_basis_state'           : '3,6,3',
            
            'always_measure'             : 'Local density,Density',
            
            'symmetry'                   : '2u1',
            'lattice_library'            : 'continuum',
            'model_library'              : 'continuum',
        }

model = {
            'LATTICE'                   : 'continuous_chain',
            'L'                         : 3,
            'Ndiscr'                    : 1,
            
            'MODEL'                     : 'optical_lattice_cons_dm',
            'Nmax'                      : 3,
            'h'                         : 1,
            'c'                         : 1,
            'u1_total_charge1'          : 1,
            'u1_total_charge2'          : -1,
            
            'MEASURE_CONTINUUM[Local density]' : 1,
        }

class mytest_2nd(apptest.DMRGTestBase):
    testname = testname + '_2nd'
    reference_file = os.path.join(reference_dir, testname+'.h5')
    
    inputs   = {
                'parms': dict( parms.items() + {'te_order': 'second'}.items() ),
                'model': model,
                  }
    observables = [
                    # apptest.observable_test.reference_file('Local density', file=reference_file,
                    #                                         load_type = 'iterations'),
                                                            
                    apptest.observable_test.reference_file('Local density', load_type = 'iterations', tolerance=0.05,
                                                           file=os.path.join(reference_dir, 'pure_bosons.diag.h5')),
                  ]

class mytest_4th(apptest.DMRGTestBase):
    testname = testname + '_4th'
    reference_file = os.path.join(reference_dir, testname+'.h5')
    
    inputs   = {
                'parms': dict( parms.items() + {'te_order': 'fourth'}.items() ),
                'model': model,
                  }
    observables = [
                    # apptest.observable_test.reference_file('Local density', file=reference_file,
                    #                                         load_type = 'iterations'),
                                                            
                    apptest.observable_test.reference_file('Local density', load_type = 'iterations', tolerance=0.05,
                                                           file=os.path.join(reference_dir, 'pure_bosons.diag.h5')),
                  ]


if __name__ == '__main__':
    apptest.main()
