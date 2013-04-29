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
            
            'measure_each'               : 10,
            'chkp_each'                  : 50,
            
            'resultfile'                 : testname+'.out.h5',
            'chkpfile'                   : testname+'.out.ckp.h5',
            
            'init_state'                 : 'basis_state',
            'init_basis_state'           : '0,1,0',
            
            'always_measure'             : 'Local density',
            
            'symmetry'                   : 'none',
            'model_library'              : 'coded',
        }

model = {
            'LATTICE'                   : 'open chain lattice',
            'L'                         : 3,
            
            'MODEL'                     : 'boson Hubbard',
            'Nmax'                      : 3,
            't'                         : 1,
            'U'                         : 1,
            
            'ENABLE_MEASURE[Local density]' : 1,
        }

class mytest_2nd(apptest.DMRGTestBase):
    testname = testname + '_2nd'
    reference_file = os.path.join(reference_dir, testname+'.h5')
    
    inputs   = {
                'parms': dict( parms.items() + {'te_order': 'second'}.items() ),
                'model': model,
                  }
    observables = [
                    apptest.observable_test.reference_file('Energy',        file=reference_file,
                                                            load_type = 'iterations'),
                    apptest.observable_test.reference_file('Local density', file=reference_file,
                                                            load_type = 'iterations'),
                                                            
                    apptest.observable_test.reference_file('Local density', load_type = 'iterations', tolerance=0.01,
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
                    apptest.observable_test.reference_file('Energy',        file=reference_file,
                                                            load_type = 'iterations'),
                    apptest.observable_test.reference_file('Local density', file=reference_file,
                                                            load_type = 'iterations'),
                                                            
                    apptest.observable_test.reference_file('Local density', load_type = 'iterations', tolerance=0.01,
                                                           file=os.path.join(reference_dir, 'pure_bosons.diag.h5')),
                  ]


if __name__ == '__main__':
    apptest.main()
