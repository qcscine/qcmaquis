#!/usr/bin/env alpspython

from pydmrg import apptest
import sys, os

testname       = os.path.splitext( os.path.basename(sys.argv[0]) )[0]
reference_dir  = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'ref/' )
reference_file = os.path.join(reference_dir, testname+'.h5')

class mytest(apptest.DMRGTestBase):
    testname = testname
    reference_file = reference_file
    
    inputs   = {
                'parms': {
                            'nsweeps'                    : 3,
                            'nmainsweeps'                : 1,
                            'ngrowsweeps'                : 1,
            
                            'max_bond_dimension'         : 200,
            
                            'truncation_initial'         : 0.0001,
                            'truncation_final'           : 1e-12,
            
                            'alpha_initial'              : 0.001,
                            'alpha_main'                 : 1e-6,
                            'alpha_final'                : 0,
                            
                            'resultfile'                 : testname+'.out.h5',
                            'chkpfile'                   : testname+'.out.ckp.h5',
                            
                            'optimization'               : 'singlesite',
                            'symmetry'                   : 'u1',
                          },
                'model': {
                            'LATTICE'                   : 'open square lattice',
                            'L'                         : 5,
                            'W'                         : 5,

                            'MODEL'                     : 'spinless fermions',
                            't'                         : 1,

                            'CONSERVED_QUANTUMNUMBERS'  : 'N',
                            'N_total'                   : 10,
                            
                            'MEASURE_LOCAL[Local density]' : 'n',
                          },
                  }
    observables = [
                    apptest.observable_test.reference_file('Energy',        file=reference_file),
                    apptest.observable_test.reference_file('Entropy',       file=reference_file),
                    apptest.observable_test.reference_file('Local density', file=reference_file),

                    apptest.observable_test.reference_file('Energy',         file=os.path.join(reference_dir, 'test_ff.diag.h5'),
                                                            tolerance=1e-4),
                    apptest.observable_test.reference_file('Local density',  file=os.path.join(reference_dir, 'test_ff.diag.h5'),
                                                            tolerance=1e-4),
                  ]


if __name__ == '__main__':
    apptest.main()
