#!/usr/bin/env alpspython

from pydmrg import apptest
import sys, os

testname = os.path.splitext( os.path.basename(sys.argv[0]) )[0]

class mytest(apptest.DMRGTestBase):
    testname = testname
    reference_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'ref/'+testname+'.h5')
    inputs   = {
                'parms': {
                            'nsweeps'                    : 3,
                            'nmainsweeps'                : 1,
                            'ngrowsweeps'                : 1,
            
                            'max_bond_dimension'         : 50,
            
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
                          },
                  }
    observables = [
                    apptest.observable_test.reference_run('Energy'),
                    apptest.observable_test.reference_run('Entropy'),
                  ]


if __name__ == '__main__':
    apptest.main()
