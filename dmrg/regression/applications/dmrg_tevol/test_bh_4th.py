#!/usr/bin/env alpspython

import sys

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
                            'nsweeps_img'                : 20,
                            'nsweeps'                    : 60,
                            
                            'max_bond_dimension'         : 100,
                            
                            'truncation_final'           : 1e-10,
                            
                            'dt'                         : 0.1,
                            'te_type'                    : 'nn',
                            'te_order'                   : 'fourth',
                            'te_optim'                   : 1,
                            
                            'measure_each'               : 20,
                            'chkp_each'                  : 100,
                            
                            'resultfile'                 : testname+'.out.h5',
                            'chkpfile'                   : testname+'.out.ckp.h5',
                            
                            'symmetry'                   : 'u1',
                          },
                'model': {
                            'LATTICE'                   : 'open chain lattice',
                            'L'                         : 20,
                            
                            'MODEL'                     : 'boson Hubbard',
                            'Nmax'                      : 2,
                            't'                         : 1,
                            'U'                         : 8,
                            
                            'CONSERVED_QUANTUMNUMBERS'  : 'N',
                            'N_total'                   : 10,
                            
                            'MEASURE_LOCAL[Local density]' : 'n',
                          },
                  }
    observables = [
                    apptest.observable_test.reference_file('Energy',        file=reference_file),
                    apptest.observable_test.reference_file('Entropy',       file=reference_file),
                    apptest.observable_test.reference_file('Local density', file=reference_file),

                    apptest.observable_test.reference_file('Energy',        file=reference_file,
                                                            load_type = 'iterations'),
                    
                  ]


if __name__ == '__main__':
    apptest.main()
