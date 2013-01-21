#!/usr/bin/env alpspython

import sys

from pydmrg import apptest
import sys, os

testname = os.path.splitext( os.path.basename(sys.argv[0]) )[0]
class mytest(apptest.DMRGTestBase):
    testname = testname
    reference_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'ref/'+testname+'.h5')
    inputs   = {
                'parms': {
                            'nsweeps'                    : 5,
                            'nmainsweeps'                : 1,
                            'ngrowsweeps'                : 1,
                            
                            'max_bond_dimension'         : 200,
                            
                            'truncation_final'           : 1e-10,
                            
                            'resultfile'                 : testname+'.out.h5',
                            'chkpfile'                   : testname+'.out.ckp.h5',
                            
                            'optimization'               : 'singlesite',
                            'symmetry'                   : 'u1',
                          },
                'model': {
                            'LATTICE'                   : 'open chain lattice',
                            'L'                         : 10,
                            
                            'MODEL'                     : 'boson Hubbard',
                            'Nmax'                      : 2,
                            't'                         : 1,
                            'U'                         : 8,
                            
                            'CONSERVED_QUANTUMNUMBERS'  : 'N',
                            'N_total'                   : 5,
                            
                            'MEASURE_LOCAL[Local density]' : 'n',
                          },
                  }
    observables = [
                    apptest.observable_test.reference_run('Energy'),
                    apptest.observable_test.reference_run('Entropy'),
                    apptest.observable_test.reference_run('Local density'),
                  ]


if __name__ == '__main__':
    apptest.main()
