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
                            'ENABLE_MEASURE[Entropy]'      : 1,
                          },
                  }
    observables = [
                    apptest.observable_test.reference_file('Energy',        file=reference_file),
                    apptest.observable_test.reference_file('Entropy',       file=reference_file),
                    apptest.observable_test.reference_file('Local density', file=reference_file),
                    
                    apptest.observable_test.reference_value('Energy',        value=-6.881090360639349,
                                                            tolerance=0.001),
                    apptest.observable_test.reference_value('Local density', value=[0.39303194698174188,
                                                                                    0.54302419722113426,
                                                                                    0.51400554735345949,
                                                                                    0.52780042423515783,
                                                                                    0.5221378842084331,
                                                                                    0.52213788420366136,
                                                                                    0.52780042424107332,
                                                                                    0.5140055473576155,
                                                                                    0.54302419721735529,
                                                                                    0.39303194698034677],
                                                            tolerance=0.001),
                  ]


if __name__ == '__main__':
    apptest.main()
