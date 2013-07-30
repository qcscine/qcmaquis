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
                            'nsweeps'                    : 2,
                            'nmainsweeps'                : 1,
                            'ngrowsweeps'                : 1,
                            
                            'max_bond_dimension'         : 400,
                            
                            'truncation_initial'         : 1e-10,
                            'truncation_final'           : 1e-10,
                            
                            'alpha_initial'              : 0.001,
                            'alpha_main'                 : 1e-6,
                            'alpha_final'                : 0,
                            
                            'resultfile'                 : testname+'.out.h5',
                            'chkpfile'                   : testname+'.out.ckp.h5',
                            
                            'optimization'               : 'twosite',
                            'symmetry'                   : '2u1',
                          },
                'model': {
                            'LATTICE'                   : 'open ladder',
                            'L'                         : 6,

                            'MODEL'                     : 'fermion Hubbard',
                            't'                         : 1,
                            "t'"                        : 1,
                            'U'                         : 8,

                            'CONSERVED_QUANTUMNUMBERS'  : 'Nup,Ndown',
                            'Nup_total'                 : 10,
                            'Ndown_total'               : 10,
                            
                            'MEASURE_LOCAL[Local density up]'                   : 'n_up',
                            'MEASURE_LOCAL[Local density down]'                 : 'n_down',
                            'MEASURE_CORRELATIONS[Onebody density matrix up]'   : 'cdag_up:c_up',
                            'MEASURE_CORRELATIONS[Onebody density matrix down]' : 'cdag_down:c_down',
                            'ENABLE_MEASURE[Entropy]'      : 1,
                          },
                  }
    observables = [
                    apptest.observable_test.reference_file('Energy',                      file=reference_file),
                    apptest.observable_test.reference_file('Entropy',                     file=reference_file),
                    apptest.observable_test.reference_file('Local density up',            file=reference_file),
                    apptest.observable_test.reference_file('Local density down',          file=reference_file),
                    apptest.observable_test.reference_file('Onebody density matrix up',   file=reference_file),
                    apptest.observable_test.reference_file('Onebody density matrix down', file=reference_file),
                    
                    apptest.observable_test.reference_value('Energy',            value=55.08693016528693,
                                                            tolerance=0.001),
                    apptest.observable_test.reference_value('Local density up',  value=[0.86015912231757885,
                                                                                       0.86015912231996861,
                                                                                       0.80381048962759327,
                                                                                       0.80381048962879098,
                                                                                       0.83603038805181895,
                                                                                       0.83603038805119045,
                                                                                       0.8360303880519383,
                                                                                       0.83603038804869678,
                                                                                       0.80381048962858803,
                                                                                       0.80381048962535662,
                                                                                       0.86015912232003822,
                                                                                       0.86015912231926772],
                                                            tolerance=0.001),
                    apptest.observable_test.reference_value('Local density down', value=[0.860159122323,
                                                                                        0.860159122321,
                                                                                        0.803810489626,
                                                                                        0.803810489623,
                                                                                        0.836030388051,
                                                                                        0.836030388051,
                                                                                        0.836030388048,
                                                                                        0.836030388052,
                                                                                        0.803810489627,
                                                                                        0.803810489629,
                                                                                        0.860159122320,
                                                                                        0.860159122320],
                                                            tolerance=0.001),
                    apptest.observable_test.reference_file('Onebody density matrix up', remove_equal_indexes = True,
                                                            file=os.path.join(reference_dir, 'test_fermiladder.diag.h5'), tolerance=0.001),
                    apptest.observable_test.reference_file('Onebody density matrix down', remove_equal_indexes = True,
                                                            file=os.path.join(reference_dir, 'test_fermiladder.diag.h5'), tolerance=0.001),
                  ]


if __name__ == '__main__':
    apptest.main()
