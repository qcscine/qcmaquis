#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

import sys,os
from os.path import expanduser
import tempfile
import shutil
from datetime import datetime as dt
import subprocess
from subprocess import check_call
import pyalps

from exception import TestFailed

def remove_if_exists(fname):
    """Remove file if exists."""
    if os.path.exists(fname):
        os.remove(fname)

def exec_with_log(cmd, logbase):
    """
    Run the command redirecting stdout and stderr to logbase.log and
    logbase.err.log, respectively. In case of failure the content of
    stderr is printed before raising CalledProcessError exception."""
    fp     = open(logbase+'.log', 'w')
    fp_err = open(logbase+'.err.log', 'w+')
    try:
        check_call(cmd, stdout=fp, stderr=fp_err)
    except subprocess.CalledProcessError as e:
        fp_err.seek(0)
        print 'The program crashed.\n',fp_err.read()
        raise e
    fp.close()
    fp_err.close()
    

class DMRGTestBase(object):
    """
    Base class for DMRG tests.
    
    List of static properties that need to be defined:
     * testname = str()
       - <parms>: application parameters
       - <model>: model parameters
     * inputs = {'parms': <parms>, 'model': <model>}
       - <parms>: application parameters
       - <model>: model parameters
     * observables = [<instance of `observable_test`>]
       list of observabled to be tested
     * (optional) reference_file = str()
       hdf5 archive with reference results
    """
    
    testname       = 'dmrgtest'
    inputs         = {'parms': {}, 'model': {}}
    observbales    = []
    reference_file = None
    
    def __init__(self):
        self.rundir_setup()
    
    def rundir_setup(self):
        if not hasattr(self, 'tmpdir'):
            self.origdir = os.getcwd()
            self.tmpdir = tempfile.mkdtemp()
            print 'Using temporary dir:', self.tmpdir
    
    def write_parameters(self):
        """Write parameters as `testname.parms/.model`."""
        pyalps.writeParameterFile(self.testname+'.parms', self.inputs['parms'])
        pyalps.writeParameterFile(self.testname+'.model', self.inputs['model'])
    
    def run(self, dmrg_app=None, meas_app=None):
        """Create parameters and execute the simulation."""
        
        self.rundir_setup()
        os.chdir(self.tmpdir)
        
        ## Write parameters
        self.write_parameters()
        
        ## Remove existing result files
        remove_if_exists( self.inputs['parms']['resultfile'] )
        remove_if_exists( self.inputs['parms']['chkpfile']   )
        remove_if_exists( self.testname+'.out.log'           )
        remove_if_exists( self.testname+'.out.meas.log'      )
        
        ## Execute DMRG app
        if dmrg_app is not None:
            cmd = [expanduser(dmrg_app), self.testname+'.parms', self.testname+'.model']
            exec_with_log(cmd, self.testname+'.out')
        
        ## Execute measure app
        if meas_app is not None:
            cmd = [expanduser(meas_app), self.testname+'.parms', self.testname+'.model']
            exec_with_log(cmd, self.testname+'.out.meas')
        
        os.chdir(self.origdir)
    
    def check(self):
        """Iterate on self.observbales to check results."""
        resfile = self.inputs['parms']['resultfile']
        
        self.rundir_setup()
        os.chdir(self.tmpdir)
        
        ## Compare observables
        passed = True
        for obstest in self.observables:
            try:
                obstest( os.path.join(self.tmpdir, resfile) )
                print obstest, 'Success!'
            except TestFailed as err:
                passed = False
                print obstest, 'Failed!'
                print 'Details:', err
                    
        os.chdir(self.origdir)
        
        if not passed:
            archive_dirname = self.testname+'_failed.%s' % dt.now().strftime('%Y%m%d%H%M%S')
            basedir = os.path.dirname(self.tmpdir)
            shutil.move( self.tmpdir, os.path.join(basedir, archive_dirname) )
            archive_file = shutil.make_archive(archive_dirname, format='gztar',
                                               root_dir=basedir, base_dir=archive_dirname)
            print 'Failed test archived in %s.' % archive_file
            
            raise TestFailed('Some tests failed!')
            
    
    def create(self, dmrg_app=None, meas_app=None):
        """Create reference file."""
        if self.reference_file is None:
            raise RuntimeError('`refence_file` has to be speficied.')
        
        self.run(dmrg_app=dmrg_app, meas_app=meas_app)
        shutil.copy(os.path.join(self.tmpdir, self.inputs['parms']['resultfile']),
                    self.reference_file)

