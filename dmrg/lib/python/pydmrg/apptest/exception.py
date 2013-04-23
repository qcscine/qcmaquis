#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

class TestFailed(Exception):
    def __init__(self, msg):
        self.msg = str(msg)
    def __str__(self):
        return self.msg

class ParametersNotMatch(TestFailed):
    def __init__(self, out_file, ref_file):
        self.out_file = str(out_file)
        self.ref_file = str(ref_file)
    def __str__(self):
        return 'Parameters in `%s` don\'t match reference file `%s`' % (self.out_file, self.ref_file)

class ObservableNotFound(TestFailed):
    def __init__(self, fname, obs):
        self.fname = str(fname)
        self.obs   = str(obs)
    def __str__(self):
        return 'Observable `%s` not found in `%s`.' % (self.obs, self.fname)

class ObservableNotMatch(TestFailed):
    def __init__(self, obs, tval, rval, tol):
        self.obs  = str(obs)
        self.tval = float(tval)
        self.rval = float(rval)
        self.tol  = float(tol)
    def __str__(self):
        return 'Observable `%s` does not match. `%s` more than %s different compared to %s.' % (self.obs, self.tval, self.tol, self.rval)
