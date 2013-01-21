#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

import os
import pyalps
import numpy as np

from exception import ObservableNotFound
from exception import ObservableNotMatch

def load_spectrum_observable(fname, observable):
    if not os.path.exists(fname):
        raise IOError('Archive `%s` not found.' % fname)
    data = pyalps.loadEigenstateMeasurements([fname], [observable])
    data = pyalps.flatten(data)
    if len(data) == 0:
        raise ObservableNotFound(fname, observbale)
    return data[0]

def load_iterations_observable(fname, observable):
    if not os.path.exists(fname):
        raise IOError('Archive `%s` not found.' % fname)
    ## TODO: load iterations (time evolution) correctly
    data = pyalps.loadEigenstateMeasurements([fname], [observable])
    data = pyalps.flatten(data)
    if len(data) == 0:
        raise ObservableNotFound(fname, observbale)
    return data[0]


class reference_run(object):
    def __init__(self, observable, tolerance=1e-6, load_type='spectrum'):
        self.observable = str(observable)
        self.tolerance  = float(tolerance)
        if load_type == 'spectrum':
            self.loader = load_spectrum_observable
        elif load_type == 'iterations':
            self.loader = load_iterations_observable
        else:
            raise RuntimeError('`%s` not a valid type.' % load_type)
    
    def __call__(self, test_file, reference_file):
        tobs = self.loader(test_file,      self.observable)
        robs = self.loader(reference_file, self.observable)
        for tval, rval in zip( np.atleast_1d(tobs.y[0]), np.atleast_1d(robs.y[0]) ):
            if np.any( np.array(abs(tval-rval)) / rval > self.tolerance ):
                raise ObservableNotMatch(self.observable, tval, rval, self.tolerance)
    
    def __str__(self):
        return 'Reference run of `%s`' % self.observable

class reference_value(object):
    def __init__(self, observable, value, tolerance=1e-6, load_type='spectrum'):
        self.observable = str(observable)
        self.value      = np.atleast_1d(value)
        self.tolerance  = float(tolerance)
        if load_type == 'spectrum':
            self.loader = load_spectrum_observable
        elif load_type == 'iterations':
            self.loader = load_iterations_observable
        else:
            raise RuntimeError('`%s` not a valid type.' % load_type)
    
    def __call__(self, test_file, reference_file=None):
        obs = self.loader(test_file, self.observable)
        for tval, rval in zip(np.atleast_1d(obs.y[0]), self.value):
            if abs(tval-rval) / rval > self.tolerance:
                raise ObservableNotMatch(self.observable, tval, rval, self.tolerance)
    
    def __str__(self):
        return 'Reference value of `%s`' % self.observable

