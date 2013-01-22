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
    if len(data) != 1:
        raise ObservableNotFound(fname, observable)
    d = data[0]
    if len(d.x) > 1 and d.props['observable'] != 'Entropy':
        x = np.array(d.x)
        if len(x.shape) > 1:
            x = x.reshape( x.shape[0], np.prod(x.shape[1:]) )
            keys = []
            for i in reversed(range(x.shape[1])):
                keys.append(x[:,i])
            ind=np.lexsort(keys)
        else:
            ind = np.argsort(x)
        d.x = d.x[ind]
        for i in range(len(d.y)):
            d.y[i] = d.y[i][ind]
    return d

def load_iterations_observable(fname, observable):
    if not os.path.exists(fname):
        raise IOError('Archive `%s` not found.' % fname)
    ## TODO: load iterations (time evolution) correctly
    data = pyalps.loadEigenstateMeasurements([fname], [observable])
    data = pyalps.flatten(data)
    if len(data) == 0:
        raise ObservableNotFound(fname, observable)
    return data[0]


class reference_file(object):
    def __init__(self, observable, file, tolerance=1e-6, load_type='spectrum'):
        self.observable     = str(observable)
        self.tolerance      = float(tolerance)
        self.reference_file = str(file)
        if load_type == 'spectrum':
            self.loader = load_spectrum_observable
        elif load_type == 'iterations':
            self.loader = load_iterations_observable
        else:
            raise RuntimeError('`%s` not a valid type.' % load_type)
    
    def __call__(self, test_file):
        tobs = self.loader(test_file,           self.observable)
        robs = self.loader(self.reference_file, self.observable)
        for ty, ry in zip(tobs.y, robs.y):
            for tval, rval in zip( np.atleast_1d(ty), np.atleast_1d(ry) ):
                if np.any( np.array(abs(tval-rval)) / rval > self.tolerance ):
                    raise ObservableNotMatch(self.observable, tval, rval, self.tolerance)
    
    def __str__(self):
        return 'Reference file test for `%s` in `%s`' % (self.observable, os.path.basename(self.reference_file))

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
    
    def __call__(self, test_file):
        obs = self.loader(test_file, self.observable)
        for tval, rval in zip(np.atleast_1d(obs.y[0]), self.value):
            if abs(tval-rval) / rval > self.tolerance:
                raise ObservableNotMatch(self.observable, tval, rval, self.tolerance)
    
    def __str__(self):
        return 'Reference value for `%s`' % self.observable

