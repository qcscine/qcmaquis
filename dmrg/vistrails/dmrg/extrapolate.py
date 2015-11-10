#############################################################################
#
# MAQUIS DMRG Project
# Vistrails package
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

import pyalps
import pydmrg
import numpy as np
from  copy import deepcopy

from core.modules.vistrails_module import Module
import core.modules.basic_modules as basic

from packages.alps.dataset import DataSets

class ExtrapolationData(Module):
    _input_ports  = []
    _output_ports = []
    
    def __init__(self, indata, x, obs, degree, num_points):
        Module.__init__(self)
        
        self.x          = str(x)
        self.obs        = str(obs)
        self.degree     = int(degree)
        self.num_points = int(num_points)
        if self.num_points == 0:
            self.num_points = None
        
        self.extrapolator = pydmrg.extrapolator(indata, x=self.x, obs=self.obs)
        

class ExtrapolateMeasurement(Module):
    
    _input_ports =  [
                        ('input',       DataSets),
                        ('x',           basic.String),
                        ('observable',  basic.String),
                        ('xval',        basic.Float),
                        ('degree',      basic.Integer, {'defaults': str([1]), 'optional': True}),
                        ('num_points',  basic.Integer, {'defaults': str([0]), 'optional': True}),
                    ]
    
    _output_ports = [('output', DataSets), ('extrapolation', ExtrapolationData)]
    
    def compute(self):
        extdata = ExtrapolationData( self.getInputFromPort('input'), x=self.getInputFromPort('x'), obs=self.getInputFromPort('observable'),
                                     degree=self.getInputFromPort('degree'), num_points=self.getInputFromPort('num_points') )
        
        xval = self.getInputFromPort('xval')
        
        ret = []
        ret.append( extdata.extrapolator.evaluate(xval=xval, degree=extdata.degree, num_points=extdata.num_points) )
        
        self.setResult('output', ret)
        self.setResult('extrapolation', extdata)


class ExtractExtrapolation(Module):
    
    _input_ports =  [
                        ('input',        ExtrapolationData),
                        ('fitting_at',   basic.Integer, {'defaults': str([0]), 'optional': True}),
                        ('fitting_xval', basic.List),
                    ]
    
    _output_ports = [
                        ('data',  DataSets),
                        ('fit',   DataSets),
                    ]
    
    def compute(self):
        extdata = self.getInputFromPort('input')
        fitting_at   = self.getInputFromPort('fitting_at')
        fitting_xval = self.getInputFromPort('fitting_xval')
        
        data = extdata.extrapolator.fit_data(xval=np.array(fitting_xval), ii=fitting_at, degree=extdata.degree, num_points=extdata.num_points)
        
        fitted = pyalps.DataSet()
        
        fitted.props = deepcopy(data.props)
        fitted.props['observable'] += ' - extrapolation(degree=%s, %s points)' % (fitted.props['fitting_degree'], fitted.props['fitting_num_points'])
        
        fitted.x = fitted.props['fitted_x']
        fitted.y = fitted.props['fitted_value']
        
        self.setResult('data', [data])
        self.setResult('fit',  [fitted])

