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

from core.modules.vistrails_module import Module
import core.modules.basic_modules as basic

from packages.alps.dataset import DataSets


class ExtrapolateMeasurement(Module):
    
    _input_ports =  [
                        ('input',       DataSets),
                        ('x',           basic.String),
                        ('observable',  basic.String),
                        ('xval',        basic.Float),
                        ('degree',      basic.Integer, {'defaults': str([1]), 'optional': True}),
                        ('num_points',  basic.Integer, {'defaults': str([0]), 'optional': True}),
                    ]
    
    _output_ports = [('output', DataSets)]
    
    def compute(self):
        ext = pydmrg.extrapolator(self.getInputFromPort('input'), x=self.getInputFromPort('x'), obs=self.getInputFromPort('observable'))
        
        xval       = self.getInputFromPort('xval')
        degree     = self.getInputFromPort('degree')
        num_points = self.getInputFromPort('num_points')
        if num_points == 0:
            num_points = None
        
        ret = []
        ret.append( ext.evaluate(xval=xval, degree=degree, num_points=num_points) )
        
        self.setResult('output', ret)
