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
import core.modules.basic_modules
basic = core.modules.basic_modules

from packages.alps.dataset import DataSets


class LoadTruncatedWeight(Module):
    """
    Load the maximum TruncatedWeight in the last optimization sweep
    from the `filename` property for each DataSet in input.
    """
    _input_ports =  [('input', DataSets)]
    _output_ports = [('output', DataSets)]
    
    def compute(self):
        sets = self.getInputFromPort('input')
        pydmrg.LoadLastTruncation(sets)
        self.setResult('output', sets)
