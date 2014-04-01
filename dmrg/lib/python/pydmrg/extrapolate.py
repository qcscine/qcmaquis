#############################################################################
#
# MAQUIS DMRG Project
#
# Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
#
#############################################################################

from copy import deepcopy
import numpy as np
from scipy import polyfit
import pyalps
import pyalps.fit_wrapper as fw

import tools

class extrapolator(object):
    
    def __init__(self, sets, x, obs):
        self.props = pyalps.dict_intersect( [d.props for d in pyalps.flatten(sets)] )
        self.props['observable'] = str(obs)
        self.xname = str(x)
        self.obsx  = []
        self.ydata = np.empty((0,0))
        self.xdata = np.empty(0)
        
        self.init_data(sets)
        self.xdata = np.array(self.xdata)
        self.ydata = np.array(self.ydata)
        
        order = np.argsort(self.xdata)
        self.xdata = self.xdata[order]
        for i in range(len(self.ydata)):
            self.ydata[i] = self.ydata[i][order]
    
    def init_data(self, sets):
        for i,q in enumerate(pyalps.flatten(sets)):
            # sorting according to q.x
            if len(q.x) > 0:
                order = tools.labels_argsort(q.x)
                qx = q.x[order]
                qy = q.y[order]
            else:
                qx = []
                qy = np.array(q.y)
            
            if i == 0:
                self.obsx = q.x
                self.ydata = np.empty((len(qy),0))
            elif not np.all(abs(self.obsx-q.x) < 1e-8):
                raise Exception("Observable `x` values don't match!")
            
            self.xdata = np.concatenate ( (self.xdata, [q.props[self.xname]]) )
            self.ydata = np.column_stack( (self.ydata, qy) )
    
    def fit_at(self, x, y, degree, xval, num_points=None):
        if num_points is None:
            fit_parms = polyfit(x, y, degree)
        else:
            fit_parms = polyfit(x[:num_points], y[:num_points], degree)
        
        val = xval*0.
        for deg in range(0,degree+1):
            val += fit_parms[deg] * xval**(degree-deg)
        
        return val
    
    def evaluate(self, xval=0, degree=1, num_points=None):
        q = pyalps.DataSet()
        q.props = deepcopy(self.props)
        q.props[self.xname] = xval
        q.x = self.obsx
        q.y = np.array([ self.fit_at(x=self.xdata, y=y, degree=degree, xval=xval, num_points=num_points) for y in self.ydata ])
        
        return q
    
    def fit_data(self, xval=0, degree=1, ii=0, num_points=None):
        q = pyalps.DataSet()
        q.props = deepcopy(self.props)
        q.x = self.xdata
        q.y = self.ydata[ii]
        q.props['ylabel'] = self.props['observable']
        q.props['xlabel'] = self.xname
        
        q.props['fitting_degree']     = degree
        q.props['fitting_num_points'] = num_points
        q.props['fitted_x']     = xval
        q.props['fitted_value'] = self.fit_at(x=q.x, y=q.y, degree=degree, xval=xval, num_points=num_points)
        
        return q
    
